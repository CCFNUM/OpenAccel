// File       : turbulentIntermittencyTransitionSSTAssemblerNodeTerms.cpp
// Created    : Wed Jan 15 2025
// Author     : Adam Fares
// Description:
// Copyright (c) 2025 CCFNUM, Lucerne University of Applied Sciences and Arts.
// SPDX-License-Identifier: BSD-3-Clause

// code
#include "turbulentIntermittencyTransitionSSTAssembler.h"

namespace accel
{

void turbulentIntermittencyTransitionSSTAssembler::
    assembleNodeTermsFusedSteady_(const domain* domain, Context* ctx)
{
    assert(phi_);

    const auto& mesh = model_->meshRef();
    Matrix& A = ctx->getAMatrix();
    Vector& b = ctx->getBVector();

    const stk::mesh::BulkData& bulkData = mesh.bulkDataRef();
    const stk::mesh::MetaData& metaData = mesh.metaDataRef();

    // space for LHS/RHS
    const label lhsSize = 1;
    const label rhsSize = 1;
    std::vector<scalar> lhs(lhsSize);
    std::vector<scalar> rhs(rhsSize);
    std::vector<label> scratchIds(rhsSize);
    std::vector<scalar> scratchVals(rhsSize);
    std::vector<stk::mesh::Entity> connectedNodes(1);

    // pointers
    scalar* p_lhs = &lhs[0];
    scalar* p_rhs = &rhs[0];

    // Get fields
    const STKScalarField* gammaSTKFieldPtr = phi_->stkFieldPtr();
    const STKScalarField* ReThetaSTKFieldPtr =
        model_->ReThetaRef().stkFieldPtr();

    const STKScalarField* omegaSTKFieldPtr = model_->omegaRef().stkFieldPtr();
    const STKScalarField* kSTKFieldPtr = model_->kRef().stkFieldPtr();

    const STKScalarField* rhoSTKFieldPtr = model_->rhoRef().stkFieldPtr();
    const STKScalarField* divSTKFieldPtr =
        model_->mDotRef().divRef().stkFieldPtr();
    const STKScalarField* muSTKFieldPtr = model_->muRef().stkFieldPtr();
    const STKScalarField* gradUSTKFieldPtr =
        model_->URef().gradRef().stkFieldPtr();

    const STKScalarField* yMinSTKFieldPtr = model_->yMinRef().stkFieldPtr();

    scalar ca1 = model_->ca1();
    scalar ca2 = model_->ca2();
    scalar ce1 = model_->ce1();
    scalar ce2 = model_->ce2();

    // Geometric fields
    const auto* volSTKFieldPtr = metaData.get_field<scalar>(
        stk::topology::NODE_RANK, this->getDualNodalVolumeID_(domain));

    // other
    scalar dt = model_->controlsRef().getPhysicalTimescale();

    // get interior parts the domain is defined on
    const stk::mesh::PartVector& partVec = domain->zonePtr()->interiorParts();

    // define some common selectors; select owned nodes
    stk::mesh::Selector selOwnedNodes =
        metaData.locally_owned_part() & stk::mesh::selectUnion(partVec);

    stk::mesh::BucketVector const& nodeBuckets =
        bulkData.get_buckets(stk::topology::NODE_RANK, selOwnedNodes);
    for (stk::mesh::BucketVector::const_iterator ib = nodeBuckets.begin();
         ib != nodeBuckets.end();
         ++ib)
    {
        stk::mesh::Bucket& nodeBucket = **ib;

        const stk::mesh::Bucket::size_type nNodesPerBucket = nodeBucket.size();

        for (stk::mesh::Bucket::size_type iNode = 0; iNode < nNodesPerBucket;
             ++iNode)
        {
            // get node
            stk::mesh::Entity node = nodeBucket[iNode];
            connectedNodes[0] = node;

            for (label i = 0; i < lhsSize; ++i)
            {
                p_lhs[i] = 0.0;
            }
            for (label i = 0; i < rhsSize; ++i)
            {
                p_rhs[i] = 0.0;
            }

            // get values of current node
            scalar gamma = *stk::mesh::field_data(*gammaSTKFieldPtr, node);
            scalar k = *stk::mesh::field_data(*kSTKFieldPtr, node);
            scalar omega = *stk::mesh::field_data(*omegaSTKFieldPtr, node);
            scalar ReTheta = *stk::mesh::field_data(*ReThetaSTKFieldPtr, node);
            scalar rho = *stk::mesh::field_data(*rhoSTKFieldPtr, node);
            scalar mu = *stk::mesh::field_data(*muSTKFieldPtr, node);
            scalar vol = *stk::mesh::field_data(*volSTKFieldPtr, node);
            scalar div = *stk::mesh::field_data(*divSTKFieldPtr, node);
            scalar y = *stk::mesh::field_data(*yMinSTKFieldPtr, node);
            scalar* dudx = stk::mesh::field_data(*gradUSTKFieldPtr, node);

            // false transient
            scalar lhsfac = rho * vol / dt;
            lhs[0] += lhsfac;

            // divergence correction
            if (div < 0)
            {
                lhs[0] += -div;
            }
            rhs[0] -= -div * gamma;

            // compute strain rate magnitude; pull pointer within the loop to
            // make it managable
            scalar sijMag = 0.0;
            for (label i = 0; i < SPATIAL_DIM; ++i)
            {
                for (label j = 0; j < SPATIAL_DIM; ++j)
                {
                    const scalar rateOfStrain =
                        0.5 *
                        (dudx[SPATIAL_DIM * i + j] + dudx[SPATIAL_DIM * j + i]);
                    sijMag += rateOfStrain * rateOfStrain;
                }
            }
            sijMag = std::sqrt(2.0 * sijMag);

            // Comute Volrticity Mag Wij
            scalar wijMag = 0.0;
            for (label i = 0; i < SPATIAL_DIM; ++i)
            {
                for (label j = 0; j < SPATIAL_DIM; ++j)
                {
                    const scalar rateOfStrain =
                        0.5 *
                        (dudx[SPATIAL_DIM * i + j] - dudx[SPATIAL_DIM * j + i]);
                    wijMag += rateOfStrain * rateOfStrain;
                }
            }
            wijMag = std::sqrt(2.0 * wijMag);

            scalar Rev = (rho * sijMag * y * y) / mu;

            scalar ReThetac;
            if (ReTheta <= 1870.0)
            {
                ReThetac =
                    (-396.035e-2) + (10120.656e-4) * ReTheta +
                    (-868.230e-6) * ReTheta * ReTheta +
                    (696.506e-9) * ReTheta * ReTheta * ReTheta +
                    (-174.105e-12) * ReTheta * ReTheta * ReTheta * ReTheta;
            }
            else if (ReTheta > 1870.0)
            {
                ReThetac = ReTheta - (593.11 + 0.482 * (ReTheta - 1870.0));
            }

            scalar fOnset1 = Rev / (2.193 * ReThetac + SMALL);
            scalar fOnset2 =
                std::min(std::max(fOnset1, std::pow(fOnset1, 4.0)), 2.0);
            scalar Rt = (rho * k) / (mu * omega + SMALL);
            scalar fOnset3 = std::max((1.0 - std::pow((Rt / 2.5), 3.0)), 0.0);
            scalar fOnset = std::max((fOnset2 - fOnset3), 0.0);
            scalar fTurb = std::exp(-1.0 * std::pow((Rt / 4.0), 4.0));

            scalar fLength1;
            if (ReTheta < 400.0)
            {
                fLength1 = 39.8189 + (-119.270e-4) * ReTheta +
                           (-132.567e-6) * ReTheta * ReTheta;
            }
            else if (ReTheta >= 400.0 && ReTheta < 596.0)
            {
                fLength1 = 263.404 + (-123.939e-2) * ReTheta +
                           (194.548e-5) * ReTheta * ReTheta +
                           (-101.695e-8) * ReTheta * ReTheta * ReTheta;
            }
            else if (ReTheta >= 596.0 && ReTheta < 1200.0)
            {
                fLength1 = 0.5 - (3.0e-4) * (ReTheta - 596.0);
            }
            else if (ReTheta >= 1200.0)
            {
                fLength1 = 0.3188;
            }

            scalar Rew = (rho * omega * y * y) / mu;
            scalar fSublayer = std::exp(-1.0 * std::pow((Rew / 200.0), 2.0));
            scalar fLength = fLength1 * (1.0 - fSublayer) + 40.0 * fSublayer;

            // Production
            const scalar Pgamma = fLength * ca1 * rho * sijMag *
                                  std::sqrt(gamma * fOnset) *
                                  (1.0 - ce1 * gamma);

            // Dissipation
            const scalar Egamma =
                ca2 * rho * wijMag * gamma * fTurb * (ce2 * gamma - 1.0);

            const scalar linearizedProductionFlux =
                fLength * ca1 * rho * sijMag * std::sqrt(gamma * fOnset) * ce1;
            const scalar linearizedDissipationFlux =
                ca2 * rho * wijMag * fTurb * ce2 * gamma;

            rhs[0] += (Pgamma - Egamma) * vol;
            lhs[0] += std::max(0.0, linearizedProductionFlux * vol) +
                      std::max(0.0, linearizedDissipationFlux * vol);

            Base::applyCoeff_(
                A, b, connectedNodes, scratchIds, scratchVals, rhs, lhs);
        }
    }
}

void turbulentIntermittencyTransitionSSTAssembler::
    assembleNodeTermsFusedFirstOrderUnsteady_(const domain* domain,
                                              Context* ctx)
{
    assert(phi_);

    const auto& mesh = model_->meshRef();
    Matrix& A = ctx->getAMatrix();
    Vector& b = ctx->getBVector();

    const stk::mesh::BulkData& bulkData = mesh.bulkDataRef();
    const stk::mesh::MetaData& metaData = mesh.metaDataRef();

    const bool meshDeforming = domain->zonePtr()->meshDeforming();

    // space for LHS/RHS
    const label lhsSize = 1;
    const label rhsSize = 1;
    std::vector<scalar> lhs(lhsSize);
    std::vector<scalar> rhs(rhsSize);
    std::vector<label> scratchIds(rhsSize);
    std::vector<scalar> scratchVals(rhsSize);
    std::vector<stk::mesh::Entity> connectedNodes(1);

    // pointers
    scalar* p_lhs = &lhs[0];
    scalar* p_rhs = &rhs[0];
    // Get fields
    const STKScalarField* gammaSTKFieldPtr = phi_->stkFieldPtr();
    const STKScalarField* gammaSTKFieldPtrOld =
        phi_->prevTimeRef().stkFieldPtr();
    const STKScalarField* ReThetaSTKFieldPtr =
        model_->ReThetaRef().stkFieldPtr();
    const STKScalarField* omegaSTKFieldPtr = model_->omegaRef().stkFieldPtr();
    const STKScalarField* kSTKFieldPtr = model_->kRef().stkFieldPtr();
    const STKScalarField* rhoSTKFieldPtr = model_->rhoRef().stkFieldPtr();
    const STKScalarField* rhoSTKFieldPtrOld =
        model_->rhoRef().prevTimeRef().stkFieldPtr();
    const STKScalarField* divSTKFieldPtr =
        model_->mDotRef().divRef().stkFieldPtr();
    const STKScalarField* muSTKFieldPtr = model_->muRef().stkFieldPtr();
    const STKScalarField* gradUSTKFieldPtr =
        model_->URef().gradRef().stkFieldPtr();
    const STKScalarField* yMinSTKFieldPtr = model_->yMinRef().stkFieldPtr();

    scalar ca1 = model_->ca1();
    scalar ca2 = model_->ca2();
    scalar ce1 = model_->ce1();
    scalar ce2 = model_->ce2();

    const STKScalarField* divUmSTKFieldPtr =
        meshDeforming ? model_->divUmRef().stkFieldPtr() : nullptr;

    // Geometric fields
    const auto* volSTKFieldPtr = metaData.get_field<scalar>(
        stk::topology::NODE_RANK, this->getDualNodalVolumeID_(domain));

    // time integrator
    const scalar dt = model_->meshRef().controlsRef().getTimestep();
    const auto c = BDF1::coeff();

    // get interior parts the domain is defined on
    const stk::mesh::PartVector& partVec = domain->zonePtr()->interiorParts();

    // define some common selectors; select owned nodes
    stk::mesh::Selector selOwnedNodes =
        metaData.locally_owned_part() & stk::mesh::selectUnion(partVec);

    stk::mesh::BucketVector const& nodeBuckets =
        bulkData.get_buckets(stk::topology::NODE_RANK, selOwnedNodes);
    for (stk::mesh::BucketVector::const_iterator ib = nodeBuckets.begin();
         ib != nodeBuckets.end();
         ++ib)
    {
        stk::mesh::Bucket& nodeBucket = **ib;

        const stk::mesh::Bucket::size_type nNodesPerBucket = nodeBucket.size();

        for (stk::mesh::Bucket::size_type iNode = 0; iNode < nNodesPerBucket;
             ++iNode)
        {
            // get node
            stk::mesh::Entity node = nodeBucket[iNode];
            connectedNodes[0] = node;

            for (label i = 0; i < lhsSize; ++i)
            {
                p_lhs[i] = 0.0;
            }
            for (label i = 0; i < rhsSize; ++i)
            {
                p_rhs[i] = 0.0;
            }

            // get values of current node
            scalar gammaOld =
                *stk::mesh::field_data(*gammaSTKFieldPtrOld, node);
            scalar gamma = *stk::mesh::field_data(*gammaSTKFieldPtr, node);
            scalar ReTheta = *stk::mesh::field_data(*ReThetaSTKFieldPtr, node);
            scalar k = *stk::mesh::field_data(*kSTKFieldPtr, node);
            scalar omega = *stk::mesh::field_data(*omegaSTKFieldPtr, node);
            scalar rho = *stk::mesh::field_data(*rhoSTKFieldPtr, node);
            scalar rhoOld = *stk::mesh::field_data(*rhoSTKFieldPtrOld, node);
            scalar mu = *stk::mesh::field_data(*muSTKFieldPtr, node);
            scalar vol = *stk::mesh::field_data(*volSTKFieldPtr, node);
            scalar div = *stk::mesh::field_data(*divSTKFieldPtr, node);
            scalar y = *stk::mesh::field_data(*yMinSTKFieldPtr, node);
            scalar* dudx = stk::mesh::field_data(*gradUSTKFieldPtr, node);

            // divergence correction
            if (div < 0)
            {
                lhs[0] += -div;
            }
            rhs[0] -= -div * gamma;

            // transient
            scalar lhsfac = c[0] * rho * vol / dt;
            scalar lhsfacOld = c[1] * rhoOld * vol / dt;

            lhs[0] += lhsfac;
            rhs[0] -= (lhsfac * gamma + lhsfacOld * gammaOld);

            // compute strain rate magnitude; pull pointer within the loop to
            // make it managable
            scalar sijMag = 0.0;
            for (label i = 0; i < SPATIAL_DIM; ++i)
            {
                for (label j = 0; j < SPATIAL_DIM; ++j)
                {
                    const scalar rateOfStrain =
                        0.5 *
                        (dudx[SPATIAL_DIM * i + j] + dudx[SPATIAL_DIM * j + i]);
                    sijMag += rateOfStrain * rateOfStrain;
                }
            }
            sijMag = std::sqrt(2.0 * sijMag);

            // compute vorticity mag
            scalar wijMag = 0.0;
            for (label i = 0; i < SPATIAL_DIM; ++i)
            {
                for (label j = 0; j < SPATIAL_DIM; ++j)
                {
                    const scalar rateOfStrain =
                        0.5 *
                        (dudx[SPATIAL_DIM * i + j] - dudx[SPATIAL_DIM * j + i]);
                    wijMag += rateOfStrain * rateOfStrain;
                }
            }
            wijMag = std::sqrt(2.0 * wijMag);

            scalar Rev = (rho * sijMag * y * y) / mu;

            scalar ReThetac;
            if (ReTheta <= 1870.0)
            {
                ReThetac =
                    (-396.035e-2) + (10120.656e-4) * ReTheta +
                    (-868.230e-6) * ReTheta * ReTheta +
                    (696.506e-9) * ReTheta * ReTheta * ReTheta +
                    (-174.105e-12) * ReTheta * ReTheta * ReTheta * ReTheta;
            }
            else if (ReTheta > 1870.0)
            {
                ReThetac = ReTheta - (593.11 + 0.482 * (ReTheta - 1870.0));
            }

            scalar fOnset1 = Rev / (2.193 * ReThetac + SMALL);
            scalar fOnset2 =
                std::min(std::max(fOnset1, std::pow(fOnset1, 4.0)), 2.0);
            scalar Rt = (rho * k) / (mu * omega + SMALL);
            scalar fOnset3 = std::max((1.0 - std::pow((Rt / 2.5), 3.0)), 0.0);
            scalar fOnset = std::max((fOnset2 - fOnset3), 0.0);
            scalar fTurb = std::exp(-1.0 * std::pow((Rt / 4.0), 4.0));

            scalar fLength1;
            if (ReTheta < 400.0)
            {
                fLength1 = 39.8189 + (-119.270e-4) * ReTheta +
                           (-132.567e-6) * ReTheta * ReTheta;
            }
            else if (ReTheta >= 400.0 && ReTheta < 596.0)
            {
                fLength1 = 263.404 + (-123.939e-2) * ReTheta +
                           (194.548e-5) * ReTheta * ReTheta +
                           (-101.695e-8) * ReTheta * ReTheta * ReTheta;
            }
            else if (ReTheta >= 596.0 && ReTheta < 1200.0)
            {
                fLength1 = 0.5 - (3.0e-4) * (ReTheta - 596.0);
            }
            else if (ReTheta >= 1200.0)
            {
                fLength1 = 0.3188;
            }

            scalar Rew = (rho * omega * y * y) / mu;
            scalar fSublayer = std::exp(-1.0 * std::pow((Rew / 200.0), 2.0));
            scalar fLength = fLength1 * (1.0 - fSublayer) + 40.0 * fSublayer;

            // Production
            const scalar Pgamma = fLength * ca1 * rho * sijMag *
                                  std::sqrt(gamma * fOnset) *
                                  (1.0 - ce1 * gamma);

            // Dissipation
            const scalar Egamma =
                ca2 * rho * wijMag * gamma * fTurb * (ce2 * gamma - 1.0);

            const scalar linearizedProductionFlux =
                fLength * ca1 * rho * sijMag * std::sqrt(gamma * fOnset) * ce1;
            const scalar linearizedDissipationFlux =
                ca2 * rho * wijMag * fTurb * ce2 * gamma;

            rhs[0] += (Pgamma - Egamma) * vol;
            lhs[0] += std::max(0.0, linearizedProductionFlux * vol) +
                      std::max(0.0, linearizedDissipationFlux * vol);

            // geometric conservative law
            if (meshDeforming)
            {
                scalar divUm = *stk::mesh::field_data(*divUmSTKFieldPtr, node);
                rhs[0] -= rho * gamma * divUm * vol;
            }

            Base::applyCoeff_(
                A, b, connectedNodes, scratchIds, scratchVals, rhs, lhs);
        }
    }
}

void turbulentIntermittencyTransitionSSTAssembler::
    assembleNodeTermsFusedSecondOrderUnsteady_(const domain* domain,
                                               Context* ctx)
{
    assert(phi_);

    const auto& mesh = model_->meshRef();
    Matrix& A = ctx->getAMatrix();
    Vector& b = ctx->getBVector();

    const stk::mesh::BulkData& bulkData = mesh.bulkDataRef();
    const stk::mesh::MetaData& metaData = mesh.metaDataRef();

    const bool meshDeforming = domain->zonePtr()->meshDeforming();

    // space for LHS/RHS
    const label lhsSize = 1;
    const label rhsSize = 1;
    std::vector<scalar> lhs(lhsSize);
    std::vector<scalar> rhs(rhsSize);
    std::vector<label> scratchIds(rhsSize);
    std::vector<scalar> scratchVals(rhsSize);
    std::vector<stk::mesh::Entity> connectedNodes(1);

    // pointers
    scalar* p_lhs = &lhs[0];
    scalar* p_rhs = &rhs[0];

    // Get fields
    const STKScalarField* gammaSTKFieldPtr = phi_->stkFieldPtr();
    const STKScalarField* gammaSTKFieldPtrOld =
        phi_->prevTimeRef().stkFieldPtr();
    const STKScalarField* gammaSTKFieldPtrOldOld =
        phi_->prevTimeRef().prevTimeRef().stkFieldPtr();
    const STKScalarField* ReThetaSTKFieldPtr =
        model_->ReThetaRef().stkFieldPtr();
    const STKScalarField* omegaSTKFieldPtr = model_->omegaRef().stkFieldPtr();
    const STKScalarField* kSTKFieldPtr = model_->kRef().stkFieldPtr();
    const STKScalarField* rhoSTKFieldPtr = model_->rhoRef().stkFieldPtr();
    const STKScalarField* rhoSTKFieldPtrOld =
        model_->rhoRef().prevTimeRef().stkFieldPtr();
    const STKScalarField* rhoSTKFieldPtrOldOld =
        model_->rhoRef().prevTimeRef().prevTimeRef().stkFieldPtr();
    const STKScalarField* divSTKFieldPtr =
        model_->mDotRef().divRef().stkFieldPtr();
    const STKScalarField* muSTKFieldPtr = model_->muRef().stkFieldPtr();
    const STKScalarField* gradUSTKFieldPtr =
        model_->URef().gradRef().stkFieldPtr();
    const STKScalarField* yMinSTKFieldPtr = model_->yMinRef().stkFieldPtr();

    const STKScalarField* divUmSTKFieldPtr =
        meshDeforming ? model_->divUmRef().stkFieldPtr() : nullptr;

    scalar ce1 = model_->ce1();
    scalar ce2 = model_->ce2();
    scalar ca1 = model_->ca1();
    scalar ca2 = model_->ca2();

    // Geometric fields
    const auto* volSTKFieldPtr = metaData.get_field<scalar>(
        stk::topology::NODE_RANK, this->getDualNodalVolumeID_(domain));

    // time integrator
    const scalar dt = model_->meshRef().controlsRef().getTimestep();
    const auto c = BDF2::coeff(dt, mesh.controlsRef().getTimestep(-1));

    // get interior parts the domain is defined on
    const stk::mesh::PartVector& partVec = domain->zonePtr()->interiorParts();

    // define some common selectors; select owned nodes
    stk::mesh::Selector selOwnedNodes =
        metaData.locally_owned_part() & stk::mesh::selectUnion(partVec);

    stk::mesh::BucketVector const& nodeBuckets =
        bulkData.get_buckets(stk::topology::NODE_RANK, selOwnedNodes);
    for (stk::mesh::BucketVector::const_iterator ib = nodeBuckets.begin();
         ib != nodeBuckets.end();
         ++ib)
    {
        stk::mesh::Bucket& nodeBucket = **ib;

        const stk::mesh::Bucket::size_type nNodesPerBucket = nodeBucket.size();

        for (stk::mesh::Bucket::size_type iNode = 0; iNode < nNodesPerBucket;
             ++iNode)
        {
            // get node
            stk::mesh::Entity node = nodeBucket[iNode];
            connectedNodes[0] = node;

            for (label i = 0; i < lhsSize; ++i)
            {
                p_lhs[i] = 0.0;
            }
            for (label i = 0; i < rhsSize; ++i)
            {
                p_rhs[i] = 0.0;
            }

            // get values of current node
            scalar gamma = *stk::mesh::field_data(*gammaSTKFieldPtr, node);
            scalar ReTheta = *stk::mesh::field_data(*ReThetaSTKFieldPtr, node);
            scalar k = *stk::mesh::field_data(*kSTKFieldPtr, node);
            scalar omega = *stk::mesh::field_data(*omegaSTKFieldPtr, node);

            scalar gammaOld =
                *stk::mesh::field_data(*gammaSTKFieldPtrOld, node);
            scalar gammaOldOld =
                *stk::mesh::field_data(*gammaSTKFieldPtrOldOld, node);

            scalar rho = *stk::mesh::field_data(*rhoSTKFieldPtr, node);
            scalar rhoOld = *stk::mesh::field_data(*rhoSTKFieldPtrOld, node);
            scalar rhoOldOld =
                *stk::mesh::field_data(*rhoSTKFieldPtrOldOld, node);
            scalar mu = *stk::mesh::field_data(*muSTKFieldPtr, node);
            scalar vol = *stk::mesh::field_data(*volSTKFieldPtr, node);
            scalar div = *stk::mesh::field_data(*divSTKFieldPtr, node);

            scalar y = *stk::mesh::field_data(*yMinSTKFieldPtr, node);
            scalar* dudx = stk::mesh::field_data(*gradUSTKFieldPtr, node);

            // divergence correction
            if (div < 0)
            {
                lhs[0] += -div;
            }
            rhs[0] -= -div * gamma;

            // transient
            scalar lhsfac = c[0] * rho * vol / dt;
            scalar lhsfacOld = c[1] * rhoOld * vol / dt;
            scalar lhsfacOldOld = c[2] * rhoOldOld * vol / dt;

            lhs[0] += lhsfac;
            rhs[0] -= (lhsfac * gamma + lhsfacOld * gammaOld +
                       lhsfacOldOld * gammaOldOld);

            // compute strain rate magnitude; pull pointer within the loop to
            // make it managable
            scalar sijMag = 0.0;
            for (label i = 0; i < SPATIAL_DIM; ++i)
            {
                for (label j = 0; j < SPATIAL_DIM; ++j)
                {
                    const scalar rateOfStrain =
                        0.5 *
                        (dudx[SPATIAL_DIM * i + j] + dudx[SPATIAL_DIM * j + i]);
                    sijMag += rateOfStrain * rateOfStrain;
                }
            }
            sijMag = std::sqrt(2.0 * sijMag);

            // compute vorticity mag
            scalar wijMag = 0.0;
            for (label i = 0; i < SPATIAL_DIM; ++i)
            {
                for (label j = 0; j < SPATIAL_DIM; ++j)
                {
                    const scalar rateOfStrain =
                        0.5 *
                        (dudx[SPATIAL_DIM * i + j] - dudx[SPATIAL_DIM * j + i]);
                    wijMag += rateOfStrain * rateOfStrain;
                }
            }
            wijMag = std::sqrt(2.0 * wijMag);

            scalar Rev = (rho * sijMag * y * y) / mu;

            scalar ReThetac;
            if (ReTheta <= 1870.0)
            {
                ReThetac =
                    (-396.035e-2) + (10120.656e-4) * ReTheta +
                    (-868.230e-6) * ReTheta * ReTheta +
                    (696.506e-9) * ReTheta * ReTheta * ReTheta +
                    (-174.105e-12) * ReTheta * ReTheta * ReTheta * ReTheta;
            }
            else if (ReTheta > 1870.0)
            {
                ReThetac = ReTheta - (593.11 + 0.482 * (ReTheta - 1870.0));
            }

            scalar fOnset1 = Rev / (2.193 * ReThetac + SMALL);
            scalar fOnset2 =
                std::min(std::max(fOnset1, std::pow(fOnset1, 4.0)), 2.0);
            scalar Rt = (rho * k) / (mu * omega + SMALL);
            scalar fOnset3 = std::max((1.0 - std::pow((Rt / 2.5), 3.0)), 0.0);
            scalar fOnset = std::max((fOnset2 - fOnset3), 0.0);
            scalar fTurb = std::exp(-1.0 * std::pow((Rt / 4.0), 4.0));

            scalar fLength1;
            if (ReTheta < 400.0)
            {
                fLength1 = 39.8189 + (-119.270e-4) * ReTheta +
                           (-132.567e-6) * ReTheta * ReTheta;
            }
            else if (ReTheta >= 400.0 && ReTheta < 596.0)
            {
                fLength1 = 263.404 + (-123.939e-2) * ReTheta +
                           (194.548e-5) * ReTheta * ReTheta +
                           (-101.695e-8) * ReTheta * ReTheta * ReTheta;
            }
            else if (ReTheta >= 596.0 && ReTheta < 1200.0)
            {
                fLength1 = 0.5 - (3.0e-4) * (ReTheta - 596.0);
            }
            else if (ReTheta >= 1200.0)
            {
                fLength1 = 0.3188;
            }

            scalar Rew = (rho * omega * y * y) / mu;
            scalar fSublayer = std::exp(-1.0 * std::pow((Rew / 200.0), 2.0));
            scalar fLength = fLength1 * (1.0 - fSublayer) + 40.0 * fSublayer;

            // Production
            const scalar Pgamma = fLength * ca1 * rho * sijMag *
                                  std::sqrt(gamma * fOnset) *
                                  (1.0 - ce1 * gamma);

            // Dissipation
            const scalar Egamma =
                ca2 * rho * wijMag * gamma * fTurb * (ce2 * gamma - 1.0);

            const scalar linearizedProductionFlux =
                fLength * ca1 * rho * sijMag * std::sqrt(gamma * fOnset) * ce1;
            const scalar linearizedDissipationFlux =
                ca2 * rho * wijMag * fTurb * ce2 * gamma;

            rhs[0] += (Pgamma - Egamma) * vol;
            lhs[0] += std::max(0.0, linearizedProductionFlux * vol) +
                      std::max(0.0, linearizedDissipationFlux * vol);

            // geometric conservative law
            if (meshDeforming)
            {
                scalar divUm = *stk::mesh::field_data(*divUmSTKFieldPtr, node);
                rhs[0] -= rho * gamma * divUm * vol;
            }

            Base::applyCoeff_(
                A, b, connectedNodes, scratchIds, scratchVals, rhs, lhs);
        }
    }
}

} /* namespace accel */
