// File       : turbulentKineticEnergyTransitionSSTAssemblerNodeTerms.cpp
// Created    : Wed Jan 15 2025
// Author     : Adam Fares
// Description:
// Copyright (c) 2025 CCFNUM, Lucerne University of Applied Sciences and Arts.
// SPDX-License-Identifier: BSD-3-Clause

#include "turbulentKineticEnergyTransitionSSTAssembler.h"

namespace accel
{

void turbulentKineticEnergyTransitionSSTAssembler::
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
    const STKScalarField* kSTKFieldPtr = phi_->stkFieldPtr();
    const STKScalarField* omegaSTKFieldPtr = model_->omegaRef().stkFieldPtr();
    const STKScalarField* ReThetaSTKFieldPtr =
        model_->ReThetaRef().stkFieldPtr();
    const STKScalarField* gammaSTKFieldPtr = model_->gammaRef().stkFieldPtr();
    const STKScalarField* rhoSTKFieldPtr = model_->rhoRef().stkFieldPtr();
    const STKScalarField* muSTKFieldPtr = model_->muRef().stkFieldPtr();
    const STKScalarField* USTKFieldPtr = model_->URef().stkFieldPtr();
    const STKScalarField* gradUSTKFieldPtr =
        model_->URef().gradRef().stkFieldPtr();
    const STKScalarField* divSTKFieldPtr =
        model_->mDotRef().divRef().stkFieldPtr();
    const STKScalarField* yMinSTKFieldPtr = model_->yMinRef().stkFieldPtr();
    const STKScalarField* PkSTKFieldPtr = model_->PkRef().stkFieldPtr();

    scalar betaStar = model_->betaStar();
    scalar tkeProdLimitRatio = model_->tkeProdLimitRatio();
    scalar s1 = model_->s1();
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
            scalar k = *stk::mesh::field_data(*kSTKFieldPtr, node);
            scalar omega = *stk::mesh::field_data(*omegaSTKFieldPtr, node);
            scalar gamma = *stk::mesh::field_data(*gammaSTKFieldPtr, node);
            scalar ReTheta = *stk::mesh::field_data(*ReThetaSTKFieldPtr, node);
            scalar rho = *stk::mesh::field_data(*rhoSTKFieldPtr, node);
            scalar mu = *stk::mesh::field_data(*muSTKFieldPtr, node);
            scalar vol = *stk::mesh::field_data(*volSTKFieldPtr, node);
            scalar* U = stk::mesh::field_data(*USTKFieldPtr, node);
            scalar div = *stk::mesh::field_data(*divSTKFieldPtr, node);
            scalar* dudx = stk::mesh::field_data(*gradUSTKFieldPtr, node);
            scalar y = *stk::mesh::field_data(*yMinSTKFieldPtr, node);
            scalar Pk = *stk::mesh::field_data(*PkSTKFieldPtr, node);

            // S MAG
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

            // false transient
            scalar lhsfac = rho * vol / dt;
            lhs[0] += lhsfac;

            // divergence correction
            if (div < 0)
            {
                lhs[0] += -div;
            }
            rhs[0] -= -div * k;

            // production/dissipation
            scalar Dk = betaStar * rho * omega * k;

            // Compute Volrticity Mag Wij
            scalar vorticityMag = 0.0;
            for (label i = 0; i < SPATIAL_DIM; ++i)
            {
                for (label j = 0; j < SPATIAL_DIM; ++j)
                {
                    const scalar rateOfStrain =
                        0.5 *
                        (dudx[SPATIAL_DIM * i + j] - dudx[SPATIAL_DIM * j + i]);
                    vorticityMag += rateOfStrain * rateOfStrain;
                }
            }
            vorticityMag = std::sqrt(2.0 * vorticityMag);

            // Compute Velocity Mag
            scalar UMag = 0;
            for (label i = 0; i < SPATIAL_DIM; ++i)
            {
                UMag += U[i] * U[i];
            }
            UMag = std::sqrt(UMag);

            scalar Rt = (rho * k) / (mu * omega);
            scalar fReattach = std::exp(-1.0 * std::pow((Rt / 20.0), 4.0));
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

            scalar ReOmega = (rho * omega * y * y) / mu;
            scalar fWake = exp(-1 * std::pow((ReOmega / 1e5), 2.0));

            scalar delta =
                (375.0 * vorticityMag * mu * ReTheta * y) / (rho * UMag * UMag);
            scalar fThetat = std::min(
                std::max(fWake * std::exp(-1.0 * std::pow((y / delta), 4.0)),
                         1.0 -
                             std::pow((ce2 * gamma - 1.0) / (ce2 - 1.0), 2.0)),
                1.0);

            scalar gammaSep =
                std::min(s1 * std::max(0.0, (Rev / (3.235 * ReThetac)) - 1.0) *
                             fReattach,
                         2.0) *
                fThetat;
            scalar gammaEff = std::max(gamma, gammaSep);

            // Modify DK for transition SST
            Dk = std::min(std::max(gammaEff, 0.1), 1.0) * Dk;

            // limit production
            if (Pk > tkeProdLimitRatio * Dk)
            {
                Pk = tkeProdLimitRatio * Dk;
            }

            // Modify Pk for Transition SST
            Pk *= gammaEff;

            rhs[0] += (Pk - Dk) * vol;
            lhs[0] += betaStar * rho * omega * vol;

            // global matrix
            Base::applyCoeff_(
                A, b, connectedNodes, scratchIds, scratchVals, rhs, lhs);
        }
    }
}

void turbulentKineticEnergyTransitionSSTAssembler::
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
    const STKScalarField* kSTKFieldPtr = phi_->stkFieldPtr();
    const STKScalarField* kSTKFieldPtrOld = phi_->prevTimeRef().stkFieldPtr();
    const STKScalarField* omegaSTKFieldPtr = model_->omegaRef().stkFieldPtr();
    const STKScalarField* gammaSTKFieldPtr = model_->gammaRef().stkFieldPtr();
    const STKScalarField* ReThetaSTKFieldPtr =
        model_->ReThetaRef().stkFieldPtr();
    const STKScalarField* rhoSTKFieldPtr = model_->rhoRef().stkFieldPtr();
    const STKScalarField* rhoSTKFieldPtrOld =
        model_->rhoRef().prevTimeRef().stkFieldPtr();
    const STKScalarField* muSTKFieldPtr = model_->muRef().stkFieldPtr();
    const STKScalarField* USTKFieldPtr = model_->URef().stkFieldPtr();
    const STKScalarField* gradUSTKFieldPtr =
        model_->URef().gradRef().stkFieldPtr();
    const STKScalarField* divSTKFieldPtr =
        model_->mDotRef().divRef().stkFieldPtr();
    const STKScalarField* yMinSTKFieldPtr = model_->yMinRef().stkFieldPtr();
    const STKScalarField* PkSTKFieldPtr = model_->PkRef().stkFieldPtr();

    scalar betaStar = model_->betaStar();
    scalar tkeProdLimitRatio = model_->tkeProdLimitRatio();
    scalar s1 = model_->s1();
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
            scalar k = *stk::mesh::field_data(*kSTKFieldPtr, node);
            scalar kOld = *stk::mesh::field_data(*kSTKFieldPtrOld, node);
            scalar omega = *stk::mesh::field_data(*omegaSTKFieldPtr, node);
            scalar gamma = *stk::mesh::field_data(*gammaSTKFieldPtr, node);
            scalar ReTheta = *stk::mesh::field_data(*ReThetaSTKFieldPtr, node);
            scalar rho = *stk::mesh::field_data(*rhoSTKFieldPtr, node);
            scalar rhoOld = *stk::mesh::field_data(*rhoSTKFieldPtrOld, node);
            scalar mu = *stk::mesh::field_data(*muSTKFieldPtr, node);
            scalar vol = *stk::mesh::field_data(*volSTKFieldPtr, node);
            scalar div = *stk::mesh::field_data(*divSTKFieldPtr, node);
            scalar Pk = *stk::mesh::field_data(*PkSTKFieldPtr, node);
            scalar* U = stk::mesh::field_data(*USTKFieldPtr, node);
            scalar* dudx = stk::mesh::field_data(*gradUSTKFieldPtr, node);
            scalar y = *stk::mesh::field_data(*yMinSTKFieldPtr, node);

            // divergence correction
            if (div < 0)
            {
                lhs[0] += -div;
            }
            rhs[0] -= -div * k;

            // S MAG
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

            // transient
            scalar lhsfac = c[0] * rho * vol / dt;
            scalar lhsfacOld = c[1] * rhoOld * vol / dt;

            lhs[0] += lhsfac;
            rhs[0] -= (lhsfac * k + lhsfacOld * kOld);

            // production/dissipation
            scalar Dk = betaStar * rho * omega * k;

            // Compute Volrticity Mag Wij
            scalar vorticityMag = 0.0;
            for (label i = 0; i < SPATIAL_DIM; ++i)
            {
                for (label j = 0; j < SPATIAL_DIM; ++j)
                {
                    const scalar rateOfStrain =
                        0.5 *
                        (dudx[SPATIAL_DIM * i + j] - dudx[SPATIAL_DIM * j + i]);
                    vorticityMag += rateOfStrain * rateOfStrain;
                }
            }
            vorticityMag = std::sqrt(2.0 * vorticityMag);

            // Compute Velocity Mag
            scalar UMag = 0;
            for (label i = 0; i < SPATIAL_DIM; ++i)
            {
                UMag += U[i] * U[i];
            }
            UMag = std::sqrt(UMag);

            scalar Rt = (rho * k) / (mu * omega);
            scalar fReattach = std::exp(-1.0 * std::pow((Rt / 20.0), 4.0));
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

            scalar ReOmega = (rho * omega * y * y) / mu;
            scalar fWake = exp(-1.0 * std::pow((ReOmega / 1E5), 2.0));

            scalar delta =
                (375.0 * vorticityMag * mu * ReTheta * y) / (rho * UMag * UMag);
            scalar fThetat = std::min(
                std::max(fWake * std::exp(-1.0 * std::pow((y / delta), 4.0)),
                         1.0 -
                             std::pow((ce2 * gamma - 1.0) / (ce2 - 1.0), 2.0)),
                1.0);

            scalar gammaSep =
                std::min(s1 * std::max(0.0, (Rev / (3.235 * ReThetac)) - 1.0) *
                             fReattach,
                         2.0) *
                fThetat;
            scalar gammaEff = std::max(gamma, gammaSep);

            // Modify DK for transition SST
            Dk = std::min(std::max(gammaEff, 0.1), 1.0) * Dk;

            // limit production
            if (Pk > tkeProdLimitRatio * Dk)
            {
                Pk = tkeProdLimitRatio * Dk;
            }

            // Modify Pk for Transition KWSST
            Pk *= gammaEff;

            rhs[0] += (Pk - Dk) * vol;
            lhs[0] += betaStar * rho * omega * vol;

            // geometric conservative law
            if (meshDeforming)
            {
                scalar divUm = *stk::mesh::field_data(*divUmSTKFieldPtr, node);
                rhs[0] -= rho * k * divUm * vol;
            }

            // global matrix
            Base::applyCoeff_(
                A, b, connectedNodes, scratchIds, scratchVals, rhs, lhs);
        }
    }
}

void turbulentKineticEnergyTransitionSSTAssembler::
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
    const STKScalarField* kSTKFieldPtr = phi_->stkFieldPtr();
    const STKScalarField* kSTKFieldPtrOld = phi_->prevTimeRef().stkFieldPtr();
    const STKScalarField* kSTKFieldPtrOldOld =
        phi_->prevTimeRef().prevTimeRef().stkFieldPtr();
    const STKScalarField* omegaSTKFieldPtr = model_->omegaRef().stkFieldPtr();
    const STKScalarField* gammaSTKFieldPtr = model_->gammaRef().stkFieldPtr();
    const STKScalarField* ReThetaSTKFieldPtr =
        model_->ReThetaRef().stkFieldPtr();
    const STKScalarField* rhoSTKFieldPtr = model_->rhoRef().stkFieldPtr();
    const STKScalarField* rhoSTKFieldPtrOld =
        model_->rhoRef().prevTimeRef().stkFieldPtr();
    const STKScalarField* rhoSTKFieldPtrOldOld =
        model_->rhoRef().prevTimeRef().prevTimeRef().stkFieldPtr();
    const STKScalarField* muSTKFieldPtr = model_->muRef().stkFieldPtr();
    const STKScalarField* divSTKFieldPtr =
        model_->mDotRef().divRef().stkFieldPtr();
    const STKScalarField* USTKFieldPtr = model_->URef().stkFieldPtr();
    const STKScalarField* gradUSTKFieldPtr =
        model_->URef().gradRef().stkFieldPtr();
    const STKScalarField* yMinSTKFieldPtr = model_->yMinRef().stkFieldPtr();
    const STKScalarField* PkSTKFieldPtr = model_->PkRef().stkFieldPtr();

    scalar betaStar = model_->betaStar();
    scalar tkeProdLimitRatio = model_->tkeProdLimitRatio();
    scalar s1 = model_->s1();
    scalar ce2 = model_->ce2();

    const STKScalarField* divUmSTKFieldPtr =
        meshDeforming ? model_->divUmRef().stkFieldPtr() : nullptr;

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
            scalar k = *stk::mesh::field_data(*kSTKFieldPtr, node);
            scalar kOld = *stk::mesh::field_data(*kSTKFieldPtrOld, node);
            scalar kOldOld = *stk::mesh::field_data(*kSTKFieldPtrOldOld, node);
            scalar omega = *stk::mesh::field_data(*omegaSTKFieldPtr, node);
            scalar gamma = *stk::mesh::field_data(*gammaSTKFieldPtr, node);
            scalar ReTheta = *stk::mesh::field_data(*ReThetaSTKFieldPtr, node);
            scalar rho = *stk::mesh::field_data(*rhoSTKFieldPtr, node);
            scalar rhoOld = *stk::mesh::field_data(*rhoSTKFieldPtrOld, node);
            scalar rhoOldOld =
                *stk::mesh::field_data(*rhoSTKFieldPtrOldOld, node);
            scalar mu = *stk::mesh::field_data(*muSTKFieldPtr, node);
            scalar vol = *stk::mesh::field_data(*volSTKFieldPtr, node);
            scalar div = *stk::mesh::field_data(*divSTKFieldPtr, node);
            scalar Pk = *stk::mesh::field_data(*PkSTKFieldPtr, node);
            scalar* U = stk::mesh::field_data(*USTKFieldPtr, node);
            scalar y = *stk::mesh::field_data(*yMinSTKFieldPtr, node);
            scalar* dudx = stk::mesh::field_data(*gradUSTKFieldPtr, node);

            // divergence correction
            if (div < 0)
            {
                lhs[0] += -div;
            }
            rhs[0] -= -div * k;

            // S MAG
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

            // transient
            scalar lhsfac = c[0] * rho * vol / dt;
            scalar lhsfacOld = c[1] * rhoOld * vol / dt;
            scalar lhsfacOldOld = c[2] * rhoOldOld * vol / dt;

            lhs[0] += lhsfac;
            rhs[0] -= (lhsfac * k + lhsfacOld * kOld + lhsfacOldOld * kOldOld);

            // production/dissipation
            scalar Dk = betaStar * rho * omega * k;

            // Compute Volrticity Mag Wij
            scalar vorticityMag = 0.0;
            for (label i = 0; i < SPATIAL_DIM; ++i)
            {
                for (label j = 0; j < SPATIAL_DIM; ++j)
                {
                    const scalar rateOfStrain =
                        0.5 *
                        (dudx[SPATIAL_DIM * i + j] - dudx[SPATIAL_DIM * j + i]);
                    vorticityMag += rateOfStrain * rateOfStrain;
                }
            }
            vorticityMag = std::sqrt(2.0 * vorticityMag);

            // Compute Velocity Mag
            scalar UMag = 0;
            for (label i = 0; i < SPATIAL_DIM; ++i)
            {
                UMag += U[i] * U[i];
            }
            UMag = std::sqrt(UMag);

            scalar Rt = (rho * k) / (mu * omega);
            scalar fReattach = std::exp(-1.0 * std::pow((Rt / 20.0), 4.0));
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

            scalar ReOmega = (rho * omega * y * y) / mu;
            scalar fWake = exp(-1.0 * std::pow((ReOmega / 1e5), 2.0));

            scalar delta =
                (375.0 * vorticityMag * mu * ReTheta * y) / (rho * UMag * UMag);
            scalar fThetat = std::min(
                std::max(fWake * std::exp(-1.0 * std::pow((y / delta), 4.0)),
                         1.0 -
                             std::pow((ce2 * gamma - 1.0) / (ce2 - 1.0), 2.0)),
                1.0);

            scalar gammaSep =
                std::min(s1 * std::max(0.0, (Rev / (3.235 * ReThetac)) - 1.0) *
                             fReattach,
                         2.0) *
                fThetat;
            scalar gammaEff = std::max(gamma, gammaSep);

            // Modify DK for transition SST
            Dk = std::min(std::max(gammaEff, 0.1), 1.0) * Dk;

            // limit production
            if (Pk > tkeProdLimitRatio * Dk)
            {
                Pk = tkeProdLimitRatio * Dk;
            }

            // Modify Pk for Transition KWSST
            Pk *= gammaEff;

            rhs[0] += (Pk - Dk) * vol;
            lhs[0] += betaStar * rho * omega * vol;

            // geometric conservative law
            if (meshDeforming)
            {
                scalar divUm = *stk::mesh::field_data(*divUmSTKFieldPtr, node);
                rhs[0] -= rho * k * divUm * vol;
            }

            // global matrix
            Base::applyCoeff_(
                A, b, connectedNodes, scratchIds, scratchVals, rhs, lhs);
        }
    }
}

} /* namespace accel */
