// File : turbulentEddyFrequencySSTAssemblerNodeTerms.cpp
// Created : Thu Feb 22 2024 13:38:51 (+0100)
// Author : Mhamad Mahdi Alloush
// Description:
// Copyright (c) 2024 CCFNUM, Lucerne University of Applied Sciences and Arts.
// SPDX-License-Identifier: BSD-3-Clause

#include "turbulentEddyFrequencySSTAssembler.h"

namespace accel
{

void turbulentEddyFrequencySSTAssembler::assembleNodeTermsFusedSteady_(
    const domain* domain,
    Context* ctx)
{
    assert(phi_);

    const auto& mesh = model_->meshRef();
    Matrix& A = ctx->getAMatrix();
    Vector& b = ctx->getBVector();

    scalar comp = domain->isMaterialCompressible() ? 1.0 : 0.0;

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
    const STKScalarField* omegaSTKFieldPtr = phi_->stkFieldPtr();
    const STKScalarField* gradOmegaSTKFieldPtr = phi_->gradRef().stkFieldPtr();
    const STKScalarField* kSTKFieldPtrPre =
        model_->kRef().prevIterRef().stkFieldPtr();
    const STKScalarField* gradKSTKFieldPtr =
        model_->kRef().gradRef().stkFieldPtr();
    const STKScalarField* rhoSTKFieldPtr = model_->rhoRef().stkFieldPtr();
    const STKScalarField* divSTKFieldPtr =
        model_->mDotRef().divRef().stkFieldPtr();
    const STKScalarField* muSTKFieldPtr = model_->muRef().stkFieldPtr();
    const STKScalarField* gradUSTKFieldPtr =
        model_->URef().gradRef().stkFieldPtr();

    const STKScalarField* fOneBlendingSTKFieldPtr =
        model_->F1Ref().stkFieldPtr();

    const STKScalarField* yMinSTKFieldPtr = model_->yMinRef().stkFieldPtr();

    scalar betaStar = model_->betaStar();
    scalar tkeProdLimitRatio = model_->tkeProdLimitRatio();

    scalar a1 = model_->aOne();
    scalar Cmu = model_->Cmu();
    scalar sigmaWTwo = model_->sigmaWTwo();
    scalar betaOne = model_->betaOne();
    scalar betaTwo = model_->betaTwo();
    scalar gammaOne = model_->gammaOne();
    scalar gammaTwo = model_->gammaTwo();

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
            scalar kPre = *stk::mesh::field_data(*kSTKFieldPtrPre, node);
            scalar omega = *stk::mesh::field_data(*omegaSTKFieldPtr, node);
            scalar rho = *stk::mesh::field_data(*rhoSTKFieldPtr, node);
            scalar mu = *stk::mesh::field_data(*muSTKFieldPtr, node);
            scalar vol = *stk::mesh::field_data(*volSTKFieldPtr, node);
            scalar div = *stk::mesh::field_data(*divSTKFieldPtr, node);
            scalar fOneBlend =
                *stk::mesh::field_data(*fOneBlendingSTKFieldPtr, node);
            scalar y = *stk::mesh::field_data(*yMinSTKFieldPtr, node);
            scalar* gradU = stk::mesh::field_data(*gradUSTKFieldPtr, node);
            scalar* gradK = stk::mesh::field_data(*gradKSTKFieldPtr, node);
            scalar* gradOmega =
                stk::mesh::field_data(*gradOmegaSTKFieldPtr, node);

            // false transient
            scalar lhsfac = rho * vol / dt;
            lhs[0] += lhsfac;

            // divergence correction
            if (div < 0)
            {
                lhs[0] += -div;
            }
            rhs[0] -= -div * omega;

            scalar PkByMut = 0;
            for (label i = 0; i < SPATIAL_DIM; ++i)
            {
                const label offSet = SPATIAL_DIM * i;
                for (label j = 0; j < SPATIAL_DIM; ++j)
                {
                    PkByMut += gradU[offSet + j] *
                               (gradU[offSet + j] + gradU[SPATIAL_DIM * j + i]);
                }
            }

            // ensure non-negative production
            PkByMut = std::max(PkByMut, 0.0);

            // Add dilation production
            scalar divergence = 0.0;
            for (label i = 0; i < SPATIAL_DIM; ++i)
            {
                divergence += gradU[SPATIAL_DIM * i + i];
            }
            scalar dilation_production = -2.0 / 3.0 * divergence * divergence;
            PkByMut += dilation_production * comp;

            // compute strain rate magnitude; pull pointer within the loop to
            // make it managable
            scalar sijMag = 0.0;
            for (label i = 0; i < SPATIAL_DIM; ++i)
            {
                for (label j = 0; j < SPATIAL_DIM; ++j)
                {
                    const scalar rateOfStrain =
                        0.5 * (gradU[SPATIAL_DIM * i + j] +
                               gradU[SPATIAL_DIM * j + i]);
                    sijMag += rateOfStrain * rateOfStrain;
                }
            }
            sijMag = std::sqrt(2.0 * sijMag);

            // FIXME: arg2 should use tke.prevIter and tef.prevIter
            scalar arg2 = std::min(
                std::max(2 * sqrt(kPre) / (Cmu * omega * y + SMALL),
                         500 * mu / (std::pow(y, 2.0) * rho * omega + SMALL)),
                scalar(100));

            scalar F2_old = tanh(pow(arg2, 2.0));

            // Production limited
            PkByMut =
                std::min(PkByMut,
                         ((tkeProdLimitRatio * betaStar) / a1) * omega *
                             std::max(a1 * omega, F2_old * sqrt(2.0) * sijMag));

            // production/dissipation
            scalar crossDiff = 0.0;
            for (label i = 0; i < SPATIAL_DIM; ++i)
            {
                crossDiff += gradK[i] * gradOmega[i];
            }

            // start the blending and constants
            const scalar om_fOneBlend = 1.0 - fOneBlend;
            const scalar beta = fOneBlend * betaOne + om_fOneBlend * betaTwo;
            const scalar gamma = fOneBlend * gammaOne + om_fOneBlend * gammaTwo;
            const scalar sigmaD = 2.0 * om_fOneBlend * sigmaWTwo;

            // Pw includes 1/tvisc scaling; tvisc may be zero at a dirichlet low
            // Re approach (clip)
            const scalar Pw = gamma * rho * PkByMut;
            const scalar Dw = beta * rho * omega * omega;
            const scalar Sw = sigmaD * rho * crossDiff / (omega + SMALL);

            rhs[0] += (Pw - Dw + Sw) * vol;
            lhs[0] += (2.0 * beta * rho * omega +
                       std::max(Sw / (omega + SMALL), 0.0)) *
                      vol;

            // global matrix
            Base::applyCoeff_(
                A, b, connectedNodes, scratchIds, scratchVals, rhs, lhs);
        }
    }
}

void turbulentEddyFrequencySSTAssembler::
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

    scalar comp = domain->isMaterialCompressible() ? 1.0 : 0.0;

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
    const STKScalarField* omegaSTKFieldPtr = phi_->stkFieldPtr();
    const STKScalarField* omegaSTKFieldPtrOld =
        phi_->prevTimeRef().stkFieldPtr();
    const STKScalarField* gradOmegaSTKFieldPtr = phi_->gradRef().stkFieldPtr();
    const STKScalarField* kPrevSTKFieldPtr =
        model_->kRef().prevIterRef().stkFieldPtr();
    const STKScalarField* gradKSTKFieldPtr =
        model_->kRef().gradRef().stkFieldPtr();
    const STKScalarField* rhoSTKFieldPtr = model_->rhoRef().stkFieldPtr();
    const STKScalarField* rhoSTKFieldPtrOld =
        model_->rhoRef().prevTimeRef().stkFieldPtr();
    const STKScalarField* divSTKFieldPtr =
        model_->mDotRef().divRef().stkFieldPtr();
    const STKScalarField* muSTKFieldPtr = model_->muRef().stkFieldPtr();
    const STKScalarField* gradUSTKFieldPtr =
        model_->URef().gradRef().stkFieldPtr();

    const STKScalarField* fOneBlendingSTKFieldPtr =
        model_->F1Ref().stkFieldPtr();

    const STKScalarField* yMinSTKFieldPtr = model_->yMinRef().stkFieldPtr();

    scalar betaStar = model_->betaStar();
    scalar tkeProdLimitRatio = model_->tkeProdLimitRatio();

    scalar a1 = model_->aOne();
    scalar Cmu = model_->Cmu();
    scalar sigmaWTwo = model_->sigmaWTwo();
    scalar betaOne = model_->betaOne();
    scalar betaTwo = model_->betaTwo();
    scalar gammaOne = model_->gammaOne();
    scalar gammaTwo = model_->gammaTwo();

    const STKScalarField* divUmSTKFieldPtr =
        meshDeforming ? model_->divUmRef().stkFieldPtr() : nullptr;

    // Geometric fields
    const auto* volSTKFieldPtr = metaData.get_field<scalar>(
        stk::topology::NODE_RANK, this->getDualNodalVolumeID_(domain));

    // time integrator
    const scalar dt = model_->controlsRef().getTimestep();
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
            scalar kPre = *stk::mesh::field_data(*kPrevSTKFieldPtr, node);
            scalar omega = *stk::mesh::field_data(*omegaSTKFieldPtr, node);
            scalar omegaOld =
                *stk::mesh::field_data(*omegaSTKFieldPtrOld, node);
            scalar rho = *stk::mesh::field_data(*rhoSTKFieldPtr, node);
            scalar rhoOld = *stk::mesh::field_data(*rhoSTKFieldPtrOld, node);
            scalar mu = *stk::mesh::field_data(*muSTKFieldPtr, node);
            scalar vol = *stk::mesh::field_data(*volSTKFieldPtr, node);
            scalar div = *stk::mesh::field_data(*divSTKFieldPtr, node);
            scalar fOneBlend =
                *stk::mesh::field_data(*fOneBlendingSTKFieldPtr, node);
            scalar y = *stk::mesh::field_data(*yMinSTKFieldPtr, node);
            scalar* gradU = stk::mesh::field_data(*gradUSTKFieldPtr, node);

            scalar* gradK = stk::mesh::field_data(*gradKSTKFieldPtr, node);
            scalar* gradOmega =
                stk::mesh::field_data(*gradOmegaSTKFieldPtr, node);

            // divergence correction
            if (div < 0)
            {
                lhs[0] += -div;
            }
            rhs[0] -= -div * omega;

            scalar PkByMut = 0;
            for (label i = 0; i < SPATIAL_DIM; ++i)
            {
                const label offSet = SPATIAL_DIM * i;
                for (label j = 0; j < SPATIAL_DIM; ++j)
                {
                    PkByMut += gradU[offSet + j] *
                               (gradU[offSet + j] + gradU[SPATIAL_DIM * j + i]);
                }
            }

            // ensure non-negative production
            PkByMut = std::max(PkByMut, 0.0);

            // Add dilation production
            scalar divergence = 0.0;
            for (label i = 0; i < SPATIAL_DIM; ++i)
            {
                divergence += gradU[SPATIAL_DIM * i + i];
            }
            scalar dilation_production = -2.0 / 3.0 * divergence * divergence;
            PkByMut += dilation_production * comp;

            // transient
            scalar lhsfac = c[0] * rho * vol / dt;
            scalar lhsfacOld = c[1] * rhoOld * vol / dt;

            lhs[0] += lhsfac;
            rhs[0] -= (lhsfac * omega + lhsfacOld * omegaOld);

            // compute strain rate magnitude; pull pointer within the loop to
            // make it managable
            scalar sijMag = 0.0;
            for (label i = 0; i < SPATIAL_DIM; ++i)
            {
                for (label j = 0; j < SPATIAL_DIM; ++j)
                {
                    const scalar rateOfStrain =
                        0.5 * (gradU[SPATIAL_DIM * i + j] +
                               gradU[SPATIAL_DIM * j + i]);
                    sijMag += rateOfStrain * rateOfStrain;
                }
            }
            sijMag = std::sqrt(2.0 * sijMag);

            // FIXME: arg2 should use tke.prevIter and tef.prevIter
            scalar arg2 = std::min(
                std::max(2 * sqrt(kPre) / (Cmu * omega * y + SMALL),
                         500 * mu / (std::pow(y, 2.0) * rho * omega + SMALL)),
                scalar(100));

            scalar F2_old = tanh(pow(arg2, 2.0));

            // Production limited
            PkByMut =
                std::min(PkByMut,
                         ((tkeProdLimitRatio * betaStar) / a1) * omega *
                             std::max(a1 * omega, F2_old * sqrt(2.0) * sijMag));

            // production/dissipation
            scalar crossDiff = 0.0;
            for (label i = 0; i < SPATIAL_DIM; ++i)
            {
                crossDiff += gradK[i] * gradOmega[i];
            }

            // start the blending and constants
            const scalar om_fOneBlend = 1.0 - fOneBlend;
            const scalar beta = fOneBlend * betaOne + om_fOneBlend * betaTwo;
            const scalar gamma = fOneBlend * gammaOne + om_fOneBlend * gammaTwo;
            const scalar sigmaD = 2.0 * om_fOneBlend * sigmaWTwo;

            // Pw includes 1/tvisc scaling; tvisc may be zero at a dirichlet low
            // Re approach (clip)
            const scalar Pw = gamma * rho * PkByMut;
            const scalar Dw = beta * rho * omega * omega;
            const scalar Sw = sigmaD * rho * crossDiff / (omega + SMALL);

            rhs[0] += (Pw - Dw + Sw) * vol;
            lhs[0] += (2.0 * beta * rho * omega +
                       std::max(Sw / (omega + SMALL), 0.0)) *
                      vol;

            // geometric conservative law
            if (meshDeforming)
            {
                scalar divUm = *stk::mesh::field_data(*divUmSTKFieldPtr, node);
                rhs[0] -= rho * omega * divUm * vol;
            }

            // global matrix
            Base::applyCoeff_(
                A, b, connectedNodes, scratchIds, scratchVals, rhs, lhs);
        }
    }
}

void turbulentEddyFrequencySSTAssembler::
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

    scalar comp = domain->isMaterialCompressible() ? 1.0 : 0.0;

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
    const STKScalarField* omegaSTKFieldPtr = phi_->stkFieldPtr();
    const STKScalarField* omegaSTKFieldPtrOld =
        phi_->prevTimeRef().stkFieldPtr();
    const STKScalarField* omegaSTKFieldPtrOldOld =
        phi_->prevTimeRef().prevTimeRef().stkFieldPtr();
    const STKScalarField* gradOmegaSTKFieldPtr = phi_->gradRef().stkFieldPtr();
    const STKScalarField* kPrevSTKFieldPtr =
        model_->kRef().prevIterRef().stkFieldPtr();
    const STKScalarField* gradKSTKFieldPtr =
        model_->kRef().gradRef().stkFieldPtr();
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

    const STKScalarField* fOneBlendingSTKFieldPtr =
        model_->F1Ref().stkFieldPtr();

    const STKScalarField* yMinSTKFieldPtr = model_->yMinRef().stkFieldPtr();

    scalar betaStar = model_->betaStar();
    scalar tkeProdLimitRatio = model_->tkeProdLimitRatio();

    scalar a1 = model_->aOne();
    scalar Cmu = model_->Cmu();
    scalar sigmaWTwo = model_->sigmaWTwo();
    scalar betaOne = model_->betaOne();
    scalar betaTwo = model_->betaTwo();
    scalar gammaOne = model_->gammaOne();
    scalar gammaTwo = model_->gammaTwo();

    const STKScalarField* divUmSTKFieldPtr =
        meshDeforming ? model_->divUmRef().stkFieldPtr() : nullptr;

    // Geometric fields
    const auto* volSTKFieldPtr = metaData.get_field<scalar>(
        stk::topology::NODE_RANK, this->getDualNodalVolumeID_(domain));

    // time integrator
    const scalar dt = model_->controlsRef().getTimestep();
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
            scalar kPre = *stk::mesh::field_data(*kPrevSTKFieldPtr, node);
            scalar omega = *stk::mesh::field_data(*omegaSTKFieldPtr, node);
            scalar omegaOld =
                *stk::mesh::field_data(*omegaSTKFieldPtrOld, node);
            scalar omegaOldOld =
                *stk::mesh::field_data(*omegaSTKFieldPtrOldOld, node);
            scalar rho = *stk::mesh::field_data(*rhoSTKFieldPtr, node);
            scalar rhoOld = *stk::mesh::field_data(*rhoSTKFieldPtrOld, node);
            scalar rhoOldOld =
                *stk::mesh::field_data(*rhoSTKFieldPtrOldOld, node);
            scalar mu = *stk::mesh::field_data(*muSTKFieldPtr, node);
            scalar vol = *stk::mesh::field_data(*volSTKFieldPtr, node);
            scalar div = *stk::mesh::field_data(*divSTKFieldPtr, node);
            scalar fOneBlend =
                *stk::mesh::field_data(*fOneBlendingSTKFieldPtr, node);
            scalar y = *stk::mesh::field_data(*yMinSTKFieldPtr, node);
            scalar* gradU = stk::mesh::field_data(*gradUSTKFieldPtr, node);

            scalar* gradK = stk::mesh::field_data(*gradKSTKFieldPtr, node);
            scalar* gradOmega =
                stk::mesh::field_data(*gradOmegaSTKFieldPtr, node);

            // divergence correction
            if (div < 0)
            {
                lhs[0] += -div;
            }
            rhs[0] -= -div * omega;

            scalar PkByMut = 0;
            for (label i = 0; i < SPATIAL_DIM; ++i)
            {
                const label offSet = SPATIAL_DIM * i;
                for (label j = 0; j < SPATIAL_DIM; ++j)
                {
                    PkByMut += gradU[offSet + j] *
                               (gradU[offSet + j] + gradU[SPATIAL_DIM * j + i]);
                }
            }

            // ensure non-negative production
            PkByMut = std::max(PkByMut, 0.0);

            // Add dilation production
            scalar divergence = 0.0;
            for (label i = 0; i < SPATIAL_DIM; ++i)
            {
                divergence += gradU[SPATIAL_DIM * i + i];
            }
            scalar dilation_production = -2.0 / 3.0 * divergence * divergence;
            PkByMut += dilation_production * comp;

            // transient
            scalar lhsfac = c[0] * rho * vol / dt;
            scalar lhsfacOld = c[1] * rhoOld * vol / dt;
            scalar lhsfacOldOld = c[2] * rhoOldOld * vol / dt;

            lhs[0] += lhsfac;
            rhs[0] -= (lhsfac * omega + lhsfacOld * omegaOld +
                       lhsfacOldOld * omegaOldOld);

            // compute strain rate magnitude; pull pointer within the loop to
            // make it managable
            scalar sijMag = 0.0;
            for (label i = 0; i < SPATIAL_DIM; ++i)
            {
                for (label j = 0; j < SPATIAL_DIM; ++j)
                {
                    const scalar rateOfStrain =
                        0.5 * (gradU[SPATIAL_DIM * i + j] +
                               gradU[SPATIAL_DIM * j + i]);
                    sijMag += rateOfStrain * rateOfStrain;
                }
            }
            sijMag = std::sqrt(2.0 * sijMag);

            // FIXME: arg2 should use tke.prevIter and tef.prevIter
            scalar arg2 = std::min(
                std::max(2 * sqrt(kPre) / (Cmu * omega * y + SMALL),
                         500 * mu / (std::pow(y, 2.0) * rho * omega + SMALL)),
                scalar(100));

            scalar F2_old = tanh(pow(arg2, 2.0));

            // Production limited
            PkByMut =
                std::min(PkByMut,
                         ((tkeProdLimitRatio * betaStar) / a1) * omega *
                             std::max(a1 * omega, F2_old * sqrt(2.0) * sijMag));

            // production/dissipation
            scalar crossDiff = 0.0;
            for (label i = 0; i < SPATIAL_DIM; ++i)
            {
                crossDiff += gradK[i] * gradOmega[i];
            }

            // start the blending and constants
            const scalar om_fOneBlend = 1.0 - fOneBlend;
            const scalar beta = fOneBlend * betaOne + om_fOneBlend * betaTwo;
            const scalar gamma = fOneBlend * gammaOne + om_fOneBlend * gammaTwo;
            const scalar sigmaD = 2.0 * om_fOneBlend * sigmaWTwo;

            // Pw includes 1/tvisc scaling; tvisc may be zero at a dirichlet low
            // Re approach (clip)
            const scalar Pw = gamma * rho * PkByMut;
            const scalar Dw = beta * rho * omega * omega;
            const scalar Sw = sigmaD * rho * crossDiff / (omega + SMALL);

            rhs[0] += (Pw - Dw + Sw) * vol;
            lhs[0] += (2.0 * beta * rho * omega +
                       std::max(Sw / (omega + SMALL), 0.0)) *
                      vol;

            // geometric conservative law
            if (meshDeforming)
            {
                scalar divUm = *stk::mesh::field_data(*divUmSTKFieldPtr, node);
                rhs[0] -= rho * omega * divUm * vol;
            }

            // global matrix
            Base::applyCoeff_(
                A, b, connectedNodes, scratchIds, scratchVals, rhs, lhs);
        }
    }
}

} /* namespace accel */
