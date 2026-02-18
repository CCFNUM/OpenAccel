// File :
// turbulentKineticEnergyCorrelationTransitionSSTAssemblerNodeTerms.cpp Created
// : Sun Dec 29 2024 Author : Adam Fares Description: Node terms for TKE
// equation in correlation-based
// transition SST model (Menter 2015)
// Copyright (c) 2024 CCFNUM, Lucerne University of Applied Sciences and Arts.
// SPDX-License-Identifier: BSD-3-Clause

#include "turbulentKineticEnergyCorrelationTransitionSSTAssembler.h"

namespace accel
{

void turbulentKineticEnergyCorrelationTransitionSSTAssembler::
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
    const STKScalarField* gammaSTKFieldPtr = model_->gammaRef().stkFieldPtr();
    const STKScalarField* rhoSTKFieldPtr = model_->rhoRef().stkFieldPtr();
    const STKScalarField* muSTKFieldPtr = model_->muRef().stkFieldPtr();
    const STKScalarField* mutSTKFieldPtr = model_->mutRef().stkFieldPtr();
    const STKScalarField* gradUSTKFieldPtr =
        model_->URef().gradRef().stkFieldPtr();
    const STKScalarField* divSTKFieldPtr =
        model_->mDotRef().divRef().stkFieldPtr();
    const STKScalarField* yMinSTKFieldPtr = model_->yMinRef().stkFieldPtr();

    // Model constants
    const scalar betaStar = model_->betaStar();
    const scalar Ck_BLT = model_->Ck_BLT();
    const scalar CSEP = model_->CSEP();
    const scalar Retclim = model_->Retclim();

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
            scalar rho = *stk::mesh::field_data(*rhoSTKFieldPtr, node);
            scalar mu = *stk::mesh::field_data(*muSTKFieldPtr, node);
            scalar mut = *stk::mesh::field_data(*mutSTKFieldPtr, node);
            scalar vol = *stk::mesh::field_data(*volSTKFieldPtr, node);
            scalar div = *stk::mesh::field_data(*divSTKFieldPtr, node);
            scalar* dudx = stk::mesh::field_data(*gradUSTKFieldPtr, node);
            scalar y = *stk::mesh::field_data(*yMinSTKFieldPtr, node);

            // Compute strain rate and vorticity magnitudes
            scalar sijMag = 1.0e-16;
            scalar vortMag = 1.0e-16;
            for (label i = 0; i < SPATIAL_DIM; ++i)
            {
                for (label j = 0; j < SPATIAL_DIM; ++j)
                {
                    const scalar rateOfStrain =
                        0.5 *
                        (dudx[SPATIAL_DIM * i + j] + dudx[SPATIAL_DIM * j + i]);
                    const scalar vortTensor = 0.5 * (dudx[SPATIAL_DIM * i + j] -
                                                     dudx[SPATIAL_DIM * j + i]);
                    sijMag += rateOfStrain * rateOfStrain;
                    vortMag += vortTensor * vortTensor;
                }
            }
            sijMag = std::sqrt(2.0 * sijMag);
            vortMag = std::sqrt(2.0 * vortMag);

            // false transient
            scalar lhsfac = rho * vol / dt;
            lhs[0] += lhsfac;

            // divergence correction
            if (div < 0)
            {
                lhs[0] += -div;
            }
            rhs[0] -= -div * k;

            // Compute Rev and Fonlim
            scalar Rev = rho * y * y * sijMag / mu;
            scalar Fonlim =
                std::min(std::max(Rev / 2.2 / Retclim - 1.0, 0.0), 3.0);

            // Production term using Kato-Launder formulation
            scalar Pk = gamma * mut * sijMag * vortMag;

            // Additional production limiter for boundary layer transition
            //
            scalar Pklim = 5.0 * Ck_BLT * std::max(gamma - 0.2, 0.0) *
                           (1.0 - gamma) * Fonlim *
                           std::max(3.0 * CSEP * mu - mut, 0.0) * sijMag *
                           vortMag;

            // Dissipation with transition modification
            scalar Dk =
                betaStar * rho * omega * k * std::max(gamma, scalar(0.1));

            rhs[0] += (Pk + Pklim - Dk) * vol;
            lhs[0] +=
                betaStar * rho * omega * std::max(gamma, scalar(0.1)) * vol;

            // global matrix
            Base::applyCoeff_(
                A, b, connectedNodes, scratchIds, scratchVals, rhs, lhs);
        }
    }
}

void turbulentKineticEnergyCorrelationTransitionSSTAssembler::
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
    const STKScalarField* rhoSTKFieldPtr = model_->rhoRef().stkFieldPtr();
    const STKScalarField* rhoSTKFieldPtrOld =
        model_->rhoRef().prevTimeRef().stkFieldPtr();
    const STKScalarField* muSTKFieldPtr = model_->muRef().stkFieldPtr();
    const STKScalarField* mutSTKFieldPtr = model_->mutRef().stkFieldPtr();
    const STKScalarField* gradUSTKFieldPtr =
        model_->URef().gradRef().stkFieldPtr();
    const STKScalarField* divSTKFieldPtr =
        model_->mDotRef().divRef().stkFieldPtr();
    const STKScalarField* yMinSTKFieldPtr = model_->yMinRef().stkFieldPtr();

    // Model constants
    const scalar betaStar = model_->betaStar();
    const scalar Ck_BLT = model_->Ck_BLT();
    const scalar CSEP = model_->CSEP();
    const scalar Retclim = model_->Retclim();

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
            scalar rho = *stk::mesh::field_data(*rhoSTKFieldPtr, node);
            scalar rhoOld = *stk::mesh::field_data(*rhoSTKFieldPtrOld, node);
            scalar mu = *stk::mesh::field_data(*muSTKFieldPtr, node);
            scalar mut = *stk::mesh::field_data(*mutSTKFieldPtr, node);
            scalar vol = *stk::mesh::field_data(*volSTKFieldPtr, node);
            scalar div = *stk::mesh::field_data(*divSTKFieldPtr, node);
            scalar* dudx = stk::mesh::field_data(*gradUSTKFieldPtr, node);
            scalar y = *stk::mesh::field_data(*yMinSTKFieldPtr, node);

            // divergence correction
            if (div < 0)
            {
                lhs[0] += -div;
            }
            rhs[0] -= -div * k;

            // transient
            scalar lhsfac = c[0] * rho * vol / dt;
            scalar lhsfacOld = c[1] * rhoOld * vol / dt;

            lhs[0] += lhsfac;
            rhs[0] -= (lhsfac * k + lhsfacOld * kOld);

            // Compute strain rate and vorticity magnitudes
            scalar sijMag = 1.0e-16;
            scalar vortMag = 1.0e-16;
            for (label i = 0; i < SPATIAL_DIM; ++i)
            {
                for (label j = 0; j < SPATIAL_DIM; ++j)
                {
                    const scalar rateOfStrain =
                        0.5 *
                        (dudx[SPATIAL_DIM * i + j] + dudx[SPATIAL_DIM * j + i]);
                    const scalar vortTensor = 0.5 * (dudx[SPATIAL_DIM * i + j] -
                                                     dudx[SPATIAL_DIM * j + i]);
                    sijMag += rateOfStrain * rateOfStrain;
                    vortMag += vortTensor * vortTensor;
                }
            }
            sijMag = std::sqrt(2.0 * sijMag);
            vortMag = std::sqrt(2.0 * vortMag);

            // Compute Rev and Fonlim
            scalar Rev = rho * y * y * sijMag / mu;
            scalar Fonlim =
                std::min(std::max(Rev / 2.2 / Retclim - 1.0, 0.0), 3.0);

            // Production term using Kato-Launder formulation
            scalar Pk = gamma * mut * sijMag * vortMag;

            // Additional production limiter for boundary layer transition
            //
            scalar Pklim = 5.0 * Ck_BLT * std::max(gamma - 0.2, 0.0) *
                           (1.0 - gamma) * Fonlim *
                           std::max(3.0 * CSEP * mu - mut, 0.0) * sijMag *
                           vortMag;

            // Dissipation with transition modification
            scalar Dk =
                betaStar * rho * omega * k * std::max(gamma, scalar(0.1));

            rhs[0] += (Pk + Pklim - Dk) * vol;
            lhs[0] +=
                betaStar * rho * omega * std::max(gamma, scalar(0.1)) * vol;

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

void turbulentKineticEnergyCorrelationTransitionSSTAssembler::
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
    const STKScalarField* rhoSTKFieldPtr = model_->rhoRef().stkFieldPtr();
    const STKScalarField* rhoSTKFieldPtrOld =
        model_->rhoRef().prevTimeRef().stkFieldPtr();
    const STKScalarField* rhoSTKFieldPtrOldOld =
        model_->rhoRef().prevTimeRef().prevTimeRef().stkFieldPtr();
    const STKScalarField* muSTKFieldPtr = model_->muRef().stkFieldPtr();
    const STKScalarField* mutSTKFieldPtr = model_->mutRef().stkFieldPtr();
    const STKScalarField* gradUSTKFieldPtr =
        model_->URef().gradRef().stkFieldPtr();
    const STKScalarField* divSTKFieldPtr =
        model_->mDotRef().divRef().stkFieldPtr();
    const STKScalarField* yMinSTKFieldPtr = model_->yMinRef().stkFieldPtr();

    // Model constants
    const scalar betaStar = model_->betaStar();
    const scalar Ck_BLT = model_->Ck_BLT();
    const scalar CSEP = model_->CSEP();
    const scalar Retclim = model_->Retclim();

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
            scalar rho = *stk::mesh::field_data(*rhoSTKFieldPtr, node);
            scalar rhoOld = *stk::mesh::field_data(*rhoSTKFieldPtrOld, node);
            scalar rhoOldOld =
                *stk::mesh::field_data(*rhoSTKFieldPtrOldOld, node);
            scalar mu = *stk::mesh::field_data(*muSTKFieldPtr, node);
            scalar mut = *stk::mesh::field_data(*mutSTKFieldPtr, node);
            scalar vol = *stk::mesh::field_data(*volSTKFieldPtr, node);
            scalar div = *stk::mesh::field_data(*divSTKFieldPtr, node);
            scalar* dudx = stk::mesh::field_data(*gradUSTKFieldPtr, node);
            scalar y = *stk::mesh::field_data(*yMinSTKFieldPtr, node);

            // divergence correction
            if (div < 0)
            {
                lhs[0] += -div;
            }
            rhs[0] -= -div * k;

            // transient
            scalar lhsfac = c[0] * rho * vol / dt;
            scalar lhsfacOld = c[1] * rhoOld * vol / dt;
            scalar lhsfacOldOld = c[2] * rhoOldOld * vol / dt;

            lhs[0] += lhsfac;
            rhs[0] -= (lhsfac * k + lhsfacOld * kOld + lhsfacOldOld * kOldOld);

            // Compute strain rate and vorticity magnitudes
            scalar sijMag = 1.0e-16;
            scalar vortMag = 1.0e-16;
            for (label i = 0; i < SPATIAL_DIM; ++i)
            {
                for (label j = 0; j < SPATIAL_DIM; ++j)
                {
                    const scalar rateOfStrain =
                        0.5 *
                        (dudx[SPATIAL_DIM * i + j] + dudx[SPATIAL_DIM * j + i]);
                    const scalar vortTensor = 0.5 * (dudx[SPATIAL_DIM * i + j] -
                                                     dudx[SPATIAL_DIM * j + i]);
                    sijMag += rateOfStrain * rateOfStrain;
                    vortMag += vortTensor * vortTensor;
                }
            }
            sijMag = std::sqrt(2.0 * sijMag);
            vortMag = std::sqrt(2.0 * vortMag);

            // Compute Rev and Fonlim
            scalar Rev = rho * y * y * sijMag / mu;
            scalar Fonlim =
                std::min(std::max(Rev / 2.2 / Retclim - 1.0, 0.0), 3.0);

            // Production term using Kato-Launder formulation
            scalar Pk = gamma * mut * sijMag * vortMag;

            // Additional production limiter for boundary layer transition
            //
            scalar Pklim = 5.0 * Ck_BLT * std::max(gamma - 0.2, 0.0) *
                           (1.0 - gamma) * Fonlim *
                           std::max(3.0 * CSEP * mu - mut, 0.0) * sijMag *
                           vortMag;

            // Dissipation with transition modification
            scalar Dk =
                betaStar * rho * omega * k * std::max(gamma, scalar(0.1));

            rhs[0] += (Pk + Pklim - Dk) * vol;
            lhs[0] +=
                betaStar * rho * omega * std::max(gamma, scalar(0.1)) * vol;

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
