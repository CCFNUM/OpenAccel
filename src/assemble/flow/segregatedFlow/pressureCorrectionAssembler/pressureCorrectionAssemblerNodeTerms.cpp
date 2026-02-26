// File       : pressureCorrectionAssemblerNodeTerms.cpp
// Created    : Wed Jan 03 2024 13:38:51 (+0100)
// Author     : Mhamad Mahdi Alloush
// Description:
// Copyright (c) 2024 CCFNUM, Lucerne University of Applied Sciences and Arts.
// SPDX-License-Identifier: BSD-3-Clause

// code
#include "flowModel.h"
#include "pressureCorrectionAssembler.h"

namespace accel
{

void pressureCorrectionAssembler::assembleNodeTermsFused_(const domain* domain,
                                                          Context* ctx)
{
    if (field_broker_->controlsRef().isTransient())
    {
        auto scheme = field_broker_->controlsRef()
                          .solverRef()
                          .solverControl_.basicSettings_.transientScheme_;
        switch (scheme)
        {
            case transientSchemeType::firstOrderBackwardEuler:
                assembleNodeTermsFusedFirstOrderUnsteady_(domain, ctx);
                break;

            case transientSchemeType::secondOrderBackwardEuler:
                assembleNodeTermsFusedSecondOrderUnsteady_(domain, ctx);
                break;

            default:
                break;
        }
    }
    else
    {
        assembleNodeTermsFusedSteady_(domain, ctx);
    }
}

void pressureCorrectionAssembler::assembleNodeTermsFusedSteady_(
    const domain* domain,
    Context* ctx)
{
    const bool compressible = domain->isMaterialCompressible();

    const bool falseMassAccumulation =
        model_->controlsRef()
            .solverRef()
            .solverControl_.expertParameters_.falseMassAccumulation_;

    if (!compressible || !falseMassAccumulation)
        return;

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
    const STKScalarField* psiSTKFieldPtr = model_->psiRef().stkFieldPtr();

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

            // false transient
            scalar psi = *stk::mesh::field_data(*psiSTKFieldPtr, node);
            scalar vol = *stk::mesh::field_data(*volSTKFieldPtr, node);

            scalar lhsfac = vol / dt * psi;
            lhs[0] += lhsfac;

            // global matrix
            Base::applyCoeff_(
                A, b, connectedNodes, scratchIds, scratchVals, rhs, lhs);
        }
    }
}

void pressureCorrectionAssembler::assembleNodeTermsFusedFirstOrderUnsteady_(
    const domain* domain,
    Context* ctx)
{
    const bool meshDeforming = domain->zonePtr()->meshDeforming();

    const bool compressible = domain->isMaterialCompressible();

    if (!meshDeforming && !compressible)
        return;

    const auto& mesh = field_broker_->meshRef();

    Matrix& A = ctx->getAMatrix();
    Vector& b = ctx->getBVector();

    const stk::mesh::BulkData& bulkData = mesh.bulkDataRef();
    const stk::mesh::MetaData& metaData = mesh.metaDataRef();

    // space for LHS/RHS
    const label lhsppSize = 1;
    const label rhsSize = 1;

    std::vector<scalar> lhs(lhsppSize);
    std::vector<scalar> rhs(rhsSize);
    std::vector<label> scratchIds(rhsSize);
    std::vector<scalar> scratchVals(rhsSize);
    std::vector<stk::mesh::Entity> connectedNodes(1);

    // pointers
    scalar* p_lhs = &lhs[0];
    scalar* p_rhs = &rhs[0];

    // Get fields
    const STKScalarField* rhoSTKFieldPtr = model_->rhoRef().stkFieldPtr();
    const STKScalarField* rhoSTKFieldPtrOld =
        compressible ? model_->rhoRef().prevTimeRef().stkFieldPtr() : nullptr;
    const STKScalarField* psiSTKFieldPtr =
        compressible ? model_->psiRef().stkFieldPtr() : nullptr;
    const STKScalarField* psiSTKFieldPtrOld =
        compressible ? model_->psiRef().prevTimeRef().stkFieldPtr() : nullptr;
    const STKScalarField* divUmSTKFieldPtr =
        meshDeforming ? field_broker_->divUmRef().stkFieldPtr() : nullptr;

#ifndef NDEBUG
    // SCL check field for debug: accumulates mesh flux to verify SCL
    auto* sclCheckSTKFieldPtr =
        meshDeforming ? metaData.get_field<scalar>(stk::topology::NODE_RANK,
                                                   mesh::scl_check_ID)
                      : nullptr;
#endif /* NDEBUG */

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

        // field chunks in bucket
        scalar* rhob = stk::mesh::field_data(*rhoSTKFieldPtr, nodeBucket);
        scalar* volb = stk::mesh::field_data(*volSTKFieldPtr, nodeBucket);

        for (stk::mesh::Bucket::size_type iNode = 0; iNode < nNodesPerBucket;
             ++iNode)
        {
            // get node
            stk::mesh::Entity node = nodeBucket[iNode];
            connectedNodes[0] = node;

            for (label i = 0; i < lhsppSize; ++i)
            {
                p_lhs[i] = 0.0;
            }
            for (label i = 0; i < rhsSize; ++i)
            {
                p_rhs[i] = 0.0;
            }

            // get values of current node
            scalar rho = rhob[iNode];
            scalar vol = volb[iNode];

            // transient term (compressible only)
            if (compressible)
            {
                // get values of current node
                scalar psi = *stk::mesh::field_data(*psiSTKFieldPtr, node);
                scalar psiOld =
                    *stk::mesh::field_data(*psiSTKFieldPtrOld, node);
                scalar rhoOld =
                    *stk::mesh::field_data(*rhoSTKFieldPtrOld, node);

                scalar lhsfac = c[0] * vol / dt * psi;
                scalar lhsfacOld = c[1] * vol / dt * psiOld;

                lhs[0] += lhsfac;
                rhs[0] -= (lhsfac * rho + lhsfacOld * rhoOld);
            }

            // geometric conservative law
            if (meshDeforming)
            {
                scalar divUm = *stk::mesh::field_data(*divUmSTKFieldPtr, node);
                rhs[0] -= rho * divUm * vol;

#ifndef NDEBUG
                // SCL check: add GCL term (opposite sign to balance IP fluxes)
                if (sclCheckSTKFieldPtr)
                {
                    scalar* scl =
                        stk::mesh::field_data(*sclCheckSTKFieldPtr, node);
                    *scl += rho * divUm * vol;
                }
#endif /* NDEBUG */
            }

            this->applyCoeff_(
                A, b, connectedNodes, scratchIds, scratchVals, rhs, lhs);
        }
    }
}

void pressureCorrectionAssembler::assembleNodeTermsFusedSecondOrderUnsteady_(
    const domain* domain,
    Context* ctx)
{
    const bool meshDeforming = domain->zonePtr()->meshDeforming();

    const bool compressible = domain->isMaterialCompressible();

    if (!meshDeforming && !compressible)
        return;

    const auto& mesh = field_broker_->meshRef();

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
    scalar* p_lhspp = &lhs[0];
    scalar* p_rhs = &rhs[0];

    // Get fields
    const STKScalarField* rhoSTKFieldPtr = model_->rhoRef().stkFieldPtr();
    const STKScalarField* rhoSTKFieldPtrOld =
        compressible ? model_->rhoRef().prevTimeRef().stkFieldPtr() : nullptr;
    const STKScalarField* rhoSTKFieldPtrOldOld =
        compressible
            ? model_->rhoRef().prevTimeRef().prevTimeRef().stkFieldPtr()
            : nullptr;
    const STKScalarField* psiSTKFieldPtr =
        compressible ? model_->psiRef().stkFieldPtr() : nullptr;
    const STKScalarField* psiSTKFieldPtrOld =
        compressible ? model_->psiRef().prevTimeRef().stkFieldPtr() : nullptr;
    const STKScalarField* psiSTKFieldPtrOldOld =
        compressible
            ? model_->psiRef().prevTimeRef().prevTimeRef().stkFieldPtr()
            : nullptr;

    const STKScalarField* divUmSTKFieldPtr =
        meshDeforming ? field_broker_->divUmRef().stkFieldPtr() : nullptr;

#ifndef NDEBUG
    // SCL check field for debug: accumulates mesh flux to verify SCL
    auto* sclCheckSTKFieldPtr =
        meshDeforming ? metaData.get_field<scalar>(stk::topology::NODE_RANK,
                                                   mesh::scl_check_ID)
                      : nullptr;
#endif /* NDEBUG */

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

        // field chunks in bucket
        scalar* rhob = stk::mesh::field_data(*rhoSTKFieldPtr, nodeBucket);
        scalar* volb = stk::mesh::field_data(*volSTKFieldPtr, nodeBucket);

        for (stk::mesh::Bucket::size_type iNode = 0; iNode < nNodesPerBucket;
             ++iNode)
        {
            // get node
            stk::mesh::Entity node = nodeBucket[iNode];
            connectedNodes[0] = node;

            for (label i = 0; i < lhsSize; ++i)
            {
                p_lhspp[i] = 0.0;
            }
            for (label i = 0; i < rhsSize; ++i)
            {
                p_rhs[i] = 0.0;
            }

            // get values of current node
            scalar rho = rhob[iNode];
            scalar vol = volb[iNode];

            // transient term (compressible only)
            if (compressible)
            {
                // get values of current node
                scalar psi = *stk::mesh::field_data(*psiSTKFieldPtr, node);
                scalar psiOld =
                    *stk::mesh::field_data(*psiSTKFieldPtrOld, node);
                scalar psiOldOld =
                    *stk::mesh::field_data(*psiSTKFieldPtrOldOld, node);
                scalar rhoOld =
                    *stk::mesh::field_data(*rhoSTKFieldPtrOld, node);
                scalar rhoOldOld =
                    *stk::mesh::field_data(*rhoSTKFieldPtrOldOld, node);

                scalar lhsfac = c[0] * vol / dt * psi;
                scalar lhsfacOld = c[1] * vol / dt * psiOld;
                scalar lhsfacOldOld = c[2] * vol / dt * psiOldOld;

                lhs[0] += lhsfac;
                rhs[0] -= (lhsfac * rho + lhsfacOld * rhoOld +
                           lhsfacOldOld * rhoOldOld);
            }

            // geometric conservative law
            if (meshDeforming)
            {
                scalar divUm = *stk::mesh::field_data(*divUmSTKFieldPtr, node);
                rhs[0] -= rho * divUm * vol;

#ifndef NDEBUG
                // SCL check: add GCL term (opposite sign to balance IP fluxes)
                if (sclCheckSTKFieldPtr)
                {
                    scalar* scl =
                        stk::mesh::field_data(*sclCheckSTKFieldPtr, node);
                    *scl += rho * divUm * vol;
                }
#endif /* NDEBUG */
            }

            this->applyCoeff_(
                A, b, connectedNodes, scratchIds, scratchVals, rhs, lhs);
        }
    }
}

} // namespace accel
