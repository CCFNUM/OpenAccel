// File : phiAssemblerNodeTerms.hpp
// Created : Thu Feb 29 2024 10:44:13 (+0100)
// Author : Mhamad Mahdi Alloush
// Description: Node based (pointwise) assembly kernel implementation
// Copyright (c) 2024 CCFNUM, Lucerne University of Applied Sciences and Arts.
// SPDX-License-Identifier: BSD-3-Clause

namespace accel
{

template <size_t N>
void phiAssembler<N>::assembleNodeTermsFused_(const domain* domain,
                                              Context* ctx)
{
    auto& mesh = field_broker_->meshRef();
    if (mesh.controlsRef().isTransient())
    {
        auto scheme = mesh.controlsRef()
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

template <size_t N>
void phiAssembler<N>::assembleNodeTermsFusedSteady_(const domain* domain,
                                                    Context* ctx)
{
    const bool includeAdv =
        (transportMode_ != diffusion) && (domain->type() == domainType::fluid);

    const auto& mesh = field_broker_->meshRef();
    Matrix& A = ctx->getAMatrix();
    Vector& b = ctx->getBVector();

    const stk::mesh::BulkData& bulkData = mesh.bulkDataRef();
    const stk::mesh::MetaData& metaData = mesh.metaDataRef();

    // space for LHS/RHS
    const label lhsSize = N * N;
    const label rhsSize = N;
    std::vector<scalar> lhs(lhsSize);
    std::vector<scalar> rhs(rhsSize);
    std::vector<label> scratchIds(rhsSize);
    std::vector<scalar> scratchVals(rhsSize);
    std::vector<stk::mesh::Entity> connectedNodes(1);

    // pointers
    scalar* p_lhs = &lhs[0];
    scalar* p_rhs = &rhs[0];

    // Get fields
    const STKScalarField* rhoSTKFieldPtr = this->rhoRef().stkFieldPtr();
    const STKScalarField* phiSTKFieldPtr = phi_->stkFieldPtr();

    // Geometric fields
    const auto* volSTKFieldPtr = metaData.get_field<scalar>(
        stk::topology::NODE_RANK, this->getDualNodalVolumeID_(domain));

    // other
    scalar dt = field_broker_->controlsRef().getPhysicalTimescale();

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
        scalar* phib = stk::mesh::field_data(*phiSTKFieldPtr, nodeBucket);

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
            scalar rho = rhob[iNode];
            scalar vol = volb[iNode];

            scalar div = includeAdv ? *stk::mesh::field_data(
                                          *divUSTKFieldPtr_, nodeBucket, iNode)
                                    : 0.0;

            // false transient: added later
            scalar lhsfac = rho * vol / dt;

            for (label i = 0; i < N; i++) // BLOCKSIZE
            {
                // divergence correction
                if (div < 0)
                {
                    lhs[i * N + i] += -div;
                }

                // divergence correction
                scalar phii = phib[N * iNode + i];
                rhs[i] -= -div * phii;

                // false transient
                lhs[i * N + i] += lhsfac;
            }

            Base::applyCoeff_(
                A, b, connectedNodes, scratchIds, scratchVals, rhs, lhs);
        }
    }
}

template <size_t N>
void phiAssembler<N>::assembleNodeTermsFusedFirstOrderUnsteady_(
    const domain* domain,
    Context* ctx)
{
    const bool includeAdv =
        (transportMode_ != diffusion) && (domain->type() == domainType::fluid);

    const auto& mesh = field_broker_->meshRef();
    Matrix& A = ctx->getAMatrix();
    Vector& b = ctx->getBVector();

    const stk::mesh::BulkData& bulkData = mesh.bulkDataRef();
    const stk::mesh::MetaData& metaData = mesh.metaDataRef();

    // space for LHS/RHS
    const label lhsSize = N * N;
    const label rhsSize = N;
    std::vector<scalar> lhs(lhsSize);
    std::vector<scalar> rhs(rhsSize);
    std::vector<label> scratchIds(rhsSize);
    std::vector<scalar> scratchVals(rhsSize);
    std::vector<stk::mesh::Entity> connectedNodes(1);

    // pointers
    scalar* p_lhs = &lhs[0];
    scalar* p_rhs = &rhs[0];

    // Get fields
    const STKScalarField* rhoSTKFieldPtr = this->rhoRef().stkFieldPtr();
    const STKScalarField* rhoSTKFieldPtrOld =
        this->rhoRef().prevTimeRef().stkFieldPtr();

    const STKScalarField* phiSTKFieldPtr = phi_->stkFieldPtr();
    const STKScalarField* phiSTKFieldPtrOld = phi_->prevTimeRef().stkFieldPtr();

    // Geometric fields
    const auto* volSTKFieldPtr = metaData.get_field<scalar>(
        stk::topology::NODE_RANK, this->getDualNodalVolumeID_(domain));

    // time integrator
    const scalar dt = mesh.controlsRef().getTimestep();
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
        scalar* rhobOld = stk::mesh::field_data(*rhoSTKFieldPtrOld, nodeBucket);
        scalar* volb = stk::mesh::field_data(*volSTKFieldPtr, nodeBucket);
        scalar* phib = stk::mesh::field_data(*phiSTKFieldPtr, nodeBucket);
        scalar* phibOld = stk::mesh::field_data(*phiSTKFieldPtrOld, nodeBucket);

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
            scalar rho = rhob[iNode];
            scalar rhoOld = rhobOld[iNode];
            scalar vol = volb[iNode];

            scalar div = includeAdv ? *stk::mesh::field_data(
                                          *divUSTKFieldPtr_, nodeBucket, iNode)
                                    : 0.0;

            // transient: added later
            scalar lhsfac = c[0] * rho * vol / dt;
            scalar lhsfacOld = c[1] * rhoOld * vol / dt;

            for (label i = 0; i < N; i++) // BLOCKSIZE
            {
                // divergence correction
                if (div < 0)
                {
                    lhs[i * N + i] += -div;
                }

                scalar phii = phib[N * iNode + i];
                scalar phiOldi = phibOld[N * iNode + i];

                // divergence correction
                rhs[i] -= -div * phii;

                // transient
                lhs[i * N + i] += lhsfac;
                rhs[i] -= (lhsfac * phii + lhsfacOld * phiOldi);
            }

            Base::applyCoeff_(
                A, b, connectedNodes, scratchIds, scratchVals, rhs, lhs);
        }
    }
}

template <size_t N>
void phiAssembler<N>::assembleNodeTermsFusedSecondOrderUnsteady_(
    const domain* domain,
    Context* ctx)
{
    const bool includeAdv =
        (transportMode_ != diffusion) && (domain->type() == domainType::fluid);

    const auto& mesh = field_broker_->meshRef();
    Matrix& A = ctx->getAMatrix();
    Vector& b = ctx->getBVector();

    const stk::mesh::BulkData& bulkData = mesh.bulkDataRef();
    const stk::mesh::MetaData& metaData = mesh.metaDataRef();

    // space for LHS/RHS
    const label lhsSize = N * N;
    const label rhsSize = N;
    std::vector<scalar> lhs(lhsSize);
    std::vector<scalar> rhs(rhsSize);
    std::vector<label> scratchIds(rhsSize);
    std::vector<scalar> scratchVals(rhsSize);
    std::vector<stk::mesh::Entity> connectedNodes(1);

    // pointers
    scalar* p_lhs = &lhs[0];
    scalar* p_rhs = &rhs[0];

    // Get fields
    const STKScalarField* rhoSTKFieldPtr = this->rhoRef().stkFieldPtr();
    const STKScalarField* rhoSTKFieldPtrOld =
        this->rhoRef().prevTimeRef().stkFieldPtr();
    const STKScalarField* rhoSTKFieldPtrOldOld =
        this->rhoRef().prevTimeRef().prevTimeRef().stkFieldPtr();

    const STKScalarField* phiSTKFieldPtr = phi_->stkFieldPtr();
    const STKScalarField* phiSTKFieldPtrOld = phi_->prevTimeRef().stkFieldPtr();
    const STKScalarField* phiSTKFieldPtrOldOld =
        phi_->prevTimeRef().prevTimeRef().stkFieldPtr();

    // Geometric fields
    const auto* volSTKFieldPtr = metaData.get_field<scalar>(
        stk::topology::NODE_RANK, this->getDualNodalVolumeID_(domain));

    // time integrator
    const scalar dt = mesh.controlsRef().getTimestep();
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
        scalar* rhobOld = stk::mesh::field_data(*rhoSTKFieldPtrOld, nodeBucket);
        scalar* rhobOldOld =
            stk::mesh::field_data(*rhoSTKFieldPtrOldOld, nodeBucket);
        scalar* volb = stk::mesh::field_data(*volSTKFieldPtr, nodeBucket);
        scalar* phib = stk::mesh::field_data(*phiSTKFieldPtr, nodeBucket);
        scalar* phibOld = stk::mesh::field_data(*phiSTKFieldPtrOld, nodeBucket);
        scalar* phibOldOld =
            stk::mesh::field_data(*phiSTKFieldPtrOldOld, nodeBucket);

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
            scalar rho = rhob[iNode];
            scalar rhoOld = rhobOld[iNode];
            scalar rhoOldOld = rhobOldOld[iNode];
            scalar vol = volb[iNode];

            scalar div = includeAdv ? *stk::mesh::field_data(
                                          *divUSTKFieldPtr_, nodeBucket, iNode)
                                    : 0.0;

            // transient: added later
            scalar lhsfac = c[0] * rho * vol / dt;
            scalar lhsfacOld = c[1] * rhoOld * vol / dt;
            scalar lhsfacOldOld = c[2] * rhoOldOld * vol / dt;

            for (label i = 0; i < N; i++) // BLOCKSIZE
            {
                // divergence correction
                if (div < 0)
                {
                    lhs[i * N + i] += -div;
                }

                scalar phii = phib[N * iNode + i];
                scalar phiOldi = phibOld[N * iNode + i];
                scalar phiOldOldi = phibOldOld[N * iNode + i];

                // divergence correction
                rhs[i] -= -div * phii;

                // transient
                lhs[i * N + i] += lhsfac;
                rhs[i] -= (lhsfac * phii + lhsfacOld * phiOldi +
                           lhsfacOldOld * phiOldOldi);
            }

            Base::applyCoeff_(
                A, b, connectedNodes, scratchIds, scratchVals, rhs, lhs);
        }
    }
}

} // namespace accel
