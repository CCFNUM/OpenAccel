// File       : totalEnergyAssemblerNodeTerms.cpp
// Created    : Mon Jun 09 2025 12:23:38 (+0100)
// Author     : Mhamad Mahdi Alloush
// Description:
// Copyright (c) 2025 CCFNUM, Lucerne University of Applied Sciences and Arts.
// SPDX-License-Identifier: BSD-3-Clause

#include "totalEnergyAssembler.h"

namespace accel
{

void totalEnergyAssembler::assembleNodeTermsFusedSteady_(const domain* domain,
                                                         Context* ctx)
{
    // retreive energy source in domain
    scalar energySourceValue = domain->energySource().value_[0];
    if (domain->energySource().option_ == sourceOption::totalSource)
    {
        energySourceValue /= domain->zonePtr()->stats().volume_;
    }

    if (domain->type() == domainType::solid)
    {
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
        scalar* p_lhs = &lhs[0];
        scalar* p_rhs = &rhs[0];

        // Geometric fields
        const auto* volSTKFieldPtr = metaData.get_field<scalar>(
            stk::topology::NODE_RANK, this->getDualNodalVolumeID_(domain));

        // other
        scalar dt = field_broker_->controlsRef().getPhysicalTimescale();

        // get interior parts the domain is defined on
        const stk::mesh::PartVector& partVec =
            domain->zonePtr()->interiorParts();

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

            const stk::mesh::Bucket::size_type nNodesPerBucket =
                nodeBucket.size();

            // field chunks in bucket
            scalar* volb = stk::mesh::field_data(*volSTKFieldPtr, nodeBucket);

            for (stk::mesh::Bucket::size_type iNode = 0;
                 iNode < nNodesPerBucket;
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
                scalar vol = volb[iNode];

                // energy source
                rhs[0] += energySourceValue * vol;

                Base::applyCoeff_(
                    A, b, connectedNodes, scratchIds, scratchVals, rhs, lhs);
            }
        }
    }
    else
    {
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
        scalar* p_lhs = &lhs[0];
        scalar* p_rhs = &rhs[0];

        // Get fields
        const STKScalarField* h0STKFieldPtr = phi_->stkFieldPtr();
        const STKScalarField* rhoSTKFieldPtr = this->rhoRef().stkFieldPtr();

        const STKScalarField* gradPSTKFieldPtr =
            model_->pRef().gradRef().stkFieldPtr();

        // Geometric fields
        const auto* volSTKFieldPtr = metaData.get_field<scalar>(
            stk::topology::NODE_RANK, this->getDualNodalVolumeID_(domain));

        const auto* coordinatesPtr = metaData.get_field<scalar>(
            stk::topology::NODE_RANK, this->getCoordinatesID_(domain));

        // required for frame motion (MFR)
        const auto coriolisMatrix = domain->zonePtr()->frameRotating()
                                        ? domain->zonePtr()
                                              ->transformationRef()
                                              .rotation()
                                              .coriolisMatrix_
                                        : utils::matrix::Zero();
        const scalar* p_mat = coriolisMatrix.data();

        const auto origin =
            domain->zonePtr()->frameRotating()
                ? domain->zonePtr()->transformationRef().rotation().origin_
                : utils::vector::Zero();
        const scalar* p_ori = origin.data();

        // Preallocate arrays for rotation work computation
        std::vector<scalar> r_vec(SPATIAL_DIM);
        std::vector<scalar> omega_cross_r(SPATIAL_DIM);

        // pointers
        scalar* p_r_vec = &r_vec[0];
        scalar* p_omega_cross_r = &omega_cross_r[0];

        // other
        scalar dt = field_broker_->controlsRef().getPhysicalTimescale();

        // get interior parts the domain is defined on
        const stk::mesh::PartVector& partVec =
            domain->zonePtr()->interiorParts();

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

            const stk::mesh::Bucket::size_type nNodesPerBucket =
                nodeBucket.size();

            // field chunks in bucket
            scalar* rhob = stk::mesh::field_data(*rhoSTKFieldPtr, nodeBucket);
            scalar* h0b = stk::mesh::field_data(*h0STKFieldPtr, nodeBucket);
            scalar* divb = stk::mesh::field_data(*divUSTKFieldPtr_, nodeBucket);
            scalar* volb = stk::mesh::field_data(*volSTKFieldPtr, nodeBucket);

            for (stk::mesh::Bucket::size_type iNode = 0;
                 iNode < nNodesPerBucket;
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
                scalar h0 = h0b[iNode];
                scalar rho = rhob[iNode];
                scalar div = divb[iNode];
                scalar vol = volb[iNode];

                // false transient
                scalar lhsfac = rho * vol / dt;
                lhs[0] += lhsfac;

                // divergence correction
                if (div < 0)
                {
                    lhs[0] += -div;
                }
                rhs[0] -= -div * h0;

                // Rotation work term: (ω×r)·∇p
                // From equation: ∇·[p(ω×r)] = (ω×r)·∇p (since ∇·(ω×r) = 0)
                const scalar* coords =
                    stk::mesh::field_data(*coordinatesPtr, node);
                const scalar* gradP =
                    stk::mesh::field_data(*gradPSTKFieldPtr, node);

                // Compute r = coords - origin
                for (label i = 0; i < SPATIAL_DIM; ++i)
                {
                    p_r_vec[i] = coords[i] - p_ori[i];
                }

                // Compute ω×r using Coriolis matrix
                for (label i = 0; i < SPATIAL_DIM; ++i)
                {
                    p_omega_cross_r[i] = 0.0;
                    for (label j = 0; j < SPATIAL_DIM; ++j)
                    {
                        p_omega_cross_r[i] +=
                            p_mat[i * SPATIAL_DIM + j] * p_r_vec[j];
                    }
                }

                // Compute (ω×r)·∇p
                scalar rotationWorkSource = 0.0;
                for (label i = 0; i < SPATIAL_DIM; ++i)
                {
                    rotationWorkSource += p_omega_cross_r[i] * gradP[i];
                }

                // Add to RHS
                rhs[0] -= rotationWorkSource * vol;

                // energy source
                rhs[0] += energySourceValue * vol;

                Base::applyCoeff_(
                    A, b, connectedNodes, scratchIds, scratchVals, rhs, lhs);
            }
        }
    }
}

void totalEnergyAssembler::assembleNodeTermsFusedFirstOrderUnsteady_(
    const domain* domain,
    Context* ctx)
{
    // retreive energy source in domain
    scalar energySourceValue = domain->energySource().value_[0];
    if (domain->energySource().option_ == sourceOption::totalSource)
    {
        energySourceValue /= domain->zonePtr()->stats().volume_;
    }

    const bool includeAdv = domain->type() == domainType::fluid;

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
    scalar* p_lhs = &lhs[0];
    scalar* p_rhs = &rhs[0];

    // Get fields
    const STKScalarField* rhoSTKFieldPtr = this->rhoRef().stkFieldPtr();
    const STKScalarField* rhoSTKFieldPtrOld =
        this->rhoRef().prevTimeRef().stkFieldPtr();
    const STKScalarField* pSTKFieldPtr = model_->pRef().stkFieldPtr();
    const STKScalarField* pSTKFieldPtrOld =
        model_->pRef().prevTimeRef().stkFieldPtr();
    const STKScalarField* h0STKFieldPtr = phi_->stkFieldPtr();
    const STKScalarField* h0STKFieldPtrOld = phi_->prevTimeRef().stkFieldPtr();

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
        scalar* h0b = stk::mesh::field_data(*h0STKFieldPtr, nodeBucket);
        scalar* h0bOld = stk::mesh::field_data(*h0STKFieldPtrOld, nodeBucket);
        scalar* pb = stk::mesh::field_data(*pSTKFieldPtr, nodeBucket);
        scalar* pOldb = stk::mesh::field_data(*pSTKFieldPtrOld, nodeBucket);
        scalar* volb = stk::mesh::field_data(*volSTKFieldPtr, nodeBucket);

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
            scalar h0 = h0b[iNode];
            scalar h0Old = h0bOld[iNode];
            scalar rho = rhob[iNode];
            scalar rhoOld = rhobOld[iNode];
            scalar vol = volb[iNode];

            scalar div = includeAdv
                             ? *stk::mesh::field_data(*divUSTKFieldPtr_, node)
                             : 0.0;

            // transient
            scalar lhsfac = c[0] * rho * vol / dt;
            scalar lhsfacOld = c[1] * rhoOld * vol / dt;
            lhs[0] += lhsfac;
            rhs[0] -= (lhsfac * h0 + lhsfacOld * h0Old);

            // divergence correction
            if (div < 0)
            {
                lhs[0] += -div;
            }
            rhs[0] -= -div * h0;

            // time-rate of change of pressure
            rhs[0] += (c[0] * pb[iNode] + c[1] * pOldb[iNode]) / dt * vol;

            // energy source
            rhs[0] += energySourceValue * vol;

            Base::applyCoeff_(
                A, b, connectedNodes, scratchIds, scratchVals, rhs, lhs);
        }
    }
}

void totalEnergyAssembler::assembleNodeTermsFusedSecondOrderUnsteady_(
    const domain* domain,
    Context* ctx)
{
    // retreive energy source in domain
    scalar energySourceValue = domain->energySource().value_[0];
    if (domain->energySource().option_ == sourceOption::totalSource)
    {
        energySourceValue /= domain->zonePtr()->stats().volume_;
    }

    const bool includeAdv = domain->type() == domainType::fluid;

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
    scalar* p_lhs = &lhs[0];
    scalar* p_rhs = &rhs[0];

    // Get fields
    const STKScalarField* rhoSTKFieldPtr = this->rhoRef().stkFieldPtr();
    const STKScalarField* rhoSTKFieldPtrOld =
        this->rhoRef().prevTimeRef().stkFieldPtr();
    const STKScalarField* rhoSTKFieldPtrOldOld =
        this->rhoRef().prevTimeRef().prevTimeRef().stkFieldPtr();
    const STKScalarField* pSTKFieldPtr = model_->pRef().stkFieldPtr();
    const STKScalarField* pSTKFieldPtrOld =
        model_->pRef().prevTimeRef().stkFieldPtr();
    const STKScalarField* pSTKFieldPtrOldOld =
        model_->pRef().prevTimeRef().prevTimeRef().stkFieldPtr();

    const STKScalarField* h0STKFieldPtr = phi_->stkFieldPtr();
    const STKScalarField* h0STKFieldPtrOld = phi_->prevTimeRef().stkFieldPtr();
    const STKScalarField* h0STKFieldPtrOldOld =
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
        scalar* h0b = stk::mesh::field_data(*h0STKFieldPtr, nodeBucket);
        scalar* h0bOld = stk::mesh::field_data(*h0STKFieldPtrOld, nodeBucket);
        scalar* h0bOldOld =
            stk::mesh::field_data(*h0STKFieldPtrOldOld, nodeBucket);
        scalar* pb = stk::mesh::field_data(*pSTKFieldPtr, nodeBucket);
        scalar* pOldb = stk::mesh::field_data(*pSTKFieldPtrOld, nodeBucket);
        scalar* pOldOldb =
            stk::mesh::field_data(*pSTKFieldPtrOldOld, nodeBucket);
        scalar* volb = stk::mesh::field_data(*volSTKFieldPtr, nodeBucket);

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
            scalar h0 = h0b[iNode];
            scalar h0Old = h0bOld[iNode];
            scalar h0OldOld = h0bOldOld[iNode];
            scalar rho = rhob[iNode];
            scalar rhoOld = rhobOld[iNode];
            scalar rhoOldOld = rhobOldOld[iNode];
            scalar vol = volb[iNode];

            scalar div = includeAdv
                             ? *stk::mesh::field_data(*divUSTKFieldPtr_, node)
                             : 0.0;

            // transient
            scalar lhsfac = c[0] * rho * vol / dt;
            scalar lhsfacOld = c[1] * rhoOld * vol / dt;
            scalar lhsfacOldOld = c[2] * rhoOldOld * vol / dt;
            lhs[0] += lhsfac;
            rhs[0] -=
                (lhsfac * h0 + lhsfacOld * h0Old + lhsfacOldOld * h0OldOld);

            // divergence correction
            if (div < 0)
            {
                lhs[0] += -div;
            }
            rhs[0] -= -div * h0;

            // time-rate of change of pressure
            rhs[0] += (c[0] * pb[iNode] + c[1] * pOldb[iNode] +
                       c[2] * pOldOldb[iNode]) /
                      dt * vol;

            // energy source
            rhs[0] += energySourceValue * vol;

            Base::applyCoeff_(
                A, b, connectedNodes, scratchIds, scratchVals, rhs, lhs);
        }
    }
}

} // namespace accel
