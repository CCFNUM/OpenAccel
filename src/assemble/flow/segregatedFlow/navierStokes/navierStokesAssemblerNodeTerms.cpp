// File       : navierStokesAssemblerNodeTerms.cpp
// Created    : Wed Jan 03 2024 13:38:51 (+0100)
// Author     : Mhamad Mahdi Alloush
// Description:
// Copyright (c) 2024 CCFNUM, Lucerne University of Applied Sciences and Arts.
// SPDX-License-Identifier: BSD-3-Clause

// code
#include "flowModel.h"
#include "navierStokesAssembler.h"
#include "zoneTransformation.h"

namespace accel
{

void navierStokesAssembler::assembleNodeTermsFusedSteady_(const domain* domain,
                                                          Context* ctx)
{
    auto& mesh = model_->meshRef();
    Matrix& A = ctx->getAMatrix();
    Vector& b = ctx->getBVector();

    stk::mesh::MetaData& metaData = mesh.metaDataRef();
    stk::mesh::BulkData& bulkData = mesh.bulkDataRef();

    // space for LHS/RHS
    const label lhsSize = (SPATIAL_DIM * SPATIAL_DIM); // BLOCKSIZE * BLOCKSIZE
    const label rhsSize = SPATIAL_DIM;                 // BLOCKSIZE
    std::vector<scalar> lhs(lhsSize);
    std::vector<scalar> rhs(rhsSize);
    std::vector<label> scratchIds(rhsSize);
    std::vector<scalar> scratchVals(rhsSize);
    std::vector<stk::mesh::Entity> connectedNodes(1);

    // pointers
    scalar* p_lhs = &lhs[0];
    scalar* p_rhs = &rhs[0];

    // Get fields
    const STKScalarField* gradPSTKFieldPtr =
        model_->pRef().gradRef().stkFieldPtr();
    const STKScalarField* rhoSTKFieldPtr = model_->rhoRef().stkFieldPtr();
    const STKScalarField* divSTKFieldPtr =
        model_->mDotRef().divRef().stkFieldPtr();
    const STKScalarField* USTKFieldPtr = model_->URef().stkFieldPtr();

    // Geometric fields
    const auto* volSTKFieldPtr = metaData.get_field<scalar>(
        stk::topology::NODE_RANK, this->getDualNodalVolumeID_(domain));

    // Body forces (use redistributed to balance pressure gradient)
    const auto* FSTKFieldPtr =
        metaData.get_field<scalar>(stk::topology::NODE_RANK, flowModel::F_ID);

    // other
    scalar dt = model_->controlsRef().getPhysicalTimescale();

    // coriolis matrix and extract pointer
    const auto coriolisMatrix =
        domain->zonePtr()->frameRotating()
            ? domain->zonePtr()->transformationRef().rotation().coriolisMatrix_
            : utils::matrix::Zero();
    const scalar* p_mat = coriolisMatrix.data();

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
        scalar* gradPb = stk::mesh::field_data(*gradPSTKFieldPtr, nodeBucket);
        scalar* rhob = stk::mesh::field_data(*rhoSTKFieldPtr, nodeBucket);
        scalar* volb = stk::mesh::field_data(*volSTKFieldPtr, nodeBucket);
        scalar* divb = stk::mesh::field_data(*divSTKFieldPtr, nodeBucket);
        scalar* Ub = stk::mesh::field_data(*USTKFieldPtr, nodeBucket);
        scalar* Fb = stk::mesh::field_data(*FSTKFieldPtr, nodeBucket);

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
            scalar div = divb[iNode];

            scalar lhsfac = rho * vol / dt;

            for (label i = 0; i < SPATIAL_DIM; i++)
            {
                // divergence correction
                if (div < 0)
                {
                    lhs[i * SPATIAL_DIM + i] += -div;
                }
                scalar Ui = Ub[SPATIAL_DIM * iNode + i];
                rhs[i] -= -div * Ui;

                // pressure gradient
                scalar gradPi = gradPb[SPATIAL_DIM * iNode + i];
                rhs[i] -= gradPi * vol;

                // false transient
                lhs[i * SPATIAL_DIM + i] += lhsfac;

                // coriolis acceleration term
                for (label j = 0; j < SPATIAL_DIM; j++)
                {
                    lhs[i * SPATIAL_DIM + j] +=
                        p_mat[i * SPATIAL_DIM + j] * vol;

                    rhs[i] -= p_mat[i * SPATIAL_DIM + j] * rho *
                              Ub[SPATIAL_DIM * iNode + j] * vol;
                }

                // body forces
                scalar Fi = Fb[SPATIAL_DIM * iNode + i];
                rhs[i] += Fi * vol;
            }

            assembler::applyCoeff_(
                A, b, connectedNodes, scratchIds, scratchVals, rhs, lhs);
        }
    }
}

void navierStokesAssembler::assembleNodeTermsFusedFirstOrderUnsteady_(
    const domain* domain,
    Context* ctx)
{
    auto& mesh = model_->meshRef();
    Matrix& A = ctx->getAMatrix();
    Vector& b = ctx->getBVector();

    stk::mesh::MetaData& metaData = mesh.metaDataRef();
    stk::mesh::BulkData& bulkData = mesh.bulkDataRef();

    // space for LHS/RHS
    const label lhsSize = SPATIAL_DIM * SPATIAL_DIM;
    const label rhsSize = SPATIAL_DIM;
    std::vector<scalar> lhs(lhsSize);
    std::vector<scalar> rhs(rhsSize);
    std::vector<label> scratchIds(rhsSize);
    std::vector<scalar> scratchVals(rhsSize);
    std::vector<stk::mesh::Entity> connectedNodes(1);

    // pointers
    scalar* p_lhs = &lhs[0];
    scalar* p_rhs = &rhs[0];

    // Get fields
    const STKScalarField* gradPSTKFieldPtr =
        model_->pRef().gradRef().stkFieldPtr();
    const STKScalarField* rhoSTKFieldPtr = model_->rhoRef().stkFieldPtr();
    const STKScalarField* rhoSTKFieldPtrOld =
        model_->rhoRef().prevTimeRef().stkFieldPtr();

    const STKScalarField* USTKFieldPtr = model_->URef().stkFieldPtr();
    const STKScalarField* USTKFieldPtrOld =
        model_->URef().prevTimeRef().stkFieldPtr();

    const STKScalarField* divSTKFieldPtr =
        model_->mDotRef().divRef().stkFieldPtr();

    // Geometric fields
    const auto* volSTKFieldPtr = metaData.get_field<scalar>(
        stk::topology::NODE_RANK, this->getDualNodalVolumeID_(domain));

    // Body forces (use redistributed to balance pressure gradient)
    const auto* FSTKFieldPtr =
        metaData.get_field<scalar>(stk::topology::NODE_RANK, flowModel::F_ID);

    // time integrator
    const scalar dt = field_broker_->meshRef().controlsRef().getTimestep();
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
        scalar* gradPb = stk::mesh::field_data(*gradPSTKFieldPtr, nodeBucket);
        scalar* rhob = stk::mesh::field_data(*rhoSTKFieldPtr, nodeBucket);
        scalar* rhobOld = stk::mesh::field_data(*rhoSTKFieldPtrOld, nodeBucket);
        scalar* volb = stk::mesh::field_data(*volSTKFieldPtr, nodeBucket);
        scalar* divb = stk::mesh::field_data(*divSTKFieldPtr, nodeBucket);
        scalar* Ub = stk::mesh::field_data(*USTKFieldPtr, nodeBucket);
        scalar* UbOld = stk::mesh::field_data(*USTKFieldPtrOld, nodeBucket);
        scalar* Fb = stk::mesh::field_data(*FSTKFieldPtr, nodeBucket);

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
            scalar div = divb[iNode];

            scalar lhsfac = c[0] * rho * vol / dt;
            scalar lhsfacOld = c[1] * rhoOld * vol / dt;

            for (label i = 0; i < SPATIAL_DIM; i++)
            {
                scalar Ui = Ub[SPATIAL_DIM * iNode + i];
                scalar UiOld = UbOld[SPATIAL_DIM * iNode + i];

                // transient
                lhs[i * SPATIAL_DIM + i] += lhsfac;
                rhs[i] -= (lhsfac * Ui + lhsfacOld * UiOld);

                // divergence correction
                if (div < 0)
                {
                    lhs[i * SPATIAL_DIM + i] += -div;
                }
                rhs[i] -= -div * Ui;

                // pressure gradient
                scalar gradPi = gradPb[SPATIAL_DIM * iNode + i];
                rhs[i] -= gradPi * vol;

                // body forces
                scalar Fi = Fb[SPATIAL_DIM * iNode + i];
                rhs[i] += Fi * vol;
            }

            this->applyCoeff_(
                A, b, connectedNodes, scratchIds, scratchVals, rhs, lhs);
        }
    }
}

void navierStokesAssembler::assembleNodeTermsFusedSecondOrderUnsteady_(
    const domain* domain,
    Context* ctx)
{
    auto& mesh = model_->meshRef();
    Matrix& A = ctx->getAMatrix();
    Vector& b = ctx->getBVector();

    stk::mesh::MetaData& metaData = mesh.metaDataRef();
    stk::mesh::BulkData& bulkData = mesh.bulkDataRef();

    // space for LHS/RHS
    const label lhsSize = SPATIAL_DIM * SPATIAL_DIM;
    const label rhsSize = SPATIAL_DIM;
    std::vector<scalar> lhs(lhsSize);
    std::vector<scalar> rhs(rhsSize);
    std::vector<label> scratchIds(rhsSize);
    std::vector<scalar> scratchVals(rhsSize);
    std::vector<stk::mesh::Entity> connectedNodes(1);

    // pointers
    scalar* p_lhs = &lhs[0];
    scalar* p_rhs = &rhs[0];

    // Get fields
    const STKScalarField* gradPSTKFieldPtr =
        model_->pRef().gradRef().stkFieldPtr();
    const STKScalarField* rhoSTKFieldPtr = model_->rhoRef().stkFieldPtr();
    const STKScalarField* rhoSTKFieldPtrOld =
        model_->rhoRef().prevTimeRef().stkFieldPtr();
    const STKScalarField* rhoSTKFieldPtrOldOld =
        model_->rhoRef().prevTimeRef().prevTimeRef().stkFieldPtr();

    const STKScalarField* USTKFieldPtr = model_->URef().stkFieldPtr();
    const STKScalarField* USTKFieldPtrOld =
        model_->URef().prevTimeRef().stkFieldPtr();
    const STKScalarField* USTKFieldPtrOldOld =
        model_->URef().prevTimeRef().prevTimeRef().stkFieldPtr();

    const STKScalarField* divSTKFieldPtr =
        model_->mDotRef().divRef().stkFieldPtr();

    // Geometric fields
    const auto* volSTKFieldPtr = metaData.get_field<scalar>(
        stk::topology::NODE_RANK, this->getDualNodalVolumeID_(domain));

    // Body forces (use redistributed to balance pressure gradient)
    const auto* FSTKFieldPtr =
        metaData.get_field<scalar>(stk::topology::NODE_RANK, flowModel::F_ID);

    // time integrator
    const scalar dt = field_broker_->meshRef().controlsRef().getTimestep();
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
        scalar* gradPb = stk::mesh::field_data(*gradPSTKFieldPtr, nodeBucket);
        scalar* rhob = stk::mesh::field_data(*rhoSTKFieldPtr, nodeBucket);
        scalar* rhobOld = stk::mesh::field_data(*rhoSTKFieldPtrOld, nodeBucket);
        scalar* rhobOldOld =
            stk::mesh::field_data(*rhoSTKFieldPtrOldOld, nodeBucket);
        scalar* volb = stk::mesh::field_data(*volSTKFieldPtr, nodeBucket);
        scalar* divb = stk::mesh::field_data(*divSTKFieldPtr, nodeBucket);
        scalar* Ub = stk::mesh::field_data(*USTKFieldPtr, nodeBucket);
        scalar* UbOld = stk::mesh::field_data(*USTKFieldPtrOld, nodeBucket);
        scalar* UbOldOld =
            stk::mesh::field_data(*USTKFieldPtrOldOld, nodeBucket);
        scalar* Fb = stk::mesh::field_data(*FSTKFieldPtr, nodeBucket);

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
            scalar div = divb[iNode];

            scalar lhsfac = c[0] * rho * vol / dt;
            scalar lhsfacOld = c[1] * rhoOld * vol / dt;
            scalar lhsfacOldOld = c[2] * rhoOldOld * vol / dt;

            for (label i = 0; i < SPATIAL_DIM; i++)
            {
                scalar Ui = Ub[SPATIAL_DIM * iNode + i];
                scalar UiOld = UbOld[SPATIAL_DIM * iNode + i];
                scalar UiOldOld = UbOldOld[SPATIAL_DIM * iNode + i];

                // transient
                lhs[i * SPATIAL_DIM + i] += lhsfac;
                rhs[i] -=
                    (lhsfac * Ui + lhsfacOld * UiOld + lhsfacOldOld * UiOldOld);

                // divergence correction
                if (div < 0)
                {
                    lhs[i * SPATIAL_DIM + i] += -div;
                }
                rhs[i] -= -div * Ui;

                // pressure gradient
                scalar gradPi = gradPb[SPATIAL_DIM * iNode + i];
                rhs[i] -= gradPi * vol;

                // body forces
                scalar Fi = Fb[SPATIAL_DIM * iNode + i];
                rhs[i] += Fi * vol;
            }

            this->applyCoeff_(
                A, b, connectedNodes, scratchIds, scratchVals, rhs, lhs);
        }
    }
}

} // namespace accel
