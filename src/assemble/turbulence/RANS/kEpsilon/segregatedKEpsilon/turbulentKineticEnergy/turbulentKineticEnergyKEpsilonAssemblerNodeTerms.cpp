// File : turbulentKineticEnergyKEpsilonAssemblerNodeTerms.cpp
// Created : Thu Feb 22 2025 13:38:51 (+0100)
// Author : Achraf Nagihi
// Description:
// Copyright (c) 2025 CCFNUM, Lucerne University of Applied Sciences and Arts.
// SPDX-License-Identifier: BSD-3-Clause

#include "turbulentKineticEnergyKEpsilonAssembler.h"

namespace accel
{

void turbulentKineticEnergyKEpsilonAssembler::assembleNodeTermsFusedSteady_(
    const domain* domain,
    Context* ctx)
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
    const STKScalarField* epsilonSTKFieldPtr =
        model_->epsilonRef().stkFieldPtr();
    const STKScalarField* rhoSTKFieldPtr = model_->rhoRef().stkFieldPtr();
    const STKScalarField* divSTKFieldPtr =
        model_->mDotRef().divRef().stkFieldPtr();
    const STKScalarField* mutSTKFieldPtr = model_->mutRef().stkFieldPtr();
    const STKScalarField* PkSTKFieldPtr = model_->PkRef().stkFieldPtr();

    // Geometric fields
    const auto* volSTKFieldPtr = metaData.get_field<scalar>(
        stk::topology::NODE_RANK, this->getDualNodalVolumeID_(domain));

    scalar tkeProdLimitRatio = model_->tkeProdLimitRatio();

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
            scalar epsilon = *stk::mesh::field_data(*epsilonSTKFieldPtr, node);
            scalar rho = *stk::mesh::field_data(*rhoSTKFieldPtr, node);
            scalar vol = *stk::mesh::field_data(*volSTKFieldPtr, node);
            scalar div = *stk::mesh::field_data(*divSTKFieldPtr, node);
            scalar Pk = *stk::mesh::field_data(*PkSTKFieldPtr, node);
            scalar mut = *stk::mesh::field_data(*mutSTKFieldPtr, node);

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

            // clip
            Pk = std::max(0.0, Pk);

            const scalar Dk = rho * epsilon;
            if (Pk > tkeProdLimitRatio * Dk)
            {
                Pk = tkeProdLimitRatio * Dk;
            }

            rhs[0] += (Pk - Dk) * vol;
            lhs[0] += rho * epsilon / (k + SMALL) * vol;

            Base::applyCoeff_(
                A, b, connectedNodes, scratchIds, scratchVals, rhs, lhs);
        }
    }
}

void turbulentKineticEnergyKEpsilonAssembler::
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
    const STKScalarField* epsilonSTKFieldPtr =
        model_->epsilonRef().stkFieldPtr();
    const STKScalarField* rhoSTKFieldPtr = model_->rhoRef().stkFieldPtr();
    const STKScalarField* rhoSTKFieldPtrOld =
        model_->rhoRef().prevTimeRef().stkFieldPtr();
    const STKScalarField* divSTKFieldPtr =
        model_->mDotRef().divRef().stkFieldPtr();
    const STKScalarField* mutSTKFieldPtr = model_->mutRef().stkFieldPtr();
    const STKScalarField* PkSTKFieldPtr = model_->PkRef().stkFieldPtr();

    const STKScalarField* divUmSTKFieldPtr =
        meshDeforming ? model_->divUmRef().stkFieldPtr() : nullptr;

    // Geometric fields
    const auto* volSTKFieldPtr = metaData.get_field<scalar>(
        stk::topology::NODE_RANK, this->getDualNodalVolumeID_(domain));

    scalar tkeProdLimitRatio = model_->tkeProdLimitRatio();

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
            scalar k = *stk::mesh::field_data(*kSTKFieldPtr, node);
            scalar kOld = *stk::mesh::field_data(*kSTKFieldPtrOld, node);
            scalar epsilon = *stk::mesh::field_data(*epsilonSTKFieldPtr, node);
            scalar rho = *stk::mesh::field_data(*rhoSTKFieldPtr, node);
            scalar rhoOld = *stk::mesh::field_data(*rhoSTKFieldPtrOld, node);
            scalar vol = *stk::mesh::field_data(*volSTKFieldPtr, node);
            scalar div = *stk::mesh::field_data(*divSTKFieldPtr, node);
            scalar Pk = *stk::mesh::field_data(*PkSTKFieldPtr, node);
            scalar mut = *stk::mesh::field_data(*mutSTKFieldPtr, node);

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

            // production/dissipation

            // clip
            Pk = std::max(0.0, Pk);

            const scalar Dk = rho * epsilon;
            if (Pk > tkeProdLimitRatio * Dk)
            {
                Pk = tkeProdLimitRatio * Dk;
            }

            rhs[0] += (Pk - Dk) * vol;
            lhs[0] += rho * epsilon / (k + SMALL) * vol;

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

void turbulentKineticEnergyKEpsilonAssembler::
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
    const STKScalarField* epsilonSTKFieldPtr =
        model_->epsilonRef().stkFieldPtr();
    const STKScalarField* rhoSTKFieldPtr = model_->rhoRef().stkFieldPtr();
    const STKScalarField* rhoSTKFieldPtrOld =
        model_->rhoRef().prevTimeRef().stkFieldPtr();
    const STKScalarField* rhoSTKFieldPtrOldOld =
        model_->rhoRef().prevTimeRef().prevTimeRef().stkFieldPtr();
    const STKScalarField* divSTKFieldPtr =
        model_->mDotRef().divRef().stkFieldPtr();
    const STKScalarField* mutSTKFieldPtr = model_->mutRef().stkFieldPtr();
    const STKScalarField* PkSTKFieldPtr = model_->PkRef().stkFieldPtr();

    const STKScalarField* divUmSTKFieldPtr =
        meshDeforming ? model_->divUmRef().stkFieldPtr() : nullptr;

    // Geometric fields
    const auto* volSTKFieldPtr = metaData.get_field<scalar>(
        stk::topology::NODE_RANK, this->getDualNodalVolumeID_(domain));

    scalar tkeProdLimitRatio = model_->tkeProdLimitRatio();

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
            scalar k = *stk::mesh::field_data(*kSTKFieldPtr, node);
            scalar kOld = *stk::mesh::field_data(*kSTKFieldPtrOld, node);
            scalar kOldOld = *stk::mesh::field_data(*kSTKFieldPtrOldOld, node);
            scalar epsilon = *stk::mesh::field_data(*epsilonSTKFieldPtr, node);
            scalar rho = *stk::mesh::field_data(*rhoSTKFieldPtr, node);
            scalar rhoOld = *stk::mesh::field_data(*rhoSTKFieldPtrOld, node);
            scalar rhoOldOld =
                *stk::mesh::field_data(*rhoSTKFieldPtrOldOld, node);
            scalar vol = *stk::mesh::field_data(*volSTKFieldPtr, node);
            scalar div = *stk::mesh::field_data(*divSTKFieldPtr, node);
            scalar Pk = *stk::mesh::field_data(*PkSTKFieldPtr, node);
            scalar mut = *stk::mesh::field_data(*mutSTKFieldPtr, node);

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

            // production/dissipation

            // clip
            Pk = std::max(0.0, Pk);

            const scalar Dk = rho * epsilon;
            if (Pk > tkeProdLimitRatio * Dk)
            {
                Pk = tkeProdLimitRatio * Dk;
            }

            rhs[0] += (Pk - Dk) * vol;
            lhs[0] += rho * epsilon / (k + SMALL) * vol;

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
