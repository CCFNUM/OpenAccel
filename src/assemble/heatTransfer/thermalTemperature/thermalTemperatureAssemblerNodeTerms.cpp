// File       : thermalTemperatureAssemblerNodeTerms.cpp
// Created    : Thu Apr 14 2024 8:36:38 (+0100)
// Author     : Fabian Wermelinger
// Description:
// Copyright (c) 2024 CCFNUM, Lucerne University of Applied Sciences and Arts.
// SPDX-License-Identifier: BSD-3-Clause

#ifdef WITH_THERMAL_TEMPERATURE

#include "thermalTemperatureAssembler.h"

namespace accel
{

void thermalTemperatureAssembler::assembleNodeTermsFusedSteady_(
    const domain* domain,
    Context* ctx)
{
    if (domain->type() == domainType::solid)
        return;

    const bool compressible = domain->isMaterialCompressible();
    const bool includeViscousWork = domain->heatTransfer_.includeViscousWork_;
    const scalar comp = domain->isMaterialCompressible() ? 1.0 : 0.0;

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
    const STKScalarField* TSTKFieldPtr = phi_->stkFieldPtr();
    const STKScalarField* rhoSTKFieldPtr = model_->rhoRef().stkFieldPtr();
    const STKScalarField* cpSTKFieldPtr = model_->cpRef().stkFieldPtr();
    const STKScalarField* gradUSTKFieldPtr =
        includeViscousWork ? model_->URef().gradRef().stkFieldPtr() : nullptr;
    const STKScalarField* muEffSTKFieldPtr =
        includeViscousWork ? model_->muEffRef().stkFieldPtr() : nullptr;
    const STKScalarField* gradPSTKFieldPtr =
        compressible ? model_->pRef().gradRef().stkFieldPtr() : nullptr;
    const STKScalarField* USTKFieldPtr =
        compressible ? model_->URef().stkFieldPtr() : nullptr;

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
        scalar* cpb = stk::mesh::field_data(*cpSTKFieldPtr, nodeBucket);
        scalar* Tb = stk::mesh::field_data(*TSTKFieldPtr, nodeBucket);
        scalar* divb = stk::mesh::field_data(*divUSTKFieldPtr_, nodeBucket);
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
            scalar T = Tb[iNode];
            scalar rho = rhob[iNode];
            scalar cp = cpb[iNode];
            scalar div = divb[iNode];
            scalar vol = volb[iNode];

            // false transient
            scalar lhsfac = rho * cp * vol / dt;
            lhs[0] += lhsfac;

            // divergence correction
            if (div < 0)
            {
                lhs[0] += -div * cp;
            }
            rhs[0] -= -div * cp * T;

            // pressure work (warning: convenient only for ideal-gas)
            if (compressible)
            {
                const scalar* U = stk::mesh::field_data(*USTKFieldPtr, node);
                const scalar* gradP =
                    stk::mesh::field_data(*gradPSTKFieldPtr, node);

                scalar uDotGp = 0.0;
                for (label i = 0; i < SPATIAL_DIM; i++)
                {
                    uDotGp += U[i] * gradP[i];
                }

                rhs[0] += uDotGp * vol;
            }

            // viscous work
            if (includeViscousWork)
            {
                scalar* gradU = stk::mesh::field_data(*gradUSTKFieldPtr, node);
                scalar muEff = *stk::mesh::field_data(*muEffSTKFieldPtr, node);

                // form divU
                scalar divU = 0.0;
                for (label j = 0; j < SPATIAL_DIM; ++j)
                {
                    const label row = j * SPATIAL_DIM;
                    divU += gradU[row + j];
                }

                scalar viscousWork = 0.0;
                for (label i = 0; i < SPATIAL_DIM; ++i)
                {
                    const label offSet = SPATIAL_DIM * i;
                    for (label j = 0; j < SPATIAL_DIM; ++j)
                    {
                        viscousWork +=
                            gradU[offSet + j] *
                            (gradU[offSet + j] + gradU[SPATIAL_DIM * j + i]);
                        if (i == j)
                            viscousWork -=
                                gradU[offSet + j] * 2.0 / 3.0 * divU * comp;
                    }
                }
                viscousWork *= muEff;

                // assemble
                rhs[0] += viscousWork * vol;
            }

            Base::applyCoeff_(
                A, b, connectedNodes, scratchIds, scratchVals, rhs, lhs);
        }
    }
}

void thermalTemperatureAssembler::assembleNodeTermsFusedFirstOrderUnsteady_(
    const domain* domain,
    Context* ctx)
{
    const bool includeAdv = domain->type() == domainType::fluid;

    const bool compressible = domain->isMaterialCompressible();
    const bool includeViscousWork =
        domain->heatTransfer_.includeViscousWork_ && includeAdv;
    const scalar comp = domain->isMaterialCompressible() ? 1.0 : 0.0;

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
    const STKScalarField* cpSTKFieldPtr = model_->cpRef().stkFieldPtr();
    const STKScalarField* cpSTKFieldPtrOld =
        model_->cpRef().prevTimeRef().stkFieldPtr();
    const STKScalarField* pSTKFieldPtr =
        compressible ? model_->pRef().stkFieldPtr() : nullptr;
    const STKScalarField* pSTKFieldPtrOld =
        compressible ? model_->pRef().prevTimeRef().stkFieldPtr() : nullptr;
    const STKScalarField* gradPSTKFieldPtr =
        compressible ? model_->pRef().gradRef().stkFieldPtr() : nullptr;
    const STKScalarField* USTKFieldPtr =
        compressible ? model_->URef().stkFieldPtr() : nullptr;

    const STKScalarField* TSTKFieldPtr = phi_->stkFieldPtr();
    const STKScalarField* TSTKFieldPtrOld = phi_->prevTimeRef().stkFieldPtr();

    const STKScalarField* gradUSTKFieldPtr =
        includeViscousWork ? model_->URef().gradRef().stkFieldPtr() : nullptr;
    const STKScalarField* muEffSTKFieldPtr =
        includeViscousWork ? model_->muEffRef().stkFieldPtr() : nullptr;

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
        scalar* cpb = stk::mesh::field_data(*cpSTKFieldPtr, nodeBucket);
        scalar* cpbOld = stk::mesh::field_data(*cpSTKFieldPtrOld, nodeBucket);
        scalar* Tb = stk::mesh::field_data(*TSTKFieldPtr, nodeBucket);
        scalar* TbOld = stk::mesh::field_data(*TSTKFieldPtrOld, nodeBucket);
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
            scalar T = Tb[iNode];
            scalar TOld = TbOld[iNode];
            scalar rho = rhob[iNode];
            scalar rhoOld = rhobOld[iNode];
            scalar cp = cpb[iNode];
            scalar cpOld = cpbOld[iNode];
            scalar vol = volb[iNode];

            scalar div = includeAdv ? *stk::mesh::field_data(
                                          *divUSTKFieldPtr_, nodeBucket, iNode)
                                    : 0.0;

            // transient
            scalar lhsfac = c[0] * rho * cp * vol / dt;
            scalar lhsfacOld = c[1] * rhoOld * cpOld * vol / dt;
            lhs[0] += lhsfac;
            rhs[0] -= (lhsfac * T + lhsfacOld * TOld);

            // divergence correction
            if (div < 0)
            {
                lhs[0] += -div * cp;
            }
            rhs[0] -= -div * cp * T;

            // pressure work (only for compressible) + pressure rate of change
            if (compressible)
            {
                scalar* U =
                    stk::mesh::field_data(*USTKFieldPtr, nodeBucket, iNode);
                scalar* gradP =
                    stk::mesh::field_data(*gradPSTKFieldPtr, nodeBucket, iNode);
                scalar p =
                    *stk::mesh::field_data(*pSTKFieldPtr, nodeBucket, iNode);
                scalar pOld =
                    *stk::mesh::field_data(*pSTKFieldPtrOld, nodeBucket, iNode);

                scalar uDotGp = 0.0;
                for (label i = 0; i < SPATIAL_DIM; i++)
                {
                    uDotGp += U[i] * gradP[i];
                }

                rhs[0] += (uDotGp + (c[0] * p + c[1] * pOld) / dt) * vol;
            }

            // viscous work
            if (includeViscousWork)
            {
                scalar* gradU = stk::mesh::field_data(*gradUSTKFieldPtr, node);
                scalar muEff = *stk::mesh::field_data(*muEffSTKFieldPtr, node);

                // form divU
                scalar divU = 0.0;
                for (label j = 0; j < SPATIAL_DIM; ++j)
                {
                    const label row = j * SPATIAL_DIM;
                    divU += gradU[row + j];
                }

                scalar viscousWork = 0.0;
                for (label i = 0; i < SPATIAL_DIM; ++i)
                {
                    const label offSet = SPATIAL_DIM * i;
                    for (label j = 0; j < SPATIAL_DIM; ++j)
                    {
                        viscousWork +=
                            gradU[offSet + j] *
                            (gradU[offSet + j] + gradU[SPATIAL_DIM * j + i]);
                        if (i == j)
                            viscousWork -=
                                gradU[offSet + j] * 2.0 / 3.0 * divU * comp;
                    }
                }
                viscousWork *= muEff;

                // assemble
                rhs[0] += viscousWork * vol;
            }

            Base::applyCoeff_(
                A, b, connectedNodes, scratchIds, scratchVals, rhs, lhs);
        }
    }
}

void thermalTemperatureAssembler::assembleNodeTermsFusedSecondOrderUnsteady_(
    const domain* domain,
    Context* ctx)
{
    const bool includeAdv = domain->type() == domainType::fluid;

    const bool compressible = domain->isMaterialCompressible();
    const bool includeViscousWork =
        domain->heatTransfer_.includeViscousWork_ && includeAdv;
    const scalar comp = domain->isMaterialCompressible() ? 1.0 : 0.0;

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
    const STKScalarField* cpSTKFieldPtr = model_->cpRef().stkFieldPtr();
    const STKScalarField* cpSTKFieldPtrOld =
        model_->cpRef().prevTimeRef().stkFieldPtr();
    const STKScalarField* cpSTKFieldPtrOldOld =
        model_->cpRef().prevTimeRef().prevTimeRef().stkFieldPtr();
    const STKScalarField* pSTKFieldPtr =
        compressible ? model_->pRef().stkFieldPtr() : nullptr;
    const STKScalarField* pSTKFieldPtrOld =
        compressible ? model_->pRef().prevTimeRef().stkFieldPtr() : nullptr;
    const STKScalarField* pSTKFieldPtrOldOld =
        compressible ? model_->pRef().prevTimeRef().prevTimeRef().stkFieldPtr()
                     : nullptr;
    const STKScalarField* gradPSTKFieldPtr =
        compressible ? model_->pRef().gradRef().stkFieldPtr() : nullptr;
    const STKScalarField* USTKFieldPtr =
        compressible ? model_->URef().stkFieldPtr() : nullptr;

    const STKScalarField* TSTKFieldPtr = phi_->stkFieldPtr();
    const STKScalarField* TSTKFieldPtrOld = phi_->prevTimeRef().stkFieldPtr();
    const STKScalarField* TSTKFieldPtrOldOld =
        phi_->prevTimeRef().prevTimeRef().stkFieldPtr();

    const STKScalarField* gradUSTKFieldPtr =
        includeViscousWork ? model_->URef().gradRef().stkFieldPtr() : nullptr;
    const STKScalarField* muEffSTKFieldPtr =
        includeViscousWork ? model_->muEffRef().stkFieldPtr() : nullptr;

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
        scalar* cpb = stk::mesh::field_data(*cpSTKFieldPtr, nodeBucket);
        scalar* cpbOld = stk::mesh::field_data(*cpSTKFieldPtrOld, nodeBucket);
        scalar* cpbOldOld =
            stk::mesh::field_data(*cpSTKFieldPtrOldOld, nodeBucket);
        scalar* Tb = stk::mesh::field_data(*TSTKFieldPtr, nodeBucket);
        scalar* TbOld = stk::mesh::field_data(*TSTKFieldPtrOld, nodeBucket);
        scalar* TbOldOld =
            stk::mesh::field_data(*TSTKFieldPtrOldOld, nodeBucket);
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
            scalar T = Tb[iNode];
            scalar TOld = TbOld[iNode];
            scalar TOldOld = TbOldOld[iNode];
            scalar rho = rhob[iNode];
            scalar rhoOld = rhobOld[iNode];
            scalar rhoOldOld = rhobOldOld[iNode];
            scalar cp = cpb[iNode];
            scalar cpOld = cpbOld[iNode];
            scalar cpOldOld = cpbOldOld[iNode];
            scalar vol = volb[iNode];

            scalar div = includeAdv ? *stk::mesh::field_data(
                                          *divUSTKFieldPtr_, nodeBucket, iNode)
                                    : 0.0;

            // transient
            scalar lhsfac = c[0] * rho * cp * vol / dt;
            scalar lhsfacOld = c[1] * rhoOld * cpOld * vol / dt;
            scalar lhsfacOldOld = c[2] * rhoOldOld * cpOldOld * vol / dt;
            lhs[0] += lhsfac;
            rhs[0] -= (lhsfac * T + lhsfacOld * TOld + lhsfacOldOld * TOldOld);

            // divergence correction
            if (div < 0)
            {
                lhs[0] += -div * cp;
            }
            rhs[0] -= -div * cp * T;

            // pressure work (only for compressible) + pressure rate of change
            if (compressible)
            {
                scalar* U =
                    stk::mesh::field_data(*USTKFieldPtr, nodeBucket, iNode);
                scalar* gradP =
                    stk::mesh::field_data(*gradPSTKFieldPtr, nodeBucket, iNode);
                scalar p =
                    *stk::mesh::field_data(*pSTKFieldPtr, nodeBucket, iNode);
                scalar pOld =
                    *stk::mesh::field_data(*pSTKFieldPtrOld, nodeBucket, iNode);
                scalar pOldOld = *stk::mesh::field_data(
                    *pSTKFieldPtrOldOld, nodeBucket, iNode);

                scalar uDotGp = 0.0;
                for (label i = 0; i < SPATIAL_DIM; i++)
                {
                    uDotGp += U[i] * gradP[i];
                }

                rhs[0] +=
                    (uDotGp + (c[0] * p + c[1] * pOld + c[2] * pOldOld) / dt) *
                    vol;
            }

            // viscous work
            if (includeViscousWork)
            {
                scalar* gradU = stk::mesh::field_data(*gradUSTKFieldPtr, node);
                scalar muEff = *stk::mesh::field_data(*muEffSTKFieldPtr, node);

                // form divU
                scalar divU = 0.0;
                for (label j = 0; j < SPATIAL_DIM; ++j)
                {
                    const label row = j * SPATIAL_DIM;
                    divU += gradU[row + j];
                }

                scalar viscousWork = 0.0;
                for (label i = 0; i < SPATIAL_DIM; ++i)
                {
                    const label offSet = SPATIAL_DIM * i;
                    for (label j = 0; j < SPATIAL_DIM; ++j)
                    {
                        viscousWork +=
                            gradU[offSet + j] *
                            (gradU[offSet + j] + gradU[SPATIAL_DIM * j + i]);
                        if (i == j)
                            viscousWork -=
                                gradU[offSet + j] * 2.0 / 3.0 * divU * comp;
                    }
                }
                viscousWork *= muEff;

                // assemble
                rhs[0] += viscousWork * vol;
            }

            Base::applyCoeff_(
                A, b, connectedNodes, scratchIds, scratchVals, rhs, lhs);
        }
    }
}

} // namespace accel

#endif
