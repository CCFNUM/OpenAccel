// File       : thermalEnergyAssemblerNodeTerms.cpp
// Created    : Mon Apr 14 2025 10:36:38 (+0100)
// Author     : Mhamad Mahdi Alloush
// Description:
// Copyright (c) 2025 CCFNUM, Lucerne University of Applied Sciences and Arts.
// SPDX-License-Identifier: BSD-3-Clause

#ifndef WITH_THERMAL_TEMPERATURE

#include "thermalEnergyAssembler.h"

namespace accel
{

void thermalEnergyAssembler::assembleNodeTermsFusedSteady_(const domain* domain,
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
        const bool materialCompressible = domain->isMaterialCompressible();
        const bool includeViscousWork =
            domain->heatTransfer_.includeViscousWork_;
        const bool includePressureWork =
            materialCompressible && domain->heatTransfer_.includePressureWork_;
        const scalar comp = materialCompressible ? 1.0 : 0.0;

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

        // Get fields (solving for enthalpy, not temperature)
        const STKScalarField* hSTKFieldPtr = phi_->stkFieldPtr();
        const STKScalarField* rhoSTKFieldPtr = model_->rhoRef().stkFieldPtr();
        const STKScalarField* gradUSTKFieldPtr =
            includeViscousWork ? model_->URef().gradRef().stkFieldPtr()
                               : nullptr;
        const STKScalarField* muEffSTKFieldPtr =
            includeViscousWork ? model_->muEffRef().stkFieldPtr() : nullptr;
        const STKScalarField* gradPSTKFieldPtr =
            includePressureWork ? model_->pRef().gradRef().stkFieldPtr()
                                : nullptr;
        const STKScalarField* USTKFieldPtr =
            includePressureWork ? model_->URef().stkFieldPtr() : nullptr;

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
            scalar* rhob = stk::mesh::field_data(*rhoSTKFieldPtr, nodeBucket);
            scalar* hb = stk::mesh::field_data(*hSTKFieldPtr, nodeBucket);
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
                scalar h = hb[iNode];
                scalar rho = rhob[iNode];
                scalar div = divb[iNode];
                scalar vol = volb[iNode];

                // false transient: d(rho*h)/dt (no cp multiplier for enthalpy)
                scalar lhsfac = rho * vol / dt;
                lhs[0] += lhsfac;

                // divergence correction
                if (div < 0)
                {
                    lhs[0] += -div;
                }
                rhs[0] -= -div * h;

                if (includePressureWork)
                {
                    const scalar* U =
                        stk::mesh::field_data(*USTKFieldPtr, node);
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
                    scalar* gradU =
                        stk::mesh::field_data(*gradUSTKFieldPtr, node);
                    scalar muEff =
                        *stk::mesh::field_data(*muEffSTKFieldPtr, node);

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
                            viscousWork += gradU[offSet + j] *
                                           (gradU[offSet + j] +
                                            gradU[SPATIAL_DIM * j + i]);
                            if (i == j)
                                viscousWork -=
                                    gradU[offSet + j] * 2.0 / 3.0 * divU * comp;
                        }
                    }
                    viscousWork *= muEff;

                    // assemble
                    rhs[0] += viscousWork * vol;
                }

                // energy source
                rhs[0] += energySourceValue * vol;

                Base::applyCoeff_(
                    A, b, connectedNodes, scratchIds, scratchVals, rhs, lhs);
            }
        }
    }
}

void thermalEnergyAssembler::assembleNodeTermsFusedFirstOrderUnsteady_(
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

    const bool materialCompressible = domain->isMaterialCompressible();
    const bool includeViscousWork =
        domain->heatTransfer_.includeViscousWork_ && includeAdv;
    const bool includePressureWork =
        materialCompressible && domain->heatTransfer_.includePressureWork_;
    const bool includeLowSpeedCompressibility =
        materialCompressible &&
        domain->heatTransfer_.includeLowSpeedCompressibility_;
    const scalar comp = materialCompressible ? 1.0 : 0.0;

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
    const STKScalarField* pSTKFieldPtr =
        includeLowSpeedCompressibility ? model_->pRef().stkFieldPtr() : nullptr;
    const STKScalarField* pSTKFieldPtrOld =
        includeLowSpeedCompressibility
            ? model_->pRef().prevTimeRef().stkFieldPtr()
            : nullptr;
    const STKScalarField* gradPSTKFieldPtr =
        includePressureWork ? model_->pRef().gradRef().stkFieldPtr() : nullptr;
    const STKScalarField* USTKFieldPtr =
        includePressureWork ? model_->URef().stkFieldPtr() : nullptr;

    const STKScalarField* hSTKFieldPtr = phi_->stkFieldPtr();
    const STKScalarField* hSTKFieldPtrOld = phi_->prevTimeRef().stkFieldPtr();

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
        scalar* hb = stk::mesh::field_data(*hSTKFieldPtr, nodeBucket);
        scalar* hbOld = stk::mesh::field_data(*hSTKFieldPtrOld, nodeBucket);
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
            scalar h = hb[iNode];
            scalar hOld = hbOld[iNode];
            scalar rho = rhob[iNode];
            scalar rhoOld = rhobOld[iNode];
            scalar vol = volb[iNode];

            scalar div = includeAdv ? *stk::mesh::field_data(
                                          *divUSTKFieldPtr_, nodeBucket, iNode)
                                    : 0.0;

            // transient: d(rho*h)/dt (no cp multiplier for enthalpy)
            scalar lhsfac = c[0] * rho * vol / dt;
            scalar lhsfacOld = c[1] * rhoOld * vol / dt;
            lhs[0] += lhsfac;
            rhs[0] -= (lhsfac * h + lhsfacOld * hOld);

            // divergence correction
            if (div < 0)
            {
                lhs[0] += -div;
            }
            rhs[0] -= -div * h;

            scalar pressureSource = 0.0;
            if (includePressureWork)
            {
                scalar* U =
                    stk::mesh::field_data(*USTKFieldPtr, nodeBucket, iNode);
                scalar* gradP =
                    stk::mesh::field_data(*gradPSTKFieldPtr, nodeBucket, iNode);

                scalar uDotGp = 0.0;
                for (label i = 0; i < SPATIAL_DIM; i++)
                {
                    uDotGp += U[i] * gradP[i];
                }

                pressureSource += uDotGp;
            }

            if (includeLowSpeedCompressibility)
            {
                const scalar p =
                    *stk::mesh::field_data(*pSTKFieldPtr, nodeBucket, iNode);
                const scalar pOld =
                    *stk::mesh::field_data(*pSTKFieldPtrOld, nodeBucket, iNode);
                pressureSource += (c[0] * p + c[1] * pOld) / dt;
            }

            if (pressureSource != 0.0)
            {
                rhs[0] += pressureSource * vol;
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

            // energy source
            rhs[0] += energySourceValue * vol;

            Base::applyCoeff_(
                A, b, connectedNodes, scratchIds, scratchVals, rhs, lhs);
        }
    }
}

void thermalEnergyAssembler::assembleNodeTermsFusedSecondOrderUnsteady_(
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

    const bool materialCompressible = domain->isMaterialCompressible();
    const bool includeViscousWork =
        domain->heatTransfer_.includeViscousWork_ && includeAdv;
    const bool includePressureWork =
        materialCompressible && domain->heatTransfer_.includePressureWork_;
    const bool includeLowSpeedCompressibility =
        materialCompressible &&
        domain->heatTransfer_.includeLowSpeedCompressibility_;
    const scalar comp = materialCompressible ? 1.0 : 0.0;

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
    const STKScalarField* pSTKFieldPtr =
        includeLowSpeedCompressibility ? model_->pRef().stkFieldPtr() : nullptr;
    const STKScalarField* pSTKFieldPtrOld =
        includeLowSpeedCompressibility
            ? model_->pRef().prevTimeRef().stkFieldPtr()
            : nullptr;
    const STKScalarField* pSTKFieldPtrOldOld =
        includeLowSpeedCompressibility
            ? model_->pRef().prevTimeRef().prevTimeRef().stkFieldPtr()
            : nullptr;
    const STKScalarField* gradPSTKFieldPtr =
        includePressureWork ? model_->pRef().gradRef().stkFieldPtr() : nullptr;
    const STKScalarField* USTKFieldPtr =
        includePressureWork ? model_->URef().stkFieldPtr() : nullptr;

    const STKScalarField* hSTKFieldPtr = phi_->stkFieldPtr();
    const STKScalarField* hSTKFieldPtrOld = phi_->prevTimeRef().stkFieldPtr();
    const STKScalarField* hSTKFieldPtrOldOld =
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
        scalar* hb = stk::mesh::field_data(*hSTKFieldPtr, nodeBucket);
        scalar* hbOld = stk::mesh::field_data(*hSTKFieldPtrOld, nodeBucket);
        scalar* hbOldOld =
            stk::mesh::field_data(*hSTKFieldPtrOldOld, nodeBucket);
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
            scalar h = hb[iNode];
            scalar hOld = hbOld[iNode];
            scalar hOldOld = hbOldOld[iNode];
            scalar rho = rhob[iNode];
            scalar rhoOld = rhobOld[iNode];
            scalar rhoOldOld = rhobOldOld[iNode];
            scalar vol = volb[iNode];

            scalar div = includeAdv ? *stk::mesh::field_data(
                                          *divUSTKFieldPtr_, nodeBucket, iNode)
                                    : 0.0;

            // transient: d(rho*h)/dt (no cp multiplier for enthalpy)
            scalar lhsfac = c[0] * rho * vol / dt;
            scalar lhsfacOld = c[1] * rhoOld * vol / dt;
            scalar lhsfacOldOld = c[2] * rhoOldOld * vol / dt;
            lhs[0] += lhsfac;
            rhs[0] -= (lhsfac * h + lhsfacOld * hOld + lhsfacOldOld * hOldOld);

            // divergence correction
            if (div < 0)
            {
                lhs[0] += -div;
            }
            rhs[0] -= -div * h;

            scalar pressureSource = 0.0;
            if (includePressureWork)
            {
                scalar* U =
                    stk::mesh::field_data(*USTKFieldPtr, nodeBucket, iNode);
                scalar* gradP =
                    stk::mesh::field_data(*gradPSTKFieldPtr, nodeBucket, iNode);

                scalar uDotGp = 0.0;
                for (label i = 0; i < SPATIAL_DIM; i++)
                {
                    uDotGp += U[i] * gradP[i];
                }

                pressureSource += uDotGp;
            }

            if (includeLowSpeedCompressibility)
            {
                const scalar p =
                    *stk::mesh::field_data(*pSTKFieldPtr, nodeBucket, iNode);
                const scalar pOld =
                    *stk::mesh::field_data(*pSTKFieldPtrOld, nodeBucket, iNode);
                const scalar pOldOld = *stk::mesh::field_data(
                    *pSTKFieldPtrOldOld, nodeBucket, iNode);
                pressureSource +=
                    (c[0] * p + c[1] * pOld + c[2] * pOldOld) / dt;
            }

            if (pressureSource != 0.0)
            {
                rhs[0] += pressureSource * vol;
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

            // energy source
            rhs[0] += energySourceValue * vol;

            Base::applyCoeff_(
                A, b, connectedNodes, scratchIds, scratchVals, rhs, lhs);
        }
    }
}

} // namespace accel

#endif
