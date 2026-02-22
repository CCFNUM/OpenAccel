// File : navierStokesAssembler.cpp
// Created : Wed Jan 03 2024 13:38:51 (+0100)
// Author : Mhamad Mahdi Alloush
// Description:
// Copyright (c) 2024 CCFNUM, Lucerne University of Applied Sciences and Arts.
// SPDX-License-Identifier: BSD-3-Clause

#include "navierStokesAssembler.h"
#include "flowModel.h"

namespace accel
{

navierStokesAssembler::navierStokesAssembler(flowModel* model)
    : phiAssembler<SPATIAL_DIM>(static_cast<fieldBroker*>(model)), model_(model)
{
}

void navierStokesAssembler::setupDUCoefficients(const domain* domain)
{
    auto& mesh = model_->meshRef();
    stk::mesh::MetaData& metaData = mesh.metaDataRef();

    stk::mesh::Field<scalar>* duSTKFieldPtr =
        metaData.get_field<scalar>(stk::topology::NODE_RANK, flowModel::du_ID);

    if (!duSTKFieldPtr)
    {
        duSTKFieldPtr = &metaData.declare_field<scalar>(
            stk::topology::NODE_RANK, flowModel::du_ID);

        // Set field output type
        stk::io::set_field_output_type(*duSTKFieldPtr, fieldType[SPATIAL_DIM]);
    }

    // Put the stk field on interior mesh parts for the current domain
    const stk::mesh::PartVector& partVec = domain->zonePtr()->interiorParts();
    for (const stk::mesh::Part* part : partVec)
    {
        // check if already defined from a previous pass
        if (!duSTKFieldPtr->defined_on(*part))
        {
            stk::mesh::put_field_on_mesh(
                *duSTKFieldPtr, *part, SPATIAL_DIM, nullptr);
        }
    }

    // consistent (SIMPLEC)?
    bool consistent = model_->controlsRef()
                          .solverRef()
                          .solverControl_.expertParameters_.consistent_;

    if (consistent)
    {
        stk::mesh::Field<scalar>* duTildeSTKFieldPtr =
            metaData.get_field<scalar>(stk::topology::NODE_RANK,
                                       flowModel::duTilde_ID);

        if (!duTildeSTKFieldPtr)
        {
            duTildeSTKFieldPtr = &metaData.declare_field<scalar>(
                stk::topology::NODE_RANK, flowModel::duTilde_ID);

            // Set field output type
            stk::io::set_field_output_type(*duTildeSTKFieldPtr,
                                           fieldType[SPATIAL_DIM]);
        }

        // Put the stk field on interior mesh parts for the current domain
        const stk::mesh::PartVector& partVec =
            domain->zonePtr()->interiorParts();
        for (const stk::mesh::Part* part : partVec)
        {
            // check if already defined from a previous pass
            if (!duTildeSTKFieldPtr->defined_on(*part))
            {
                stk::mesh::put_field_on_mesh(
                    *duTildeSTKFieldPtr, *part, SPATIAL_DIM, nullptr);
            }
        }
    }
}

void navierStokesAssembler::computeDUCoefficients(const domain* domain,
                                                  Context* ctx)
{
    auto& mesh = model_->meshRef();
    stk::mesh::BulkData& bulkData = mesh.bulkDataRef();
    stk::mesh::MetaData& metaData = mesh.metaDataRef();

    // other
    const bool is_transient = model_->controlsRef().isTransient();
    const scalar dt = is_transient
                          ? model_->controlsRef().getTimestep()
                          : model_->controlsRef().getPhysicalTimescale();

    // fractional-step method
    scalar fsmflag =
        model_->controlsRef()
                .solverRef()
                .solverControl_.expertParameters_.fractionalStepMethod_
            ? 1.0
            : 0.0;

    auto scheme = model_->controlsRef()
                      .solverRef()
                      .solverControl_.basicSettings_.transientScheme_;

    scalar gamma1 = 1.0; // steady
    if (is_transient && scheme == transientSchemeType::firstOrderBackwardEuler)
    {
        const auto c = BDF1::coeff();
        gamma1 = c[0];
    }
    else if (is_transient &&
             scheme == transientSchemeType::secondOrderBackwardEuler)
    {
        const auto c = BDF2::coeff(dt, model_->controlsRef().getTimestep(-1));
        gamma1 = c[0];
    }

    // consistent (SIMPLEC)?
    bool consistent = model_->controlsRef()
                          .solverRef()
                          .solverControl_.expertParameters_.consistent_;

    // get matrix
    Matrix& A = ctx->getAMatrix();

    // get DU field
    auto* duSTKFieldPtr =
        metaData.get_field<scalar>(stk::topology::NODE_RANK, flowModel::du_ID);

    stk::mesh::Field<scalar>* duTildeSTKFieldPtr = nullptr;
    if (consistent)
    {
        duTildeSTKFieldPtr = metaData.get_field<scalar>(
            stk::topology::NODE_RANK, flowModel::duTilde_ID);
    }

    // get density field
    const auto* rhoSTKFieldPtr = model_->rhoRef().stkFieldPtr();

    // Get volumes
    const auto* volSTKFieldPtr = metaData.get_field<scalar>(
        stk::topology::NODE_RANK, this->getDualNodalVolumeID_(domain));

    // get interior parts of the domain it is defined on
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
        scalar* dub = stk::mesh::field_data(*duSTKFieldPtr, nodeBucket);
        const scalar* volb = stk::mesh::field_data(*volSTKFieldPtr, nodeBucket);
        const scalar* rhob = stk::mesh::field_data(*rhoSTKFieldPtr, nodeBucket);

        for (stk::mesh::Bucket::size_type iNode = 0; iNode < nNodesPerBucket;
             ++iNode)
        {
            // get node
            stk::mesh::Entity node = nodeBucket[iNode];

            label localID = bulkData.local_id(node);

            auto rowVals = A.rowVals(localID);
            auto rowCols = A.rowCols(localID);

            scalar vol = volb[iNode];
            scalar rho = rhob[iNode];

            for (label i = 0; i < SPATIAL_DIM; i++)
            {
                const scalar& di = A.dofDiag(localID, i);

                dub[SPATIAL_DIM * iNode + i] =
                    fsmflag * dt / gamma1 / rho +
                    (1.0 - fsmflag) * vol / (di + SMALL);
            }

            if (consistent)
            {
                scalar* duTilde =
                    stk::mesh::field_data(*duTildeSTKFieldPtr, node);

                const auto& diagOffset = A.diagOffsetRef();

                for (label i = 0; i < SPATIAL_DIM; i++)
                {
                    const scalar& di = A.dofDiag(localID, i);

                    scalar sumOffDiagi = 0.0;
                    for (label icol = 0;
                         icol < static_cast<label>(rowCols.size());
                         icol++)
                    {
                        if (diagOffset[localID] == icol)
                            continue;

                        sumOffDiagi += rowVals[BLOCKSIZE * BLOCKSIZE * icol +
                                               BLOCKSIZE * i + i];
                    }

                    duTilde[i] = vol / (di + sumOffDiagi + SMALL);
                }
            }
        }
    }

    // synchronize DU coefficients
    if (messager::parallel())
    {
        stk::mesh::communicate_field_data(bulkData, {duSTKFieldPtr});

        if (consistent)
        {
            stk::mesh::communicate_field_data(bulkData, {duTildeSTKFieldPtr});
        }
    }
}

void navierStokesAssembler::postAssemble_(const domain* domain, Context* ctx)
{
    phiAssembler::postAssemble_(domain, ctx);
    assembleBoundaryRelaxation_(domain, ctx->getBVector(), 0.75);
    applySymmetryConditions_(domain, ctx->getBVector());
}

void navierStokesAssembler::assembleBoundaryRelaxation_(const domain* domain,
                                                        Vector& b,
                                                        const scalar urf)
{
    using Bucket = stk::mesh::Bucket;
    using BucketVec = stk::mesh::BucketVector;

    // select all locally owned nodes for this domain
    const auto& mesh = field_broker_->meshRef();
    const stk::mesh::MetaData& metaData = mesh.metaDataRef();
    const stk::mesh::BulkData& bulkData = mesh.bulkDataRef();

    const zone* zonePtr = domain->zonePtr();

    // Collect all boundary parts: It is ensured in this way that no
    // duplicate relaxation or more are made to corner nodes (belong to two or
    // more boundaries)
    stk::mesh::PartVector partVec;

    for (label iBoundary = 0; iBoundary < zonePtr->nBoundaries(); iBoundary++)
    {
        const auto& boundaryRef = zonePtr->boundaryRef(iBoundary);
        const stk::mesh::PartVector& parts = boundaryRef.parts();

        boundaryPhysicalType type = boundaryRef.type();
        switch (type)
        {
            case boundaryPhysicalType::wall:
                {
                    if (domain->type() == domainType::fluid)
                    {
                        const boundaryConditionType UBCType =
                            field_broker_->URef()
                                .boundaryConditionRef(domain->index(),
                                                      iBoundary)
                                .type();

                        switch (UBCType)
                        {
                            case boundaryConditionType::noSlip:
                                {
                                    for (auto part : parts)
                                    {
                                        partVec.push_back(part);
                                    }
                                }
                                break;

                            case boundaryConditionType::slip:
                                {
                                    if (mesh.controlsRef().isTransient())
                                    {
                                        // only no-slip walls in case of
                                        // transient are relaxed
                                    }
                                    else
                                    {
                                        for (auto part : parts)
                                        {
                                            partVec.push_back(part);
                                        }
                                    }
                                }
                                break;

                            default:
                                errorMsg("invalid velocity boundary "
                                         "condition at wall");
                        }
                    }
                    else
                    {
                        for (auto part : parts)
                        {
                            partVec.push_back(part);
                        }
                    }
                }
                break;

            default:
                {
                    if (mesh.controlsRef().isTransient())
                    {
                        // only no-slip walls in case of transient are relaxed
                    }
                    else
                    {
                        for (auto part : parts)
                        {
                            partVec.push_back(part);
                        }
                    }
                }
                break;
        }
    }

#ifdef HAS_INTERFACE
    // add interface parts for relaxation under certain circumstances:
    // 1) if the interface is a fluid-solid interface
    // 2) if the interface whether inter-domain are connecting multiple
    //    domains, we only consider nodes nearest to exposed ip's
    for (const interface* interf : domain->interfacesRef())
    {
        if (interf->isFluidSolidType())
        {
            for (auto part :
                 interf->interfaceSideInfoPtr(domain->index())->currentPartVec_)
            {
                partVec.push_back(part);
            }
        }
    }
#endif /* HAS_INTERFACE */

    // Apply relaxation
    {
        stk::mesh::Selector selOwnedNodes =
            metaData.locally_owned_part() & stk::mesh::selectUnion(partVec);

        const BucketVec& nodeBuckets =
            bulkData.get_buckets(stk::topology::NODE_RANK, selOwnedNodes);

        // loop over local nodes and relax associated matrix rows
        for (size_t ib = 0; ib < nodeBuckets.size(); ib++)
        {
            const Bucket& bucket = *nodeBuckets[ib];
            const Bucket::size_type n_entities = bucket.size();
            for (Bucket::size_type i = 0; i < n_entities; ++i)
            {
                const stk::mesh::Entity entity = bucket[i];
                const auto lid = bulkData.local_id(entity);

                // relax diagonal
                assert(lid < b.size() / BLOCKSIZE);
                scalar* rhs_val = &b[BLOCKSIZE * lid];
                for (int k = 0; k < BLOCKSIZE; k++)
                {
                    rhs_val[k] *= urf;
                }
            }
        }
    }
}

void navierStokesAssembler::applySymmetryConditions_(const domain* domain,
                                                     Vector& b)
{
    // select all locally owned nodes for this domain
    const auto& mesh = model_->meshRef();
    const stk::mesh::MetaData& metaData = mesh.metaDataRef();
    const stk::mesh::BulkData& bulkData = mesh.bulkDataRef();

    const zone* zonePtr = domain->zonePtr();

    // Collect symmetry boundary parts
    stk::mesh::PartVector partVec;

    for (label iBoundary = 0; iBoundary < zonePtr->nBoundaries(); iBoundary++)
    {
        const auto& boundaryRef = zonePtr->boundaryRef(iBoundary);
        const stk::mesh::PartVector& parts = boundaryRef.parts();

        boundaryPhysicalType type = boundaryRef.type();
        switch (type)
        {
            case boundaryPhysicalType::symmetry:
                {
                    for (auto part : parts)
                    {
                        partVec.push_back(part);
                    }
                }
                break;

            default:
                break;
        }
    }

    // remove symmetry normal component from residual vector
    // select all locally owned nodes for this domain
    const auto& assembledSymmSTKFieldRef = *metaData.template get_field<scalar>(
        stk::topology::NODE_RANK, mesh::assembled_symm_area_ID);

    stk::mesh::Selector selOwnedNodes =
        metaData.locally_owned_part() & stk::mesh::selectUnion(partVec);
    const auto& sideNodeBuckets =
        bulkData.get_buckets(stk::topology::NODE_RANK, selOwnedNodes);

    for (const stk::mesh::Bucket* bucket : sideNodeBuckets)
    {
        const stk::mesh::Bucket& sideNodeBucket = *bucket;
        const auto nSideNodesPerBucket = sideNodeBucket.size();

        for (size_t iNode = 0; iNode < nSideNodesPerBucket; ++iNode)
        {
            stk::mesh::Entity node = sideNodeBucket[iNode];

            const auto lid = bulkData.local_id(node);

            scalar* rhs_val = &b[BLOCKSIZE * lid];

            const scalar* aarea =
                stk::mesh::field_data(assembledSymmSTKFieldRef, node);

            scalar asq = 0.0;
            for (label j = 0; j < SPATIAL_DIM; ++j)
            {
                const scalar axj = aarea[j];
                asq += axj * axj;
            }
            const scalar amag = std::sqrt(asq);

            scalar dot = 0.0;
            for (label j = 0; j < SPATIAL_DIM; ++j)
            {
                dot += rhs_val[j] * aarea[j] / amag;
            }

            for (int j = 0; j < SPATIAL_DIM; j++)
            {
                rhs_val[j] -= dot * aarea[j] / amag;
            }
        }
    }
}

void navierStokesAssembler::assembleNormalRelaxation(const domain* domain,
                                                     Context* ctx,
                                                     const scalar urf)
{
    using Bucket = stk::mesh::Bucket;
    using BucketVec = stk::mesh::BucketVector;

    // select all locally owned nodes for this domain
    const auto& mesh = model_->meshRef();
    const stk::mesh::MetaData& metaData = mesh.metaDataRef();
    const stk::mesh::BulkData& bulkData = mesh.bulkDataRef();

    const zone* zonePtr = domain->zonePtr();

    // no-slip walls and specified-velocity inlets
    {
        // Collect all boundary parts: It is ensured in this way that no
        // duplicate relaxation or more are made to corner nodes (belong to two
        // or more boundaries)
        stk::mesh::PartVector partVec;

        for (label iBoundary = 0; iBoundary < zonePtr->nBoundaries();
             iBoundary++)
        {
            const auto& boundaryRef = zonePtr->boundaryRef(iBoundary);
            const stk::mesh::PartVector& parts = boundaryRef.parts();

            boundaryConditionType UBCType =
                model_->URef()
                    .boundaryConditionRef(domain->index(), iBoundary)
                    .type();

            boundaryPhysicalType type = boundaryRef.type();
            switch (type)
            {
                case boundaryPhysicalType::inlet:
                    {
                        switch (UBCType)
                        {
                            case boundaryConditionType::specifiedValue:
                            case boundaryConditionType::normalSpeed:
                            case boundaryConditionType::massFlowRate:
                                {
                                    for (auto part : parts)
                                    {
                                        partVec.push_back(part);
                                    }
                                }
                                break;

                            default:
                                break;
                        }
                    }
                    break;

                case boundaryPhysicalType::wall:
                    {
                        switch (UBCType)
                        {
                            case boundaryConditionType::noSlip:
                                {
                                    for (auto part : parts)
                                    {
                                        partVec.push_back(part);
                                    }
                                }
                                break;

                            default:
                                break;
                        }
                    }
                    break;

                default:
                    break;
            }
        }

#ifdef HAS_INTERFACE
        // add interface parts for relaxation if fluid-solid interface
        for (const interface* interf : domain->interfacesRef())
        {
            if (interf->isFluidSolidType())
            {
                auto parts = interf->interfaceSideInfoPtr(domain->index())
                                 ->currentPartVec_;
                for (auto part : parts)
                {
                    partVec.push_back(part);
                }
            }
        }
#endif /* HAS_INTERFACE */

        // Apply relaxation
        {
            // not implemented yet
        }
    }

    // symmetry planes
    {
        // Collect all boundary parts: It is ensured in this way that no
        // duplicate relaxation or more are made to corner nodes (belong to two
        // or more boundaries)
        stk::mesh::PartVector partVec;

        for (label iBoundary = 0; iBoundary < zonePtr->nBoundaries();
             iBoundary++)
        {
            const auto& boundaryRef = zonePtr->boundaryRef(iBoundary);
            const stk::mesh::PartVector& parts = boundaryRef.parts();

            boundaryPhysicalType type = boundaryRef.type();
            switch (type)
            {
                case boundaryPhysicalType::symmetry:
                    {
                        for (auto part : parts)
                        {
                            partVec.push_back(part);
                        }
                    }
                    break;

                default:
                    break;
            }
        }

        // Apply relaxation
        {
            Matrix& A = ctx->getAMatrix();

            stk::mesh::Selector selOwnedNodes =
                metaData.locally_owned_part() & stk::mesh::selectUnion(partVec);

            const BucketVec& sideNodeBuckets =
                bulkData.get_buckets(stk::topology::NODE_RANK, selOwnedNodes);

            const auto& assembledSymmSTKFieldRef =
                *metaData.template get_field<scalar>(
                    stk::topology::NODE_RANK, mesh::assembled_symm_area_ID);

            std::vector<scalar> tran(SPATIAL_DIM * SPATIAL_DIM, 0.0);
            std::vector<scalar> modf(SPATIAL_DIM * SPATIAL_DIM, 0.0);

            // loop over local nodes and relax associated matrix rows
            for (size_t ib = 0; ib < sideNodeBuckets.size(); ib++)
            {
                const Bucket& sideNodeBucket = *sideNodeBuckets[ib];
                const Bucket::size_type nSideNodesPerBucket =
                    sideNodeBucket.size();
                for (Bucket::size_type iNode = 0; iNode < nSideNodesPerBucket;
                     ++iNode)
                {
                    const stk::mesh::Entity node = sideNodeBucket[iNode];

                    const auto lid = bulkData.local_id(node);

                    assert(lid < A.nRows());
                    scalar* diag = A.diag(lid);

                    const scalar* aarea =
                        stk::mesh::field_data(assembledSymmSTKFieldRef, node);

                    scalar amag = 0.0;
                    for (label i = 0; i < SPATIAL_DIM; i++)
                    {
                        amag += aarea[i] * aarea[i];
                    }
                    amag = std::sqrt(amag);

                    // zero out transformation matrix
                    for (label i = 0; i < SPATIAL_DIM; i++)
                    {
                        for (label j = 0; j < SPATIAL_DIM; j++)
                        {
                            tran[SPATIAL_DIM * i + j] = 0;
                            modf[SPATIAL_DIM * i + j] = 0;
                        }
                    }

                    // calculate transformation matrix = ninj
                    for (label i = 0; i < SPATIAL_DIM; i++)
                    {
                        scalar ni = aarea[i] / amag;
                        for (label j = 0; j < SPATIAL_DIM; j++)
                        {
                            scalar nj = aarea[j] / amag;
                            tran[SPATIAL_DIM * i + j] += ni * nj;
                        }
                    }

                    // multiply matrices
                    for (label i = 0; i < SPATIAL_DIM; i++)
                    {
                        for (label j = 0; j < SPATIAL_DIM; j++)
                        {
                            for (label k = 0; k < SPATIAL_DIM; ++k)
                            {
                                modf[SPATIAL_DIM * i + j] +=
                                    diag[SPATIAL_DIM * i + k] *
                                    tran[SPATIAL_DIM * k + j];
                            }
                        }
                    }

                    // modify momentum diagonal coefficient
                    for (label i = 0; i < SPATIAL_DIM; i++)
                    {
                        for (label j = 0; j < SPATIAL_DIM; j++)
                        {
                            diag[SPATIAL_DIM * i + j] +=
                                (1.0 - urf) * modf[SPATIAL_DIM * i + j];
                        }
                    }
                }
            }
        }
    }
}

} // namespace accel
