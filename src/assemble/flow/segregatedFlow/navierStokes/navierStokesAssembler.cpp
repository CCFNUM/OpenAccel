// File       : navierStokesAssembler.cpp
// Created    : Wed Jan 03 2024 13:38:51 (+0100)
// Author     : Mhamad Mahdi Alloush
// Description:
// Copyright 2024 CCFNUM HSLU T&A. All Rights Reserved.

#include "navierStokesAssembler.h"
#include "flowModel.h"
#include "interface.h"

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

            const int64_t localID =
                A.getGraph()->localToRow(bulkData.local_id(node));
            if (localID < 0) // node not in this (subset) graph
                continue;

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

    // conformal interface: give master and slave the effective (vol1+vol2)
    // volume in their D coefficients, mirroring the coupled assembler
    for (const interface* interf : domain->interfacesRef())
    {
        if (!interf->isConformalTreatment() ||
            !interf->isMasterZone(domain->index()))
            continue;

        std::vector<std::pair<std::string, std::string>> fields = {
            {flowModel::du_ID, flowModel::du_ID}};
        if (consistent)
        {
            fields.push_back({flowModel::duTilde_ID, flowModel::duTilde_ID});
        }

        conformalDataTransfer transfer(const_cast<interface*>(interf),
                                       std::string(flowModel::du_ID) +
                                           "_xfer_" + interf->name(),
                                       fields,
                                       dataTransferType::customized);
        transfer.setup();
        transfer.initialize();
        transfer.update();
    }
}

void navierStokesAssembler::postAssemble(const domain* domain, Context* ctx)
{
    phiAssembler::postAssemble(domain, ctx);
    computeDUCoefficients(domain, ctx);
    assembleBoundaryRelaxation_(domain, ctx, 0.75);
    applySymmetryConditions_(domain, ctx);
}

void navierStokesAssembler::assembleBoundaryRelaxation_(const domain* domain,
                                                        Context* ctx,
                                                        const scalar urf)
{
    using Bucket = stk::mesh::Bucket;
    using BucketVec = stk::mesh::BucketVector;

    Vector& b = ctx->getBVector();
    const auto* graph = ctx->getAMatrix().getGraph();

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
                const int64_t row =
                    graph->localToRow(bulkData.local_id(entity));
                if (row < 0) // node not part of this (subset) system
                    continue;

                // relax diagonal
                assert(row < static_cast<int64_t>(b.size() / BLOCKSIZE));
                scalar* rhs_val = &b[BLOCKSIZE * row];
                for (int k = 0; k < BLOCKSIZE; k++)
                {
                    rhs_val[k] *= urf;
                }
            }
        }
    }
}

void navierStokesAssembler::applySymmetryConditions_(const domain* domain,
                                                     Context* ctx)
{
    Matrix& A = ctx->getAMatrix();
    Vector& b = ctx->getBVector();

    const auto& mesh = model_->meshRef();
    const stk::mesh::MetaData& metaData = mesh.metaDataRef();
    const stk::mesh::BulkData& bulkData = mesh.bulkDataRef();
    const zone* zonePtr = domain->zonePtr();

    // collect symmetry boundary parts
    stk::mesh::PartVector partVec;
    for (label iB = 0; iB < zonePtr->nBoundaries(); ++iB)
        if (zonePtr->boundaryRef(iB).type() == boundaryPhysicalType::symmetry)
            for (auto* pp : zonePtr->boundaryRef(iB).parts())
                partVec.push_back(pp);
    if (partVec.empty())
        return;

    const auto& symmArea = *metaData.template get_field<scalar>(
        stk::topology::NODE_RANK, mesh::assembled_symm_area_ID);

    // check if any interface slave nodes shared by any symmetry boundary.
    // Construct a set of these salve nodes to avoid applying any relaxation or
    // normal-residual cancellation
    std::unordered_set<std::size_t> slaveIfcNodes;
    for (const interface* interf : domain->interfacesRef())
    {
        if (!interf->isConformalTreatment())
            continue;

        for (const auto& np : interf->matchingNodePairVector())
            slaveIfcNodes.insert(np.second.local_offset());
    }

    // fixed-size containers
    scalar n[SPATIAL_DIM];

    // we require diagonal offsets
    const auto& diagOffsets = A.diagOffsetRef();

    stk::mesh::Selector sel =
        metaData.locally_owned_part() & stk::mesh::selectUnion(partVec);
    for (const stk::mesh::Bucket* bkt :
         bulkData.get_buckets(stk::topology::NODE_RANK, sel))
    {
        for (size_t iN = 0; iN < bkt->size(); ++iN)
        {
            const stk::mesh::Entity node = (*bkt)[iN];
            if (slaveIfcNodes.count(node.local_offset()))
                continue;

            const int64_t row = ctx->getAMatrix().getGraph()->localToRow(
                bulkData.local_id(node));
            if (row < 0) // node not part of this (subset) system
                continue;

            // unit symmetry normal
            const scalar* aarea = stk::mesh::field_data(symmArea, node);
            scalar asq = 0.0;
            for (label i = 0; i < SPATIAL_DIM; ++i)
            {
                n[i] = aarea[i];
                asq += n[i] * n[i];
            }
            if (asq < SMALL)
                continue;
            const scalar amag = std::sqrt(asq);
            for (label i = 0; i < SPATIAL_DIM; ++i)
                n[i] /= amag;

            auto vals = A.rowVals(row);
            const label nBlocks =
                static_cast<label>(vals.size()) / (SPATIAL_DIM * SPATIAL_DIM);
            const label diagBk = diagOffsets[row];

            // normal stiffness scale = n^T D n (keeps du physical after decoup)
            scalar* Dblk = &vals[SPATIAL_DIM * SPATIAL_DIM * diagBk];
            scalar scale = 0.0;
            for (label i = 0; i < SPATIAL_DIM; ++i)
                for (label j = 0; j < SPATIAL_DIM; ++j)
                    scale += n[i] * Dblk[i * SPATIAL_DIM + j] * n[j];
            if (std::abs(scale) < SMALL)
                scale = amag;

            // For every block B in the row: B <- T B (zero its normal row), and
            // for the diagonal block restore the normal coefficient (+ scale
            // N). Result: the normal equation is scale * (n.U) = 0 with no
            // tangential or neighbour coupling, i.e. U.n = 0 fully decoupled.
            for (label bk = 0; bk < nBlocks; ++bk)
            {
                scalar* B = &vals[SPATIAL_DIM * SPATIAL_DIM * bk];
                scalar colProj[SPATIAL_DIM];
                for (label j = 0; j < SPATIAL_DIM; ++j)
                {
                    scalar s = 0.0;
                    for (label k = 0; k < SPATIAL_DIM; ++k)
                        s += n[k] * B[k * SPATIAL_DIM + j];
                    colProj[j] = s;
                }
                for (label i = 0; i < SPATIAL_DIM; ++i)
                    for (label j = 0; j < SPATIAL_DIM; ++j)
                        B[i * SPATIAL_DIM + j] -= n[i] * colProj[j];
                if (bk == diagBk)
                    for (label i = 0; i < SPATIAL_DIM; ++i)
                        for (label j = 0; j < SPATIAL_DIM; ++j)
                            B[i * SPATIAL_DIM + j] += scale * n[i] * n[j];
            }

            // RHS: remove the normal component (normal equation rhs = 0)
            scalar* rhs = &b[BLOCKSIZE * row];
            scalar bproj = 0.0;
            for (label i = 0; i < SPATIAL_DIM; ++i)
                bproj += n[i] * rhs[i];
            for (label i = 0; i < SPATIAL_DIM; ++i)
                rhs[i] -= n[i] * bproj;
        }
    }
}

} // namespace accel
