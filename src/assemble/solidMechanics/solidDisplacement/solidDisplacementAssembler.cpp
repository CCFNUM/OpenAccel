// File       : solidDisplacementAssembler.cpp
// Created    : Sun Feb 01 2026 02:30:10 (+0100)
// Author     : Mhamad Mahdi Alloush
// Description:
// Copyright 2026 CCFNUM HSLU T&A. All Rights Reserved.

#include "solidDisplacementAssembler.h"

namespace accel
{

void solidDisplacementAssembler::postAssemble(const domain* domain,
                                              Context* ctx)
{
    if (field_broker_->controlsRef().isCvfemSolidMechanics())
    {
        phiAssembler<SPATIAL_DIM>::postAssemble(domain, ctx);
    }
    else
    {
        // SFEM tangent is consistent; diagonal relaxation would break it
        this->applyConstraints(domain, ctx);
    }
    applySymmetryConditions_(domain, ctx);
}

void solidDisplacementAssembler::applySymmetryConditions_(const domain* domain,
                                                          Context* ctx)
{
    if (field_broker_->controlsRef().isCvfemSolidMechanics())
    {
        const auto& mesh = model_->meshRef();
        const stk::mesh::MetaData& metaData = mesh.metaDataRef();
        const stk::mesh::BulkData& bulkData = mesh.bulkDataRef();
        const zone* zonePtr = domain->zonePtr();

        const auto& assembledSymmSTKFieldRef =
            *metaData.template get_field<scalar>(stk::topology::NODE_RANK,
                                                 mesh::assembled_symm_area_ID);

        Matrix& A = ctx->getAMatrix();
        Vector& b = ctx->getBVector();

        stk::mesh::PartVector partVec;
        for (label iBoundary = 0; iBoundary < zonePtr->nBoundaries();
             iBoundary++)
        {
            const auto& bcType =
                model_->DRef()
                    .boundaryConditionRef(domain->index(), iBoundary)
                    .type();

            if (bcType != boundaryConditionType::symmetry)
                continue;

            for (auto* part : zonePtr->boundaryRef(iBoundary).parts())
            {
                partVec.push_back(part);
            }
        }

        if (partVec.empty())
            return;

        // fixed-size containers
        scalar n[SPATIAL_DIM];

        // we require diagonal offsets
        const auto& diagOffsets = A.diagOffsetRef();

        stk::mesh::Selector selOwnedNodes =
            metaData.locally_owned_part() & stk::mesh::selectUnion(partVec);
        const auto& nodeBuckets =
            bulkData.get_buckets(stk::topology::NODE_RANK, selOwnedNodes);

        for (const stk::mesh::Bucket* bucket : nodeBuckets)
        {
            for (size_t iNode = 0; iNode < bucket->size(); ++iNode)
            {
                stk::mesh::Entity node = (*bucket)[iNode];
                const int64_t lid =
                    A.getGraph()->localToRow(bulkData.local_id(node));
                if (lid < 0) // node not part of this (subset) system
                    continue;

                const scalar* aarea =
                    stk::mesh::field_data(assembledSymmSTKFieldRef, node);
                scalar amagSq = 0.0;
                for (label j = 0; j < SPATIAL_DIM; ++j)
                    amagSq += aarea[j] * aarea[j];

                if (amagSq < 1.0e-30)
                    continue;
                const scalar amag = std::sqrt(amagSq);

                // unit symmetry normal
                for (label j = 0; j < SPATIAL_DIM; ++j)
                    n[j] = aarea[j] / amag;

                auto vals = A.rowVals(lid);
                const label nBlocks = static_cast<label>(vals.size()) /
                                      (SPATIAL_DIM * SPATIAL_DIM);
                const label diagBk = diagOffsets[lid];

                // normal stiffness scale = n^T D n
                scalar* Dblk = &vals[SPATIAL_DIM * SPATIAL_DIM * diagBk];
                scalar scale = 0.0;
                for (label i = 0; i < SPATIAL_DIM; ++i)
                    for (label j = 0; j < SPATIAL_DIM; ++j)
                        scale += n[i] * Dblk[i * SPATIAL_DIM + j] * n[j];
                if (std::abs(scale) < SMALL)
                    scale = amag;

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

                // RHS: remove the normal component
                scalar bproj = 0.0;
                for (label i = 0; i < SPATIAL_DIM; ++i)
                    bproj += n[i] * b[BLOCKSIZE * lid + i];
                for (label i = 0; i < SPATIAL_DIM; ++i)
                    b[BLOCKSIZE * lid + i] -= n[i] * bproj;
            }
        }
    }
    else
    {
        const auto& mesh = model_->meshRef();
        const stk::mesh::MetaData& metaData = mesh.metaDataRef();
        const stk::mesh::BulkData& bulkData = mesh.bulkDataRef();
        const zone* zonePtr = domain->zonePtr();

        Matrix& A = ctx->getAMatrix();
        Vector& b = ctx->getBVector();

        stk::mesh::PartVector partVec;
        for (label iBoundary = 0; iBoundary < zonePtr->nBoundaries();
             iBoundary++)
        {
            const auto& bcType =
                model_->DRef()
                    .boundaryConditionRef(domain->index(), iBoundary)
                    .type();

            if (bcType != boundaryConditionType::symmetry)
                continue;

            for (auto* part : zonePtr->boundaryRef(iBoundary).parts())
            {
                partVec.push_back(part);
            }
        }

        if (partVec.empty())
            return;

        const auto& exposedAreaVectorField = *metaData.get_field<scalar>(
            metaData.side_rank(), this->getExposedAreaVectorID_(domain));
        const stk::mesh::Selector selectedSides =
            metaData.universal_part() & stk::mesh::selectUnion(partVec);
        const auto& sideBuckets =
            bulkData.get_buckets(metaData.side_rank(), selectedSides);

        std::unordered_map<label, std::array<scalar, SPATIAL_DIM>>
            assembledArea;

        for (const stk::mesh::Bucket* bucket : sideBuckets)
        {
            MasterElement* faceMasterElement =
                MasterElementRepo::get_surface_master_element(
                    bucket->topology());
            const label integrationPoints = faceMasterElement->numIntPoints_;

            for (stk::mesh::Entity side : *bucket)
            {
                const scalar* areaVector =
                    stk::mesh::field_data(exposedAreaVectorField, side);
                std::array<scalar, SPATIAL_DIM> area{};
                for (label ip = 0; ip < integrationPoints; ++ip)
                {
                    for (label dim = 0; dim < SPATIAL_DIM; ++dim)
                        area[dim] += areaVector[ip * SPATIAL_DIM + dim];
                }

                const stk::mesh::Entity* sideNodes = bulkData.begin_nodes(side);
                const unsigned nodesPerSide = bulkData.num_nodes(side);
                for (unsigned node = 0; node < nodesPerSide; ++node)
                {
                    const label lid = bulkData.local_id(sideNodes[node]);
                    auto& nodalArea = assembledArea[lid];
                    for (label dim = 0; dim < SPATIAL_DIM; ++dim)
                        nodalArea[dim] += area[dim];
                }
            }
        }

        stk::mesh::Selector selOwnedNodes =
            metaData.locally_owned_part() & stk::mesh::selectUnion(partVec);
        const auto& nodeBuckets =
            bulkData.get_buckets(stk::topology::NODE_RANK, selOwnedNodes);

        for (const stk::mesh::Bucket* bucket : nodeBuckets)
        {
            for (size_t iNode = 0; iNode < bucket->size(); ++iNode)
            {
                stk::mesh::Entity node = (*bucket)[iNode];
                const auto lid = bulkData.local_id(node);
                const auto assembled = assembledArea.find(lid);
                if (assembled == assembledArea.end())
                    continue;

                scalar amagSq = 0.0;
                for (label dim = 0; dim < SPATIAL_DIM; ++dim)
                    amagSq += assembled->second[dim] * assembled->second[dim];
                if (amagSq < 1.0e-30)
                    continue;

                std::array<scalar, SPATIAL_DIM> normal{};
                const scalar amag = std::sqrt(amagSq);
                for (label dim = 0; dim < SPATIAL_DIM; ++dim)
                    normal[dim] = assembled->second[dim] / amag;

                scalar* const diag = A.diag(lid);
                scalar scale = 0.0;
                for (label j = 0; j < SPATIAL_DIM; ++j)
                    scale = std::max(scale, std::abs(diag[BLOCKSIZE * j + j]));
                if (scale == 0.0)
                    scale = 1.0;

                scalar bDotN = 0.0;
                for (label j = 0; j < SPATIAL_DIM; ++j)
                    bDotN += b[BLOCKSIZE * lid + j] * normal[j];

                for (label j = 0; j < SPATIAL_DIM; ++j)
                    b[BLOCKSIZE * lid + j] -= bDotN * normal[j];

                for (label i = 0; i < SPATIAL_DIM; ++i)
                {
                    for (label j = 0; j < SPATIAL_DIM; ++j)
                    {
                        diag[BLOCKSIZE * i + j] +=
                            scale * normal[i] * normal[j];
                    }
                }
            }
        }
    }
}

} /* namespace accel */
