// File       : displacementDiffusionAssembler.cpp
// Created    : Sun Feb 01 2026
// Author     : Mhamad Mahdi Alloush
// Description:
// Copyright 2026 CCFNUM HSLU T&A. All Rights Reserved.

#include "displacementDiffusionAssembler.h"

namespace accel
{

void displacementDiffusionAssembler::postAssemble(const domain* domain,
                                                  Context* ctx)
{
    phiAssembler<SPATIAL_DIM>::postAssemble(domain, ctx);
    applySymmetryConditions_(domain, ctx);
}

void displacementDiffusionAssembler::applySymmetryConditions_(
    const domain* domain,
    Context* ctx)
{
    // get system matrix and rhs vector
    Matrix& A = ctx->getAMatrix();
    Vector& b = ctx->getBVector();

    // select all locally owned nodes for this domain
    const auto& mesh = field_broker_->meshRef();
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

    if (partVec.empty())
        return;

    // remove symmetry normal component from residual vector
    // select all locally owned nodes for this domain
    const auto& assembledSymmSTKFieldRef = *metaData.template get_field<scalar>(
        stk::topology::NODE_RANK, mesh::assembled_symm_area_ID);

    // fixed-size containers
    scalar n[SPATIAL_DIM];

    // we require diagonal offsets
    const auto& diagOffsets = A.diagOffsetRef();

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

            const int64_t row = ctx->getAMatrix().getGraph()->localToRow(
                bulkData.local_id(node));
            if (row < 0) // node not part of this (subset) system
                continue;

            // unit symmetry normal
            const scalar* aarea =
                stk::mesh::field_data(assembledSymmSTKFieldRef, node);
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

} /* namespace accel */
