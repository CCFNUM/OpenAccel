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
#ifdef USE_CVFEM_SOLID_MECHANICS
    phiAssembler<SPATIAL_DIM>::postAssemble(domain, ctx);
    applySymmetryConditions_(domain, ctx);
#else
    errorMsg("FEM solid mechanics not implemented yet");
#endif
}

void solidDisplacementAssembler::applySymmetryConditions_(const domain* domain,
                                                          Context* ctx)
{
#ifdef USE_CVFEM_SOLID_MECHANICS
    const auto& mesh = model_->meshRef();
    const stk::mesh::MetaData& metaData = mesh.metaDataRef();
    const stk::mesh::BulkData& bulkData = mesh.bulkDataRef();
    const zone* zonePtr = domain->zonePtr();

    const auto& assembledSymmSTKFieldRef = *metaData.template get_field<scalar>(
        stk::topology::NODE_RANK, mesh::assembled_symm_area_ID);

    Matrix& A = ctx->getAMatrix();
    Vector& b = ctx->getBVector();

    stk::mesh::PartVector partVec;
    for (label iBoundary = 0; iBoundary < zonePtr->nBoundaries(); iBoundary++)
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
            const label nBlocks =
                static_cast<label>(vals.size()) / (SPATIAL_DIM * SPATIAL_DIM);
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
#else
    errorMsg("FEM solid mechanics not implemented yet");
#endif
}

} /* namespace accel */
