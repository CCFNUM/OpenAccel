// File : solidDisplacementAssembler.cpp
// Created : Sun Feb 01 2026 02:30:10 (+0100)
// Author : Mhamad Mahdi Alloush
// Description:
// Copyright (c) 2026 CCFNUM, Lucerne University of Applied Sciences and Arts.
// SPDX-License-Identifier: BSD-3-Clause

#include "solidDisplacementAssembler.h"

namespace accel
{

void solidDisplacementAssembler::postAssemble_(const domain* domain,
                                               Context* ctx)
{
#ifdef USE_CVFEM_SOLID_MECHANICS
    phiAssembler<SPATIAL_DIM>::postAssemble_(domain, ctx);
    applySymmetryConditions_(domain, ctx->getBVector());
#else
    errorMsg("FEM solid mechanics not implemented yet");
#endif
}

void solidDisplacementAssembler::applySymmetryConditions_(const domain* domain,
                                                          Vector& b)
{
#ifdef USE_CVFEM_SOLID_MECHANICS
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
#else
    errorMsg("FEM solid mechanics not implemented yet");
#endif
}

} /* namespace accel */
