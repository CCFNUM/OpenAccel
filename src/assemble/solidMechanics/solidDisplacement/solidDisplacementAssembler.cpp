// File       : solidDisplacementAssembler.cpp
// Created    : Sun Feb 01 2026 02:30:10 (+0100)
// Author     : Mhamad Mahdi Alloush
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
    std::vector<scalar> n(SPATIAL_DIM);

    // pointers ..
    scalar* p_n = &n[0];

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

            const scalar* aarea =
                stk::mesh::field_data(assembledSymmSTKFieldRef, node);
            scalar amagSq = 0.0;
            for (label j = 0; j < SPATIAL_DIM; ++j)
            {
                amagSq += aarea[j] * aarea[j];
            }

            if (amagSq < 1.0e-30)
                continue;
            const scalar amag = std::sqrt(amagSq);

            // 1. Compute Unit Normal n
            scalar n[SPATIAL_DIM];
            for (label j = 0; j < SPATIAL_DIM; ++j)
            {
                p_n[j] = aarea[j] / amag;
            }

            // 2. Identify the characteristic scale (stiffness) for this node
            scalar* const diag = A.diag(lid);
            scalar scale = 0.0;
            for (label j = 0; j < SPATIAL_DIM; ++j)
            {
                scale = std::max(scale, std::abs(diag[BLOCKSIZE * j + j]));
            }

            // 3. Project the RHS: b = b - (b \cdot n)n
            // This removes any "force" component pushing out of the plane.
            scalar b_dot_n = 0.0;
            for (label j = 0; j < SPATIAL_DIM; ++j)
            {
                b_dot_n += b[BLOCKSIZE * lid + j] * p_n[j];
            }

            for (label j = 0; j < SPATIAL_DIM; ++j)
            {
                b[BLOCKSIZE * lid + j] -= b_dot_n * p_n[j];
            }

            // 4. Modify the Diagonal Block to enforce u \cdot n = 0
            // We use a "Penalty-lite" approach:
            // We strip the diagonal's normal stiffness and replace it with
            // 'scale'. This is equivalent to u_n = 0 without biasing a specific
            // Cartesian axis.
            for (label i = 0; i < SPATIAL_DIM; ++i)
            {
                for (label j = 0; j < SPATIAL_DIM; ++j)
                {
                    // This creates a contribution proportional to n_i * n_j
                    // effectively adding stiffness only in the normal
                    // direction.
                    diag[BLOCKSIZE * i + j] += scale * p_n[i] * p_n[j];
                }
            }
        }
    }
#else
    errorMsg("FEM solid mechanics not implemented yet");
#endif
}

} /* namespace accel */
