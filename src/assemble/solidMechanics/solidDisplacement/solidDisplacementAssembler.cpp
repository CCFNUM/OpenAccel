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

    // Collect all symmetry boundary parts
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
            partVec.push_back(part);
    }

    if (partVec.empty())
        return;

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

            // Compute unit normal from the assembled symmetry area vector
            const scalar* aarea =
                stk::mesh::field_data(assembledSymmSTKFieldRef, node);

            scalar asq = 0.0;
            for (label j = 0; j < SPATIAL_DIM; ++j)
                asq += aarea[j] * aarea[j];
            const scalar amag = std::sqrt(asq);

            if (amag < 1.0e-30)
                continue;

            scalar n_hat[SPATIAL_DIM];
            for (label j = 0; j < SPATIAL_DIM; ++j)
                n_hat[j] = aarea[j] / amag;

            // Choose the row to replace: the Cartesian DOF with the largest
            // normal component. This maximises diagonal dominance.
            label constrained_dof = 0;
            scalar maxAbs = std::abs(n_hat[0]);
            for (label j = 1; j < SPATIAL_DIM; ++j)
            {
                if (std::abs(n_hat[j]) > maxAbs)
                {
                    maxAbs = std::abs(n_hat[j]);
                    constrained_dof = j;
                }
            }

            // Zero row constrained_dof across all blocks (rowVals includes
            // the diagonal block), then write n_hat into the diagonal entry.
            // The diagonal write must come after the zero sweep because
            // rowVals and diag alias the same storage.
            auto rowVals = A.rowVals(lid);
            auto rowCols = A.rowCols(lid);
            for (label icol = 0; icol < static_cast<label>(rowCols.size());
                 ++icol)
            {
                for (label j = 0; j < BLOCKSIZE; ++j)
                    rowVals[BLOCKSIZE * BLOCKSIZE * icol +
                            BLOCKSIZE * constrained_dof + j] = 0.0;
            }

            // Replace row constrained_dof of the diagonal block with n_hat.
            // For axis-aligned normals this reduces to the standard Dirichlet
            // row (1 on the diagonal, 0 elsewhere). For oblique normals it
            // enforces the general constraint n_hat · u = 0.
            auto* diag = A.diag(lid);
            for (label j = 0; j < SPATIAL_DIM; ++j)
                diag[BLOCKSIZE * constrained_dof + j] = n_hat[j];

            // Zero the RHS entry for the constrained DOF
            b[BLOCKSIZE * lid + constrained_dof] = 0.0;
        }
    }
#else
    errorMsg("FEM solid mechanics not implemented yet");
#endif
}

} /* namespace accel */
