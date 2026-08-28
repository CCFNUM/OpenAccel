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
        // FEM solid mechanics assembles a consistent Newton tangent. Relaxing
        // only its diagonal destroys rigid-body consistency and, together
        // with the relaxed field correction, applies the displacement URF
        // twice.
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
                const int64_t lid =
                    A.getGraph()->localToRow(bulkData.local_id(node));
                if (lid < 0) // node not part of this (subset) system
                    continue;

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

                // 2. Identify the characteristic scale (stiffness) for this
                // node
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
                // 'scale'. This is equivalent to u_n = 0 without biasing a
                // specific Cartesian axis.
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
