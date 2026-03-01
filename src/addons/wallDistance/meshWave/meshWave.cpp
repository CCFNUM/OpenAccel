// File       : meshWave.cpp
// Created    : Sat Mar 01 2026
// Author     : Mhamad Mahdi Alloush
// Description: Wall distance computation via direct Euclidean distance to
// nearest wall node
// Copyright (c) 2026 CCFNUM, Lucerne University of Applied Sciences and Arts.
// SPDX-License-Identifier: BSD-3-Clause

// code
#include "meshWave.h"
#include "messager.h"
#include "realm.h"
#include "simulation.h"
#include "wallDistance.h"

namespace accel
{

meshWave::meshWave(wallDistance* wallDistancePtr)
    : wallDistancePtr_(wallDistancePtr)
{
}

void meshWave::setup()
{
    // no equation to set up
}

void meshWave::initialize()
{
    computeWallDistance_();
}

void meshWave::update()
{
    // only recompute if mesh is deforming
    for (auto domain :
         wallDistancePtr_->realmRef().simulationRef().domainVector())
    {
        if (domain->isWallDistanceRequired() &&
            domain->zonePtr()->deformationRef().specification() !=
                meshDeformationSpecificationType::none)
        {
            computeWallDistance_();
            return;
        }
    }
}

void meshWave::computeWallDistance_()
{
    auto& mesh = wallDistancePtr_->meshRef();
    stk::mesh::BulkData& bulkData = mesh.bulkDataRef();
    stk::mesh::MetaData& metaData = mesh.metaDataRef();

    // get coordinate field
    STKScalarField* coordinates = metaData.get_field<scalar>(
        stk::topology::NODE_RANK, mesh::coordinates_ID);

    // collect locally-owned wall node coordinates (owned only to avoid
    // duplicates when gathering across MPI ranks)
    stk::mesh::Selector selWallNodes =
        metaData.locally_owned_part() &
        stk::mesh::selectUnion(mesh.wallBoundaryActiveParts());

    const label nLocalWallNodes = stk::mesh::count_entities(
        bulkData, stk::topology::NODE_RANK, selWallNodes);

    std::vector<scalar> localWallCoords(nLocalWallNodes * SPATIAL_DIM);
    label idx = 0;

    stk::mesh::BucketVector const& wallBuckets =
        bulkData.get_buckets(stk::topology::NODE_RANK, selWallNodes);

    for (auto ib = wallBuckets.begin(); ib != wallBuckets.end(); ++ib)
    {
        stk::mesh::Bucket& b = **ib;
        const stk::mesh::Bucket::size_type length = b.size();

        for (stk::mesh::Bucket::size_type k = 0; k < length; ++k)
        {
            stk::mesh::Entity node = b[k];
            const scalar* coords = stk::mesh::field_data(*coordinates, node);
            for (label d = 0; d < SPATIAL_DIM; ++d)
            {
                localWallCoords[idx++] = coords[d];
            }
        }
    }

    // gather all wall node coordinates across MPI ranks
    // gatherVector requires the output vector to be pre-sized
    int localSize = localWallCoords.size();
    int globalSize = localSize;
    if (messager::parallel())
    {
        MPI_Allreduce(
            &localSize, &globalSize, 1, MPI_INT, MPI_SUM, messager::comm());
    }
    std::vector<scalar> globalWallCoords(globalSize);
    messager::gatherVector(localWallCoords, globalWallCoords);

    const label nGlobalWallNodes = globalSize / SPATIAL_DIM;

    // compute minimum distance for each interior node
    STKScalarField* yMinField = wallDistancePtr_->yMinRef().stkFieldPtr();

    for (auto domain :
         wallDistancePtr_->realmRef().simulationRef().domainVector())
    {
        if (domain->isWallDistanceRequired())
        {
            stk::mesh::Selector selInterior =
                metaData.locally_owned_part() &
                stk::mesh::selectUnion(domain->zonePtr()->interiorParts());

            stk::mesh::BucketVector const& nodeBuckets =
                bulkData.get_buckets(stk::topology::NODE_RANK, selInterior);

            for (auto ib = nodeBuckets.begin(); ib != nodeBuckets.end(); ++ib)
            {
                stk::mesh::Bucket& b = **ib;
                const stk::mesh::Bucket::size_type length = b.size();

                for (stk::mesh::Bucket::size_type k = 0; k < length; ++k)
                {
                    stk::mesh::Entity node = b[k];
                    const scalar* nodeCoords =
                        stk::mesh::field_data(*coordinates, node);
                    scalar* yMin = stk::mesh::field_data(*yMinField, node);

                    scalar minDist = std::numeric_limits<scalar>::max();

                    for (label w = 0; w < nGlobalWallNodes; ++w)
                    {
                        scalar distSq = 0.0;
                        for (label d = 0; d < SPATIAL_DIM; ++d)
                        {
                            const scalar diff =
                                nodeCoords[d] -
                                globalWallCoords[w * SPATIAL_DIM + d];
                            distSq += diff * diff;
                        }
                        minDist = std::min(minDist, distSq);
                    }

                    *yMin = std::sqrt(minDist);
                }
            }
        }
    }

    // synchronize for ghosted/shared nodes
    if (messager::parallel())
    {
        stk::mesh::communicate_field_data(bulkData, {yMinField});
    }
}

} // namespace accel
