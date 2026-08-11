// File       : signedDistanceFunction.cpp
// Description: Wall distance computation using the SSDF library

// code
#include "signedDistanceFunction.h"

#include "messager.h"
#include "realm.h"
#include "simulation.h"
#include "wallDistance.h"

// paths are those exported by the in-tree ssdf target, not by an installed SSDF
#include <ssdf.hpp>
#ifdef SSDF_ENABLE_MPI
#include <distributed/psdf.hpp>
#endif

#include <cstddef>
#include <limits>
#include <unordered_map>
#include <unordered_set>
#include <vector>

namespace accel
{

signedDistanceFunction::signedDistanceFunction(wallDistance* wallDistancePtr)
    : wallDistancePtr_(wallDistancePtr)
{
}

void signedDistanceFunction::setup()
{
    // SSDF is configured and linked by CMake; no equation setup is required.
}

void signedDistanceFunction::initialize()
{
    computeWallDistance_();
}

void signedDistanceFunction::update()
{
    // Recompute only when a contributing mesh zone is deforming.
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

void signedDistanceFunction::computeWallDistance_()
{
    auto& mesh = wallDistancePtr_->meshRef();
    stk::mesh::BulkData& bulkData = mesh.bulkDataRef();
    stk::mesh::MetaData& metaData = mesh.metaDataRef();

    STKScalarField* coordinates = metaData.get_field<scalar>(
        stk::topology::NODE_RANK, mesh::coordinates_ID);
    STKScalarField* yMinField = wallDistancePtr_->yMinRef().stkFieldPtr();

    if (!coordinates || !yMinField)
    {
        errorMsg("SSDF wall distance requires coordinates and yMin fields");
    }

    // SSDF uses structure-of-arrays storage for query points.
    std::vector<scalar> x;
    std::vector<scalar> y;
    std::vector<scalar> z;
    std::vector<stk::mesh::Entity> queryNodes;
    std::unordered_set<stk::mesh::EntityId> insertedQueryNodes;

    for (auto domain :
         wallDistancePtr_->realmRef().simulationRef().domainVector())
    {
        if (!domain->isWallDistanceRequired())
        {
            continue;
        }

        const stk::mesh::Selector selInterior =
            metaData.locally_owned_part() &
            stk::mesh::selectUnion(domain->zonePtr()->interiorParts());

        const stk::mesh::BucketVector& nodeBuckets =
            bulkData.get_buckets(stk::topology::NODE_RANK, selInterior);

        for (const stk::mesh::Bucket* bucketPtr : nodeBuckets)
        {
            const stk::mesh::Bucket& bucket = *bucketPtr;
            for (stk::mesh::Bucket::size_type k = 0; k < bucket.size(); ++k)
            {
                const stk::mesh::Entity node = bucket[k];
                const stk::mesh::EntityId id = bulkData.identifier(node);
                if (!insertedQueryNodes.insert(id).second)
                {
                    continue;
                }

                const scalar* xyz = stk::mesh::field_data(*coordinates, node);
                queryNodes.push_back(node);
                x.push_back(xyz[0]);
                y.push_back(xyz[1]);
#if SPATIAL_DIM == 3
                z.push_back(xyz[2]);
#else
                z.push_back(0.0);
#endif
            }
        }
    }

    // Extract the locally owned wall surface. Each rank supplies its own
    // portion to pedf; pedf handles communication between ranks.
    std::vector<scalar> sx;
    std::vector<scalar> sy;
    std::vector<scalar> sz;
    std::vector<label> s0;
    std::vector<label> s1;
    std::vector<label> s2;
    std::unordered_map<stk::mesh::EntityId, label> surfaceNodeToIndex;

    const stk::mesh::Selector selWallSides =
        metaData.locally_owned_part() &
        stk::mesh::selectUnion(mesh.wallBoundaryActiveParts());

    const stk::mesh::BucketVector& sideBuckets =
        bulkData.get_buckets(metaData.side_rank(), selWallSides);

    auto surfaceIndex = [&](const stk::mesh::Entity node) -> label
    {
        const stk::mesh::EntityId id = bulkData.identifier(node);
        const auto found = surfaceNodeToIndex.find(id);
        if (found != surfaceNodeToIndex.end())
        {
            return found->second;
        }

        const label index = static_cast<label>(sx.size());
        surfaceNodeToIndex.emplace(id, index);

        const scalar* xyz = stk::mesh::field_data(*coordinates, node);
        sx.push_back(xyz[0]);
        sy.push_back(xyz[1]);
#if SPATIAL_DIM == 3
        sz.push_back(xyz[2]);
#else
        sz.push_back(0.0);
#endif
        return index;
    };

    for (const stk::mesh::Bucket* bucketPtr : sideBuckets)
    {
        const stk::mesh::Bucket& bucket = *bucketPtr;
        for (stk::mesh::Bucket::size_type k = 0; k < bucket.size(); ++k)
        {
            const stk::mesh::Entity side = bucket[k];
            const stk::mesh::Entity* nodes = bulkData.begin_nodes(side);
            const unsigned numNodes = bulkData.num_nodes(side);

            if (numNodes < 2)
            {
                continue;
            }

            const label first = surfaceIndex(nodes[0]);

            if (numNodes == 2)
            {
                // In a 2D OpenAccel build, a wall side is an edge. SSDF's
                // point-triangle kernel handles this degenerate triangle as a
                // segment.
                const label second = surfaceIndex(nodes[1]);
                s0.push_back(first);
                s1.push_back(second);
                s2.push_back(second);
                continue;
            }

            // Fan-triangulate triangles, quads, and higher-order polygonal
            // sides while preserving the STK side-node order.
            label previous = surfaceIndex(nodes[1]);
            for (unsigned n = 2; n < numNodes; ++n)
            {
                const label current = surfaceIndex(nodes[n]);
                s0.push_back(first);
                s1.push_back(previous);
                s2.push_back(current);
                previous = current;
            }
        }
    }

    // an maskedRegion is a wall the mesh does not carry: only one rank adds
    // its facets, since all ranks hold the same copy of them
    const auto* overridesPtr =
        wallDistancePtr_->realmRef().simulationRef().overridesPtr();
    const masking* solidsPtr =
        overridesPtr ? overridesPtr->maskingPtr() : nullptr;

    if (solidsPtr && messager::master())
    {
        const label offset = static_cast<label>(sx.size());

        sx.insert(sx.end(),
                  solidsPtr->surfaceX().begin(),
                  solidsPtr->surfaceX().end());
        sy.insert(sy.end(),
                  solidsPtr->surfaceY().begin(),
                  solidsPtr->surfaceY().end());
        sz.insert(sz.end(),
                  solidsPtr->surfaceZ().begin(),
                  solidsPtr->surfaceZ().end());

        for (label f = 0; f < solidsPtr->facetCount(); ++f)
        {
            s0.push_back(offset + solidsPtr->facetNode0()[f]);
            s1.push_back(offset + solidsPtr->facetNode1()[f]);
            s2.push_back(offset + solidsPtr->facetNode2()[f]);
        }
    }

    long long globalTriangleCount = static_cast<long long>(s0.size());
    MPI_Allreduce(MPI_IN_PLACE,
                  &globalTriangleCount,
                  1,
                  MPI_LONG_LONG,
                  MPI_SUM,
                  messager::comm());

    if (globalTriangleCount == 0)
    {
        errorMsg("SSDF wall distance found no active wall surface sides");
    }

    std::vector<scalar> distance(queryNodes.size(),
                                 std::numeric_limits<scalar>::max());

    int result = 0;
#ifdef SSDF_ENABLE_MPI
    result = ssdf::pedf<scalar, scalar, label>(
        messager::comm(),
        static_cast<std::ptrdiff_t>(queryNodes.size()),
        x.data(),
        y.data(),
        z.data(),
        static_cast<std::ptrdiff_t>(s0.size()),
        s0.data(),
        s1.data(),
        s2.data(),
        static_cast<std::ptrdiff_t>(sx.size()),
        sx.data(),
        sy.data(),
        sz.data(),
        distance.data());
#else
    // A non-MPI SSDF build is valid only for a single-rank OpenAccel run.
    if (messager::parallel())
    {
        errorMsg("Parallel OpenAccel requires SSDF_ENABLE_MPI");
    }

    result = ssdf::edf_select<scalar, scalar, label>(
        static_cast<std::ptrdiff_t>(queryNodes.size()),
        x.data(),
        y.data(),
        z.data(),
        static_cast<std::ptrdiff_t>(s0.size()),
        s0.data(),
        s1.data(),
        s2.data(),
        static_cast<std::ptrdiff_t>(sx.size()),
        sx.data(),
        sy.data(),
        sz.data(),
        distance.data());
#endif

    if (result != 0)
    {
        errorMsg("SSDF wall distance computation failed");
    }

    for (std::size_t i = 0; i < queryNodes.size(); ++i)
    {
        scalar* yMin = stk::mesh::field_data(*yMinField, queryNodes[i]);
        *yMin = distance[i];
    }

    if (messager::parallel())
    {
        stk::mesh::communicate_field_data(bulkData, {yMinField});
    }

    if (solidsPtr)
    {
        computeImmersedSolidDistance_(*solidsPtr, queryNodes, x, y, z);
    }
}

void signedDistanceFunction::computeImmersedSolidDistance_(
    const masking& masks,
    const std::vector<stk::mesh::Entity>& queryNodes,
    const std::vector<scalar>& x,
    const std::vector<scalar>& y,
    const std::vector<scalar>& z)
{
    auto& mesh = wallDistancePtr_->meshRef();
    stk::mesh::BulkData& bulkData = mesh.bulkDataRef();
    stk::mesh::MetaData& metaData = mesh.metaDataRef();

    STKScalarField* solidDistanceField = metaData.get_field<scalar>(
        stk::topology::NODE_RANK, masking::distance_ID);

    if (!solidDistanceField)
    {
        return;
    }

    // the body surface is replicated on every rank, so the serial kernel is
    // enough and no communication is needed
    std::vector<scalar> distance(queryNodes.size(),
                                 std::numeric_limits<scalar>::max());

    const int result = ssdf::edf_select<scalar, scalar, label>(
        static_cast<std::ptrdiff_t>(queryNodes.size()),
        x.data(),
        y.data(),
        z.data(),
        static_cast<std::ptrdiff_t>(masks.facetCount()),
        masks.facetNode0().data(),
        masks.facetNode1().data(),
        masks.facetNode2().data(),
        static_cast<std::ptrdiff_t>(masks.surfaceNodeCount()),
        masks.surfaceX().data(),
        masks.surfaceY().data(),
        masks.surfaceZ().data(),
        distance.data());

    if (result != 0)
    {
        errorMsg("SSDF maskedRegion distance computation failed");
    }

    for (std::size_t i = 0; i < queryNodes.size(); ++i)
    {
        scalar* solidDistance =
            stk::mesh::field_data(*solidDistanceField, queryNodes[i]);
        *solidDistance = distance[i];
    }

    if (messager::parallel())
    {
        stk::mesh::communicate_field_data(bulkData, {solidDistanceField});
    }
}

} // namespace accel
