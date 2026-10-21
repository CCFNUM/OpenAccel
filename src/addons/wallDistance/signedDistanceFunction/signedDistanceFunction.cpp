// File       : signedDistanceFunction.cpp
// Description: Wall distance computation using the SSDF library

// code
#include "signedDistanceFunction.h"
#include "immersedGeometry.h"

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
#include <stk_io/IossBridge.hpp>
#include <algorithm>
#include <cmath>
#include <cstdint>
#include <map>
#include <utility>

namespace accel
{

signedDistanceFunction::signedDistanceFunction(wallDistance* wallDistancePtr)
    : wallDistancePtr_(wallDistancePtr)
{
    declareImmersedDistanceField_();
}

void signedDistanceFunction::declareImmersedDistanceField_()
{
    const auto analysis = wallDistancePtr_->realmRef().simulationRef()
                              .getYAMLPhysicalAnalysisNode();
    const auto overrides = analysis["overrides"];
    if (!overrides || !overrides["masking"]) return;

    auto& mesh = wallDistancePtr_->meshRef();
    auto& meta = mesh.metaDataRef();
    auto& field = meta.declare_field<scalar>(
        stk::topology::NODE_RANK, immersedDistanceFieldName_);
    const scalar initial = std::numeric_limits<scalar>::quiet_NaN();
    stk::mesh::put_field_on_mesh(
        field, stk::mesh::selectUnion(mesh.interiorActiveParts()), 1, &initial);
    stk::io::set_field_output_type(field, stk::io::FieldOutputType::SCALAR);
}

void signedDistanceFunction::setup()
{
#if SPATIAL_DIM == 3 && !defined(SSDF_ENABLE_CUBIQL)
    const auto* overrides =
        wallDistancePtr_->realmRef().simulationRef().overridesPtr();
    if (overrides && overrides->maskingPtr())
        errorMsg("3D masked signed distance requires SSDF_ENABLE_CUBIQL=ON");
#endif
}

void signedDistanceFunction::initialize()
{
    computeWallDistance_();
    computeMaskedImmersedDistance_();
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
            computeMaskedImmersedDistance_();
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

}

void signedDistanceFunction::computeMaskedImmersedDistance_()
{
    const auto* overrides =
        wallDistancePtr_->realmRef().simulationRef().overridesPtr();
    const auto* masks = overrides ? overrides->maskingPtr() : nullptr;
    if (!masks) return;

    namespace geometry = immersedGeometry;
    auto& mesh = wallDistancePtr_->meshRef();
    auto& bulk = mesh.bulkDataRef();
    auto& meta = mesh.metaDataRef();
    const auto comm = messager::comm();
    auto* coordinates = meta.get_field<scalar>(
        stk::topology::NODE_RANK, mesh::coordinates_ID);
    auto* output = meta.get_field<scalar>(
        stk::topology::NODE_RANK, immersedDistanceFieldName_);
    auto* maskDistance = meta.get_field<scalar>(
        stk::topology::NODE_RANK, masking::distance_ID);

    try
    {
        geometry::collectiveCheck(coordinates && output && maskDistance, comm,
                                  "Missing masked distance fields");
        stk::mesh::PartVector parts;
        for (label r = 0; r < masks->regionCount(); ++r)
        {
            const auto regionParts = masks->regionRef(r).domainParts();
            parts.insert(parts.end(), regionParts.begin(), regionParts.end());
        }

        std::vector<stk::mesh::Entity> nodes;
        std::vector<geometry::Point> points;
        const auto selector =
            meta.locally_owned_part() & stk::mesh::selectUnion(parts);
        for (const auto* bucket :
             bulk.get_buckets(stk::topology::NODE_RANK, selector))
            for (auto node : *bucket)
            {
                const auto* xyz = stk::mesh::field_data(*coordinates, node);
                nodes.push_back(node);
#if SPATIAL_DIM == 3
                points.push_back({xyz[0], xyz[1], xyz[2]});
#else
                points.push_back({xyz[0], xyz[1], 0.0});
#endif
            }

        geometry::Combined combined(nodes.size());
        for (label r = 0; r < masks->regionCount(); ++r)
        {
            const auto& region = masks->regionRef(r);
            const auto regionSelector =
                stk::mesh::selectUnion(region.domainParts());
            bool okay = true;
            std::string detail;
            try
            {
                // Masking already holds the complete surface on every rank.
                geometry::Surface surface;
                surface.dimension = SPATIAL_DIM;
                for (label v = 0; v < region.surfaceNodeCount(); ++v)
                    surface.points.push_back({region.surfaceX()[v],
                                              region.surfaceY()[v],
                                              region.surfaceZ()[v]});
                for (label f = 0; f < region.facetCount(); ++f)
                    surface.faces.push_back({region.facetNode0()[f],
                                             region.facetNode1()[f],
                                             region.facetNode2()[f]});

                std::vector<std::size_t> indices;
                std::vector<geometry::Point> queries;
                for (std::size_t i = 0; i < nodes.size(); ++i)
                    if (regionSelector(bulk.bucket(nodes[i])))
                    {
                        indices.push_back(i);
                        queries.push_back(points[i]);
                    }
                const auto result = geometry::evaluate(surface, queries);
                for (std::size_t j = 0; j < indices.size(); ++j)
                {
                    const auto i = indices[j];
                    if (!std::isfinite(result.distance[j]) ||
                        result.distance[j] < 0 ||
                        result.distance[j] == std::numeric_limits<double>::max())
                        throw std::runtime_error("Invalid immersed distance");
                    combined.distance[i] =
                        std::min(combined.distance[i], result.distance[j]);
                    combined.inside[i] |= result.inside[j];
                }
            }
            catch (const std::exception& e)
            {
                okay = false;
                detail = e.what();
            }
            geometry::collectiveCheck(
                okay, comm, "Masked distance failed for " + region.name() +
                            (detail.empty() ? "" : ": " + detail));
        }

        // Keep masking's internal unsigned distance, and expose the signed field.
        for (std::size_t i = 0; i < nodes.size(); ++i)
        {
            *stk::mesh::field_data(*output, nodes[i]) =
                combined.signedValue(i, true);
            *stk::mesh::field_data(*maskDistance, nodes[i]) =
                combined.distance[i];
        }
        stk::mesh::copy_owned_to_shared(bulk, {output, maskDistance});
        stk::mesh::communicate_field_data(bulk, {output, maskDistance});
    }
    catch (const std::exception& e)
    {
        errorMsg(std::string("Masked immersed distance: ") + e.what());
    }
}

} // namespace accel
