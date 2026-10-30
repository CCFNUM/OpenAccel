#include "immersedDistance.h"
#include "immersedGeometry.h"
#include "masking.h"
#include "messager.h"
#include "realm.h"

#include <stk_io/IossBridge.hpp>
#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>
#include <string>
#include <vector>

namespace accel
{

void declareImmersedDistance(realm& realm)
{
    auto& mesh = realm.meshRef();
    auto& meta = mesh.metaDataRef();

    auto& field = meta.declare_field<scalar>(
        stk::topology::NODE_RANK, "signed_distance");

    const scalar initial =
        std::numeric_limits<scalar>::quiet_NaN();

    stk::mesh::put_field_on_mesh(
        field,
        stk::mesh::selectUnion(mesh.interiorActiveParts()),
        1,
        &initial);

    stk::io::set_field_output_type(
        field, stk::io::FieldOutputType::SCALAR);
}

void computeImmersedDistance(realm& realm, const masking& maskingRef,
                             ImmersedDistanceTimings* timings)
{
    if (timings) *timings = {};
#if SPATIAL_DIM == 3 && !defined(SSDF_ENABLE_CUBIQL)
    errorMsg(
        "3D immersed signed distance requires SSDF_ENABLE_CUBIQL=ON");
#endif

    const auto* masks = &maskingRef;

    namespace geometry = immersedGeometry;

    auto& mesh = realm.meshRef();
    auto& bulk = mesh.bulkDataRef();
    auto& meta = mesh.metaDataRef();
    const auto comm = messager::comm();

    auto* coordinates = meta.get_field<scalar>(
        stk::topology::NODE_RANK, mesh::coordinates_ID);

    auto* output = meta.get_field<scalar>(
        stk::topology::NODE_RANK, "signed_distance");

    try
    {
        geometry::collectiveCheck(coordinates && output, comm,
                "Missing immersed signed-distance fields");

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
                if (timings)
                {
                    timings->distance += result.distanceSeconds;
                    timings->sign += result.signSeconds;
                    timings->combined2D += result.combined2DSeconds;
                }
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

        const double fieldStart = MPI_Wtime();
        // Write signed-distance field without modifying masking's fields.
        for (std::size_t i = 0; i < nodes.size(); ++i)
        {
            *stk::mesh::field_data(*output, nodes[i]) =
                combined.signedValue(i, true);
        }

        stk::mesh::copy_owned_to_shared(bulk, {output});
        stk::mesh::communicate_field_data(bulk, {output});
        if (timings) timings->field = MPI_Wtime() - fieldStart;
    }
    catch (const std::exception& e)
    {
        errorMsg(std::string("Masked immersed distance: ") + e.what());
    }
}

}