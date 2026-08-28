// File       : maskedRegion.cpp
// Created    : Tue Aug 11 2026
// Author     : Mhamad Mahdi Alloush
// Description: A solid body immersed in the fluid mesh
// Copyright 2026 CCFNUM HSLU T&A. All Rights Reserved.

// code
#include "maskedRegion.h"
#include "domain.h"
#include "macros.h"
#include "masking.h"
#include "messager.h"
#include "realm.h"
#include "simulation.h"

#include <ssdf.hpp>

#include <unordered_map>

namespace accel
{

scalar maskFrictionVelocity(const scalar uTangential,
                            const scalar distance,
                            const scalar nu,
                            const scalar kappa,
                            const scalar E)
{
    const scalar yPlusLimit = 11.06;

    // viscous branch, which is also the starting guess
    scalar uTau = std::sqrt(nu * uTangential / std::max(distance, SMALL));

    if (uTau * distance / nu <= yPlusLimit)
    {
        return uTau;
    }

    // Newton on u = uTau ln(E yPlus) / kappa, with yPlus = uTau d / nu
    for (label iteration = 0; iteration < 20; ++iteration)
    {
        const scalar yPlus = std::max(uTau * distance / nu, yPlusLimit);
        const scalar uPlus = std::log(E * yPlus) / kappa;

        const scalar increment =
            (uTau * uPlus - uTangential) / (uPlus + 1.0 / kappa);
        uTau -= increment;

        if (uTau <= 0.0)
        {
            return std::sqrt(nu * uTangential / std::max(distance, SMALL));
        }

        if (std::abs(increment) < 1.0e-10 * std::max(uTau, SMALL))
        {
            break;
        }
    }

    return uTau;
}

maskedRegion::maskedRegion(realm* realmPtr, label id)
    : realmPtr_(realmPtr), id_(id)
{
}

void maskedRegion::read(const YAML::Node& regionNode)
{
    if (!regionNode["name"])
    {
        errorMsg("maskedRegion requires a name");
    }

    name_ = regionNode["name"].template as<std::string>();

    if (!regionNode["shape"])
    {
        errorMsg("maskedRegion `" + name_ + "` requires a shape");
    }

    const YAML::Node& shapeNode = regionNode["shape"];

    if (!shapeNode["option"])
    {
        errorMsg("maskedRegion `" + name_ +
                 "`: shape requires option: mesh, box, sphere or cylinder");
    }

    shape_ = convertMaskShapeFromString(
        shapeNode["option"].template as<std::string>());

    auto readPoint = [&](const char* key, scalar* target)
    {
        if (!shapeNode[key])
        {
            errorMsg("maskedRegion `" + name_ + "`: shape requires " + key);
        }

        label component = 0;
        for (const auto& value : shapeNode[key])
        {
            target[component++] = value.template as<scalar>();
        }
    };

    auto readScalar = [&](const char* key, scalar& target)
    {
        if (!shapeNode[key])
        {
            errorMsg("maskedRegion `" + name_ + "`: shape requires " + key);
        }

        target = shapeNode[key].template as<scalar>();
    };

    switch (shape_)
    {
        case maskShape::mesh:
            {
                if (!shapeNode["path"])
                {
                    errorMsg("maskedRegion `" + name_ +
                             "`: shape mesh requires path");
                }

                surfaceFilePath_ = shapeNode["path"].template as<std::string>();

                if (shapeNode["scale_factor"])
                {
                    scaleFactor_ =
                        shapeNode["scale_factor"].template as<scalar>();
                }
                break;
            }

        case maskShape::box:
            {
                readPoint("min", boxMin_);
                readPoint("max", boxMax_);
                break;
            }

        case maskShape::sphere:
            {
                readPoint("centre", centre_);
                readScalar("radius", radius_);
                break;
            }

        case maskShape::cylinder:
            {
                readPoint("centre", centre_);
                readScalar("radius", radius_);
#if SPATIAL_DIM == 3
                readPoint("axis", axis_);
                readScalar("height", height_);
#endif
                break;
            }
    }

    if (shapeNode["resolution"])
    {
        resolution_ = shapeNode["resolution"].template as<label>();
    }

    if (regionNode["domains"])
    {
        for (const auto& domainName : regionNode["domains"])
        {
            domainNames_.push_back(domainName.template as<std::string>());
        }
    }

    if (regionNode["forced_layers"])
    {
        forcedLayers_ = regionNode["forced_layers"].template as<label>();
        referenceLayer_ = forcedLayers_ + 1;
    }

    if (regionNode["reference_layer"])
    {
        referenceLayer_ = regionNode["reference_layer"].template as<label>();
    }

    if (regionNode["wall_treatment"])
    {
        wallTreatment_ = convertMaskWallTreatmentFromString(
            regionNode["wall_treatment"].template as<std::string>());
    }
}

stk::mesh::Selector maskedRegion::selectDomains_() const
{
    return stk::mesh::selectUnion(domainParts());
}

stk::mesh::PartVector maskedRegion::domainParts() const
{
    stk::mesh::PartVector parts;

    for (const auto& domain : realmPtr_->simulationRef().domainVector())
    {
        if (!domainNames_.empty() &&
            std::find(domainNames_.begin(),
                      domainNames_.end(),
                      domain->name()) == domainNames_.end())
        {
            continue;
        }

        if (domain->type() != domainType::fluid)
        {
            continue;
        }

        for (auto part : domain->zonePtr()->interiorParts())
        {
            parts.push_back(part);
        }
    }

    if (parts.empty())
    {
        errorMsg("masked region `" + name_ +
                 "` acts on no fluid domain: check its domains entry");
    }

    return parts;
}

void maskedRegion::initialize()
{
    buildSurface_();
}

void maskedRegion::buildSurface_()
{
    switch (shape_)
    {
        case maskShape::mesh:
            readSurfaceMesh_();
            break;
        case maskShape::box:
            tessellateBox_();
            break;
        case maskShape::sphere:
            tessellateSphere_();
            break;
        case maskShape::cylinder:
            tessellateCylinder_();
            break;
    }

    computeBoundingBox_();

    if (messager::master())
    {
        std::cout << "Masked region `" << name_ << "`: " << surfaceNodeCount()
                  << " surface nodes, " << facetCount() << " facets"
                  << std::endl;
    }
}

label maskedRegion::addSurfacePoint_(scalar x, scalar y, scalar z)
{
    const label index = static_cast<label>(sx_.size());
    sx_.push_back(x);
    sy_.push_back(y);
#if SPATIAL_DIM == 3
    sz_.push_back(z);
#else
    (void)z;
    sz_.push_back(0.0);
#endif
    return index;
}

void maskedRegion::addFacet_(label a, label b, label c)
{
    s0_.push_back(a);
    s1_.push_back(b);
    s2_.push_back(c);
}

void maskedRegion::tessellateBox_()
{
#if SPATIAL_DIM == 2
    const label c0 = addSurfacePoint_(boxMin_[0], boxMin_[1], 0.0);
    const label c1 = addSurfacePoint_(boxMax_[0], boxMin_[1], 0.0);
    const label c2 = addSurfacePoint_(boxMax_[0], boxMax_[1], 0.0);
    const label c3 = addSurfacePoint_(boxMin_[0], boxMax_[1], 0.0);

    // a 2D facet is an edge, kept as a degenerate triangle
    addFacet_(c0, c1, c1);
    addFacet_(c1, c2, c2);
    addFacet_(c2, c3, c3);
    addFacet_(c3, c0, c0);
#else
    label corner[8];
    label n = 0;
    for (label k = 0; k < 2; ++k)
    {
        for (label j = 0; j < 2; ++j)
        {
            for (label i = 0; i < 2; ++i)
            {
                corner[n++] = addSurfacePoint_(i ? boxMax_[0] : boxMin_[0],
                                               j ? boxMax_[1] : boxMin_[1],
                                               k ? boxMax_[2] : boxMin_[2]);
            }
        }
    }

    // two triangles per face, corner index bits are (k j i)
    const label faces[6][4] = {{0, 1, 3, 2},
                               {4, 5, 7, 6},
                               {0, 1, 5, 4},
                               {2, 3, 7, 6},
                               {0, 2, 6, 4},
                               {1, 3, 7, 5}};

    for (const auto& f : faces)
    {
        addFacet_(corner[f[0]], corner[f[1]], corner[f[2]]);
        addFacet_(corner[f[0]], corner[f[2]], corner[f[3]]);
    }
#endif
}

void maskedRegion::tessellateSphere_()
{
    const scalar pi = std::acos(scalar(-1));

#if SPATIAL_DIM == 2
    // a circle, closed with degenerate triangles
    std::vector<label> ring(resolution_);
    for (label i = 0; i < resolution_; ++i)
    {
        const scalar angle = 2.0 * pi * scalar(i) / scalar(resolution_);
        ring[i] = addSurfacePoint_(centre_[0] + radius_ * std::cos(angle),
                                   centre_[1] + radius_ * std::sin(angle),
                                   0.0);
    }

    for (label i = 0; i < resolution_; ++i)
    {
        const label next = ring[(i + 1) % resolution_];
        addFacet_(ring[i], next, next);
    }
#else
    const label nPhi = resolution_;
    const label nTheta = std::max<label>(resolution_ / 2, 2);

    std::vector<std::vector<label>> grid(nTheta + 1);
    for (label t = 0; t <= nTheta; ++t)
    {
        const scalar theta = pi * scalar(t) / scalar(nTheta);
        grid[t].resize(nPhi);
        for (label p = 0; p < nPhi; ++p)
        {
            const scalar phi = 2.0 * pi * scalar(p) / scalar(nPhi);
            grid[t][p] = addSurfacePoint_(
                centre_[0] + radius_ * std::sin(theta) * std::cos(phi),
                centre_[1] + radius_ * std::sin(theta) * std::sin(phi),
                centre_[2] + radius_ * std::cos(theta));
        }
    }

    for (label t = 0; t < nTheta; ++t)
    {
        for (label p = 0; p < nPhi; ++p)
        {
            const label q = (p + 1) % nPhi;
            addFacet_(grid[t][p], grid[t + 1][p], grid[t + 1][q]);
            addFacet_(grid[t][p], grid[t + 1][q], grid[t][q]);
        }
    }
#endif
}

void maskedRegion::tessellateCylinder_()
{
#if SPATIAL_DIM == 2
    // in 2D a cylinder is the circle of its cross-section
    tessellateSphere_();
#else
    const scalar pi = std::acos(scalar(-1));

    // orthonormal frame around the axis
    scalar e3[3] = {axis_[0], axis_[1], axis_[2]};
    const scalar norm =
        std::sqrt(e3[0] * e3[0] + e3[1] * e3[1] + e3[2] * e3[2]);
    if (norm < SMALL)
    {
        errorMsg("maskedRegion `" + name_ + "`: cylinder axis is degenerate");
    }
    for (label d = 0; d < 3; ++d)
    {
        e3[d] /= norm;
    }

    scalar helper[3] = {1.0, 0.0, 0.0};
    if (std::abs(e3[0]) > 0.9)
    {
        helper[0] = 0.0;
        helper[1] = 1.0;
    }

    scalar e1[3] = {helper[1] * e3[2] - helper[2] * e3[1],
                    helper[2] * e3[0] - helper[0] * e3[2],
                    helper[0] * e3[1] - helper[1] * e3[0]};
    const scalar e1Norm =
        std::sqrt(e1[0] * e1[0] + e1[1] * e1[1] + e1[2] * e1[2]);
    for (label d = 0; d < 3; ++d)
    {
        e1[d] /= e1Norm;
    }

    const scalar e2[3] = {e3[1] * e1[2] - e3[2] * e1[1],
                          e3[2] * e1[0] - e3[0] * e1[2],
                          e3[0] * e1[1] - e3[1] * e1[0]};

    std::vector<label> bottom(resolution_), top(resolution_);
    for (label i = 0; i < resolution_; ++i)
    {
        const scalar angle = 2.0 * pi * scalar(i) / scalar(resolution_);
        const scalar cosine = radius_ * std::cos(angle);
        const scalar sine = radius_ * std::sin(angle);

        scalar low[3], high[3];
        for (label d = 0; d < 3; ++d)
        {
            const scalar offset = cosine * e1[d] + sine * e2[d];
            low[d] = centre_[d] + offset - 0.5 * height_ * e3[d];
            high[d] = centre_[d] + offset + 0.5 * height_ * e3[d];
        }

        bottom[i] = addSurfacePoint_(low[0], low[1], low[2]);
        top[i] = addSurfacePoint_(high[0], high[1], high[2]);
    }

    const label bottomCentre =
        addSurfacePoint_(centre_[0] - 0.5 * height_ * e3[0],
                         centre_[1] - 0.5 * height_ * e3[1],
                         centre_[2] - 0.5 * height_ * e3[2]);
    const label topCentre =
        addSurfacePoint_(centre_[0] + 0.5 * height_ * e3[0],
                         centre_[1] + 0.5 * height_ * e3[1],
                         centre_[2] + 0.5 * height_ * e3[2]);

    for (label i = 0; i < resolution_; ++i)
    {
        const label j = (i + 1) % resolution_;

        addFacet_(bottom[i], top[i], top[j]);
        addFacet_(bottom[i], top[j], bottom[j]);

        addFacet_(bottomCentre, bottom[j], bottom[i]);
        addFacet_(topCentre, top[i], top[j]);
    }
#endif
}

void maskedRegion::readSurfaceMesh_()
{
    // every rank reads the whole body: the surface is small and a local copy
    // removes any communication from the covering and closest-facet queries
    stk::io::StkMeshIoBroker broker(MPI_COMM_SELF);
    broker.add_mesh_database(surfaceFilePath_, stk::io::READ_MESH);
    broker.create_input_mesh();

    stk::mesh::MetaData& meta = broker.meta_data();

    if (static_cast<label>(meta.spatial_dimension()) != SPATIAL_DIM)
    {
        errorMsg("maskedRegion `" + name_ + "`: mesh `" + surfaceFilePath_ +
                 "` has a spatial dimension that does not match the build");
    }

    stk::mesh::Part& skinPart =
        meta.declare_part("masked_region_skin", meta.side_rank());

    broker.populate_bulk_data();

    stk::mesh::BulkData& bulk = broker.bulk_data();

    // a body given as a volume (2D: area) mesh is skinned, a body already given
    // as a surface mesh is taken as is
    bool isVolumeMesh = false;
    for (const stk::mesh::Bucket* bucketPtr :
         bulk.get_buckets(stk::topology::ELEMENT_RANK, meta.universal_part()))
    {
        if (static_cast<label>(bucketPtr->topology().dimension()) ==
            SPATIAL_DIM)
        {
            isVolumeMesh = true;
            break;
        }
    }

    stk::mesh::EntityRank facetRank = stk::topology::ELEMENT_RANK;
    stk::mesh::Selector selFacets = meta.universal_part();

    if (isVolumeMesh)
    {
        stk::mesh::create_exposed_block_boundary_sides(
            bulk, meta.universal_part(), {&skinPart});
        facetRank = meta.side_rank();
        selFacets = skinPart;
    }

    const stk::mesh::FieldBase* coordinates = meta.coordinate_field();

    std::unordered_map<stk::mesh::EntityId, label> surfaceNodeToIndex;

    auto surfaceIndex = [&](const stk::mesh::Entity node) -> label
    {
        const stk::mesh::EntityId id = bulk.identifier(node);
        const auto found = surfaceNodeToIndex.find(id);
        if (found != surfaceNodeToIndex.end())
        {
            return found->second;
        }

        const label index = static_cast<label>(sx_.size());
        surfaceNodeToIndex.emplace(id, index);

        const scalar* xyz = static_cast<const scalar*>(
            stk::mesh::field_data(*coordinates, node));
        sx_.push_back(scaleFactor_ * xyz[0]);
        sy_.push_back(scaleFactor_ * xyz[1]);
#if SPATIAL_DIM == 3
        sz_.push_back(scaleFactor_ * xyz[2]);
#else
        sz_.push_back(0.0);
#endif
        return index;
    };

    for (const stk::mesh::Bucket* bucketPtr :
         bulk.get_buckets(facetRank, selFacets))
    {
        const stk::mesh::Bucket& bucket = *bucketPtr;
        for (stk::mesh::Bucket::size_type k = 0; k < bucket.size(); ++k)
        {
            const stk::mesh::Entity facet = bucket[k];
            const stk::mesh::Entity* nodes = bulk.begin_nodes(facet);
            const unsigned numNodes = bulk.num_nodes(facet);

            if (numNodes < 2)
            {
                continue;
            }

            const label first = surfaceIndex(nodes[0]);

            if (numNodes == 2)
            {
                // a 2D facet is an edge, kept as a degenerate triangle
                const label second = surfaceIndex(nodes[1]);
                s0_.push_back(first);
                s1_.push_back(second);
                s2_.push_back(second);
                continue;
            }

            label previous = surfaceIndex(nodes[1]);
            for (unsigned n = 2; n < numNodes; ++n)
            {
                const label current = surfaceIndex(nodes[n]);
                s0_.push_back(first);
                s1_.push_back(previous);
                s2_.push_back(current);
                previous = current;
            }
        }
    }

    if (s0_.empty())
    {
        errorMsg("maskedRegion `" + name_ + "`: mesh `" + surfaceFilePath_ +
                 "` produced no surface facets");
    }
}

void maskedRegion::computeBoundingBox_()
{
    bBoxMin_[0] = *std::min_element(sx_.begin(), sx_.end());
    bBoxMin_[1] = *std::min_element(sy_.begin(), sy_.end());
    bBoxMin_[2] = *std::min_element(sz_.begin(), sz_.end());
    bBoxMax_[0] = *std::max_element(sx_.begin(), sx_.end());
    bBoxMax_[1] = *std::max_element(sy_.begin(), sy_.end());
    bBoxMax_[2] = *std::max_element(sz_.begin(), sz_.end());

    // a curved region is tessellated inside its true surface, so grow the box
    // back to the analytic extent before the margin is added
    if (shape_ == maskShape::sphere || shape_ == maskShape::cylinder)
    {
        for (label d = 0; d < 3; ++d)
        {
            bBoxMin_[d] = std::min(bBoxMin_[d], centre_[d] - radius_);
            bBoxMax_[d] = std::max(bBoxMax_[d], centre_[d] + radius_);
        }
    }

    // a margin keeps a point exactly on the box from being rejected
    for (label d = 0; d < 3; ++d)
    {
        const scalar margin =
            SMALL * (1.0 + std::abs(bBoxMax_[d] - bBoxMin_[d]));
        bBoxMin_[d] -= margin;
        bBoxMax_[d] += margin;
    }
}

bool maskedRegion::isCovered_(const scalar* point) const
{
    // outside the bounding box nothing can be covered, whatever the shape
    const scalar px = point[0];
    const scalar py = point[1];
#if SPATIAL_DIM == 3
    const scalar pz = point[2];
#else
    const scalar pz = 0.0;
#endif

    if (px < bBoxMin_[0] || px > bBoxMax_[0] || py < bBoxMin_[1] ||
        py > bBoxMax_[1] || pz < bBoxMin_[2] || pz > bBoxMax_[2])
    {
        return false;
    }

    // an analytic region is tested exactly, a meshed body by ray casting on
    // its facets
    return (shape_ == maskShape::mesh) ? isCoveredByFacets_(point)
                                       : isCoveredByShape_(point);
}

bool maskedRegion::isCoveredByShape_(const scalar* point) const
{
    const scalar px = point[0];
    const scalar py = point[1];
#if SPATIAL_DIM == 3
    const scalar pz = point[2];
#else
    const scalar pz = 0.0;
#endif

    switch (shape_)
    {
        case maskShape::box:
            {
                const scalar p[3] = {px, py, pz};
                for (label d = 0; d < SPATIAL_DIM; ++d)
                {
                    if (p[d] < boxMin_[d] || p[d] > boxMax_[d])
                    {
                        return false;
                    }
                }
                return true;
            }

        case maskShape::sphere:
            {
                scalar distanceSq = (px - centre_[0]) * (px - centre_[0]) +
                                    (py - centre_[1]) * (py - centre_[1]);
#if SPATIAL_DIM == 3
                distanceSq += (pz - centre_[2]) * (pz - centre_[2]);
#endif
                return distanceSq <= radius_ * radius_;
            }

        case maskShape::cylinder:
            {
#if SPATIAL_DIM == 2
                const scalar distanceSq =
                    (px - centre_[0]) * (px - centre_[0]) +
                    (py - centre_[1]) * (py - centre_[1]);
                return distanceSq <= radius_ * radius_;
#else
                scalar e3[3] = {axis_[0], axis_[1], axis_[2]};
                const scalar norm =
                    std::sqrt(e3[0] * e3[0] + e3[1] * e3[1] + e3[2] * e3[2]);
                for (label d = 0; d < 3; ++d)
                {
                    e3[d] /= std::max(norm, SMALL);
                }

                const scalar r[3] = {
                    px - centre_[0], py - centre_[1], pz - centre_[2]};
                const scalar along = r[0] * e3[0] + r[1] * e3[1] + r[2] * e3[2];

                if (std::abs(along) > 0.5 * height_)
                {
                    return false;
                }

                scalar radialSq = 0.0;
                for (label d = 0; d < 3; ++d)
                {
                    const scalar radial = r[d] - along * e3[d];
                    radialSq += radial * radial;
                }
                return radialSq <= radius_ * radius_;
#endif
            }

        case maskShape::mesh:
            return isCoveredByFacets_(point);
    }

    return false;
}

bool maskedRegion::isCoveredByFacets_(const scalar* point) const
{
    const scalar px = point[0];
    const scalar py = point[1];
#if SPATIAL_DIM == 3
    const scalar pz = point[2];
#else
    const scalar pz = 0.0;
#endif

    if (px < bBoxMin_[0] || px > bBoxMax_[0] || py < bBoxMin_[1] ||
        py > bBoxMax_[1] || pz < bBoxMin_[2] || pz > bBoxMax_[2])
    {
        return false;
    }

    // odd number of crossings of a ray cast along +x means the point is inside
    label crossings = 0;

#if SPATIAL_DIM == 2
    for (size_t f = 0; f < s0_.size(); ++f)
    {
        const scalar ax = sx_[s0_[f]], ay = sy_[s0_[f]];
        const scalar bx = sx_[s1_[f]], by = sy_[s1_[f]];

        if ((ay > py) == (by > py))
        {
            continue;
        }

        const scalar t = (py - ay) / (by - ay);
        if (ax + t * (bx - ax) > px)
        {
            ++crossings;
        }
    }
#else
    // slightly tilted ray so that vertex and edge hits stay unlikely
    const scalar dx = 1.0;
    const scalar dy = 1.0e-3;
    const scalar dz = 1.7e-3;

    for (size_t f = 0; f < s0_.size(); ++f)
    {
        const scalar ax = sx_[s0_[f]], ay = sy_[s0_[f]], az = sz_[s0_[f]];
        const scalar bx = sx_[s1_[f]], by = sy_[s1_[f]], bz = sz_[s1_[f]];
        const scalar cx = sx_[s2_[f]], cy = sy_[s2_[f]], cz = sz_[s2_[f]];

        // Moeller-Trumbore
        const scalar e1x = bx - ax, e1y = by - ay, e1z = bz - az;
        const scalar e2x = cx - ax, e2y = cy - ay, e2z = cz - az;

        const scalar hx = dy * e2z - dz * e2y;
        const scalar hy = dz * e2x - dx * e2z;
        const scalar hz = dx * e2y - dy * e2x;

        const scalar det = e1x * hx + e1y * hy + e1z * hz;
        if (std::abs(det) < SMALL)
        {
            continue;
        }

        const scalar invDet = 1.0 / det;
        const scalar spx = px - ax, spy = py - ay, spz = pz - az;

        const scalar u = invDet * (spx * hx + spy * hy + spz * hz);
        if (u < 0.0 || u > 1.0)
        {
            continue;
        }

        const scalar qx = spy * e1z - spz * e1y;
        const scalar qy = spz * e1x - spx * e1z;
        const scalar qz = spx * e1y - spy * e1x;

        const scalar v = invDet * (dx * qx + dy * qy + dz * qz);
        if (v < 0.0 || u + v > 1.0)
        {
            continue;
        }

        const scalar t = invDet * (e2x * qx + e2y * qy + e2z * qz);
        if (t > 0.0)
        {
            ++crossings;
        }
    }
#endif

    return (crossings % 2) == 1;
}

void maskedRegion::closestFacet(const scalar* point,
                                scalar& distance,
                                scalar* normal) const
{
    const scalar px = point[0];
    const scalar py = point[1];
#if SPATIAL_DIM == 3
    const scalar pz = point[2];
#else
    const scalar pz = 0.0;
#endif

    scalar bestDistanceSq = std::numeric_limits<scalar>::max();
    scalar bestPoint[3] = {0.0, 0.0, 0.0};

    for (size_t f = 0; f < s0_.size(); ++f)
    {
        scalar qx = 0.0, qy = 0.0, qz = 0.0;
        ssdf::point_triangle_closest_point<scalar>(px,
                                                   py,
                                                   pz,
                                                   sx_[s0_[f]],
                                                   sy_[s0_[f]],
                                                   sz_[s0_[f]],
                                                   sx_[s1_[f]],
                                                   sy_[s1_[f]],
                                                   sz_[s1_[f]],
                                                   sx_[s2_[f]],
                                                   sy_[s2_[f]],
                                                   sz_[s2_[f]],
                                                   &qx,
                                                   &qy,
                                                   &qz);

        const scalar distanceSq = (px - qx) * (px - qx) +
                                  (py - qy) * (py - qy) + (pz - qz) * (pz - qz);

        if (distanceSq < bestDistanceSq)
        {
            bestDistanceSq = distanceSq;
            bestPoint[0] = qx;
            bestPoint[1] = qy;
            bestPoint[2] = qz;
        }
    }

    distance = std::sqrt(bestDistanceSq);

    // direction away from the body, taken from the closest point itself so that
    // no facet orientation is required
    const scalar norm = std::max(distance, SMALL);
    normal[0] = (px - bestPoint[0]) / norm;
    normal[1] = (py - bestPoint[1]) / norm;
#if SPATIAL_DIM == 3
    normal[2] = (pz - bestPoint[2]) / norm;
#endif
}

void maskedRegion::markCovered_(std::vector<uint8_t>& seeds)
{
    auto& mesh = realmPtr_->meshRef();
    const auto& bulkData = mesh.bulkDataRef();
    auto& metaData = mesh.metaDataRef();

    const STKScalarField* coordinates = metaData.get_field<scalar>(
        stk::topology::NODE_RANK, mesh::coordinates_ID);

    const auto* masks = realmPtr_->simulationRef().overridesPtr()->maskingPtr();

    // covering is a purely geometric test, so every rank agrees on it for the
    // nodes it owns as well as for the ones it only shares
    for (label i = 0; i < masks->graphNodeCount(); ++i)
    {
        const stk::mesh::Entity node = masks->graphNode(i);
        const scalar* xyz = stk::mesh::field_data(*coordinates, node);

        seeds[i] = isCovered_(xyz) ? 0 : 255;
    }
}

void maskedRegion::takeLayers_(const std::vector<uint8_t>& layers)
{
    auto& mesh = realmPtr_->meshRef();
    const auto& bulkData = mesh.bulkDataRef();
    auto& metaData = mesh.metaDataRef();

    const STKScalarField* coordinates = metaData.get_field<scalar>(
        stk::topology::NODE_RANK, mesh::coordinates_ID);

    const auto* masks = realmPtr_->simulationRef().overridesPtr()->maskingPtr();

    coveredNodes_.clear();
    ringNodes_.clear();
    ringNormals_.clear();
    ringDistances_.clear();
    ringReferenceNodes_.clear();
    ringReferenceDistances_.clear();
    ringVelocities_.clear();

    const bool owned = true;

    for (label i = 0; i < masks->graphNodeCount(); ++i)
    {
        const stk::mesh::Entity node = masks->graphNode(i);

        // only the equations' own rows are constrained, so ghosts are skipped
        if (owned && !bulkData.bucket(node).owned())
        {
            continue;
        }

        const label layer = static_cast<label>(layers[i]);

        if (layer == 0)
        {
            coveredNodes_.push_back(node);
            continue;
        }

        if (layer > forcedLayers_)
        {
            continue;
        }

        const scalar* xyz = stk::mesh::field_data(*coordinates, node);

        scalar distance = 0.0;
        scalar normal[3] = {0.0, 0.0, 0.0};
        closestFacet(xyz, distance, normal);

        ringNodes_.push_back(node);
        ringDistances_.push_back(distance);
        for (label d = 0; d < SPATIAL_DIM; ++d)
        {
            ringNormals_.push_back(normal[d]);
        }

        // anchor the wall law on the reference layer, taking the neighbour of
        // that layer which sits furthest from the region
        stk::mesh::Entity reference = node;
        scalar referenceDistance = distance;

        const stk::mesh::Entity* elements = bulkData.begin_elements(node);
        const unsigned numElements = bulkData.num_elements(node);

        for (unsigned e = 0; e < numElements; ++e)
        {
            const stk::mesh::Entity* nodes = bulkData.begin_nodes(elements[e]);
            const unsigned numNodes = bulkData.num_nodes(elements[e]);

            for (unsigned n = 0; n < numNodes; ++n)
            {
                const label index = masks->graphIndexOf(nodes[n]);
                if (index < 0)
                {
                    continue;
                }

                if (static_cast<label>(layers[index]) != referenceLayer_)
                {
                    continue;
                }

                const scalar* neighbourCoordinates =
                    stk::mesh::field_data(*coordinates, nodes[n]);

                scalar neighbourDistance = 0.0;
                scalar neighbourNormal[3] = {0.0, 0.0, 0.0};
                closestFacet(
                    neighbourCoordinates, neighbourDistance, neighbourNormal);

                if (neighbourDistance > referenceDistance)
                {
                    referenceDistance = neighbourDistance;
                    reference = nodes[n];
                }
            }
        }

        ringReferenceNodes_.push_back(reference);
        ringReferenceDistances_.push_back(
            std::max(referenceDistance, distance));
    }

    ringVelocities_.assign(SPATIAL_DIM * ringNodes_.size(), 0.0);

    label counts[2] = {static_cast<label>(coveredNodes_.size()),
                       static_cast<label>(ringNodes_.size())};
    if (messager::parallel())
    {
        MPI_Allreduce(
            MPI_IN_PLACE, counts, 2, MPI_INT, MPI_SUM, messager::comm());
    }

    if (messager::master())
    {
        std::cout << "Masked region `" << name_ << "`: " << counts[0]
                  << " covered nodes, " << counts[1] << " nodes in "
                  << forcedLayers_ << " forced layer(s) around it" << std::endl;

        if (counts[0] == 0)
        {
            warningMsg("masked region `" + name_ +
                       "` covers no node: it is smaller than the local cell "
                       "spacing and masks nothing");
        }
    }
}

void maskedRegion::setAtCovered(STKScalarField& field, scalar value) const
{
    const auto& bulkData = realmPtr_->meshRef().bulkDataRef();

    for (const stk::mesh::Entity node : coveredNodes_)
    {
        scalar* data = stk::mesh::field_data(field, node);
        const unsigned numComponents =
            stk::mesh::field_scalars_per_entity(field, bulkData.bucket(node));

        for (unsigned c = 0; c < numComponents; ++c)
        {
            data[c] = value;
        }
    }
}

void maskedRegion::copyAtCovered(STKScalarField& target,
                                 const STKScalarField& source) const
{
    const auto& bulkData = realmPtr_->meshRef().bulkDataRef();

    for (const stk::mesh::Entity node : coveredNodes_)
    {
        scalar* to = stk::mesh::field_data(target, node);
        const scalar* from = stk::mesh::field_data(source, node);
        const unsigned numComponents =
            stk::mesh::field_scalars_per_entity(target, bulkData.bucket(node));

        for (unsigned c = 0; c < numComponents; ++c)
        {
            to[c] = from[c];
        }
    }
}

scalar maskedRegion::frictionVelocity(size_t i,
                                      const STKScalarField& U,
                                      const STKScalarField& rho,
                                      const STKScalarField& mu,
                                      scalar kappa,
                                      scalar B) const
{
    const stk::mesh::Entity node = ringNodes_[i];
    const stk::mesh::Entity reference = ringReferenceNodes_[i];

    const scalar density = *stk::mesh::field_data(rho, node);
    const scalar viscosity = *stk::mesh::field_data(mu, node);
    const scalar nu = viscosity / std::max(density, SMALL);

    const scalar* referenceVelocity = stk::mesh::field_data(U, reference);
    const scalar* normal = ringNormal(i);

    scalar normalComponent = 0.0;
    for (label d = 0; d < SPATIAL_DIM; ++d)
    {
        normalComponent += referenceVelocity[d] * normal[d];
    }

    scalar tangentialMagnitude = 0.0;
    for (label d = 0; d < SPATIAL_DIM; ++d)
    {
        const scalar t = referenceVelocity[d] - normalComponent * normal[d];
        tangentialMagnitude += t * t;
    }
    tangentialMagnitude = std::sqrt(tangentialMagnitude);

    if (tangentialMagnitude < SMALL)
    {
        return 0.0;
    }

    return maskFrictionVelocity(tangentialMagnitude,
                                ringReferenceDistances_[i],
                                nu,
                                kappa,
                                std::exp(kappa * B));
}

void maskedRegion::computeRingVelocities(const STKScalarField& U,
                                         const STKScalarField& rho,
                                         const STKScalarField& mu,
                                         scalar kappa,
                                         scalar B)
{
    const scalar E = std::exp(kappa * B);
    const scalar yPlusLimit = 11.06;

    ringVelocities_.assign(SPATIAL_DIM * ringNodes_.size(), 0.0);

    for (size_t i = 0; i < ringNodes_.size(); ++i)
    {
        // no slip holds the ring at rest, the wall law puts it on the profile
        // the outer flow implies
        if (wallTreatment_ != maskWallTreatment::wallFunction)
        {
            continue;
        }

        const stk::mesh::Entity node = ringNodes_[i];
        const stk::mesh::Entity reference = ringReferenceNodes_[i];

        const scalar* referenceVelocity = stk::mesh::field_data(U, reference);
        const scalar* normal = ringNormal(i);

        scalar normalComponent = 0.0;
        for (label d = 0; d < SPATIAL_DIM; ++d)
        {
            normalComponent += referenceVelocity[d] * normal[d];
        }

        scalar tangential[SPATIAL_DIM];
        scalar tangentialMagnitude = 0.0;
        for (label d = 0; d < SPATIAL_DIM; ++d)
        {
            tangential[d] = referenceVelocity[d] - normalComponent * normal[d];
            tangentialMagnitude += tangential[d] * tangential[d];
        }
        tangentialMagnitude = std::sqrt(tangentialMagnitude);

        if (tangentialMagnitude < SMALL)
        {
            continue;
        }

        const scalar density = *stk::mesh::field_data(rho, node);
        const scalar viscosity = *stk::mesh::field_data(mu, node);
        const scalar nu = viscosity / std::max(density, SMALL);

        const scalar uTau = maskFrictionVelocity(
            tangentialMagnitude, ringReferenceDistances_[i], nu, kappa, E);

        const scalar yPlus = uTau * std::max(ringDistances_[i], SMALL) / nu;
        const scalar uPlus =
            (yPlus <= yPlusLimit) ? yPlus : std::log(E * yPlus) / kappa;

        for (label d = 0; d < SPATIAL_DIM; ++d)
        {
            ringVelocities_[SPATIAL_DIM * i + d] =
                uTau * uPlus * tangential[d] / tangentialMagnitude;
        }
    }
}

void maskedRegion::reapplyVelocity(STKScalarField& U) const
{
    for (const stk::mesh::Entity node : coveredNodes_)
    {
        scalar* velocity = stk::mesh::field_data(U, node);
        for (label d = 0; d < SPATIAL_DIM; ++d)
        {
            velocity[d] = 0.0;
        }
    }

    for (size_t i = 0; i < ringNodes_.size(); ++i)
    {
        scalar* velocity = stk::mesh::field_data(U, ringNodes_[i]);
        for (label d = 0; d < SPATIAL_DIM; ++d)
        {
            velocity[d] = ringVelocities_[SPATIAL_DIM * i + d];
        }
    }
}

void maskedRegion::resetForce()
{
    for (label d = 0; d < 3; ++d)
    {
        force_[d] = 0.0;
    }
}

void maskedRegion::accumulateForce(const scalar* nodeForce)
{
    for (label d = 0; d < SPATIAL_DIM; ++d)
    {
        force_[d] += nodeForce[d];
    }
}

void maskedRegion::reduceForce()
{
    if (messager::parallel())
    {
        MPI_Allreduce(
            MPI_IN_PLACE, force_, 3, MPI_DOUBLE, MPI_SUM, messager::comm());
    }
}

} // namespace accel
