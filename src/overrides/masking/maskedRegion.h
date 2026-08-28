// File       : maskedRegion.h
// Created    : Tue Aug 11 2026
// Author     : Mhamad Mahdi Alloush
// Description: A solid body immersed in the fluid mesh
// Copyright 2026 CCFNUM HSLU T&A. All Rights Reserved.

#ifndef MASKEDREGION_H
#define MASKEDREGION_H

// code
#include "types.h"

namespace accel
{

class realm;

// Classification of a fluid node with respect to an maskedRegion
enum class maskState : int
{
    fluid = 0,
    covered = 1,
    ring = 2
};

// friction velocity that reproduces a tangential velocity at a given distance
scalar maskFrictionVelocity(const scalar uTangential,
                            const scalar distance,
                            const scalar nu,
                            const scalar kappa,
                            const scalar E);

// A solid body that is not part of the fluid mesh. It masks the fluid nodes it
// covers, contributes its surface to the wall distance, and imposes a wall
// treatment on the first ring of fluid nodes around it. The body is static:
// motion (rigid or from a structural solver) enters through update().
class maskedRegion
{
private:
    realm* realmPtr_;

    std::string name_;
    label id_ = 1;

    // geometry: either a mesh of the body or an analytic region
    maskShape shape_ = maskShape::mesh;

    std::string surfaceFilePath_;
    scalar scaleFactor_ = 1.0;

    scalar boxMin_[3] = {0.0, 0.0, 0.0};
    scalar boxMax_[3] = {0.0, 0.0, 0.0};

    scalar centre_[3] = {0.0, 0.0, 0.0};
    scalar axis_[3] = {0.0, 0.0, 1.0};
    scalar radius_ = 0.0;
    scalar height_ = 0.0;

    // facets used to represent a curved analytic region
    label resolution_ = 64;

    // domains the body acts on, empty means every fluid domain
    std::vector<std::string> domainNames_;

    // layers of nodes around the region the wall treatment is imposed on, and
    // the layer the wall law is anchored on
    label forcedLayers_ = 1;
    label referenceLayer_ = 2;

    maskWallTreatment wallTreatment_ = maskWallTreatment::wallFunction;

    // body surface as an SSDF triangle soup (a 2D side is a degenerate
    // triangle)
    std::vector<scalar> sx_;
    std::vector<scalar> sy_;
    std::vector<scalar> sz_;
    std::vector<label> s0_;
    std::vector<label> s1_;
    std::vector<label> s2_;

    // surface bounding box, used to skip nodes that cannot be covered
    scalar bBoxMin_[3] = {0.0, 0.0, 0.0};
    scalar bBoxMax_[3] = {0.0, 0.0, 0.0};

    // locally owned fluid nodes covered by the body and in the first ring
    // outside it
    std::vector<stk::mesh::Entity> coveredNodes_;
    std::vector<stk::mesh::Entity> ringNodes_;

    // per ring node: unit normal into the fluid, distance to the surface, the
    // node further out that anchors the wall law and its distance
    std::vector<scalar> ringNormals_;
    std::vector<scalar> ringDistances_;
    std::vector<stk::mesh::Entity> ringReferenceNodes_;
    std::vector<scalar> ringReferenceDistances_;

    // velocity the wall treatment imposes, per ring node
    std::vector<scalar> ringVelocities_;

    // force the body exerts on the fluid, from the momentum the constraint
    // absorbs at the masked nodes
    scalar force_[3] = {0.0, 0.0, 0.0};

    void buildSurface_();

    void readSurfaceMesh_();

    void tessellateBox_();

    void tessellateSphere_();

    void tessellateCylinder_();

    // add one facet from three surface points, degenerate in 2D
    label addSurfacePoint_(scalar x, scalar y, scalar z);

    void addFacet_(label a, label b, label c);

    void computeBoundingBox_();

    void markCovered_(std::vector<uint8_t>& seeds);

    void takeLayers_(const std::vector<uint8_t>& layers);

    bool isCovered_(const scalar* point) const;

    bool isCoveredByShape_(const scalar* point) const;

    bool isCoveredByFacets_(const scalar* point) const;

    stk::mesh::Selector selectDomains_() const;

public:
    maskedRegion(realm* realmPtr, label id);

    // operations

    void read(const YAML::Node& regionNode);

    void initialize();

    // seed the covering test, then keep the nodes the layer walk assigned
    void classify(std::vector<uint8_t>& seeds, bool marking);

    void takeLayers(const std::vector<uint8_t>& layers)
    {
        takeLayers_(layers);
    }

    void markCovered(std::vector<uint8_t>& seeds)
    {
        markCovered_(seeds);
    }

    // parts of the domains this region acts on
    stk::mesh::PartVector domainParts() const;

    label forcedLayers() const
    {
        return forcedLayers_;
    }

    label referenceLayer() const
    {
        return referenceLayer_;
    }

    // closest point on the body surface: distance and unit normal into the
    // fluid
    void
    closestFacet(const scalar* point, scalar& distance, scalar* normal) const;

    // hold a node field at a constant value inside the body
    void setAtCovered(STKScalarField& field, scalar value) const;

    // copy a node field into another one inside the body
    void copyAtCovered(STKScalarField& target,
                       const STKScalarField& source) const;

    // wall-law velocity the first ring must be driven to: no penetration plus
    // the law-of-the-wall magnitude built on the reference node further out
    void computeRingVelocities(const STKScalarField& U,
                               const STKScalarField& rho,
                               const STKScalarField& mu,
                               scalar kappa,
                               scalar B);

    const std::vector<scalar>& ringVelocities() const
    {
        return ringVelocities_;
    }

    // re-impose the velocity the body owns, after any explicit correction
    void reapplyVelocity(STKScalarField& U) const;

    // friction velocity at a ring node, from its reference node
    scalar frictionVelocity(size_t i,
                            const STKScalarField& U,
                            const STKScalarField& rho,
                            const STKScalarField& mu,
                            scalar kappa,
                            scalar B) const;

    void accumulateForce(const scalar* nodeForce);

    void resetForce();

    void reduceForce();

    // access

    const std::string& name() const
    {
        return name_;
    }

    label id() const
    {
        return id_;
    }

    maskShape shape() const
    {
        return shape_;
    }

    maskWallTreatment wallTreatment() const
    {
        return wallTreatment_;
    }

    label facetCount() const
    {
        return static_cast<label>(s0_.size());
    }

    label surfaceNodeCount() const
    {
        return static_cast<label>(sx_.size());
    }

    const std::vector<scalar>& surfaceX() const
    {
        return sx_;
    }

    const std::vector<scalar>& surfaceY() const
    {
        return sy_;
    }

    const std::vector<scalar>& surfaceZ() const
    {
        return sz_;
    }

    const std::vector<label>& facetNode0() const
    {
        return s0_;
    }

    const std::vector<label>& facetNode1() const
    {
        return s1_;
    }

    const std::vector<label>& facetNode2() const
    {
        return s2_;
    }

    const std::vector<stk::mesh::Entity>& coveredNodes() const
    {
        return coveredNodes_;
    }

    const std::vector<stk::mesh::Entity>& ringNodes() const
    {
        return ringNodes_;
    }

    const scalar* ringNormal(size_t i) const
    {
        return &ringNormals_[SPATIAL_DIM * i];
    }

    scalar ringDistance(size_t i) const
    {
        return ringDistances_[i];
    }

    const scalar* force() const
    {
        return force_;
    }
};

} // namespace accel

#endif // MASKEDREGION_H
