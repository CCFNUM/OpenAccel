// File       : masking.h
// Created    : Tue Aug 11 2026
// Author     : Mhamad Mahdi Alloush
// Description: Manager of the solid bodies immersed in the fluid mesh
// Copyright 2026 CCFNUM HSLU T&A. All Rights Reserved.

#ifndef MASKING_H
#define MASKING_H

// code
#include "fieldBroker.h"
#include "maskedRegion.h"

#include "types.h"
#include <cstdint>

namespace accel
{

// Owns the masking and is the single place that imposes what they mask
// on the fluid state. Equations and models never touch a body directly: they
// mask through here, so a new explicit field update cannot silently escape the
// masking.
class masking : public fieldBroker
{
private:
    std::vector<std::unique_ptr<maskedRegion>> regionVector_;

    // log-law constants of the wall treatment: they belong to the treatment,
    // not to a turbulence model, so masking works for a laminar run too
    scalar kappa_ = 0.41;
    scalar B_ = 5.2;
    scalar Cmu_ = 0.09;

    // union of the bodies, in the order the equations consume them
    std::vector<stk::mesh::Entity> coveredNodes_;
    std::vector<stk::mesh::Entity> ringNodes_;

    // union of the body surfaces as one SSDF triangle soup
    std::vector<scalar> sx_;
    std::vector<scalar> sy_;
    std::vector<scalar> sz_;
    std::vector<label> s0_;
    std::vector<label> s1_;
    std::vector<label> s2_;

    // node-to-element graph of the masked domains, in the layout ssdf::layers
    // expects: one array per element-node slot, indexed by element
    std::vector<stk::mesh::Entity> graphNodes_;
    std::vector<label> graphIndex_;
    std::vector<std::vector<label>> elementNodes_;
    std::vector<label*> elementNodeColumns_;
    label nodesPerElement_ = 0;

    void buildGraph_();

    // multi-source BFS from the seeded nodes, layer 0 being the seeds. Runs
    // ssdf::layers locally and, in parallel, keeps exchanging and continuing
    // with ssdf::layers_iterative until no rank changes any more
    void computeLayers_(const std::vector<uint8_t>& seeds,
                        std::vector<uint8_t>& layers,
                        label maxLayers);

    void collect_();

    void classify_();

public:
    static constexpr char layer_ID[] = "mask_layer";
    static constexpr char id_ID[] = "mask_id";
    static constexpr char distance_ID[] = "mask_distance";

    masking(realm* realmPtr);

    // operations

    void read(const YAML::Node& maskingArrayNode);

    // declare the fields, must run before the mesh is populated
    void setup();

    void initialize();

    // re-classify, for a deforming mesh or a body that moved
    void update();

    // hold a node field at a value inside every masked region
    void setAtCovered(STKScalarField& field, scalar value);

    // copy a node field into another one inside every masked region
    void copyAtCovered(STKScalarField& target, const STKScalarField& source);

    // re-impose the velocity the regions own, after an explicit correction that
    // ran over every node
    void reapplyVelocity();

    scalar kappa() const
    {
        return kappa_;
    }

    scalar B() const
    {
        return B_;
    }

    scalar Cmu() const
    {
        return Cmu_;
    }

    void reportForces() const;

    // access

    label graphNodeCount() const
    {
        return static_cast<label>(graphNodes_.size());
    }

    stk::mesh::Entity graphNode(label i) const
    {
        return graphNodes_[i];
    }

    label graphIndexOf(const stk::mesh::Entity node) const;

    label regionCount() const
    {
        return static_cast<label>(regionVector_.size());
    }

    label size() const
    {
        return static_cast<label>(regionVector_.size());
    }

    maskedRegion& regionRef(label i)
    {
        return *regionVector_[i];
    }

    const maskedRegion& regionRef(label i) const
    {
        return *regionVector_[i];
    }

    // nodes every equation has to hold, over all bodies
    const std::vector<stk::mesh::Entity>& coveredNodes() const
    {
        return coveredNodes_;
    }

    const std::vector<stk::mesh::Entity>& ringNodes() const
    {
        return ringNodes_;
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
};

} // namespace accel

#endif // MASKING_H
