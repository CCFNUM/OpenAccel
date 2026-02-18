// File : zone.h
// Created : Fri Aug 25 2023 12:55:24 (+0100)
// Author : Mhamad Mahdi Alloush
// Description: Mesh zone with interior parts, boundaries, and motion
// capabilities
// Copyright (c) 2023 CCFNUM, Lucerne University of Applied
// Sciences and Arts. SPDX-License-Identifier: BSD-3-Clause

#ifndef ZONE_H
#define ZONE_H

// code
#include "mesh.h"
#include "zoneDeformation.h"
#include "zoneTransformation.h"

namespace accel
{

struct zoneStats
{
    // current
    scalar volume_ = 0.0;

    // original: required in case of deforming mesh
    scalar volume0_ = 0.0;
};

class boundary;

class zone
{
private:
    // Stats

    zoneStats stats_;

    void computeStats0_();

    void computeStats_();

protected:
    // A pointer to the mesh necessary to access bulk and meta data
    mesh* meshPtr_;

    // all zone parts
    stk::mesh::PartVector interiorParts_;

    // subset of the zone parts that will be ignored for moving
    stk::mesh::PartVector stationaryParts_;

    std::vector<std::unique_ptr<boundary>> boundaryVector_;

    // Zone index (id)
    label index_ = -1;

    std::string name_;

    // domain motion

    bool frameRotating_ = false;

    bool meshTransforming_ = false;

    mutable std::unique_ptr<zoneTransformation> transformationPtr_ = nullptr;

    // domain deformation

    bool meshDeforming_ = false;

    mutable std::unique_ptr<zoneDeformation> deformationPtr_ = nullptr;

public:
    zone(mesh* meshPtr, label index, std::string name);

    // IO

    void read(const YAML::Node& inputNode);

    // Methods

    void setup();

    void initialize();

    void update();

    // Access

    // primary attributes

    label index() const
    {
        return index_;
    }

    std::string name() const
    {
        return name_;
    }

    // interior parts of the zone
    stk::mesh::PartVector interiorParts() const
    {
        return interiorParts_;
    }

    // stationary parts of the zone (subset of interior)
    stk::mesh::PartVector stationaryParts() const
    {
        return stationaryParts_;
    }

    zoneStats stats() const
    {
        return stats_;
    }

    // For zone motion, a transient motion implies a real
    // mesh motion, but if steady-state, then a frame motion
    // is utilized (MFR)

    // domain motion
    zoneTransformation& transformationRef();

    const zoneTransformation& transformationRef() const;

    // mesh deformation
    zoneDeformation& deformationRef();

    const zoneDeformation& deformationRef() const;

    // frame rotation: steady-state only
    bool frameRotating() const
    {
        return frameRotating_;
    }

    // mesh transformation: transient only
    bool meshTransforming() const
    {
        return meshTransforming_;
    }

    void setMeshTransforming(bool state)
    {
        meshTransforming_ = state;
    }

    // mesh deformation: transient only
    bool meshDeforming() const
    {
        return meshDeforming_;
    }

    void setMeshDeforming(bool state)
    {
        meshDeforming_ = state;
    }

    // mesh moving: transient only
    bool meshMoving() const
    {
        return meshTransforming_ || meshDeforming_;
    }

    // boundaries: it contains physical type of boundary and the corresponding
    // part

    label nBoundaries() const;

    boundary* boundaryPtr(label iBoundary);

    const boundary* boundaryPtr(label iBoundary) const;

    boundary& boundaryRef(label iBoundary);

    const boundary& boundaryRef(label iBoundary) const;

    // access the mesh: bulk and meta data

    mesh* meshPtr();

    const mesh* meshPtr() const;

    mesh& meshRef();

    const mesh& meshRef() const;
};

} // namespace accel

#endif // ZONE_H
