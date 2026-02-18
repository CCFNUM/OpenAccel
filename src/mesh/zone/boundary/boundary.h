// File : boundary.h
// Created : Fri Aug 25 2023 12:55:24 (+0100)
// Author : Mhamad Mahdi Alloush
// Description: Boundary face representation with physical type and reference
// frame
// Copyright (c) 2023 CCFNUM, Lucerne University of Applied Sciences and
// Arts. SPDX-License-Identifier: BSD-3-Clause

#ifndef BOUNDARY_H
#define BOUNDARY_H

// code
#include "zone.h"

namespace accel
{

struct boundaryStats
{
    // current
    scalar area_ = 0.0;
    utils::vector centroid_;

    // original: required in case of moving mesh (deforming/transforming)
    scalar area0_ = 0.0;
    utils::vector centroid0_;
};

class boundary
{
private:
    zone* zonePtr_;

    std::string name_;

    // Local index of boundary in the containing zone
    label index_;

    boundaryPhysicalType physicalType_ = boundaryPhysicalType::symmetry;

    stk::mesh::PartVector parts_;

    // In case of a rotating zone (rigid rotational motion), the frame type of
    // the boundary can be either stationary or rotating. If stationary, this
    // means that quantities assigned for this boundary are relative to the
    // absolute (stationary) frame of reference of the zone containing this
    // boundary. If the frame type is rotating, then the quantities on the
    // boundary are relative to the rotating frame of reference of the zone
    // containing the boundary. Note: a translating frame of reference does not
    // require any transformation
    boundaryRelativeFrameType frameType_ = boundaryRelativeFrameType::absolute;

    // Stats

    boundaryStats stats_;

    void computeStats0_();

    void computeStats_();

public:
    // Constructors/Destructors

    boundary(zone* zonePtr, std::string name, label index);

    virtual ~boundary()
    {
    }

    // IO

    void read(const YAML::Node& node);

    // Methods

    void setup();

    void initialize();

    void update();

    // Access

    std::string name() const
    {
        return name_;
    };

    label index() const
    {
        return index_;
    };

    boundaryPhysicalType type() const
    {
        return physicalType_;
    };

    stk::mesh::PartVector parts() const
    {
        return parts_;
    };

    boundaryRelativeFrameType frameType() const
    {
        return frameType_;
    };

    boundaryStats stats() const
    {
        return stats_;
    }

    const zone* zonePtr() const
    {
        return zonePtr_;
    }
};

class wall : public boundary
{
private:
    bool counterRotating_ = false;

    bool wallVelocityRelativeToMeshMotion_ = true;

public:
    // Constructors
    using boundary::boundary;

    // Getters and Setters

    bool counterRotating() const
    {
        return counterRotating_;
    };

    void setCounterRotating(bool state)
    {
        counterRotating_ = state;
    };

    bool wallVelocityRelativeToMeshMotion() const
    {
        return wallVelocityRelativeToMeshMotion_;
    };

    void setWallVelocityRelativeToMeshMotion(bool state)
    {
        wallVelocityRelativeToMeshMotion_ = state;
    };
};

} // namespace accel

#endif // BOUNDARY_H
