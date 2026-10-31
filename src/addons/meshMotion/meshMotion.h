// File       : meshMotion.h
// Created    : Fri Dec 13 2024 12:55:24 (+0100)
// Author     : Mhamad Mahdi Alloush
// Description: Mesh motion manager combining deformation and transformation
// Copyright 2024 CCFNUM HSLU T&A. All Rights Reserved.

#ifndef MESHMOTION_H
#define MESHMOTION_H

// code
#include "deformation.h"
#include "fieldBroker.h"
#include "transformation.h"

namespace accel
{

class meshMotion : public fieldBroker
{
private:
    std::unique_ptr<deformation> deformationPtr_ = nullptr;

    std::unique_ptr<transformation> transformationPtr_ = nullptr;

    // aWallBc holds a previous a_n for rn_accel_relaxation
    bool hasAWallBc_ = false;

    // a solved boundary displacement is newer than the last update()
    bool boundaryMoved_ = false;

public:
    // Constructors
    meshMotion(realm* realm);

    using fieldBroker::DRef;
    using fieldBroker::DtRef;

    // Methods

    void reset();

    void setup();

    void initialize();

    void update();

    // true if fluid-solid interfaces move the mesh: it follows the solution
    bool followsSolution() const;

    bool boundaryMoved() const
    {
        return boundaryMoved_;
    }

    void markBoundaryMoved()
    {
        boundaryMoved_ = true;
    }

protected:
    // protected update methods
    void updateCoordinates_();

    void updateMeshVelocityField_();

    // exact w x r on a rigidly rotating zone (the BDF chord velocity is not
    // divergence-free)
    void updateRigidRotationMeshVelocity_(label iZone);

    void updateMeshVelocityDivergenceField_();

    // no-slip wall velocity from the Dt history at the fluid-solid BIPs
    void updateFsiWallVelocity_();

    // RN wall normal acceleration (solid d2dt2 operator on Dt) at the BIPs
    void updateFsiWallAcceleration_();
};

} // namespace accel

#endif // MESHMOTION_H
