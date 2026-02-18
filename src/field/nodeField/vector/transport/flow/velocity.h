// File : velocity.h
// Created : Tue Apr 20 2024 12:55:24 (+0100)
// Author : Mhamad Mahdi Alloush
// Description: Velocity vector field with flow reversal detection and direction
// fields
// Copyright (c) 2024 CCFNUM, Lucerne University of Applied Sciences and
// Arts. SPDX-License-Identifier: BSD-3-Clause

#ifndef VELOCITY_H
#define VELOCITY_H

#include "nodeVectorField.h"

namespace accel
{

class velocity : public nodeVectorField
{
private:
    // acts as a label field to detect flow reversal
    std::unique_ptr<sideField<label, 1>> reversalFlagPtr_ = nullptr;

    std::unique_ptr<nodeSideField<scalar, SPATIAL_DIM>>
        nodeSideFlowDirectionFieldPtr_ = nullptr;

    std::unique_ptr<sideField<scalar, SPATIAL_DIM>> sideFlowDirectionFieldPtr_ =
        nullptr;

public:
    // Constructors

    velocity(realm* realmPtr,
             const std::string name,
             unsigned numberOfStates,
             bool highResolution);

    // Update

    void updateBoundarySideField(label iZone, label iBoundary) override;

    void updateBoundarySideFieldSpecifiedValue(label iZone,
                                               label iBoundary) override;

    void updateBoundarySideFieldNoSlipWall(label iZone, label iBoundary);

    void updateBoundarySideFieldNormalSpeed(label iZone, label iBoundary);

    void updateBoundarySideDirectionFields(label iZone, label iBoundary);

    void registerSideFlowDirectionFields(label iZone, label iBoundary);

    // Access

    sideField<label, 1>& reversalFlagRef();

    const sideField<label, 1>& reversalFlagRef() const;

    nodeSideField<scalar, SPATIAL_DIM>& nodeSideFlowDirectionFieldRef();

    const nodeSideField<scalar, SPATIAL_DIM>&
    nodeSideFlowDirectionFieldRef() const;

    sideField<scalar, SPATIAL_DIM>& sideFlowDirectionFieldRef();

    const sideField<scalar, SPATIAL_DIM>& sideFlowDirectionFieldRef() const;
};

} // namespace accel

#endif // VELOCITY_H
