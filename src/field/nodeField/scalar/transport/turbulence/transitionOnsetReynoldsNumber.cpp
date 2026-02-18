// File : transitionOnsetReynoldsNumber.cpp
// Created : Tue Jan 14 2025
// Author : Adam Fares
// Description:
// Copyright (c) 2025 CCFNUM, Lucerne University of Applied Sciences and Arts.
// SPDX-License-Identifier: BSD-3-Clause

#include "transitionOnsetReynoldsNumber.h"
#include "controls.h"
#include "realm.h"
#include "simulation.h"

namespace accel
{

transitionOnsetReynoldsNumber::transitionOnsetReynoldsNumber(
    realm* realmPtr,
    const std::string name,
    unsigned numberOfStates,
    bool highResolution)
    : nodeScalarField(realmPtr,
                      name,
                      numberOfStates,
                      true,
                      highResolution,
                      true,
                      false)
{
    realmPtr->registerRestartField(name);

    interpolationScheme_ =
        meshRef()
            .controlsRef()
            .solverRef()
            .solverControl_.basicSettings_.interpolationSchemeType_
            .transitionOnsetReynoldsNumberInterpolationType_;

    gradientInterpolationScheme_ =
        meshRef()
            .controlsRef()
            .solverRef()
            .solverControl_.basicSettings_.interpolationSchemeType_
            .transitionOnsetReynoldsNumberGradientInterpolationType_;

    // turbulent dissipation rate may only apply to fluid domains
    mediumIndependent_ = false;
}

void transitionOnsetReynoldsNumber::updateBoundarySideField(label iZone,
                                                            label iBoundary)
{
    assert(this->isZoneSet(iZone));

    boundaryPhysicalType physicalType =
        this->meshPtr()->zonePtr(iZone)->boundaryRef(iBoundary).type();
    const auto& bcType = this->boundaryConditionRef(iZone, iBoundary).type();

    switch (physicalType)
    {
        case boundaryPhysicalType::inlet:
        case boundaryPhysicalType::opening:
            {
                // model-based update
            }
            break;

        default:
            nodeScalarField::updateBoundarySideField(iZone, iBoundary);
            break;
    }
}

} // namespace accel
