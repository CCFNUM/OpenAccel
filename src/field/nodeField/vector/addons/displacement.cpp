// File : displacement.cpp
// Created : Sat Dec 06 2025 10:22:38 (+0100)
// Author : Mhamad Mahdi Alloush
// Description:
// Copyright (c) 2025 CCFNUM, Lucerne University of Applied Sciences and Arts.
// SPDX-License-Identifier: BSD-3-Clause

#include "displacement.h"
#include "controls.h"
#include "realm.h"
#include "simulation.h"

namespace accel
{

displacement::displacement(realm* realmPtr,
                           const std::string name,
                           unsigned numberOfStates)
    : nodeVectorField(realmPtr,
                      name,
                      realmPtr->meshPtr()->controlsRef().isTransient() &&
                              numberOfStates == 2
                          ? numberOfStates + 1 // force 1 more state
                          : numberOfStates,
                      true,
                      false,
                      true,
                      true)
{
    realmPtr->registerRestartField(name);

    interpolationScheme_ =
        meshRef()
            .controlsRef()
            .solverRef()
            .solverControl_.basicSettings_.interpolationSchemeType_
            .displacementInterpolationType_;

    gradientInterpolationScheme_ =
        meshRef()
            .controlsRef()
            .solverRef()
            .solverControl_.basicSettings_.interpolationSchemeType_
            .displacementGradientInterpolationType_;
}

void displacement::updateBoundarySideField(label iZone, label iBoundary)
{
    assert(this->isZoneSet(iZone));

    boundaryPhysicalType physicalType =
        this->meshPtr()->zonePtr(iZone)->boundaryRef(iBoundary).type();
    const auto& bcType = this->boundaryConditionRef(iZone, iBoundary).type();

    switch (physicalType)
    {
        case boundaryPhysicalType::wall:
            {
                switch (bcType)
                {
                    case boundaryConditionType::periodicDisplacement:
                    case boundaryConditionType::rigidBodySolution:
                    case boundaryConditionType::specifiedFlux:
                        {
                            // model-related
                        }
                        break;

                    default:
                        nodeVectorField::updateBoundarySideField(iZone,
                                                                 iBoundary);
                        break;
                }
            }
            break;

        default:
            {
                switch (bcType)
                {
                    case boundaryConditionType::periodicDisplacement:
                    case boundaryConditionType::rigidBodySolution:
                        {
                            // model-related
                        }
                        break;

                    default:
                        nodeVectorField::updateBoundarySideField(iZone,
                                                                 iBoundary);
                        break;
                }
            }
            break;
    }
}

} // namespace accel
