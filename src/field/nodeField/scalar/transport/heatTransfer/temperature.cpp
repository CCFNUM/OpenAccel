// File       : temperature.cpp
// Created    : Tue Apr 20 2024 12:55:24 (+0100)
// Author     : Mhamad Mahdi Alloush
// Description:
// Copyright (c) 2024 CCFNUM, Lucerne University of Applied Sciences and Arts.
// SPDX-License-Identifier: BSD-3-Clause

#include "temperature.h"
#include "controls.h"
#include "realm.h"

namespace accel
{

temperature::temperature(realm* realmPtr,
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
            .temperatureInterpolationType_;

    gradientInterpolationScheme_ =
        meshRef()
            .controlsRef()
            .solverRef()
            .solverControl_.basicSettings_.interpolationSchemeType_
            .temperatureGradientInterpolationType_;
}

void temperature::updateBoundarySideField(label iZone, label iBoundary)
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
                switch (bcType)
                {
                    case boundaryConditionType::staticTemperature:
                        {
                            nodeScalarField::
                                updateBoundarySideFieldSpecifiedValue(
                                    iZone, iBoundary);
                        }
                        break;

                    case boundaryConditionType::totalTemperature:
                        {
                            // model-related: done in HT model

                            assert(nodeSideFieldPtr_);
                            assert(sideFieldPtr_);
                        }
                        break;

                    default:
                        nodeScalarField::updateBoundarySideField(iZone,
                                                                 iBoundary);
                        break;
                }
            }
            break;

        default:
            nodeScalarField::updateBoundarySideField(iZone, iBoundary);
            break;
    }
}

} // namespace accel
