// File : pressure.cpp
// Created : Tue Apr 20 2024 12:55:24 (+0100)
// Author : Mhamad Mahdi Alloush
// Description:
// Copyright (c) 2024 CCFNUM, Lucerne University of Applied Sciences and Arts.
// SPDX-License-Identifier: BSD-3-Clause

#include "pressure.h"
#include "controls.h"
#include "realm.h"

namespace accel
{

pressure::pressure(realm* realmPtr,
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
            .pressureInterpolationType_;

    gradientInterpolationScheme_ =
        meshRef()
            .controlsRef()
            .solverRef()
            .solverControl_.basicSettings_.interpolationSchemeType_
            .pressureGradientInterpolationType_;

    // pressure may only apply to fluid domains
    mediumIndependent_ = false;

    // force correct gradient for pressure: remove symmetric contributions to
    // gradient
    correctGradient_ = true;
}

void pressure::updateBoundarySideField(label iZone, label iBoundary)
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
                    case boundaryConditionType::staticPressure:
                        {
                            nodeScalarField::
                                updateBoundarySideFieldSpecifiedValue(
                                    iZone, iBoundary);
                        }
                        break;

                    case boundaryConditionType::totalPressure:
                        {
                            // model-related: done in NS model

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

        case boundaryPhysicalType::outlet:
            {
                switch (bcType)
                {
                    case boundaryConditionType::staticPressure:
                        {
                            nodeScalarField::
                                updateBoundarySideFieldSpecifiedValue(
                                    iZone, iBoundary);
                        }
                        break;

                    case boundaryConditionType::averageStaticPressure:
                        {
                            // model-related
                        }
                        break;

                    case boundaryConditionType::massFlowRate:
                        {
                            // model-related
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
