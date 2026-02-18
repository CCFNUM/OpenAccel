// File : turbulentDissipationRate.cpp
// Created : Tue Apr 20 2024 12:55:24 (+0100)
// Author : Achraf Nagihi
// Description:
// Copyright (c) 2024 CCFNUM, Lucerne University of Applied Sciences and Arts.
// SPDX-License-Identifier: BSD-3-Clause

#include "turbulentDissipationRate.h"
#include "controls.h"
#include "realm.h"
#include "simulation.h"

namespace accel
{

turbulentDissipationRate::turbulentDissipationRate(realm* realmPtr,
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
            .turbulentDissipationRateInterpolationType_;

    gradientInterpolationScheme_ =
        meshRef()
            .controlsRef()
            .solverRef()
            .solverControl_.basicSettings_.interpolationSchemeType_
            .turbulentDissipationRateGradientInterpolationType_;

    // turbulent dissipation rate may only apply to fluid domains
    mediumIndependent_ = false;
}

void turbulentDissipationRate::updateBoundarySideField(label iZone,
                                                       label iBoundary)
{
    assert(this->isZoneSet(iZone));

    boundaryPhysicalType physicalType =
        this->meshPtr()->zonePtr(iZone)->boundaryRef(iBoundary).type();
    const auto& bcType = this->boundaryConditionRef(iZone, iBoundary).type();

    switch (physicalType)
    {
        case boundaryPhysicalType::inlet:
            {
                switch (bcType)
                {
                    case boundaryConditionType::intensityAndLengthScale:
                    case boundaryConditionType::intensityAndEddyViscosityRatio:
                        {
                            // model-related: done in RANS model

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
