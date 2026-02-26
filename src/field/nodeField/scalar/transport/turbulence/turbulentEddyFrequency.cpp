// File       : turbulentEddyFrequency.cpp
// Created    : Tue Apr 20 2024 12:55:24 (+0100)
// Author     : Mhamad Mahdi Alloush
// Description:
// Copyright (c) 2024 CCFNUM, Lucerne University of Applied Sciences and Arts.
// SPDX-License-Identifier: BSD-3-Clause

#include "turbulentEddyFrequency.h"
#include "controls.h"
#include "realm.h"
#include "simulation.h"

namespace accel
{

turbulentEddyFrequency::turbulentEddyFrequency(realm* realmPtr,
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
            .turbulentEddyFrequencyInterpolationType_;

    gradientInterpolationScheme_ =
        meshRef()
            .controlsRef()
            .solverRef()
            .solverControl_.basicSettings_.interpolationSchemeType_
            .turbulentEddyFrequencyGradientInterpolationType_;

    // turbulent eddy frequency may only apply to fluid domains
    mediumIndependent_ = false;
}

void turbulentEddyFrequency::updateBoundarySideField(label iZone,
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
