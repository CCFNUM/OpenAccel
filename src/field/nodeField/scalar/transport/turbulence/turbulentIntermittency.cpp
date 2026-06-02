// File       : turbulentIntermittency.cpp
// Created    : Tue Jan 14 2025
// Author     : Adam Fares
// Description:
// Copyright 2025 CCFNUM HSLU T&A. All Rights Reserved.

#include "turbulentIntermittency.h"
#include "controls.h"
#include "realm.h"
#include "simulation.h"

namespace accel
{

turbulentIntermittency::turbulentIntermittency(realm* realmPtr,
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
            .turbulentIntermittencyInterpolationType_;

    gradientInterpolationScheme_ =
        meshRef()
            .controlsRef()
            .solverRef()
            .solverControl_.basicSettings_.interpolationSchemeType_
            .turbulentIntermittencyGradientInterpolationType_;

    // turbulent dissipation rate may only apply to fluid domains
    mediumIndependent_ = false;
}

} // namespace accel
