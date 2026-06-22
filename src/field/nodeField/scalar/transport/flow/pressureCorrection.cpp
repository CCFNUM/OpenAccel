// File       : pressureCorrection.cpp
// Created    : Mon May 26 2026
// Author     : Mhamad Mahdi Alloush
// Description:
// Copyright 2026 CCFNUM HSLU T&A. All Rights Reserved.

#include "pressureCorrection.h"
#include "controls.h"
#include "realm.h"

namespace accel
{

pressureCorrection::pressureCorrection(realm* realmPtr,
                                       const std::string name,
                                       unsigned numberOfStates,
                                       bool highResolution)
    : nodeScalarField(realmPtr,
                      name,
                      numberOfStates,
                      /* prevIter                   */ false,
                      /* highResolution             */ false,
                      /* computeGradient            */ true,
                      /* correctedBoundaryNodeValues */ false)
{
    // Mirror the pressure field's interpolation so that grad(p') is computed
    // with the very same Green-Gauss discretization as grad(p). This keeps the
    // segregated velocity correction U -= D*grad(p') consistent with the
    // pressure gradient used as the momentum source.
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

    // force correct gradient: remove symmetric contributions to the gradient
    correctGradient_ = true;
}

} // namespace accel
