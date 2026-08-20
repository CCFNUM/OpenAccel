// File       : fluidEquations.h
// Created    : Wed Feb 21 2024 12:34:15 (+0100)
// Author     : Fabian Wermelinger
// Description: Convenience header that includes all fluid physical equations
// Copyright 2024 CCFNUM HSLU T&A. All Rights Reserved.

#ifndef FLUIDEQUATIONSINCLUDES_H
#define FLUIDEQUATIONSINCLUDES_H
// clang-format off

// scalar equations
#include "pressureCorrectionEquation.h"
#ifdef WITH_THERMAL_TEMPERATURE
#include "thermalTemperatureEquation.h"
#else
#include "thermalEnergyEquation.h"
#endif
#include "totalEnergyEquation.h"
#include "solidDisplacementEquation.h"

// vector equations
#include "navierStokesEquation.h" // coupled momentum

// higher-level equations
#include "segregatedFlowEquations.h" // Navier-Stokes (segregated: coupled momentum, scalar pressure correction)

// turbulence
#include "segregatedShearStressTransportEquations.h"
#include "segregatedKEpsilonEquations.h"
#include "segregatedTransitionShearStressTransportEquations.h"
#include "segregatedCorrelationTransitionShearStressTransportEquations.h"

// multiphase
#include "segregatedEulerEulerFlowEquations.h"
#include "segregatedFreeSurfaceFlowEquations.h"

// clang-format on
#endif /* FLUIDEQUATIONSINCLUDES_H */
