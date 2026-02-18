// File : heatFlowRate.cpp
// Created : Sat Jul 05 2025 11:10:00 (+0100)
// Author : Mhamad Mahdi Alloush
// Description:
// Copyright (c) 2025 CCFNUM, Lucerne University of Applied Sciences and Arts.
// SPDX-License-Identifier: BSD-3-Clause

#include "heatFlowRate.h"
#include "boundary.h"
#include "domain.h"
#include "mesh.h"
#include "realm.h" // required for initialization
#include "simulation.h"
#include "zone.h"

namespace accel
{

heatFlowRate::heatFlowRate(realm* realmPtr,
                           const std::string name,
                           unsigned numberOfStates)
    : elementScalarField(realmPtr, name, numberOfStates)
{
}

// Methods

void heatFlowRate::putFieldOnRegisteredParts_()
{
    // skip ... no heat flow rate data on interior ip's are of any meaning; let
    // us not allocate any memory for that
}

} // namespace accel
