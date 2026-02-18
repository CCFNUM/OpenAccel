// File : turbulentThermalConductivity.cpp
// Created : Tue Apr 20 2024 12:55:24 (+0100)
// Author : Mhamad Mahdi Alloush
// Description:
// Copyright (c) 2024 CCFNUM, Lucerne University of Applied Sciences and Arts.
// SPDX-License-Identifier: BSD-3-Clause

#include "turbulentThermalConductivity.h"
#include "realm.h"

namespace accel
{

turbulentThermalConductivity::turbulentThermalConductivity(
    realm* realmPtr,
    const std::string name,
    unsigned numberOfStates,
    bool highResolution,
    bool computeGradient)
    : property(realmPtr,
               name,
               numberOfStates,
               true,
               highResolution,
               computeGradient)
{
    realmPtr->registerRestartField(name);
}

} // namespace accel
