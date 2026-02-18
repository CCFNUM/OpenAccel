// File : thermalExpansivity.cpp
// Created : Mon Apr 21 2025 09:18:24 (+0100)
// Author : Mhamad Mahdi Alloush
// Description:
// Copyright (c) 2025 CCFNUM, Lucerne University of Applied Sciences and Arts.
// SPDX-License-Identifier: BSD-3-Clause

#include "thermalExpansivity.h"
#include "realm.h"

namespace accel
{

thermalExpansivity::thermalExpansivity(realm* realmPtr,
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
