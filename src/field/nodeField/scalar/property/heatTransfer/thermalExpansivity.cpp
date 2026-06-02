// File       : thermalExpansivity.cpp
// Created    : Mon Apr 21 2025 09:18:24 (+0100)
// Author     : Mhamad Mahdi Alloush
// Description:
// Copyright 2025 CCFNUM HSLU T&A. All Rights Reserved.

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
