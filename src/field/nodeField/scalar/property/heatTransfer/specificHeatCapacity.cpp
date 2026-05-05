// File       : specificHeatCapacity.cpp
// Created    : Tue Apr 20 2024 12:55:24 (+0100)
// Author     : Mhamad Mahdi Alloush
// Description:
// Copyright 2024 CCFNUM HSLU T&A. All Rights Reserved.

#include "specificHeatCapacity.h"
#include "realm.h"

namespace accel
{

specificHeatCapacity::specificHeatCapacity(realm* realmPtr,
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
