// File       : thermalExpansionCoefficient.cpp
// Created    : Sun Nov 08 2026 00:18:57 (+0200)
// Author     : Mohammad Nasreddine
// Description:
// Copyright 2026 CCFNUM HSLU T&A. All Rights Reserved.

#include "thermalExpansionCoefficient.h"
#include "realm.h"

namespace accel
{

thermalExpansionCoefficient::thermalExpansionCoefficient(
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
