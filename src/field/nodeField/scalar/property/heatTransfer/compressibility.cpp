// File       : compressibility.cpp
// Created    : Fri Apr 11 2025 18:10:24 (+0100)
// Author     : Mhamad Mahdi Alloush
// Description:
// Copyright 2025 CCFNUM HSLU T&A. All Rights Reserved.

#include "compressibility.h"
#include "realm.h"

namespace accel
{

compressibility::compressibility(realm* realmPtr,
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
