// File       : youngModulus.cpp
// Created    : Thu Dec 04 2025 08:42:10 (+0100)
// Author     : Mhamad Mahdi Alloush
// Description:
// Copyright 2025 CCFNUM HSLU T&A. All Rights Reserved.

#include "youngModulus.h"
#include "realm.h"

namespace accel
{

youngModulus::youngModulus(realm* realmPtr,
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
