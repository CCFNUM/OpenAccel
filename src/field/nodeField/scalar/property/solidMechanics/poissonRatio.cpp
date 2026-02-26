// File       : poissonRatio.cpp
// Created    : Thu Dec 04 2025 08:42:10 (+0100)
// Author     : Mhamad Mahdi Alloush
// Description:
// Copyright (c) 2025 CCFNUM, Lucerne University of Applied Sciences and Arts.
// SPDX-License-Identifier: BSD-3-Clause

#include "poissonRatio.h"
#include "realm.h"

namespace accel
{

poissonRatio::poissonRatio(realm* realmPtr,
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
