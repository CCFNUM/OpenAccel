// File : compressibility.cpp
// Created : Fri Apr 11 2025 18:10:24 (+0100)
// Author : Mhamad Mahdi Alloush
// Description:
// Copyright (c) 2025 CCFNUM, Lucerne University of Applied Sciences and Arts.
// SPDX-License-Identifier: BSD-3-Clause

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
