// File : property.cpp
// Created : Tue Apr 20 2024 12:55:24 (+0100)
// Author : Mhamad Mahdi Alloush
// Description:
// Copyright (c) 2024 CCFNUM, Lucerne University of Applied Sciences and Arts.
// SPDX-License-Identifier: BSD-3-Clause

#include "property.h"
#include "realm.h"

namespace accel
{

property::property(realm* realmPtr,
                   std::string name,
                   unsigned numberOfStates,
                   bool prevIter,
                   bool highResolution,
                   bool computeGradient)
    : nodeScalarField(realmPtr,
                      name,
                      numberOfStates,
                      prevIter,
                      highResolution,
                      computeGradient,
                      false)
{
}

} // namespace accel
