// File       : property.cpp
// Created    : Tue Apr 20 2024 12:55:24 (+0100)
// Author     : Mhamad Mahdi Alloush
// Description:
// Copyright 2024 CCFNUM HSLU T&A. All Rights Reserved.

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
