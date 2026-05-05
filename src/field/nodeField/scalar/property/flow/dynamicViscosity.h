// File       : dynamicViscosity.h
// Created    : Tue Apr 20 2024 12:55:24 (+0100)
// Author     : Mhamad Mahdi Alloush
// Description: Node scalar property field for fluid dynamic viscosity
// Copyright 2024 CCFNUM HSLU T&A. All Rights Reserved.

#ifndef DYNAMICVISCOSITY_H
#define DYNAMICVISCOSITY_H

#include "property.h"

namespace accel
{

class dynamicViscosity : public property
{
public:
    // Constructors

    dynamicViscosity(realm* realmPtr,
                     const std::string name,
                     unsigned numberOfStates,
                     bool highResolution,
                     bool computeGradient);
};

} // namespace accel

#endif // DYNAMICVISCOSITY_H
