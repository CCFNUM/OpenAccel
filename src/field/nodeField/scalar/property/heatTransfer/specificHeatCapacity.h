// File       : specificHeatCapacity.h
// Created    : Tue Apr 20 2024 12:55:24 (+0100)
// Author     : Mhamad Mahdi Alloush
// Description: Node scalar property field for specific heat capacity
// Copyright 2024 CCFNUM HSLU T&A. All Rights Reserved.

#ifndef SPECIFICHEATCAPACITY_H
#define SPECIFICHEATCAPACITY_H

#include "property.h"

namespace accel
{

class specificHeatCapacity : public property
{
public:
    // Constructors

    specificHeatCapacity(realm* realmPtr,
                         const std::string name,
                         unsigned numberOfStates,
                         bool highResolution,
                         bool computeGradient);
};

} // namespace accel

#endif // SPECIFICHEATCAPACITY_H
