// File : specificHeatCapacity.h
// Created : Tue Apr 20 2024 12:55:24 (+0100)
// Author : Mhamad Mahdi Alloush
// Description: Node scalar property field for specific heat capacity
// Copyright (c) 2024 CCFNUM, Lucerne University of Applied Sciences and Arts.
// SPDX-License-Identifier: BSD-3-Clause

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
