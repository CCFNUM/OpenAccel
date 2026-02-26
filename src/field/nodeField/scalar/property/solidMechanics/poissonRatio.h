// File       : poissonRatio.h
// Created    : Thu Dec 04 2025 08:42:10 (+0100)
// Author     : Mhamad Mahdi Alloush
// Description: Node scalar property field for Poisson's ratio in solid
// mechanics
// Copyright (c) 2025 CCFNUM, Lucerne University of Applied Sciences
// and Arts.
// SPDX-License-Identifier: BSD-3-Clause

#ifndef POISSONRATIO_H
#define POISSONRATIO_H

#include "property.h"

namespace accel
{

class poissonRatio : public property
{
public:
    // Constructors

    poissonRatio(realm* realmPtr,
                 const std::string name,
                 unsigned numberOfStates,
                 bool highResolution,
                 bool computeGradient);
};

} // namespace accel

#endif // POISSONRATIO_H
