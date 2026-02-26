// File       : density.h
// Created    : Tue Apr 20 2024 12:55:24 (+0100)
// Author     : Mhamad Mahdi Alloush
// Description: Node scalar property field for fluid density
// Copyright (c) 2024 CCFNUM, Lucerne University of Applied Sciences and Arts.
// SPDX-License-Identifier: BSD-3-Clause

#ifndef DENSITY_H
#define DENSITY_H

#include "property.h"

#include "types.h"
#include <cassert>
#include <string>
#include <vector>

namespace accel
{

class density : public property
{
public:
    // Constructors

    density(realm* realmPtr,
            const std::string name,
            unsigned numberOfStates,
            bool highResolution,
            bool computeGradient);
};

} // namespace accel

#endif // DENSITY_H
