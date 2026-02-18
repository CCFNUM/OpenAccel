// File : thermalExpansivity.h
// Created : Mon Apr 21 2025 09:18:24 (+0100)
// Author : Mhamad Mahdi Alloush
// Description: Material property field for thermal expansivity
// Copyright (c) 2025 CCFNUM, Lucerne University of Applied Sciences and Arts.
// SPDX-License-Identifier: BSD-3-Clause

#ifndef THERMALEXPANSIVITY_H
#define THERMALEXPANSIVITY_H

#include "property.h"

namespace accel
{

class thermalExpansivity : public property
{
public:
    // Constructors

    thermalExpansivity(realm* realmPtr,
                       const std::string name,
                       unsigned numberOfStates,
                       bool highResolution,
                       bool computeGradient);
};

} // namespace accel

#endif // THERMALEXPANSIVITY_H
