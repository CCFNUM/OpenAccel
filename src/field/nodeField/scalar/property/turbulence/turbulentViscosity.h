// File : turbulentViscosity.h
// Created : Tue Apr 20 2024 12:55:24 (+0100)
// Author : Mhamad Mahdi Alloush
// Description: Node scalar property field for turbulent eddy viscosity
// Copyright (c) 2024 CCFNUM, Lucerne University of Applied Sciences and Arts.
// SPDX-License-Identifier: BSD-3-Clause

#ifndef TURBULENTVISCOSITY_H
#define TURBULENTVISCOSITY_H

#include "property.h"

namespace accel
{

class turbulentViscosity : public property
{
public:
    // Constructors

    turbulentViscosity(realm* realmPtr,
                       const std::string name,
                       unsigned numberOfStates,
                       bool highResolution,
                       bool computeGradient);
};

} // namespace accel

#endif // TURBULENTVISCOSITY_H
