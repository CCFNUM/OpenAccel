// File : turbulentThermalConductivity.h
// Created : Tue Apr 20 2024 12:55:24 (+0100)
// Author : Mhamad Mahdi Alloush
// Description: Node scalar property field for turbulent thermal conductivity
// Copyright (c) 2024 CCFNUM, Lucerne University of Applied Sciences and Arts.
// SPDX-License-Identifier: BSD-3-Clause

#ifndef TURBULENTTHERMALCONDUCTIVITY_H
#define TURBULENTTHERMALCONDUCTIVITY_H

#include "property.h"

namespace accel
{

class turbulentThermalConductivity : public property
{
public:
    // Constructors

    turbulentThermalConductivity(realm* realmPtr,
                                 const std::string name,
                                 unsigned numberOfStates,
                                 bool highResolution,
                                 bool computeGradient);
};

} // namespace accel

#endif // TURBULENTTHERMALCONDUCTIVITY_H
