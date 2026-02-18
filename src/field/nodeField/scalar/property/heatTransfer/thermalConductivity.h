// File : thermalConductivity.h
// Created : Tue Apr 20 2024 12:55:24 (+0100)
// Author : Mhamad Mahdi Alloush
// Description: Material property field for thermal conductivity
// Copyright (c) 2024 CCFNUM, Lucerne University of Applied Sciences and Arts.
// SPDX-License-Identifier: BSD-3-Clause

#ifndef THERMALCONDUCTIVITY_H
#define THERMALCONDUCTIVITY_H

#include "property.h"

namespace accel
{

class thermalConductivity : public property
{
public:
    // Constructors

    thermalConductivity(realm* realmPtr,
                        const std::string name,
                        unsigned numberOfStates,
                        bool highResolution,
                        bool computeGradient);
};

} // namespace accel

#endif // THERMALCONDUCTIVITY_H
