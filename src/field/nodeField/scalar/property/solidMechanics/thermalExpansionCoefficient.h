// File       : thermalExpansionCoefficient.h
// Created    : Sun Nov 08 2026 00:18:57 (+0200)
// Author     : Mohammad Nasreddine
// Description: Node scalar property field for the linear thermal expansion
//              coefficient in solid mechanics (thermal stress)
// Copyright 2026 CCFNUM HSLU T&A. All Rights Reserved.

#ifndef THERMALEXPANSIONCOEFFICIENT_H
#define THERMALEXPANSIONCOEFFICIENT_H

#include "property.h"

namespace accel
{

class thermalExpansionCoefficient : public property
{
public:
    // Constructors

    thermalExpansionCoefficient(realm* realmPtr,
                                const std::string name,
                                unsigned numberOfStates,
                                bool highResolution,
                                bool computeGradient);
};

} // namespace accel

#endif // THERMALEXPANSIONCOEFFICIENT_H
