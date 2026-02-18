// File : property.h
// Created : Tue Apr 20 2024 12:55:24 (+0100)
// Author : Mhamad Mahdi Alloush
// Description: Base class for material property node scalar fields
// Copyright (c) 2024 CCFNUM, Lucerne University of Applied Sciences and Arts.
// SPDX-License-Identifier: BSD-3-Clause

#ifndef PROPERTY_H
#define PROPERTY_H

#include "nodeScalarField.h"

namespace accel
{

class property : public nodeScalarField
{
public:
    // Constructors

    property(realm* realmPtr,
             std::string name,
             unsigned numberOfStates,
             bool prevIter,
             bool highResolution,
             bool computeGradient);
};

} // namespace accel

#endif // PROPERTY_H
