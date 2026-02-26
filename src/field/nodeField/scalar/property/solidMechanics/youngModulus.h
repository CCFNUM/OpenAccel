// File       : youngModulus.h
// Created    : Thu Dec 04 2025 08:42:10 (+0100)
// Author     : Mhamad Mahdi Alloush
// Description: Node scalar property field for Young's modulus in solid
// mechanics
// Copyright (c) 2025 CCFNUM, Lucerne University of Applied Sciences
// and Arts.
// SPDX-License-Identifier: BSD-3-Clause

#ifndef YOUNGMODULUS_H
#define YOUNGMODULUS_H

#include "property.h"

namespace accel
{

class youngModulus : public property
{
public:
    // Constructors

    youngModulus(realm* realmPtr,
                 const std::string name,
                 unsigned numberOfStates,
                 bool highResolution,
                 bool computeGradient);
};

} // namespace accel

#endif // YOUNGMODULUS_H
