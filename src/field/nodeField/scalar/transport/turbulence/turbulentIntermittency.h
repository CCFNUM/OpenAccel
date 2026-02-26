// File       : turbulentIntermittency.h
// Created    : Tue Jan 14 2025
// Author     : Adam Fares
// Description: Node scalar transport field for turbulent intermittency (gamma)
// Copyright (c) 2025 CCFNUM, Lucerne University of Applied Sciences and Arts.
// SPDX-License-Identifier: BSD-3-Clause

#ifndef TURBULENTINTERMITTENCY_H
#define TURBULENTINTERMITTENCY_H

#include "nodeScalarField.h"

namespace accel
{

class turbulentIntermittency : public nodeScalarField
{
public:
    // Constructors

    turbulentIntermittency(realm* realmPtr,
                           const std::string name,
                           unsigned numberOfStates,
                           bool highResolution);
};

} // namespace accel

#endif // TURBULENTINTERMITTENCY_H
