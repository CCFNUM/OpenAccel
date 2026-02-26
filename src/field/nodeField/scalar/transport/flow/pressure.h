// File       : pressure.h
// Created    : Tue Apr 20 2024 12:55:24 (+0100)
// Author     : Mhamad Mahdi Alloush
// Description: Node scalar transport field for pressure
// Copyright (c) 2024 CCFNUM, Lucerne University of Applied Sciences and Arts.
// SPDX-License-Identifier: BSD-3-Clause

#ifndef PRESSURE_H
#define PRESSURE_H

#include "nodeScalarField.h"

namespace accel
{

class pressure : public nodeScalarField
{
public:
    // Constructors

    pressure(realm* realmPtr,
             const std::string name,
             unsigned numberOfStates,
             bool highResolution);

    // Update

    void updateBoundarySideField(label iZone, label iBoundary) override;
};

} // namespace accel

#endif // PRESSURE_H
