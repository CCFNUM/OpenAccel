// File       : temperature.h
// Created    : Tue Apr 20 2024 12:55:24 (+0100)
// Author     : Mhamad Mahdi Alloush
// Description: Node scalar transport field for temperature
// Copyright (c) 2024 CCFNUM, Lucerne University of Applied Sciences and Arts.
// SPDX-License-Identifier: BSD-3-Clause

#ifndef TEMPERATURE_H
#define TEMPERATURE_H

#include "nodeScalarField.h"

namespace accel
{

class temperature : public nodeScalarField
{
public:
    // Constructors

    temperature(realm* realmPtr,
                const std::string name,
                unsigned numberOfStates,
                bool highResolution);

    // Update

    void updateBoundarySideField(label iZone, label iBoundary) override;
};

} // namespace accel

#endif // TEMPERATURE_H
