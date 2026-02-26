// File       : transitionOnsetReynoldsNumber.h
// Created    : Tue Jan 14 2025
// Author     : Adam Fares
// Description: Node scalar transport field for transition onset Reynolds number
// Copyright (c) 2025 CCFNUM, Lucerne University of Applied Sciences and Arts.
// SPDX-License-Identifier: BSD-3-Clause

#ifndef TRANSITIONONSETREYNOLDSNUMBER_H
#define TRANSITIONONSETREYNOLDSNUMBER_H

#include "nodeScalarField.h"

namespace accel
{

class transitionOnsetReynoldsNumber : public nodeScalarField
{
public:
    // Constructors

    transitionOnsetReynoldsNumber(realm* realmPtr,
                                  const std::string name,
                                  unsigned numberOfStates,
                                  bool highResolution);

    // Update

    void updateBoundarySideField(label iZone, label iBoundary) override;
};

} // namespace accel

#endif // TRANSITIONONSETREYNOLDSNUMBER_H
