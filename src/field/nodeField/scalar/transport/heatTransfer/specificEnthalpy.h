// File : specificEnthalpy.h
// Created : Tue Apr 20 2024 12:55:24 (+0100)
// Author : Mhamad Mahdi Alloush
// Description: Node scalar transport field for specific enthalpy
// Copyright (c) 2024 CCFNUM, Lucerne University of Applied Sciences and Arts.
// SPDX-License-Identifier: BSD-3-Clause

#ifndef SPECIFICENTHALPY_H
#define SPECIFICENTHALPY_H

#include "nodeScalarField.h"

namespace accel
{

class specificEnthalpy : public nodeScalarField
{
public:
    // Constructors

    specificEnthalpy(realm* realmPtr,
                     const std::string name,
                     unsigned numberOfStates,
                     bool highResolution);

    // Update

    void updateBoundarySideField(label iZone, label iBoundary) override;
};

} // namespace accel

#endif // SPECIFICENTHALPY_H
