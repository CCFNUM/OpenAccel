// File : turbulentDissipationRate.h
// Created : Tue Apr 20 2024 12:55:24 (+0100)
// Author : Achraf Nagihi
// Description: Node scalar transport field for turbulent dissipation rate
// (epsilon)
// Copyright (c) 2024 CCFNUM, Lucerne University of Applied Sciences
// and Arts. SPDX-License-Identifier: BSD-3-Clause

#ifndef TURBULENTDISSIPATIONRATE_H
#define TURBULENTDISSIPATIONRATE_H

#include "nodeScalarField.h"

namespace accel
{

class turbulentDissipationRate : public nodeScalarField
{
public:
    // Constructors

    turbulentDissipationRate(realm* realmPtr,
                             const std::string name,
                             unsigned numberOfStates,
                             bool highResolution);

    // Update

    void updateBoundarySideField(label iZone, label iBoundary) override;
};

} // namespace accel

#endif // TURBULENTDISSIPATIONRATE_H
