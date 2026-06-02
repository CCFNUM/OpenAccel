// File       : turbulentDissipationRate.h
// Created    : Tue Apr 20 2024 12:55:24 (+0100)
// Author     : Achraf Nagihi
// Description: Node scalar transport field for turbulent dissipation rate
//              (epsilon)
// Copyright 2024 CCFNUM HSLU T&A. All Rights Reserved.

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
