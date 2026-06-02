// File       : turbulentKineticEnergy.h
// Created    : Tue Apr 20 2024 12:55:24 (+0100)
// Author     : Mhamad Mahdi Alloush
// Description: Node scalar transport field for turbulent kinetic energy (k)
// Copyright 2024 CCFNUM HSLU T&A. All Rights Reserved.

#ifndef TURBULENTKINETICENERGY_H
#define TURBULENTKINETICENERGY_H

#include "nodeScalarField.h"

namespace accel
{

class turbulentKineticEnergy : public nodeScalarField
{
public:
    // Constructors

    turbulentKineticEnergy(realm* realmPtr,
                           const std::string name,
                           unsigned numberOfStates,
                           bool highResolution);

    // Update

    void updateBoundarySideField(label iZone, label iBoundary) override;
};

} // namespace accel

#endif // TURBULENTKINETICENERGY_H
