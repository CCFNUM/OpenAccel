// File       : turbulentEddyFrequency.h
// Created    : Tue Apr 20 2024 12:55:24 (+0100)
// Author     : Mhamad Mahdi Alloush
// Description: Node scalar transport field for turbulent eddy frequency (omega)
// Copyright 2024 CCFNUM HSLU T&A. All Rights Reserved.

#ifndef TURBULENTEDDYFREQUENCY_H
#define TURBULENTEDDYFREQUENCY_H

#include "nodeScalarField.h"

namespace accel
{

class turbulentEddyFrequency : public nodeScalarField
{
public:
    // Constructors

    turbulentEddyFrequency(realm* realmPtr,
                           const std::string name,
                           unsigned numberOfStates,
                           bool highResolution);

    // Update

    void updateBoundarySideField(label iZone, label iBoundary) override;
};

} // namespace accel

#endif // TURBULENTEDDYFREQUENCY_H
