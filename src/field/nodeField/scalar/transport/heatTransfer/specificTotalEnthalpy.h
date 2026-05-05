// File       : specificTotalEnthalpy.h
// Created    : Tue Apr 20 2024 12:55:24 (+0100)
// Author     : Mhamad Mahdi Alloush
// Description: Node scalar transport field for specific total enthalpy
// Copyright 2024 CCFNUM HSLU T&A. All Rights Reserved.

#ifndef SPECIFICTOTALENTHALPY_H
#define SPECIFICTOTALENTHALPY_H

#include "nodeScalarField.h"

namespace accel
{

class specificTotalEnthalpy : public nodeScalarField
{
public:
    // Constructors

    specificTotalEnthalpy(realm* realmPtr,
                          const std::string name,
                          unsigned numberOfStates,
                          bool highResolution);

    // Update

    void updateBoundarySideField(label iZone, label iBoundary) override;
};

} // namespace accel

#endif // SPECIFICTOTALENTHALPY_H
