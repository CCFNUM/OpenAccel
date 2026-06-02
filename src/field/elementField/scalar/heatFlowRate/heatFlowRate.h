// File       : heatFlowRate.h
// Created    : Sat Jul 05 2025 11:10:00 (+0100)
// Author     : Mhamad Mahdi Alloush
// Description: Element scalar field for heat flow rate on mesh parts
// Copyright 2025 CCFNUM HSLU T&A. All Rights Reserved.

#ifndef HEATFLOWRATE_H
#define HEATFLOWRATE_H

#include "elementScalarField.h"

namespace accel
{

class heatFlowRate : public elementScalarField
{
public:
    // Constructors

    heatFlowRate(realm* realmPtr,
                 const std::string name,
                 unsigned numberOfStates);

    // Methods

    void putFieldOnRegisteredParts_() override;
};

} // namespace accel

#endif // HEATFLOWRATE_H
