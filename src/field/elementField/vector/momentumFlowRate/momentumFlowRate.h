// File       : momentumFlowRate.h
// Created    : Sun May 25 2026
// Author     : Mhamad Mahdi Alloush
// Description: Element vector field for momentum flow rate on mesh parts
//              (interface/boundary side field; vector analogue of heatFlowRate)
// Copyright 2026 CCFNUM HSLU T&A. All Rights Reserved.

#ifndef MOMENTUMFLOWRATE_H
#define MOMENTUMFLOWRATE_H

#include "elementField.h"

namespace accel
{

class realm;

class momentumFlowRate : public elementField<scalar, SPATIAL_DIM>
{
public:
    // Constructors

    momentumFlowRate(realm* realmPtr,
                     const std::string name,
                     unsigned numberOfStates);

    // Methods

    void putFieldOnRegisteredParts_() override;
};

} // namespace accel

#endif // MOMENTUMFLOWRATE_H
