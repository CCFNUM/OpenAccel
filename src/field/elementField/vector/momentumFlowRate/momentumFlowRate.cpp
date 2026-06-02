// File       : momentumFlowRate.cpp
// Created    : Sun May 25 2026
// Author     : Mhamad Mahdi Alloush
// Description:
// Copyright 2026 CCFNUM HSLU T&A. All Rights Reserved.

#include "momentumFlowRate.h"
#include "domain.h"
#include "realm.h" // required for initialization

namespace accel
{

momentumFlowRate::momentumFlowRate(realm* realmPtr,
                                   const std::string name,
                                   unsigned numberOfStates)
    : elementField<scalar, SPATIAL_DIM>(realmPtr->meshPtr(),
                                        name,
                                        numberOfStates)
{
}

// Methods

void momentumFlowRate::putFieldOnRegisteredParts_()
{
    // skip ... no momentum flow rate data on interior ip's are of any meaning;
    // only the interface/boundary side field is allocated
}

} // namespace accel
