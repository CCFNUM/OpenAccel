// File       : interfaceScalarField.cpp
// Created    : Wed Oct 07 2026 13:34:00 (+0300)
// Author     : Mohammad Nasreddine
// Description: Element scalar field that lives on interface/boundary sides
//              only (no interior allocation)
// Copyright 2026 CCFNUM HSLU T&A. All Rights Reserved.

#include "interfaceScalarField.h"
#include "boundary.h"
#include "domain.h"
#include "mesh.h"
#include "realm.h" // required for initialization
#include "simulation.h"
#include "zone.h"

namespace accel
{

interfaceScalarField::interfaceScalarField(realm* realmPtr,
                                           const std::string name,
                                           unsigned numberOfStates)
    : elementScalarField(realmPtr, name, numberOfStates)
{
}

// Methods

void interfaceScalarField::putFieldOnRegisteredParts_()
{
    // interface-side field only: no interior allocation
}

} // namespace accel
