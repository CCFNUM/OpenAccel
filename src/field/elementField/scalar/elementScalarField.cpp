// File       : elementScalarField.cpp
// Created    : Fri Aug 25 2023 12:55:24 (+0100)
// Author     : Mhamad Mahdi Alloush
// Description:
// Copyright (c) 2023 CCFNUM, Lucerne University of Applied Sciences and Arts.
// SPDX-License-Identifier: BSD-3-Clause

#include "elementScalarField.h"
#include "boundary.h"
#include "domain.h"
#include "realm.h"
#include "simulation.h"
#include "zone.h"

namespace accel
{

elementScalarField::elementScalarField(realm* realmPtr,
                                       std::string name,
                                       unsigned numberOfStates)
    : elementField<scalar, 1>(realmPtr->meshPtr(), name, numberOfStates),
      realmPtr_(realmPtr)
{
}

// Access

simulation* elementScalarField::simulationPtr()
{
    return realmRef().simulationPtr();
}

const simulation* elementScalarField::simulationPtr() const
{
    return realmRef().simulationPtr();
}

simulation& elementScalarField::simulationRef()
{
    return realmRef().simulationRef();
}

const simulation& elementScalarField::simulationRef() const
{
    return realmRef().simulationRef();
}

realm* elementScalarField::realmPtr()
{
    return realmPtr_;
}

const realm* elementScalarField::realmPtr() const
{
    return realmPtr_;
}

realm& elementScalarField::realmRef()
{
    return *realmPtr_;
}

const realm& elementScalarField::realmRef() const
{
    return *realmPtr_;
}

} // namespace accel
