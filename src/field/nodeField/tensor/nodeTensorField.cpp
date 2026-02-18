// File : nodeTensorField.cpp
// Created : Tue Apr 20 2024 12:55:24 (+0100)
// Author : Mhamad Mahdi Alloush
// Description:
// Copyright (c) 2024 CCFNUM, Lucerne University of Applied Sciences and Arts.
// SPDX-License-Identifier: BSD-3-Clause

#include "nodeTensorField.h"
#include "boundary.h"
#include "domain.h"
#include "realm.h"
#include "simulation.h"
#include "zone.h"

namespace accel
{

nodeTensorField::nodeTensorField(realm* realmPtr,
                                 std::string name,
                                 unsigned numberOfStates,
                                 bool prevIter)
    : nodeField<SPATIAL_DIM * SPATIAL_DIM,
                SPATIAL_DIM * SPATIAL_DIM * SPATIAL_DIM>(realmPtr->meshPtr(),
                                                         name,
                                                         numberOfStates,
                                                         prevIter),
      realmPtr_(realmPtr)
{
}

// Access

simulation* nodeTensorField::simulationPtr()
{
    return realmRef().simulationPtr();
}

const simulation* nodeTensorField::simulationPtr() const
{
    return realmRef().simulationPtr();
}

simulation& nodeTensorField::simulationRef()
{
    return realmRef().simulationRef();
}

const simulation& nodeTensorField::simulationRef() const
{
    return realmRef().simulationRef();
}

realm* nodeTensorField::realmPtr()
{
    return realmPtr_;
}

const realm* nodeTensorField::realmPtr() const
{
    return realmPtr_;
}

realm& nodeTensorField::realmRef()
{
    return *realmPtr_;
}

const realm& nodeTensorField::realmRef() const
{
    return *realmPtr_;
}

} // namespace accel
