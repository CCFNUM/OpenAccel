// File       : nodeTensorField.h
// Created    : Tue Apr 20 2024 12:55:24 (+0100)
// Author     : Mhamad Mahdi Alloush
// Description: Node-based tensor field for second-order tensor quantities
// Copyright (c) 2024 CCFNUM, Lucerne University of Applied Sciences and Arts.
// SPDX-License-Identifier: BSD-3-Clause

#ifndef NODETENSORFIELD_H
#define NODETENSORFIELD_H

#include "nodeField.h"

namespace accel
{

class realm;
class simulation;

class nodeTensorField
    : public nodeField<SPATIAL_DIM * SPATIAL_DIM,
                       SPATIAL_DIM * SPATIAL_DIM * SPATIAL_DIM>
{
protected:
    realm* realmPtr_;

public:
    // Constructors

    nodeTensorField(realm* realmPtr,
                    std::string name,
                    unsigned numberOfStates,
                    bool prevIter);

    // Access

    simulation* simulationPtr();

    const simulation* simulationPtr() const;

    simulation& simulationRef();

    const simulation& simulationRef() const;

    realm* realmPtr();

    const realm* realmPtr() const;

    realm& realmRef();

    const realm& realmRef() const;
};

} // namespace accel

#endif // NODETENSORFIELD_H
