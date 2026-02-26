// File       : elementScalarField.h
// Created    : Fri Aug 25 2023 12:55:24 (+0100)
// Author     : Mhamad Mahdi Alloush
// Description: Base class for scalar-valued element fields with realm access
// Copyright (c) 2023 CCFNUM, Lucerne University of Applied Sciences and Arts.
// SPDX-License-Identifier: BSD-3-Clause

#ifndef ELEMENTSCALARFIELD_H
#define ELEMENTSCALARFIELD_H

#include "elementField.h"

namespace accel
{

class simulation;
class realm;
class zone;

class elementScalarField : public elementField<scalar, 1>
{
protected:
    realm* realmPtr_;

public:
    // Constructors

    elementScalarField(realm* realmPtr,
                       std::string name,
                       unsigned numberOfStates);

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

#endif // ELEMENTSCALARFIELD_H
