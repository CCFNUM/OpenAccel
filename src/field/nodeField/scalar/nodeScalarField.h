// File : nodeScalarField.h
// Created : Tue Apr 20 2024 12:55:24 (+0100)
// Author : Mhamad Mahdi Alloush
// Description: Base class for scalar-valued node fields with gradient
// correction
// Copyright (c) 2024 CCFNUM, Lucerne University of Applied Sciences
// and Arts.
// SPDX-License-Identifier: BSD-3-Clause

#ifndef NODESCALARFIELD_H
#define NODESCALARFIELD_H

#include "nodeField.h"

namespace accel
{

class realm;
class simulation;

class nodeScalarField : public nodeField<1, SPATIAL_DIM>
{
protected:
    realm* realmPtr_;

    // Methods

    void correctGradientField_(label iZone) override;

public:
    // Constructors

    // main constructor
    nodeScalarField(realm* realmPtr,
                    std::string name,
                    unsigned numberOfStates,
                    bool prevIter,
                    bool highResolution,
                    bool computeGradient,
                    bool correctedBoundaryNodeValues);

    // special constructor to use the class functionalities for an existing stk
    // field (i.e. grad computation)
    nodeScalarField(realm* realmPtr, STKScalarField* stkField_ptr);

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

#endif // NODESCALARFIELD_H
