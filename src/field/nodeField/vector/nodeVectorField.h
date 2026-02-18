// File : nodeVectorField.h
// Created : Tue Apr 20 2024 12:55:24 (+0100)
// Author : Mhamad Mahdi Alloush
// Description: Base class for vector-valued node fields with gradient
// correction
// Copyright (c) 2024 CCFNUM, Lucerne University of Applied Sciences
// and Arts.
// SPDX-License-Identifier: BSD-3-Clause

#ifndef NODEVECTORFIELD_H
#define NODEVECTORFIELD_H

#include "nodeField.h"

namespace accel
{

class realm;
class simulation;

class nodeVectorField : public nodeField<SPATIAL_DIM, SPATIAL_DIM * SPATIAL_DIM>
{
protected:
    realm* realmPtr_;

    // Methods

    void correctGradientField_(label iZone) override;

public:
    // Constructors

    // main constructor
    nodeVectorField(realm* realmPtr,
                    std::string name,
                    unsigned numberOfStates,
                    bool prevIter,
                    bool highResolution,
                    bool computeGradient,
                    bool correctedBoundaryNodeValues);

    // special constructor to use the class functionalities for an existing stk
    // field (i.e. grad computation)
    nodeVectorField(realm* realmPtr, STKScalarField* stkField_ptr);

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

#endif // NODEVECTORFIELD_H
