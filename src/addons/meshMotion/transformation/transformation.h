// File : transformation.h
// Created : Fri Feb 14 2025 12:55:24 (+0100)
// Author : Mhamad Mahdi Alloush
// Description: Rigid mesh transformation for prescribed mesh motion
// Copyright (c) 2025 CCFNUM, Lucerne University of Applied Sciences and Arts.
// SPDX-License-Identifier: BSD-3-Clause

#ifndef TRANSFORMATION_H
#define TRANSFORMATION_H

// code
#include "mesh.h"
#include "types.h"

namespace accel
{

class meshMotion;

class transformation
{
private:
    meshMotion* meshMotionPtr_;

public:
    // Constructors
    transformation(meshMotion* meshMotionPtr);

    // Methods

    void setup();

    void initialize();

    void update();
};

} /* namespace accel */

#endif // TRANSFORMATION_H
