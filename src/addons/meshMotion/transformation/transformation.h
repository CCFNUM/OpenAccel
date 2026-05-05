// File       : transformation.h
// Created    : Fri Feb 14 2025 12:55:24 (+0100)
// Author     : Mhamad Mahdi Alloush
// Description: Rigid mesh transformation for prescribed mesh motion
// Copyright 2025 CCFNUM HSLU T&A. All Rights Reserved.

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
