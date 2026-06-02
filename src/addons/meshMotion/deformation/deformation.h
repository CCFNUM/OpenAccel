// File       : deformation.h
// Created    : Fri Feb 14 2025 12:55:24 (+0100)
// Author     : Mhamad Mahdi Alloush
// Description: Mesh deformation via displacement diffusion equation
// Copyright 2025 CCFNUM HSLU T&A. All Rights Reserved.

#ifndef DEFORMATION_H
#define DEFORMATION_H

// code
#include "equation.h"
#include "types.h"

namespace accel
{

class meshMotion;

class deformation
{
private:
    meshMotion* meshMotionPtr_;

    std::unique_ptr<equation> displacementDiffusionEquation_ = nullptr;

public:
    // Constructors
    deformation(meshMotion* meshMotionPtr);

    // Methods

    void setup();

    void initialize();

    void update();
};

} /* namespace accel */

#endif // DEFORMATION_H
