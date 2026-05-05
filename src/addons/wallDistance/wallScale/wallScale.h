// File       : wallScale.h
// Created    : Tue Feb 10 2026 12:50:04 (+0100)
// Author     : Mhamad Mahdi Alloush
// Description: Wall scale computation via diffusion equation for wall distance
// Copyright 2026 CCFNUM HSLU T&A. All Rights Reserved.

#ifndef WALLSCALE_H
#define WALLSCALE_H

// code
#include "equation.h"
#include "types.h"

namespace accel
{

class wallDistance;

class wallScale
{
private:
    wallDistance* wallDistancePtr_;

    std::unique_ptr<equation> wallScaleDiffusionEquation_ = nullptr;

public:
    // Constructors
    wallScale(wallDistance* wallDistancePtr);

    // Methods

    void setup();

    void initialize();

    void update();
};

} // namespace accel

#endif // WALLSCALE_H
