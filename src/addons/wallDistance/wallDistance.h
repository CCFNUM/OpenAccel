// File       : wallDistance.h
// Created    : Tue Feb 10 2026 12:50:04 (+0100)
// Author     : Mhamad Mahdi Alloush
// Description: Computes and manages wall distance fields using wall scale
//              methods
// Copyright 2026 CCFNUM HSLU T&A. All Rights Reserved.

#ifndef WALLDISTANCE_H
#define WALLDISTANCE_H

// code
#include "equation.h"
#include "fieldBroker.h"
#include "meshWave.h"
#include "wallScale.h"

namespace accel
{

class wallDistance : public fieldBroker
{
private:
    std::unique_ptr<wallScale> wallScalePtr_ = nullptr;
    std::unique_ptr<meshWave> meshWavePtr_ = nullptr;

public:
    // Constructors
    wallDistance(realm* realm);

    using fieldBroker::yMinRef;
    using fieldBroker::yScaleRef;

    // Methods

    void reset();

    void setup();

    void initialize();

    void update();
};

} // namespace accel

#endif // WALLDISTANCE_H
