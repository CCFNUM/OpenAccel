// File       : meshMotion.h
// Created    : Fri Dec 13 2024 12:55:24 (+0100)
// Author     : Mhamad Mahdi Alloush
// Description: Mesh motion manager combining deformation and transformation
// Copyright (c) 2024 CCFNUM, Lucerne University of Applied Sciences and Arts.
// SPDX-License-Identifier: BSD-3-Clause

#ifndef MESHMOTION_H
#define MESHMOTION_H

// code
#include "deformation.h"
#include "fieldBroker.h"
#include "transformation.h"

namespace accel
{

class meshMotion : public fieldBroker
{
private:
    std::unique_ptr<deformation> deformationPtr_ = nullptr;

    std::unique_ptr<transformation> transformationPtr_ = nullptr;

public:
    // Constructors
    meshMotion(realm* realm);

    using fieldBroker::DRef;
    using fieldBroker::DtRef;

    // Methods

    void reset();

    void setup();

    void initialize();

    void update();

protected:
    // protected update methods
    void updateCoordinates_();

    void updateMeshVelocityField_();

    void updateMeshVelocityDivergenceField_();
};

} // namespace accel

#endif // MESHMOTION_H
