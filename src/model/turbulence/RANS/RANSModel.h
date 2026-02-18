// File : RANSModel.h
// Created : Fri Mar 15 2024 15:06:38 (+0100)
// Author : Mhamad Mahdi Alloush
// Description: Base RANS turbulence model with common closure coefficients
// Copyright (c) 2024 CCFNUM, Lucerne University of Applied Sciences and Arts.
// SPDX-License-Identifier: BSD-3-Clause

#ifndef RANSMODEL_H
#define RANSMODEL_H

// code
#include "turbulenceModel.h"

namespace accel
{

class RANSModel : public turbulenceModel
{
private:
    scalar betaStar_ = 0.09;

    scalar tkeProdLimitRatio_ = 10.0;

public:
    RANSModel(realm* realm);

    // Access

    scalar betaStar() const
    {
        return betaStar_;
    }

    scalar tkeProdLimitRatio() const
    {
        return tkeProdLimitRatio_;
    }
};

} /* namespace accel */

#endif // RANSMODEL_H
