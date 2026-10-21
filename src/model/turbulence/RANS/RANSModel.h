// File       : RANSModel.h
// Created    : Fri Mar 15 2024 15:06:38 (+0100)
// Author     : Mhamad Mahdi Alloush
// Description: Base RANS turbulence model with common closure coefficients
// Copyright 2024 CCFNUM HSLU T&A. All Rights Reserved.

#ifndef RANSMODEL_H
#define RANSMODEL_H

// code
#include "nodeTensorField.h"
#include "turbulenceModel.h"

namespace accel
{

class RANSModel : public turbulenceModel
{
private:
    scalar betaStar_ = 0.09;

    scalar tkeProdLimitRatio_ = 10.0;

    // Spalart-Shur / Smirnov-Menter production correction.
    scalar cR1_ = 1.0;

    scalar cR2_ = 2.0;

    scalar cR3_ = 1.0;

    std::unique_ptr<nodeTensorField> SijPtr_ = nullptr;

public:
    RANSModel(realm* realm);

    // Shared two-equation production correction (SST and k-epsilon).
    void setupCurvatureCorrection(const std::shared_ptr<domain> domain);

    void updateCurvatureCorrection(const std::shared_ptr<domain> domain);

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
