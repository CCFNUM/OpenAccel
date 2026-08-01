// File       : SmagorinskyLESModel.h
// Created    : Fri Jul 31 2026
// Author     : GitHub Copilot
// Description: LES Smagorinsky turbulence model
// Copyright 2026 CCFNUM HSLU T&A. All Rights Reserved.

#ifndef SMAGORINSKYLESMODEL_H
#define SMAGORINSKYLESMODEL_H

#include "turbulenceModel.h"

namespace accel
{

class SmagorinskyLESModel : public turbulenceModel
{
private:
    scalar Cs_ = 0.17;

public:
    SmagorinskyLESModel(realm* realm);

    void updateTurbulentDynamicViscosity(
        const std::shared_ptr<domain> domain) override;

    scalar Cs() const
    {
        return Cs_;
    }
};

} /* namespace accel */

#endif // SMAGORINSKYLESMODEL_H
