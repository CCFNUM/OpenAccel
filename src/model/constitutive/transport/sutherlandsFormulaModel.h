// File       : sutherlandsFormulaModel.h
// Created    : Thu Apr 03 2025 17:05:11 (+0100)
// Author     : Mhamad Mahdi Alloush
// Description: Sutherland's formula for temperature-dependent viscosity and
//              conductivity
// Copyright 2025 CCFNUM HSLU T&A. All Rights Reserved.

#ifndef SUTHERLANDSFORMULAMODEL_H
#define SUTHERLANDSFORMULAMODEL_H

// code
#include "transportModel.h"

namespace accel
{

class sutherlandsFormulaModel : public transportModel
{
public:
    sutherlandsFormulaModel(realm* realm);

    virtual ~sutherlandsFormulaModel()
    {
    }

    virtual void
    initializeDynamicViscosity(const std::shared_ptr<domain> domain) override;

    virtual void
    initializeDynamicViscosity(const std::shared_ptr<domain> domain,
                               label iPhase) override;

    virtual void initializeThermalConductivity(
        const std::shared_ptr<domain> domain) override;

    virtual void
    initializeThermalConductivity(const std::shared_ptr<domain> domain,
                                  label iPhase) override;

    virtual void
    updateDynamicViscosity(const std::shared_ptr<domain> domain) override;

    virtual void updateDynamicViscosity(const std::shared_ptr<domain> domain,
                                        label iPhase) override;

    virtual void
    updateThermalConductivity(const std::shared_ptr<domain> domain) override;

    virtual void updateThermalConductivity(const std::shared_ptr<domain> domain,
                                           label iPhase) override;
};

} /* namespace accel */

#endif // SUTHERLANDSFORMULAMODEL_H
