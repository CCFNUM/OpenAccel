// File : idealGasModel.h
// Created : Thu Apr 03 2025 17:05:11 (+0100)
// Author : Mhamad Mahdi Alloush
// Description: Ideal gas thermodynamic model for compressible flow simulations
// Copyright (c) 2025 CCFNUM, Lucerne University of Applied Sciences and Arts.
// SPDX-License-Identifier: BSD-3-Clause

#ifndef IDEALGASMODEL_H
#define IDEALGASMODEL_H

// code
#include "thermoModel.h"

namespace accel
{

class idealGasModel : public thermoModel
{
public:
    idealGasModel(realm* realm);

    virtual ~idealGasModel()
    {
    }

    using fieldBroker::h0Ref;
    using fieldBroker::hRef;
    using fieldBroker::pRef;
    using fieldBroker::psiRef;

    virtual void
    initializeDensity(const std::shared_ptr<domain> domain) override;

    virtual void initializeDensity(const std::shared_ptr<domain> domain,
                                   label iPhase) override;

    virtual void
    initializeCompressibility(const std::shared_ptr<domain> domain) override;

    virtual void initializeCompressibility(const std::shared_ptr<domain> domain,
                                           label iPhase) override;

    virtual void initializeSpecificHeatCapacity(
        const std::shared_ptr<domain> domain) override;

    virtual void
    initializeSpecificHeatCapacity(const std::shared_ptr<domain> domain,
                                   label iPhase) override;

    virtual void
    initializeSpecificEnthalpy(const std::shared_ptr<domain> domain) override;

    virtual void initializeSpecificTotalEnthalpy(
        const std::shared_ptr<domain> domain) override;

    virtual void
    initializeTemperature(const std::shared_ptr<domain> domain) override;

    virtual void updateDensity(const std::shared_ptr<domain> domain) override;

    virtual void updateDensity(const std::shared_ptr<domain> domain,
                               label iPhase) override;

    virtual void
    updateCompressibility(const std::shared_ptr<domain> domain) override;

    virtual void updateCompressibility(const std::shared_ptr<domain> domain,
                                       label iPhase) override;

    virtual void
    updateSpecificHeatCapacity(const std::shared_ptr<domain> domain) override;

    virtual void
    updateSpecificHeatCapacity(const std::shared_ptr<domain> domain,
                               label iPhase) override;

    virtual void
    updateSpecificEnthalpy(const std::shared_ptr<domain> domain) override;

    virtual void
    updateSpecificTotalEnthalpy(const std::shared_ptr<domain> domain) override;

    virtual void
    updateTemperature(const std::shared_ptr<domain> domain) override;
};

} /* namespace accel */

#endif // IDEALGASMODEL_H
