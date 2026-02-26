// File       : thermoModel.h
// Created    : Thu Apr 03 2025 17:05:11 (+0100)
// Author     : Mhamad Mahdi Alloush
// Description: Base thermodynamic model for density, enthalpy, and temperature
// properties
// Copyright (c) 2025 CCFNUM, Lucerne University of Applied Sciences
// and Arts.
// SPDX-License-Identifier: BSD-3-Clause

#ifndef THERMOMODEL_H
#define THERMOMODEL_H

// code
#include "model.h"

namespace accel
{

class thermoModel : public model
{
public:
    thermoModel(realm* realm);

    virtual ~thermoModel()
    {
    }

    using fieldBroker::cpRef;
    using fieldBroker::lambdaRef;
    using fieldBroker::TRef;

    virtual void
    initializeDensity(const std::shared_ptr<domain> domain) override
    {
    }

    virtual void initializeDensity(const std::shared_ptr<domain> domain,
                                   label iPhase) override
    {
    }

    virtual void
    initializeCompressibility(const std::shared_ptr<domain> domain) override
    {
    }

    virtual void initializeCompressibility(const std::shared_ptr<domain> domain,
                                           label iPhase) override
    {
    }

    virtual void initializeSpecificHeatCapacity(
        const std::shared_ptr<domain> domain) override
    {
    }

    virtual void
    initializeSpecificHeatCapacity(const std::shared_ptr<domain> domain,
                                   label iPhase) override
    {
    }

    virtual void
    initializeSpecificEnthalpy(const std::shared_ptr<domain> domain) override
    {
    }

    virtual void initializeSpecificTotalEnthalpy(
        const std::shared_ptr<domain> domain) override
    {
    }

    virtual void
    initializeTemperature(const std::shared_ptr<domain> domain) override
    {
    }

    virtual void updateDensity(const std::shared_ptr<domain> domain) override
    {
    }

    virtual void updateDensity(const std::shared_ptr<domain> domain,
                               label iPhase) override
    {
    }

    virtual void
    updateCompressibility(const std::shared_ptr<domain> domain) override
    {
    }

    virtual void updateCompressibility(const std::shared_ptr<domain> domain,
                                       label iPhase) override
    {
    }

    virtual void
    updateSpecificHeatCapacity(const std::shared_ptr<domain> domain) override
    {
    }

    virtual void
    updateSpecificHeatCapacity(const std::shared_ptr<domain> domain,
                               label iPhase) override
    {
    }

    virtual void
    updateSpecificEnthalpy(const std::shared_ptr<domain> domain) override
    {
    }

    virtual void
    updateSpecificTotalEnthalpy(const std::shared_ptr<domain> domain) override
    {
    }

    virtual void
    updateTemperature(const std::shared_ptr<domain> domain) override
    {
    }

    static const scalar universalGasConstant_;
};

} /* namespace accel */

#endif // THERMOMODEL_H
