// File : transportModel.h
// Created : Thu Apr 03 2025 17:05:11 (+0100)
// Author : Mhamad Mahdi Alloush
// Description: Base transport model interface for viscosity and thermal
// conductivity
// Copyright (c) 2025 CCFNUM, Lucerne University of Applied
// Sciences and Arts.
// SPDX-License-Identifier: BSD-3-Clause

#ifndef TRANSPORTMODEL_H
#define TRANSPORTMODEL_H

// code
#include "model.h"

namespace accel
{

class transportModel : public model
{
public:
    transportModel(realm* realm);

    virtual ~transportModel()
    {
    }

    using fieldBroker::lambdaRef;
    using fieldBroker::muRef;
    using fieldBroker::TRef;

    virtual void
    initializeDynamicViscosity(const std::shared_ptr<domain> domain) override {
    };

    virtual void
    initializeDynamicViscosity(const std::shared_ptr<domain> domain,
                               label iPhase) override {};

    virtual void initializeThermalConductivity(
        const std::shared_ptr<domain> domain) override {};

    virtual void
    initializeThermalConductivity(const std::shared_ptr<domain> domain,
                                  label iPhase) override {};

    virtual void
    updateDynamicViscosity(const std::shared_ptr<domain> domain) override {};

    virtual void updateDynamicViscosity(const std::shared_ptr<domain> domain,
                                        label iPhase) override {};

    virtual void
    updateThermalConductivity(const std::shared_ptr<domain> domain) override {};

    virtual void updateThermalConductivity(const std::shared_ptr<domain> domain,
                                           label iPhase) override {};
};

} /* namespace accel */

#endif // TRANSPORTMODEL_H
