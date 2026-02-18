// File : kineticTheoryModel.h
// Created : Tue Jun 17 2025 17:05:11 (+0100)
// Author : Mhamad Mahdi Alloush
// Description: Kinetic theory model for gas viscosity and thermal conductivity
// Copyright (c) 2025 CCFNUM, Lucerne University of Applied Sciences and Arts.
// SPDX-License-Identifier: BSD-3-Clause

#ifndef KINETICTHEORYMODEL_H
#define KINETICTHEORYMODEL_H

// code
#include "transportModel.h"

namespace accel
{

class kineticTheoryModel : public transportModel
{
public:
    kineticTheoryModel(realm* realm);

    virtual ~kineticTheoryModel()
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

#endif // KINETICTHEORYMODEL_H
