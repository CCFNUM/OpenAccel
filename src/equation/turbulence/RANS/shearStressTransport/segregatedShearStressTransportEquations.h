// File : segregatedShearStressTransportEquations.h
// Created : Fri Mar 15 2024 15:06:38 (+0100)
// Author : Mhamad Mahdi Alloush
// Description: Segregated solution of the SST k-omega turbulence equations
// Copyright (c) 2024 CCFNUM, Lucerne University of Applied Sciences and Arts.
// SPDX-License-Identifier: BSD-3-Clause

#ifndef SEGREGATEDSHEARSTRESSTRANSPORTEQUATIONS_H
#define SEGREGATEDSHEARSTRESSTRANSPORTEQUATIONS_H

// code
#include "shearStressTransportModel.h"
#include "turbulentEddyFrequencySSTEquation.h"
#include "turbulentKineticEnergySSTEquation.h"

namespace accel
{

class segregatedShearStressTransportEquations : public equation,
                                                public shearStressTransportModel
{
public:
    static constexpr equationID ID = equationID::segregatedShearStressTransport;

    segregatedShearStressTransportEquations(realm* realm);

    void addDomain(std::shared_ptr<domain> domain) override;

    bool isConverged() const override;

    void setup() override;

    void initialize() override;

    void postInitialize() override;

    void preSolve() override;

    void solve() override;

    void postSolve() override;

    void preTimeStep() override;

    void printScales() override;

    equationID getID() override
    {
        return ID;
    }

private:
    // instances for segregated equations

    std::unique_ptr<turbulentKineticEnergySSTEquation> tke_eq_;

    std::unique_ptr<turbulentEddyFrequencySSTEquation> tef_eq_;
};

} /* namespace accel */

#endif // SEGREGATEDSHEARSTRESSTRANSPORTEQUATIONS_H
