// File       : segregatedTransitionShearStressTransportEquations.h
// Created    : Mon Jan 14 2025
// Author     : Adam Fares
// Description: Segregated solution of the transition SST turbulence equations
// Copyright (c) 2025 CCFNUM, Lucerne University of Applied Sciences and Arts.
// SPDX-License-Identifier: BSD-3-Clause

#ifndef SEGREGATEDTRANSITIONSHEARSTRESSTRANSPORTEQUATIONS_H
#define SEGREGATEDTRANSITIONSHEARSTRESSTRANSPORTEQUATIONS_H

#include "transitionOnsetReynoldsNumberTransitionSSTEquation.h"
#include "transitionShearStressTransportModel.h"
#include "turbulentEddyFrequencySSTEquation.h"
#include "turbulentIntermittencyTransitionSSTEquation.h"
#include "turbulentKineticEnergyTransitionSSTEquation.h"

namespace accel
{

class segregatedTransitionShearStressTransportEquations
    : public equation,
      public transitionShearStressTransportModel
{
public:
    static constexpr equationID ID =
        equationID::segregatedTransitionShearStressTransport;

    segregatedTransitionShearStressTransportEquations(realm* realm);

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

    std::unique_ptr<turbulentKineticEnergyTransitionSSTEquation> tke_eq_;

    std::unique_ptr<turbulentEddyFrequencySSTEquation> tef_eq_;

    std::unique_ptr<transitionOnsetReynoldsNumberTransitionSSTEquation> tor_eq_;

    std::unique_ptr<turbulentIntermittencyTransitionSSTEquation> ti_eq_;
};

} /* namespace accel */

#endif // SEGREGATEDTRANSITIONSHEARSTRESSTRANSPORTEQUATIONS_H
