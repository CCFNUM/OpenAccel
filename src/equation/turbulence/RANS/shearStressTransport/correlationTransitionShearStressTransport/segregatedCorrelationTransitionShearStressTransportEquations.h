// File : segregatedCorrelationTransitionShearStressTransportEquations.h
// Created : Sun Dec 29 2024
// Author : Adam Fares
// Description: Segregated equations for correlation-based transition SST
// model (Menter 2015). Solves k, omega, and gamma equations.
// Unlike the full transition model, ReTheta is not solved but
// computed locally using correlations.
// Copyright (c) 2024 CCFNUM, Lucerne University of Applied Sciences and Arts.
// SPDX-License-Identifier: BSD-3-Clause

#ifndef SEGREGATEDCORRELATIONTRANSITIONSHEARSTRESSTRANSPORTEQUATIONS_H
#define SEGREGATEDCORRELATIONTRANSITIONSHEARSTRESSTRANSPORTEQUATIONS_H

#include "correlationTransitionShearStressTransportModel.h"
#include "turbulentEddyFrequencySSTEquation.h"
#include "turbulentIntermittencyCorrelationTransitionSSTEquation.h"
#include "turbulentKineticEnergyCorrelationTransitionSSTEquation.h"

namespace accel
{

class segregatedCorrelationTransitionShearStressTransportEquations
    : public equation,
      public correlationTransitionShearStressTransportModel
{
public:
    static constexpr equationID ID =
        equationID::segregatedCorrelationTransitionShearStressTransport;

    segregatedCorrelationTransitionShearStressTransportEquations(realm* realm);

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

    std::unique_ptr<turbulentKineticEnergyCorrelationTransitionSSTEquation>
        tke_eq_;

    std::unique_ptr<turbulentEddyFrequencySSTEquation> tef_eq_;

    // Note: No ReTheta equation - uses local correlations instead
    std::unique_ptr<turbulentIntermittencyCorrelationTransitionSSTEquation>
        ti_eq_;
};

} /* namespace accel */

#endif // SEGREGATEDCORRELATIONTRANSITIONSHEARSTRESSTRANSPORTEQUATIONS_H
