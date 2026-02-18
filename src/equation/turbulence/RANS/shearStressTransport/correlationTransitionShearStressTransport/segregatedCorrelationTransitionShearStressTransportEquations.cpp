// File : segregatedCorrelationTransitionShearStressTransportEquations.cpp
// Created : Sun Dec 29 2024
// Author : Adam Fares
// Description: Segregated equations for correlation-based transition SST
// model (Menter 2015)
// Copyright (c) 2024 CCFNUM, Lucerne University of Applied Sciences and Arts.
// SPDX-License-Identifier: BSD-3-Clause

#include "segregatedCorrelationTransitionShearStressTransportEquations.h"
#include "realm.h"

namespace accel
{

segregatedCorrelationTransitionShearStressTransportEquations::
    segregatedCorrelationTransitionShearStressTransportEquations(realm* realm)
    : equation("Segregated Correlation-Transition-SST"),
      correlationTransitionShearStressTransportModel(realm)
{
    tke_eq_ = std::make_unique<
        turbulentKineticEnergyCorrelationTransitionSSTEquation>(realm, this);
    tef_eq_ = std::make_unique<turbulentEddyFrequencySSTEquation>(realm, this);
    // Note: No ReTheta equation - correlation-based model uses local
    // correlations
    ti_eq_ = std::make_unique<
        turbulentIntermittencyCorrelationTransitionSSTEquation>(realm, this);
}

void segregatedCorrelationTransitionShearStressTransportEquations::addDomain(
    std::shared_ptr<domain> domain)
{
    equation::addDomain(domain);

    tke_eq_->addDomain(domain);
    tef_eq_->addDomain(domain);
    ti_eq_->addDomain(domain);
}

bool segregatedCorrelationTransitionShearStressTransportEquations::isConverged()
    const
{
    return tke_eq_->isConverged() && tef_eq_->isConverged() &&
           ti_eq_->isConverged();
}

void segregatedCorrelationTransitionShearStressTransportEquations::setup()
{
    tke_eq_->setup();
    tef_eq_->setup();
    ti_eq_->setup();

    equation::isCreated_ =
        tke_eq_->isCreated() && tef_eq_->isCreated() && ti_eq_->isCreated();
}

void segregatedCorrelationTransitionShearStressTransportEquations::initialize()
{
    tke_eq_->initialize();
    tef_eq_->initialize();
    ti_eq_->initialize();

    FOREACH_DOMAIN(clipValues);
    FOREACH_DOMAIN(clipMinDistToWall);

    equation::isInitialized_ = tke_eq_->isInitialized() &&
                               tef_eq_->isInitialized() &&
                               ti_eq_->isInitialized();
}

void segregatedCorrelationTransitionShearStressTransportEquations::
    postInitialize()
{
}

void segregatedCorrelationTransitionShearStressTransportEquations::preSolve()
{
    // Wall-functions
    FOREACH_DOMAIN(updateUStar);
    FOREACH_DOMAIN(updateYStar);
    FOREACH_DOMAIN(updateYPlus);
    FOREACH_DOMAIN(updateUTau);
    FOREACH_DOMAIN(updateUPlus);
    FOREACH_DOMAIN(updateUWallCoeffs);
    FOREACH_DOMAIN(updateWallShearStress);
    FOREACH_DOMAIN(updateDuPlusDyPlus);
    FOREACH_DOMAIN(updateTPlus);
    FOREACH_DOMAIN(updateTWallCoeffs);

    // Other
    FOREACH_DOMAIN(updateTurbulentDynamicViscosity);
    FOREACH_DOMAIN(updateEffectiveDynamicViscosity);
    FOREACH_DOMAIN(updateTurbulentThermalConductivity);
    FOREACH_DOMAIN(updateEffectiveThermalConductivity);

    // omega depends on the previous being executed first
    tke_eq_->preSolve();
    tef_eq_->preSolve();
    ti_eq_->preSolve();
}

void segregatedCorrelationTransitionShearStressTransportEquations::solve()
{
    // compute blending for SST model
    FOREACH_DOMAIN(updateFOneBlending);
    FOREACH_DOMAIN(updateTurbulentProduction);

    tke_eq_->solve();
    tef_eq_->solve();
    ti_eq_->solve();

    // correction in solve() done using old gradient values -> update gradients
    // of both fields here post-correction
    FOREACH_DOMAIN(updateTurbulentKineticEnergyGradientField);
    FOREACH_DOMAIN(updateTurbulentEddyFrequencyGradientField);
    FOREACH_DOMAIN(updateTurbulentIntermittencyGradientField);
}

void segregatedCorrelationTransitionShearStressTransportEquations::postSolve()
{
    tke_eq_->postSolve();
    tef_eq_->postSolve();
    ti_eq_->postSolve();
}

void segregatedCorrelationTransitionShearStressTransportEquations::preTimeStep()
{
    tke_eq_->preTimeStep();
    tef_eq_->preTimeStep();
    ti_eq_->preTimeStep();
}

void segregatedCorrelationTransitionShearStressTransportEquations::printScales()
{
    tke_eq_->printScales();
    tef_eq_->printScales();
    ti_eq_->printScales();
}

} /* namespace accel */
