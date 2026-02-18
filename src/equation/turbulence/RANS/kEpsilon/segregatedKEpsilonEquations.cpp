// File : segregatedKEpsilonEquations.cpp
// Created : Thu Feb 22 2025 13:38:51 (+0100)
// Author : Achraf Nagihi
// Description:
// Copyright (c) 2025 CCFNUM, Lucerne University of Applied Sciences and Arts.
// SPDX-License-Identifier: BSD-3-Clause

#include "segregatedKEpsilonEquations.h"
#include "realm.h"

namespace accel
{

segregatedKEpsilonEquations::segregatedKEpsilonEquations(realm* realm)
    : equation("Segregated k-epsilon"), kEpsilonModel(realm)
{
    tke_eq_ =
        std::make_unique<turbulentKineticEnergyKEpsilonEquation>(realm, this);
    tdr_eq_ = std::make_unique<turbulentDissipationRateEquation>(realm, this);
}

void segregatedKEpsilonEquations::addDomain(std::shared_ptr<domain> domain)
{
    equation::addDomain(domain);

    tke_eq_->addDomain(domain);
    tdr_eq_->addDomain(domain);
}

bool segregatedKEpsilonEquations::isConverged() const
{
    return tke_eq_->isConverged() && tdr_eq_->isConverged();
}

void segregatedKEpsilonEquations::setup()
{
    tke_eq_->setup();
    tdr_eq_->setup();

    equation::isCreated_ = tke_eq_->isCreated() && tdr_eq_->isCreated();
}

void segregatedKEpsilonEquations::initialize()
{
    tke_eq_->initialize();
    tdr_eq_->initialize();

    FOREACH_DOMAIN(clipValues);
    FOREACH_DOMAIN(clipMinDistToWall);

    equation::isInitialized_ =
        tke_eq_->isInitialized() && tdr_eq_->isInitialized();
}

void segregatedKEpsilonEquations::postInitialize()
{
}

void segregatedKEpsilonEquations::preSolve()
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
    tdr_eq_->preSolve();
}

void segregatedKEpsilonEquations::solve()
{
    FOREACH_DOMAIN(updateTurbulentProduction);

    tke_eq_->solve();
    tdr_eq_->solve();

    // correction in solve() done using old gradient values -> update gradients
    // of both fields here post-correction
    FOREACH_DOMAIN(updateTurbulentKineticEnergyGradientField);
    FOREACH_DOMAIN(updateTurbulentDissipationRateGradientField);
}

void segregatedKEpsilonEquations::postSolve()
{
    tke_eq_->postSolve();
    tdr_eq_->postSolve();
}

void segregatedKEpsilonEquations::preTimeStep()
{
    tke_eq_->preTimeStep();
    tdr_eq_->preTimeStep();
}

} /* namespace accel */
