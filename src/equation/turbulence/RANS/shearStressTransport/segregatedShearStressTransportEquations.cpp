// File       : segregatedShearStressTransportEquations.cpp
// Created    : Fri Mar 15 2024 15:06:38 (+0100)
// Author     : Mhamad Mahdi Alloush
// Description:
// Copyright (c) 2024 CCFNUM, Lucerne University of Applied Sciences and Arts.
// SPDX-License-Identifier: BSD-3-Clause

#include "segregatedShearStressTransportEquations.h"
#include "realm.h"

namespace accel
{

segregatedShearStressTransportEquations::
    segregatedShearStressTransportEquations(realm* realm)
    : equation("Segregated Shear-Stress-Transport"),
      shearStressTransportModel(realm)
{
    tke_eq_ = std::make_unique<turbulentKineticEnergySSTEquation>(realm, this);
    tef_eq_ = std::make_unique<turbulentEddyFrequencySSTEquation>(realm, this);
}

void segregatedShearStressTransportEquations::addDomain(
    std::shared_ptr<domain> domain)
{
    equation::addDomain(domain);

    tke_eq_->addDomain(domain);
    tef_eq_->addDomain(domain);
}

bool segregatedShearStressTransportEquations::isConverged() const
{
    return tke_eq_->isConverged() && tef_eq_->isConverged();
}

void segregatedShearStressTransportEquations::setup()
{
    tke_eq_->setup();
    tef_eq_->setup();

    equation::isCreated_ = tke_eq_->isCreated() && tef_eq_->isCreated();
}

void segregatedShearStressTransportEquations::initialize()
{
    tke_eq_->initialize();
    tef_eq_->initialize();

    FOREACH_DOMAIN(clipValues);
    FOREACH_DOMAIN(clipMinDistToWall);

    equation::isInitialized_ =
        tke_eq_->isInitialized() && tef_eq_->isInitialized();
}

void segregatedShearStressTransportEquations::postInitialize()
{
}

void segregatedShearStressTransportEquations::preSolve()
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
}

void segregatedShearStressTransportEquations::solve()
{
    // compute blending for SST model
    FOREACH_DOMAIN(updateFOneBlending);
    FOREACH_DOMAIN(updateTurbulentProduction);

    tke_eq_->solve();
    tef_eq_->solve();

    // correction in solve() done using old gradient values -> update gradients
    // of both fields here post-correction
    FOREACH_DOMAIN(updateTurbulentKineticEnergyGradientField);
    FOREACH_DOMAIN(updateTurbulentEddyFrequencyGradientField);
}

void segregatedShearStressTransportEquations::postSolve()
{
    tke_eq_->postSolve();
    tef_eq_->postSolve();
}

void segregatedShearStressTransportEquations::preTimeStep()
{
    tke_eq_->preTimeStep();
    tef_eq_->preTimeStep();
}

void segregatedShearStressTransportEquations::printScales()
{
    tke_eq_->printScales();
    tef_eq_->printScales();
}

} /* namespace accel */
