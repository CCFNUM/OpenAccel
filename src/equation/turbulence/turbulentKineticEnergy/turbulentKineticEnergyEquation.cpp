// File       : turbulentKineticEnergyEquation.cpp
// Created    : Thu Mar 14 2024 12:50:04 (+0100)
// Author     : Mhamad Mahdi Alloush
// Description:
// Copyright (c) 2024 CCFNUM, Lucerne University of Applied Sciences and Arts.
// SPDX-License-Identifier: BSD-3-Clause

#include "turbulentKineticEnergyEquation.h"
#include "realm.h"

namespace accel
{

turbulentKineticEnergyEquation::turbulentKineticEnergyEquation(
    realm* realm,
    turbulenceModel* model)
    : equation("Turbulent Kinetic Energy", true),
      linearSystem(realm->simulationRef()), model_(model)
{
    this->setEquationName({"TKE"});

    // set relaxation for k
    model_->kRef().setURF(
        model_->meshRef()
            .controlsRef()
            .solverRef()
            .solverControl_.basicSettings_.convergenceControl_
            .relaxationParameters_.turbulenceRelaxationFactor_);
}

void turbulentKineticEnergyEquation::checkDomain(
    const std::shared_ptr<domain> domain)
{
    assert(domain->type() == domainType::fluid);
}

bool turbulentKineticEnergyEquation::isConverged() const
{
    return linearSystem::isConverged();
}

void turbulentKineticEnergyEquation::setup()
{
    // setup of fields (initial values and boundary conditions)
    FOREACH_DOMAIN(model_->setupTurbulentKineticEnergy);

    // linear solver
    linearSystem::setupSolver(this->name(), model_->meshRef());

    equation::isCreated_ = true;
}

void turbulentKineticEnergyEquation::initialize()
{
    // initialization of fields (initial values and boundary conditions)
    FOREACH_DOMAIN(model_->initializeTurbulentKineticEnergy);

    // post initialization

    // 1) update gradient
    FOREACH_DOMAIN(model_->updateTurbulentKineticEnergyGradientField);

    // 2) update high-res fields
    FOREACH_DOMAIN(model_->updateTurbulentKineticEnergyBlendingFactorField);

    // 3) update scale
    model_->kRef().updateScale();

    equation::isInitialized_ = true;
}

void turbulentKineticEnergyEquation::postInitialize()
{
}

void turbulentKineticEnergyEquation::preSolve()
{
    // raw update
    FOREACH_DOMAIN(model_->updateTurbulentKineticEnergyPrevIterField);
    FOREACH_DOMAIN(model_->updateTurbulentKineticEnergy);
}

void turbulentKineticEnergyEquation::preTimeStep()
{
    FOREACH_DOMAIN(model_->updateTurbulentKineticEnergyPrevTimeField);
}

void turbulentKineticEnergyEquation::setResidualScales_()
{
    auto& residual_scales = linearSystem::getContext()->getResidualScales();
    residual_scales = {1.0 / (model_->kRef().scale() + ::accel::SMALL)};
}

void turbulentKineticEnergyEquation::printScales()
{
    if (messager::master())
    {
        std::cout << model_->kRef().name() << " scales:" << std::endl;
        std::cout << "\tmin: " << std::scientific << std::setprecision(8)
                  << model_->kRef().minScale() << std::endl;
        std::cout << "\tmax: " << std::scientific << std::setprecision(8)
                  << model_->kRef().maxScale() << std::endl;
        std::cout << "\tscale: " << std::scientific << std::setprecision(8)
                  << model_->kRef().scale() << std::endl
                  << std::endl;
    }
}

} /* namespace accel */
