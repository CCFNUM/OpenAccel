// File       : turbulentEddyFrequencyEquation.cpp
// Created    : Thu Mar 14 2024 12:50:04 (+0100)
// Author     : Mhamad Mahdi Alloush
// Description:
// Copyright (c) 2024 CCFNUM, Lucerne University of Applied Sciences and Arts.
// SPDX-License-Identifier: BSD-3-Clause

#include "turbulentEddyFrequencyEquation.h"
#include "realm.h"

namespace accel
{

turbulentEddyFrequencyEquation::turbulentEddyFrequencyEquation(
    realm* realm,
    turbulenceModel* model)
    : equation("Turbulent Eddy Frequency", true),
      linearSystem(realm->simulationRef()), model_(model)
{
    this->setEquationName({"TEF"});

    // set relaxation for omega
    model_->omegaRef().setURF(
        model_->meshRef()
            .controlsRef()
            .solverRef()
            .solverControl_.basicSettings_.convergenceControl_
            .relaxationParameters_.turbulenceRelaxationFactor_);
}

void turbulentEddyFrequencyEquation::checkDomain(
    const std::shared_ptr<domain> domain)
{
    assert(domain->type() == domainType::fluid);
}

bool turbulentEddyFrequencyEquation::isConverged() const
{
    return linearSystem::isConverged();
}

void turbulentEddyFrequencyEquation::setup()
{
    // setup of fields (initial values and boundary conditions)
    FOREACH_DOMAIN(model_->setupTurbulentEddyFrequency);

    // linear solver
    linearSystem::setupSolver(this->name(), model_->meshRef());

    equation::isCreated_ = true;
}

void turbulentEddyFrequencyEquation::initialize()
{
    // initialization of fields (initial values and boundary conditions)
    FOREACH_DOMAIN(model_->initializeTurbulentEddyFrequency);

    // post initialization

    // 1) update gradient
    FOREACH_DOMAIN(model_->updateTurbulentEddyFrequencyGradientField);

    // 2) update high-res fields
    FOREACH_DOMAIN(model_->updateTurbulentEddyFrequencyBlendingFactorField);

    // 3) update scale
    model_->omegaRef().updateScale();

    equation::isInitialized_ = true;
}

void turbulentEddyFrequencyEquation::postInitialize()
{
}

void turbulentEddyFrequencyEquation::preSolve()
{
    // raw update
    FOREACH_DOMAIN(model_->updateTurbulentEddyFrequencyPrevIterField);
    FOREACH_DOMAIN(model_->updateTurbulentEddyFrequency);
}

void turbulentEddyFrequencyEquation::preTimeStep()
{
    FOREACH_DOMAIN(model_->updateTurbulentEddyFrequencyPrevTimeField);
}

void turbulentEddyFrequencyEquation::setResidualScales_()
{
    auto& residual_scales = linearSystem::getContext()->getResidualScales();
    residual_scales = {1.0 / (model_->omegaRef().scale() + ::accel::SMALL)};
}

void turbulentEddyFrequencyEquation::printScales()
{
    if (messager::master())
    {
        std::cout << model_->omegaRef().name() << " scales:" << std::endl;
        std::cout << "\tmin: " << std::scientific << std::setprecision(8)
                  << model_->omegaRef().minScale() << std::endl;
        std::cout << "\tmax: " << std::scientific << std::setprecision(8)
                  << model_->omegaRef().maxScale() << std::endl;
        std::cout << "\tscale: " << std::scientific << std::setprecision(8)
                  << model_->omegaRef().scale() << std::endl
                  << std::endl;
    }
}

} /* namespace accel */
