// File : pressureCorrectionEquation.h
// Created : Thu Mar 14 2024 12:50:04 (+0100)
// Author : Fabian Wermelinger
// Description: Pressure correction (continuity) equation implementation details
// Copyright (c) 2024 CCFNUM, Lucerne University of Applied Sciences and Arts.
// SPDX-License-Identifier: BSD-3-Clause

#include "pressureCorrectionEquation.h"
#include "realm.h"

namespace accel
{

pressureCorrectionEquation::pressureCorrectionEquation(realm* realm,
                                                       flowModel* model)
    : equation("Pressure Correction", true), model_(model),
      linearSystem(realm->simulationRef()),
      assembler_(std::make_unique<Assembler>(model))
{
    this->setEquationName({"P-Mass"});

    // set relaxation for pressure: will be used as explicit
    model_->pRef().setURF(model_->controlsRef()
                              .solverRef()
                              .solverControl_.basicSettings_.convergenceControl_
                              .relaxationParameters_.pressureRelaxationFactor_);

    // set sub-iterations
    subIters_ = model_->controlsRef()
                    .solverRef()
                    .solverControl_.advancedOptions_.equationControls_
                    .subIterations_.pressureCorrection_;
}

void pressureCorrectionEquation::checkDomain(
    const std::shared_ptr<domain> domain)
{
    assert(domain->type() == domainType::fluid);
}

bool pressureCorrectionEquation::isConverged() const
{
    return linearSystem::isConverged();
}

void pressureCorrectionEquation::setup()
{
    // setup pressure
    FOREACH_DOMAIN(model_->setupPressure);

    // setup assembler
    assembler_->setup(&model_->pRef(), null, domainVector_, nullptr);

    // linear solver
    linearSystem::setupSolver(this->name(), model_->meshRef());

    equation::isCreated_ = true;
}

void pressureCorrectionEquation::initialize()
{
    // initialization
    FOREACH_DOMAIN(model_->initializePressure);

    // update gradient
    FOREACH_DOMAIN(model_->updatePressureGradientField);

    // update scale: model-based
    model_->updatePressureScale();

    equation::isInitialized_ = true;
}

void pressureCorrectionEquation::postInitialize()
{
}

void pressureCorrectionEquation::preSolve()
{
    FOREACH_DOMAIN(model_->updatePressurePrevIterField);
    FOREACH_DOMAIN(model_->updatePressure);

#ifndef NDEBUG
    // Zero SCL check field before assembly
    FOREACH_DOMAIN(model_->zeroSCLCheckField);
#endif /* NDEBUG */
}

void pressureCorrectionEquation::solve()
{
    auto ctx = linearSystem::getContext();
    ctx->zeroSystemStorage();

    // assembly
    FOREACH_DOMAIN_PTR(assembler_->assemble, ctx.get());

    // assemble diagonal-domainance for pressure correction equation in
    // compressible domains
    FOREACH_DOMAIN_PTR_IF(assembler_->assembleDiagonalDominance,
                          domain->isMaterialCompressible(),
                          ctx.get());

    // set reference pressure in closed zones
    FOREACH_DOMAIN_PTR_IF(assembler_->adjustMatrixForPressureReference,
                          domain->pressureLevelRequired(),
                          ctx.get());

    // fix system in domains where the model is not active
    assembler_->fix(
        this->collectInactiveInteriorParts(), {}, ctx.get(), {}, true);

    // solve linear system
    linearSystem::solve();

    // correction
    auto& p = model_->pRef();
    for (const auto& domain : domainVector_)
    {
        scalar effectiveRelaxationFactor = p.urf();

        // find the ramp value (explicit relaxation) for compressible domains
        if (domain->isMaterialCompressible())
        {
            scalar iter = model_->controlsRef().iter;
            label rampIter = 20;
            if (iter <= rampIter)
            {
                effectiveRelaxationFactor *=
                    std::max(scalar(iter) / scalar(rampIter), 0.1);
            }
        }

        this->template correctField_<linearSystem::BLOCKSIZE>(
            domain.get(),
            ctx->getXVector(),
            stk::topology::NODE_RANK,
            p.stkFieldRef(),
            effectiveRelaxationFactor);

        // synchronize
        p.synchronizeGhostedEntities(domain->index());
    }

    // post correction

    // pressure gradient must be only be updated after velocity correction ...

    // update scale (gradient postponed to after mass correction)
    p.updateScale();
}

void pressureCorrectionEquation::postSolve()
{
#ifndef NDEBUG
    // Sync SCL check field after assembly
    FOREACH_DOMAIN(model_->syncSCLCheckField);
#endif /* NDEBUG */
}

void pressureCorrectionEquation::preTimeStep()
{
    FOREACH_DOMAIN(model_->updatePressurePrevTimeField);
}

void pressureCorrectionEquation::setResidualScales_()
{
    auto& residual_scales = linearSystem::getContext()->getResidualScales();
    residual_scales = {1.0 / (model_->pRef().scale() + ::accel::SMALL)};
}

void pressureCorrectionEquation::printScales()
{
    if (messager::master())
    {
        std::cout << model_->pRef().name() << " scales:" << std::endl;
        std::cout << "\tmin: " << std::scientific << std::setprecision(8)
                  << model_->pRef().minScale() << std::endl;
        std::cout << "\tmax: " << std::scientific << std::setprecision(8)
                  << model_->pRef().maxScale() << std::endl;
        std::cout << "\tscale: " << std::scientific << std::setprecision(8)
                  << model_->pRef().scale() << std::endl
                  << std::endl;
    }
}

} /* namespace accel */
