// File : bulkPressureCorrectionEquation.cpp
// Created : Sun Feb 02 2025 20:44:56 (+0100)
// Author : Mhamad Mahdi Alloush
// Description:
// Copyright (c) 2025 CCFNUM, Lucerne University of Applied Sciences and Arts.
// SPDX-License-Identifier: BSD-3-Clause

#include "bulkPressureCorrectionEquation.h"
#include "realm.h"

namespace accel
{

bulkPressureCorrectionEquation::bulkPressureCorrectionEquation(
    realm* realm,
    freeSurfaceFlowModel* model)
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

void bulkPressureCorrectionEquation::checkDomain(
    const std::shared_ptr<domain> domain)
{
    assert(domain->type() == domainType::fluid);
}

bool bulkPressureCorrectionEquation::isConverged() const
{
    return linearSystem::isConverged();
}

void bulkPressureCorrectionEquation::setup()
{
    // setup of fields (initial values and boundary conditions)
    FOREACH_DOMAIN(model_->setupPressure);

    // setup assembler
    assembler_->setup(&model_->pRef(), null, domainVector_, nullptr);

    // linear solver
    linearSystem::setupSolver(this->name(), model_->meshRef());

    equation::isCreated_ = true;
}

void bulkPressureCorrectionEquation::initialize()
{
    // initialization
    FOREACH_DOMAIN(model_->initializePressure);

    // update gradient
    FOREACH_DOMAIN(model_->updatePressureGradientField);

    // update scale: model-based
    model_->updatePressureScale();

    equation::isInitialized_ = true;
}

void bulkPressureCorrectionEquation::postInitialize()
{
}

void bulkPressureCorrectionEquation::preSolve()
{
    // raw update
    FOREACH_DOMAIN(model_->updatePressurePrevIterField);
    FOREACH_DOMAIN(model_->updatePressure);

#ifndef NDEBUG
    // Zero SCL check field before assembly
    FOREACH_DOMAIN(model_->zeroSCLCheckField);
#endif /* NDEBUG */
}

void bulkPressureCorrectionEquation::solve()
{
    auto ctx = linearSystem::getContext();
    ctx->zeroSystemStorage();

    // assembly
    for (const auto& domain : domainVector_)
    {
        for (label iPhase = 0; iPhase < domain->nMaterials(); iPhase++)
        {
            label phaseIndex = domain->localToGlobalMaterialIndex(iPhase);

            assembler_->setPhaseIndex(phaseIndex);
            assembler_->assemble(domain.get(), ctx.get());
        }
    }

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

    // pressure gradient must be only updated after velocity correction ...

    // Update (gradient postponed to after mass correction)
    p.updateScale();
}

void bulkPressureCorrectionEquation::postSolve()
{
#ifndef NDEBUG
    // Sync SCL check field after assembly
    FOREACH_DOMAIN(model_->syncSCLCheckField);
#endif /* NDEBUG */
}

void bulkPressureCorrectionEquation::preTimeStep()
{
    FOREACH_DOMAIN(model_->updatePressurePrevTimeField);
}

void bulkPressureCorrectionEquation::setResidualScales_()
{
    auto& residual_scales = linearSystem::getContext()->getResidualScales();
    residual_scales = {1.0 / (model_->pRef().scale() + ::accel::SMALL)};
}

void bulkPressureCorrectionEquation::printScales()
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
