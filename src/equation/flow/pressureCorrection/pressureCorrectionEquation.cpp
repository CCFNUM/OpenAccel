// File       : pressureCorrectionEquation.h
// Created    : Thu Mar 14 2024 12:50:04 (+0100)
// Author     : Fabian Wermelinger
// Description: Pressure correction (continuity) equation implementation details
// Copyright 2024 CCFNUM HSLU T&A. All Rights Reserved.

#include "pressureCorrectionEquation.h"
#include "realm.h"
#include "scaling.h"
#include "simulation.h"

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

    // setup the dedicated pressure correction field
    FOREACH_DOMAIN(model_->setupPressureCorrection);

    // Robin-Neumann wall-acceleration side field (RN FSI only)
    if (anyDomainHasRobinNeumannFsi_())
    {
        FOREACH_DOMAIN(model_->setupRnWallAcceleration);

        // Robin-Neumann stored wall flux side field (Tukovic Eq. 30)
        FOREACH_DOMAIN(model_->setupRnWallFlux);

        // Robin-Neumann frozen explicit Robin term side field
        FOREACH_DOMAIN(model_->setupRnWallFluxE);
    }

    // setup assembler
    assembler_->setup(&model_->pRef(), null, domainVector_, nullptr);

    // linear solver
    linearSystem::setupSolver(this->name(),
                              model_->meshRef(),
                              this->domainZones_(),
                              this->fallbackName());

    equation::isCreated_ = true;
}

void pressureCorrectionEquation::initialize()
{
    // initialization
    FOREACH_DOMAIN(model_->initializePressure);

    // update gradient
    FOREACH_DOMAIN(model_->updatePressureGradientField);

    equation::isInitialized_ = true;
}

void pressureCorrectionEquation::postInitialize()
{
    // update scale: depends on rhoRef().scale() and URef().scale().  This call
    // must happen after the postInitialize() call that sets these scales!
    model_->updatePressureScale();
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

    // assemble
    FOREACH_DOMAIN_PTR(assembler_->preAssemble, ctx.get());
    FOREACH_DOMAIN_PTR(assembler_->assemble, ctx.get());
    FOREACH_DOMAIN_PTR(assembler_->postAssemble, ctx.get());

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
    assembler_->fix(this->collectInactiveInteriorParts(), {}, ctx.get(), {});

    // solve linear system
    if (ctx->getGraph()->isGraphMember())
    {
        linearSystem::solve();
    }
    messager::barrier();

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
            ctx->coeffs().get(),
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

void pressureCorrectionEquation::startOfOuterIteration()
{
    if (!anyRobinNeumannFsi(domainVector_))
    {
        return;
    }

    // reset Q_vol; E is recomputed and frozen at the first corrector
    FOREACH_DOMAIN(model_->zeroRnWallFlux);
    FOREACH_DOMAIN(model_->zeroRnWallFluxE);
    assembler_->requestFreezeE();
}

void pressureCorrectionEquation::postSolve()
{
#ifndef NDEBUG
    // Sync SCL check field after assembly
    FOREACH_DOMAIN(model_->syncSCLCheckField);
#endif /* NDEBUG */

    // Tukovic Eq. 30: accumulate Q_vol from the p' just solved
    if (anyRobinNeumannFsi(domainVector_))
    {
        for (const auto& domain : domainVector_)
        {
            assembler_->accumulateRnWallFlux(domain.get());
        }
    }

    // E stored by the first corrector; correctors 2..N reuse it
    assembler_->clearFreezeE();

    // TODO: Tukovic Eq. 31 (normal-velocity reconstruction)
}

bool pressureCorrectionEquation::anyDomainHasRobinNeumannFsi_() const
{
    return anyRobinNeumannFsi(domainVector_);
}

void pressureCorrectionEquation::preTimeStep()
{
    FOREACH_DOMAIN(model_->updatePressurePrevTimeField);

    // RN parameters (h_s, rho_s) once per timestep; no-slip fallback if invalid
    if (!domainVector_.empty())
    {
        const scalar dt = model_->controlsRef().getTimestep();
        const auto rnp = computeRobinNeumannParams(
            domainVector_[0]->simulationRef(), domainVector_, dt);
        assembler_->resetC1Stats();
        assembler_->invalidateRobinParams();

        // reset the RN stored wall flux (Tukovic Eq. 30) at step start
        if (anyRobinNeumannFsi(domainVector_))
        {
            FOREACH_DOMAIN(model_->zeroRnWallFlux);
        }

        if (rnp.valid)
        {
            assembler_->setRobinParams(rnp.hs, rnp.rhoSolid);
        }
    }
}

void pressureCorrectionEquation::postTimeStep()
{
    // robin_alpha.out: one row per physical timestep (after all outer FSI
    // iterations and pressure sub-iterations)
    writeRobinC1Diagnostics(*assembler_,
                            domainVector_,
                            robinAlphaStream_,
                            model_->controlsRef().time);
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
