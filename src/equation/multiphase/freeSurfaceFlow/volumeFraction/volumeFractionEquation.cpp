// File : volumeFractionEquation.cpp
// Created : Sun Jan 26 2025 22:57:38 (+0100)
// Author : Mhamad Mahdi Alloush
// Description:
// Copyright (c) 2025 CCFNUM, Lucerne University of Applied Sciences and Arts.
// SPDX-License-Identifier: BSD-3-Clause

#include "volumeFractionEquation.h"

namespace accel
{

volumeFractionEquation::volumeFractionEquation(realm* realm,
                                               freeSurfaceFlowModel* model,
                                               label phaseIndex)
    : equation("Volume Fraction - " +
                   realm->simulationRef().materialName(phaseIndex),
               true),
      model_(model), linearSystem(realm->simulationRef()),
      assembler_(std::make_unique<Assembler>(model, phaseIndex)),
      phaseIndex_(phaseIndex)
{
    this->setEquationName(
        {"alpha." + realm->simulationRef().materialName(phaseIndex)});

    // set sub-iterations
    subIters_ = model_->controlsRef()
                    .solverRef()
                    .solverControl_.advancedOptions_.equationControls_
                    .subIterations_.volumeFraction_;
}

void volumeFractionEquation::checkDomain(const std::shared_ptr<domain> domain)
{
    assert(domain->type() == domainType::fluid);
}

bool volumeFractionEquation::isConverged() const
{
    return linearSystem::isConverged();
}

void volumeFractionEquation::setup()
{
    // setup of fields on defined domains: rho field
    // might have been already initialized over other domains
    FOREACH_DOMAIN(model_->setupVolumeFraction, phaseIndex_);

    // setup FCT fields
    FOREACH_DOMAIN_IF(
        model_->setupFCTFields,
        domain->multiphase_.freeSurfaceModel_.fluxCorrectedTransport_,
        phaseIndex_);

    using Bucket = typename Assembler::Bucket;
    using BucketVector = typename Assembler::BucketVector;

    // setup assembler
    assembler_->setup(&model_->alphaRef(phaseIndex_), advection, domainVector_);

    // setup linear solver
    // FIXME: Consider passing mesh argument or
    // connectivity arrays passed to initialize() directly is more flexible
    // rather than this->meshRef() which is set through simulation object
    // obtained via realm in fieldBroker
    linearSystem::setupSolver(this->name(), model_->meshRef());

    equation::isCreated_ = true;
}

void volumeFractionEquation::initialize()
{
    // initialization of fields on defined domains: rho field
    // might have been already initialized over other domains
    FOREACH_DOMAIN(model_->initializeVolumeFraction, phaseIndex_);

    // post initialization

    // 1) update gradient
    FOREACH_DOMAIN(model_->updateVolumeFractionGradientField, phaseIndex_);

    // 2) update high-res fields
    FOREACH_DOMAIN(model_->updateVolumeFractionBlendingFactorField,
                   phaseIndex_);

    // 3) update scale
    model_->alphaRef(phaseIndex_).updateScale();

    equation::isInitialized_ = true;
}

void volumeFractionEquation::postInitialize()
{
}

void volumeFractionEquation::preSolve()
{
    FOREACH_DOMAIN(model_->updateVolumeFractionPrevIterField, phaseIndex_);
    FOREACH_DOMAIN(model_->updateVolumeFraction, phaseIndex_);

    // Note: smoothing iters for alpha and interface normal
    FOREACH_DOMAIN(model_->applyVOFSmoothing, phaseIndex_);
    FOREACH_DOMAIN(model_->updateInterfaceNormal, phaseIndex_);
}

void volumeFractionEquation::solve()
{
    auto ctx = linearSystem::getContext();
    ctx->zeroSystemStorage();

    // assembly
    linearSystem::simulationRef().getProfiler().push("linear_system_assembly");

    FOREACH_DOMAIN_PTR(assembler_->assemble, ctx.get());

    // fix system in domains where the model is not active
    assembler_->fix(
        this->collectInactiveInteriorParts(), {}, ctx.get(), {}, true);

    linearSystem::simulationRef().getProfiler().pop();

    // solve linear system
    linearSystem::solve();

    // correction
    // clip values in source field to `lower_clip_value` and `upper_clip_value`
    static constexpr int CLIP = 1; // true
    const scalar relax_value = 1.0;
    const scalar lower_clip_value = 0;
    const scalar upper_clip_value = 1;
    for (const auto& domain : domainVector_)
    {
        this->template correctField_<linearSystem::BLOCKSIZE, 1, 0, CLIP>(
            domain.get(),
            ctx->getXVector(),
            stk::topology::NODE_RANK,
            model_->alphaRef(phaseIndex_).stkFieldRef(),
            relax_value,
            lower_clip_value,
            upper_clip_value);

        // synchronize
        model_->alphaRef(phaseIndex_)
            .synchronizeGhostedEntities(domain->index());
    }

    // post correction

    // 1) update gradient
    FOREACH_DOMAIN(model_->updateVolumeFractionGradientField, phaseIndex_);

    // 2) update high-res fields
    FOREACH_DOMAIN(model_->updateVolumeFractionBlendingFactorField,
                   phaseIndex_);

    // 3) correct alpha using FCT
    FOREACH_DOMAIN_IF(
        model_->correctFCT,
        domain->multiphase_.freeSurfaceModel_.fluxCorrectedTransport_,
        phaseIndex_);

    // 4) update scale
    model_->alphaRef(phaseIndex_).updateScale();
}

void volumeFractionEquation::preTimeStep()
{
    FOREACH_DOMAIN(model_->updateVolumeFractionPrevTimeField, phaseIndex_);
}

void volumeFractionEquation::setResidualScales_()
{
    auto& residual_scales = linearSystem::getContext()->getResidualScales();
    residual_scales = {
        1.0 / (model_->alphaRef(phaseIndex_).scale() + ::accel::SMALL)};
}

void volumeFractionEquation::printScales()
{
    if (messager::master())
    {
        std::cout << model_->alphaRef(phaseIndex_).name()
                  << " scales:" << std::endl;
        std::cout << "\tmin: " << std::scientific << std::setprecision(8)
                  << model_->alphaRef(phaseIndex_).minScale() << std::endl;
        std::cout << "\tmax: " << std::scientific << std::setprecision(8)
                  << model_->alphaRef(phaseIndex_).maxScale() << std::endl;
        std::cout << "\tscale: " << std::scientific << std::setprecision(8)
                  << model_->alphaRef(phaseIndex_).scale() << std::endl
                  << std::endl;
    }
}

} /* namespace accel */
