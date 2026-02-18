// File : navierStokesEquation.cpp
// Created : Wed Mar 13 2024 13:41:00 (+0100)
// Author : Fabian Wermelinger
// Description: Coupled momentum equation implementation details
// Copyright (c) 2024 CCFNUM, Lucerne University of Applied Sciences and Arts.
// SPDX-License-Identifier: BSD-3-Clause

#include "navierStokesEquation.h"
#include "realm.h"
#include "turbulenceModel.h"

namespace accel
{

navierStokesEquation::navierStokesEquation(realm* realm, flowModel* model)
    : equation("Coupled Navier-Stokes", true), model_(model),
      linearSystem(realm->simulationRef()),
      assembler_(std::make_unique<Assembler>(model))
{
    // Collect equation names
    std::vector<std::string> eqNames;

#if SPATIAL_DIM == 2
    eqNames.push_back("U-Mom");
    eqNames.push_back("V-Mom");
#elif SPATIAL_DIM == 3
    eqNames.push_back("U-Mom");
    eqNames.push_back("V-Mom");
    eqNames.push_back("W-Mom");
#endif

    // set
    this->setEquationName(eqNames);

    // set relaxation for velocity
    model_->URef().setURF(model_->controlsRef()
                              .solverRef()
                              .solverControl_.basicSettings_.convergenceControl_
                              .relaxationParameters_.velocityRelaxationFactor_);
}

void navierStokesEquation::checkDomain(const std::shared_ptr<domain> domain)
{
    assert(domain->type() == domainType::fluid);
}

bool navierStokesEquation::isConverged() const
{
    return linearSystem::isConverged();
}

void navierStokesEquation::setup()
{
    // setup DU coefficients and body forces
    FOREACH_DOMAIN_PTR(assembler_->setupDUCoefficients);
    FOREACH_DOMAIN(model_->setupBodyForces);

    // force creation of beta field for U if not yet
    if (model_->URef().blendingFactorPtr() == nullptr)
    {
        model_->URef().setupBlendingFactorField(
            /*enable only a dummy field*/ true);
    }

    if (model_->rhoRef().blendingFactorPtr() == nullptr)
    {
        model_->rhoRef().setupBlendingFactorField(
            /*enable only a dummy field*/ true);
    }

    // setup required fields
    FOREACH_DOMAIN(model_->setupDensity);
    FOREACH_DOMAIN(model_->setupVelocity);
    FOREACH_DOMAIN(model_->setupMassFlowRate);
    FOREACH_DOMAIN(model_->setupDynamicViscosity);

    // setup assembler
    assembler_->setup(&model_->URef(), null, domainVector_, nullptr);

    // setup linear solver
    linearSystem::setupSolver(this->name(), model_->meshRef());

    equation::isCreated_ = true;
}

void navierStokesEquation::initialize()
{
    // raw initialization of velocity
    FOREACH_DOMAIN(model_->initializeVelocity);

    // update velocity gradient
    FOREACH_DOMAIN(model_->updateVelocityGradientField);

    // update high-res fields
    FOREACH_DOMAIN(model_->updateVelocityBlendingFactorField);

    // update velocity scale
    model_->URef().updateScale();
}

void navierStokesEquation::postInitialize()
{
    // initialization of properties
    FOREACH_DOMAIN(model_->initializeDensity);
    FOREACH_DOMAIN_IF(model_->updateDensityGradientField,
                      domain->isMaterialCompressible());
    FOREACH_DOMAIN_IF(model_->updateDensityBlendingFactorField,
                      domain->isMaterialCompressible());
    FOREACH_DOMAIN(model_->initializeDynamicViscosity);

    // post initialization
    model_->rhoRef().updateScale();
    model_->muRef().updateScale();

    // raw initialization of mass flux: can be safely done now after
    // initialization of velocity and density at all active domains
    FOREACH_DOMAIN(model_->initializeMassFlowRate);

    // transform mDot to relative: availability of a moving frame
    // will be checked inside. Then, update divergence fields
    FOREACH_DOMAIN(model_->transformMassFlowRateToRelative);

    // in case of transient, old density must be updated before div update
    if (model_->controlsRef().isTransient())
    {
        FOREACH_DOMAIN(model_->updateDensityPrevTimeField);
    }

    // now update div field
    FOREACH_DOMAIN(model_->updateMassDivergenceField);

    equation::isInitialized_ = true;
}

void navierStokesEquation::preSolve()
{
    // velocity update
    FOREACH_DOMAIN(model_->updateVelocityPrevIterField);
    FOREACH_DOMAIN(model_->updateVelocity);

    // properties update
    FOREACH_DOMAIN(model_->updateDensity);
    FOREACH_DOMAIN(model_->updateDynamicViscosity);
    FOREACH_DOMAIN(model_->updateEffectiveDynamicViscosity); // laminar

    // update density-related fields
    FOREACH_DOMAIN_IF(model_->updateDensityGradientField,
                      domain->isMaterialCompressible());
    FOREACH_DOMAIN_IF(model_->updateDensityBlendingFactorField,
                      domain->isMaterialCompressible());

    // update body forces
    FOREACH_DOMAIN(model_->computeBodyForces);
    FOREACH_DOMAIN(model_->redistributeBodyForces);
}

void navierStokesEquation::solve()
{
    auto ctx = linearSystem::getContext();
    ctx->zeroSystemStorage();

    // assembly
    FOREACH_DOMAIN_PTR(assembler_->assemble, ctx.get());

    // compute pressure-diffusivity
    FOREACH_DOMAIN_PTR(assembler_->computeDUCoefficients, ctx.get());

    // relax momentum at mass specified (or nearly specified) boundaries
    FOREACH_DOMAIN_PTR(assembler_->assembleNormalRelaxation, ctx.get(), 0.75);

    // fix system in domains where the model is not active
    assembler_->fix(
        this->collectInactiveInteriorParts(), {}, ctx.get(), {}, true);

    if (!model_->controlsRef()
             .solverRef()
             .solverControl_.expertParameters_.disableMomentumPredictor_)
    {
        // solve linear system
        linearSystem::solve();

        // correction
        auto& U = model_->URef();
        for (const auto& domain : domainVector_)
        {
            scalar effectiveRelaxationFactor = 1.0;

            // find the ramp value (explicit relaxation) for compressible
            // domains
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

            this->template correctField_<linearSystem::BLOCKSIZE, SPATIAL_DIM>(
                domain.get(),
                ctx->getXVector(),
                stk::topology::NODE_RANK,
                U.stkFieldRef(),
                effectiveRelaxationFactor);

            // synchronize
            U.synchronizeGhostedEntities(domain->index());
        }

        // post correction

        // necessary to have a momentum-satisfying mass flux
        FOREACH_DOMAIN(model_->updateMassFlowRate);

        // velocity gradient must be only updated after corrector step ...

        // update scale
        U.updateScale();
    }
}

void navierStokesEquation::preTimeStep()
{
    FOREACH_DOMAIN(model_->updateVelocityPrevTimeField);
    FOREACH_DOMAIN(model_->updateDensityPrevTimeField);
}

void navierStokesEquation::setResidualScales_()
{
    auto& residual_scales = linearSystem::getContext()->getResidualScales();
    const scalar U_scale_inv = 1.0 / (model_->URef().scale() + ::accel::SMALL);
    for (int i = 0; i < SPATIAL_DIM; i++)
    {
        residual_scales[i] = U_scale_inv;
    }
}

void navierStokesEquation::printScales()
{
    if (messager::master())
    {
        std::cout << model_->URef().name() << " scales:" << std::endl;
        std::cout << "\tmin: " << std::scientific << std::setprecision(8)
                  << model_->URef().minScale() << std::endl;
        std::cout << "\tmax: " << std::scientific << std::setprecision(8)
                  << model_->URef().maxScale() << std::endl;
        std::cout << "\tscale: " << std::scientific << std::setprecision(8)
                  << model_->URef().scale() << std::endl
                  << std::endl;
    }
}

} /* namespace accel */
