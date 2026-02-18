// File : solidDisplacementEquation.cpp
// Created : Thu Dec 04 2025 08:42:10 (+0100)
// Author : Mhamad Mahdi Alloush
// Description:
// Copyright (c) 2025 CCFNUM, Lucerne University of Applied Sciences and Arts.
// SPDX-License-Identifier: BSD-3-Clause

#include "solidDisplacementEquation.h"

namespace accel
{

solidDisplacementEquation::solidDisplacementEquation(realm* realm)
    : equation("Solid Displacement"), solidMechanicsModel(realm),
      linearSystem(realm->simulationRef()),
      assembler_(std::make_unique<Assembler>(this))
{
    // Collect equation names
    std::vector<std::string> eqNames;

#if SPATIAL_DIM == 2
    eqNames.push_back("x-Disp");
    eqNames.push_back("y-Disp");
#elif SPATIAL_DIM == 3
    eqNames.push_back("x-Disp");
    eqNames.push_back("y-Disp");
    eqNames.push_back("z-Disp");
#endif

    // set
    this->setEquationName(eqNames);

    // set relaxation factor for solid displacement field. Default is 1
    this->DRef().setURF(
        controlsRef()
            .solverRef()
            .solverControl_.basicSettings_.convergenceControl_
            .relaxationParameters_.solidDisplacementRelaxationFactor_);

    // set sub-iterations
    subIters_ = this->controlsRef()
                    .solverRef()
                    .solverControl_.advancedOptions_.equationControls_
                    .subIterations_.solidDisplacement_;

    // Initialize Aitken acceleration parameters from config
    const auto& accel = controlsRef()
                            .solverRef()
                            .solverControl_.advancedOptions_.equationControls_
                            .acceleration_.solidDisplacement_;
    useAitken_ = accel.aitkenEnabled_;
    aitkenOmegaInit_ = accel.aitkenInitialOmega_;
    aitkenOmegaMin_ = accel.aitkenOmegaMin_;
    aitkenOmegaMax_ = accel.aitkenOmegaMax_;
    aitkenOmega_ = aitkenOmegaInit_;
}

void solidDisplacementEquation::checkDomain(
    const std::shared_ptr<domain> domain)
{
    assert(domain->type() == domainType::solid);
}

bool solidDisplacementEquation::isConverged() const
{
    return linearSystem::isConverged();
}

void solidDisplacementEquation::setup()
{
    // setup of fields on defined domains
    FOREACH_DOMAIN(setupDensity);
    FOREACH_DOMAIN(setupYoungModulus);
    FOREACH_DOMAIN(setupPoissonRatio);
    FOREACH_DOMAIN(setupDisplacement);

    using Bucket = typename Assembler::Bucket;
    using BucketVector = typename Assembler::BucketVector;

    // setup assembler: Gamma is dummy
    assembler_->setup(&DRef(), diffusion, domainVector_, 0.0);

    // setup linear solver
    // FIXME: [2024-03-13] Consider passing mesh argument or
    // connectivity arrays passed to initialize() directly is more flexible
    // rather than this->meshRef() which is set through simulation object
    // obtained via realm in fieldBroker
    linearSystem::setupSolver(this->name(), fieldBroker::meshRef());

    equation::isCreated_ = true;
}

void solidDisplacementEquation::initialize()
{
    // raw initialization of solid displacement
    FOREACH_DOMAIN(initializeDisplacement);

    // raw initialization of properties
    FOREACH_DOMAIN(initializeDensity);
    FOREACH_DOMAIN(initializeYoungModulus);
    FOREACH_DOMAIN(initializePoissonRatio);

    // update displacement gradient
    FOREACH_DOMAIN(updateDisplacementGradientField);

    // calculate initial stress and strain
    FOREACH_DOMAIN(updateStressAndStrain_);

    // update enthalpy scale
    DRef().updateScale();
    rhoRef().updateScale();
    ERef().updateScale();
    nuRef().updateScale();
}

void solidDisplacementEquation::postInitialize()
{
    equation::isInitialized_ = true;
}

void solidDisplacementEquation::preSolve()
{
    // solid displacement updates
    FOREACH_DOMAIN(updateDisplacementPrevIterField);
    FOREACH_DOMAIN(updateDisplacement);

    // other property updates
    FOREACH_DOMAIN(updateDensity);
    FOREACH_DOMAIN(updateYoungModulus);
    FOREACH_DOMAIN(updatePoissonRatio);
}

void solidDisplacementEquation::solve()
{
    auto ctx = linearSystem::getContext();
    ctx->zeroSystemStorage();

    // assembly
    linearSystem::simulationRef().getProfiler().push("linear_system_assembly");

    FOREACH_DOMAIN_PTR(assembler_->assemble, ctx.get());

    // fix system in domains where the model is not active
    assembler_->fix(this->collectInactiveInteriorParts(), {}, ctx.get());

    // fix system on all dirichlet boundaries
    assembler_->fix(this->collectDirichletBoundaryParts_(), {}, ctx.get());

    // fix system on all mixed boundaries: pass ignored dof's
    for (const auto& domain : domainVector_)
    {
        for (label iBoundary = 0; iBoundary < domain->zonePtr()->nBoundaries();
             iBoundary++)
        {
            boundaryConditionType bcType =
                this->DRef()
                    .boundaryConditionRef(domain->index(), iBoundary)
                    .type();

            switch (bcType)
            {
                case boundaryConditionType::mixed:
                    {
                        const scalar* fixedValueFlag =
                            (this->DRef()
                                 .boundaryConditionRef(domain->index(),
                                                       iBoundary)
                                 .template data<SPATIAL_DIM>(
                                     "fixed_value_flag"))
                                .value();

                        std::vector<label> ignoreDofs;
                        for (label i = 0; i < SPATIAL_DIM; i++)
                        {
                            if (fixedValueFlag[i] < 0.5)
                            {
                                ignoreDofs.push_back(i);
                            }
                        }
                        assembler_->fix(
                            domain->zonePtr()->boundaryPtr(iBoundary)->parts(),
                            {},
                            ctx.get(),
                            ignoreDofs);
                    }
                    break;

                default:
                    break;
            }
        }
    }

    linearSystem::simulationRef().getProfiler().pop();

    // solve linear system
    linearSystem::solve();

    // Compute relaxation factor (Aitken or static URF)
    const Vector& correction = ctx->getXVector();
    scalar relaxationFactor =
        useAitken_ ? computeAitkenOmega_(correction) : DRef().urf();

    // correction
    for (const auto& domain : domainVector_)
    {
        this->template correctField_<linearSystem::BLOCKSIZE, SPATIAL_DIM>(
            domain.get(),
            correction,
            stk::topology::NODE_RANK,
            DRef().stkFieldRef(),
            relaxationFactor);

        // synchronize
        DRef().synchronizeGhostedEntities(domain->index());
    }

    // post correction

    // 1) update gradient
    FOREACH_DOMAIN(updateDisplacementGradientField);

    // 2) update scales
    DRef().updateScale();
    rhoRef().updateScale();
    ERef().updateScale();
    nuRef().updateScale();
}

void solidDisplacementEquation::preTimeStep()
{
    FOREACH_DOMAIN(updateDisplacementPrevTimeField);

    // Reset Aitken state for new timestep
    if (useAitken_)
    {
        aitkenIter_ = 0;
        aitkenOmega_ = aitkenOmegaInit_;
        aitkenResidualPrev_.clear();
    }
}

void solidDisplacementEquation::postSolve()
{
    // calculate post-processed quantities
    FOREACH_DOMAIN(updateStressAndStrain_);
}

stk::mesh::PartVector
solidDisplacementEquation::collectDirichletBoundaryParts_()
{
    stk::mesh::PartVector incPartVec;
    for (const auto& domain : domainVector_)
    {
        for (label iBoundary = 0; iBoundary < domain->zonePtr()->nBoundaries();
             iBoundary++)
        {
            boundaryConditionType bcType =
                this->DRef()
                    .boundaryConditionRef(domain->index(), iBoundary)
                    .type();

            switch (bcType)
            {
                case boundaryConditionType::specifiedValue:
                    {
                        const auto& boundaryRef =
                            domain->zonePtr()->boundaryRef(iBoundary);

                        for (auto part : boundaryRef.parts())
                        {
                            incPartVec.push_back(part);
                        }
                    }
                    break;

                default:
                    break;
            }
        }
    }

    return incPartVec;
}

void solidDisplacementEquation::setResidualScales_()
{
    auto& residual_scales = linearSystem::getContext()->getResidualScales();
    const scalar disp_scale_inv = 1.0 / (this->DRef().scale() + ::accel::SMALL);
    for (int i = 0; i < SPATIAL_DIM; i++)
    {
        residual_scales[i] = disp_scale_inv;
    }
}

void solidDisplacementEquation::applyDependencyUpdates_(
    const domain* domain,
    const stk::mesh::EntityRank entityRank,
    STKScalarField& /*stk_dst*/)
{
    assert(entityRank == stk::mesh::EntityRank::NODE_RANK);

    // if the zone deforms, and a fluid-solid interface exists, transfer
    // displacement from solid side to fluid side
    if (domain->zonePtr()->meshDeforming())
    {
    }
}

void solidDisplacementEquation::printScales()
{
    if (messager::master())
    {
        std::cout << this->DRef().name() << " scales:" << std::endl;
        std::cout << "\tmin: " << std::scientific << std::setprecision(8)
                  << DRef().minScale() << std::endl;
        std::cout << "\tmax: " << std::scientific << std::setprecision(8)
                  << DRef().maxScale() << std::endl;
        std::cout << "\tscale: " << std::scientific << std::setprecision(8)
                  << DRef().scale() << std::endl
                  << std::endl;

        std::cout << rhoRef().name() << " scales:" << std::endl;
        std::cout << "\tmin: " << std::scientific << std::setprecision(8)
                  << rhoRef().minScale() << std::endl;
        std::cout << "\tmax: " << std::scientific << std::setprecision(8)
                  << rhoRef().maxScale() << std::endl;
        std::cout << "\tscale: " << std::scientific << std::setprecision(8)
                  << rhoRef().scale() << std::endl
                  << std::endl;

        std::cout << ERef().name() << " scales:" << std::endl;
        std::cout << "\tmin: " << std::scientific << std::setprecision(8)
                  << ERef().minScale() << std::endl;
        std::cout << "\tmax: " << std::scientific << std::setprecision(8)
                  << ERef().maxScale() << std::endl;
        std::cout << "\tscale: " << std::scientific << std::setprecision(8)
                  << ERef().scale() << std::endl
                  << std::endl;

        std::cout << nuRef().name() << " scales:" << std::endl;
        std::cout << "\tmin: " << std::scientific << std::setprecision(8)
                  << nuRef().minScale() << std::endl;
        std::cout << "\tmax: " << std::scientific << std::setprecision(8)
                  << nuRef().maxScale() << std::endl;
        std::cout << "\tscale: " << std::scientific << std::setprecision(8)
                  << nuRef().scale() << std::endl
                  << std::endl;
    }
}

scalar solidDisplacementEquation::computeAitkenOmega_(const Vector& correction)
{
    const size_t n = correction.size();

    if (aitkenIter_ == 0)
    {
        // First iteration: store correction, use initial omega
        aitkenResidualPrev_ = correction;
        aitkenIter_++;
        return aitkenOmegaInit_;
    }

    // Compute norms and dot products for Aitken formula
    // r_prev = previous correction, r_curr = current correction
    scalar normPrevSq = 0.0;  // |r_{k-1}|²
    scalar normCurrSq = 0.0;  // |r_k|²
    scalar dotPrevCurr = 0.0; // r_{k-1} · r_k

    for (size_t i = 0; i < n; ++i)
    {
        normPrevSq += aitkenResidualPrev_[i] * aitkenResidualPrev_[i];
        normCurrSq += correction[i] * correction[i];
        dotPrevCurr += aitkenResidualPrev_[i] * correction[i];
    }

    // MPI reduction for parallel runs
    messager::sumReduce(normPrevSq);
    messager::sumReduce(normCurrSq);
    messager::sumReduce(dotPrevCurr);

    // Compute Aitken relaxation factor:
    // ω_k = ω_{k-1} * (|r_{k-1}|² - r_{k-1}·r_k) / |Δr|²
    // where |Δr|² = |r_k - r_{k-1}|² = |r_k|² - 2*r_{k-1}·r_k + |r_{k-1}|²
    scalar normDrSq = normCurrSq - 2.0 * dotPrevCurr + normPrevSq;
    scalar numerator = normPrevSq - dotPrevCurr;

    if (normDrSq > SMALL)
    {
        // Compute Aitken omega using standard formula (preCICE-style)
        scalar omegaNew = aitkenOmega_ * (numerator / normDrSq);
        aitkenOmega_ = std::clamp(omegaNew, aitkenOmegaMin_, aitkenOmegaMax_);
    }

    // Store current correction for next iteration
    aitkenResidualPrev_ = correction;
    aitkenIter_++;

#ifndef NDEBUG
    // Print Aitken omega for monitoring
    if (messager::master())
    {
        std::cout << "Aitken omega: " << std::scientific << std::setprecision(4)
                  << aitkenOmega_ << " (iter " << aitkenIter_ << ")"
                  << std::endl;
    }
#endif

    return aitkenOmega_;
}

} /* namespace accel */
