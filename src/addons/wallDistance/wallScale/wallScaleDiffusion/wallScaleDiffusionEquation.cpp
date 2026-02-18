// File : wallScaleDiffusionEquation.cpp
// Created : Thu Mar 14 2024 12:50:04 (+0100)
// Author : Mhamad Mahdi Alloush
// Description:
// Copyright (c) 2024 CCFNUM, Lucerne University of Applied Sciences and Arts.
// SPDX-License-Identifier: BSD-3-Clause

#include "wallScaleDiffusionEquation.h"
#include "realm.h"

namespace accel
{

wallScaleDiffusionEquation::wallScaleDiffusionEquation(realm* realm)
    : equation("Wall Scale Diffusion", true), wallScaleDiffusionModel(realm),
      linearSystem(realm->simulationRef()),
      assembler_(std::make_unique<Assembler>(this))
{
    this->setEquationName({"YSCALE"});

    // set settings for wall scale field
    this->yScaleRef().setURF(
        this->meshRef()
            .controlsRef()
            .solverRef()
            .solverControl_.basicSettings_.convergenceControl_
            .relaxationParameters_.wallScaleRelaxationFactor_);
    yScaleRef().setMediumIndependent(false);
    yScaleRef().setInterpolationScheme(
        meshRef()
            .controlsRef()
            .solverRef()
            .solverControl_.basicSettings_.interpolationSchemeType_
            .wallScaleInterpolationType_);
    yScaleRef().setGradientInterpolationScheme(
        meshRef()
            .controlsRef()
            .solverRef()
            .solverControl_.basicSettings_.interpolationSchemeType_
            .wallScaleGradientInterpolationType_);
}

void wallScaleDiffusionEquation::checkDomain(
    const std::shared_ptr<domain> domain)
{
    assert(domain->type() == domainType::fluid);
}

bool wallScaleDiffusionEquation::isConverged() const
{
    return linearSystem::isConverged();
}

void wallScaleDiffusionEquation::setup()
{
    // setup of fields (initial values and boundary conditions)
    FOREACH_DOMAIN(this->setupWallScale);

    using Bucket = typename Assembler::Bucket;
    using BucketVector = typename Assembler::BucketVector;

    assembler_->setup(&this->yScaleRef(), diffusion, domainVector_, 1.0);

    // linear solver
    linearSystem::setupSolver(this->name(), this->meshRef());

    equation::isCreated_ = true;
}

void wallScaleDiffusionEquation::initialize()
{
    // initialization of fields (initial values and boundary conditions)
    FOREACH_DOMAIN(this->initializeWallScale);

    // post initialization

    // 1) update gradient
    FOREACH_DOMAIN(this->updateWallScaleGradientField);

    // 2) update scale
    this->yScaleRef().updateScale();

    equation::isInitialized_ = true;
}

void wallScaleDiffusionEquation::postInitialize()
{
}

void wallScaleDiffusionEquation::preSolve()
{
    // raw updates
    FOREACH_DOMAIN(this->updateWallScalePrevIterField);
    FOREACH_DOMAIN(this->updateWallScale);
}

void wallScaleDiffusionEquation::solve()
{
    auto ctx = linearSystem::getContext();
    ctx->zeroSystemStorage();

    // assembly
    FOREACH_DOMAIN_PTR(assembler_->assemble, ctx.get());

    // fix system in domains where the model is not active
    assembler_->fix(
        this->collectInactiveInteriorParts(), {}, ctx.get(), {}, true);

    // fix system on all dirichlet boundaries if required
    if (this->yScaleRef().correctedBoundaryNodeValues())
    {
        assembler_->fix(this->collectDirichletBoundaryParts_(), {}, ctx.get());
    }

    // solve linear system
    linearSystem::solve();

    // correction
    const scalar relax_value = 0.75;
    auto& y_scale = this->yScaleRef();
    for (const auto& domain : domainVector_)
    {
        this->template correctField_<linearSystem::BLOCKSIZE>(
            domain.get(),
            ctx->getXVector(),
            stk::topology::NODE_RANK,
            y_scale.stkFieldRef(),
            relax_value);

        // synchronize
        y_scale.synchronizeGhostedEntities(domain->index());
    }

    // post correction

    // update gradient
    FOREACH_DOMAIN(this->updateWallScaleGradientField);

    // update scale
    y_scale.updateScale();
}

void wallScaleDiffusionEquation::postSolve()
{
}

stk::mesh::PartVector
wallScaleDiffusionEquation::collectDirichletBoundaryParts_()
{
    stk::mesh::PartVector incPartVec;
    for (const auto& domain : domainVector_)
    {
#ifdef HAS_INTERFACE
        for (const interface* interf : domain->interfacesRef())
        {
            if (interf->isFluidSolidType())
            {
                for (auto part : interf->interfaceSideInfoPtr(domain->index())
                                     ->currentPartVec_)
                {
                    incPartVec.push_back(part);
                }
            }
        }
#endif /* HAS_INTERFACE */

        for (label iBoundary = 0; iBoundary < domain->zonePtr()->nBoundaries();
             iBoundary++)
        {
            boundaryConditionType bcType =
                this->yScaleRef()
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

void wallScaleDiffusionEquation::preTimeStep()
{
}

void wallScaleDiffusionEquation::setResidualScales_()
{
    auto& residual_scales = linearSystem::getContext()->getResidualScales();
    residual_scales = {1.0 / (this->yScaleRef().scale() + ::accel::SMALL)};
}

void wallScaleDiffusionEquation::printScales()
{
    if (messager::master())
    {
        std::cout << this->yScaleRef().name() << " scales:" << std::endl;
        std::cout << "\tmin: " << std::scientific << std::setprecision(8)
                  << this->yScaleRef().minScale() << std::endl;
        std::cout << "\tmax: " << std::scientific << std::setprecision(8)
                  << this->yScaleRef().maxScale() << std::endl;
        std::cout << "\tscale: " << std::scientific << std::setprecision(8)
                  << this->yScaleRef().scale() << std::endl
                  << std::endl;
    }
}

} /* namespace accel */
