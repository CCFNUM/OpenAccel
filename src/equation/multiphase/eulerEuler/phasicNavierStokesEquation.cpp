// File       : phasicNavierStokesEquation.cpp
// Description: One conservative momentum equation for an Euler-Euler phase

#include "phasicNavierStokesEquation.h"

namespace accel
{

phasicNavierStokesEquation::phasicNavierStokesEquation(
    realm* realm,
    eulerEulerModel* model,
    label phaseIndex)
    : equation("Coupled Navier-Stokes - " +
                   realm->simulationRef().materialName(phaseIndex),
               true),
      linearSystem(realm->simulationRef()), model_(model),
      phaseIndex_(phaseIndex),
      assembler_(
          std::make_unique<phasicNavierStokesAssembler>(model, phaseIndex))
{
    std::vector<std::string> names;
    const std::string suffix =
        "." + realm->simulationRef().materialName(phaseIndex);
    names.push_back("U-Mom" + suffix);
    names.push_back("V-Mom" + suffix);
#if SPATIAL_DIM == 3
    names.push_back("W-Mom" + suffix);
#endif
    setEquationName(names);
    model_->URef(phaseIndex)
        .setURF(model_->controlsRef()
                    .solverRef()
                    .solverControl_.basicSettings_.convergenceControl_
                    .relaxationParameters_.velocityRelaxationFactor_);
}

void phasicNavierStokesEquation::checkDomain(
    const std::shared_ptr<domain> domain)
{
    assert(domain->type() == domainType::fluid);
    assert(model_->isPhaseActive(domain, phaseIndex_));
}

bool phasicNavierStokesEquation::isConverged() const
{
    return linearSystem::isConverged();
}

void phasicNavierStokesEquation::setup()
{
    assembler_->setup(domainVector_);
    linearSystem::setupSolver(name(),
                              model_->meshRef(),
                              domainZones_(),
                              fallbackName());
    // let the nodal field take its boundary values from the node-side field
    model_->URef(phaseIndex_).setCorrectedBoundaryNodeValues(true);

    isCreated_ = true;
}

void phasicNavierStokesEquation::initialize()
{
    FOREACH_DOMAIN(model_->updateVelocityGradientField, phaseIndex_);
    FOREACH_DOMAIN(model_->updateVelocityBlendingFactorField, phaseIndex_);
    model_->updateVelocityScale(phaseIndex_);
}

void phasicNavierStokesEquation::postInitialize()
{
    isInitialized_ = true;
}

void phasicNavierStokesEquation::preSolve()
{
    FOREACH_DOMAIN(model_->updateVelocityPrevIterField, phaseIndex_);
    FOREACH_DOMAIN(model_->updateVelocity, phaseIndex_);
    FOREACH_DOMAIN(model_->applyPhaseVelocityDirichlet, phaseIndex_);
}

void phasicNavierStokesEquation::solve()
{
    auto ctx = linearSystem::getContext();
    ctx->zeroSystemStorage();
    FOREACH_DOMAIN_PTR(assembler_->preAssemble, ctx.get());
    FOREACH_DOMAIN_PTR(assembler_->assemble, ctx.get());
    FOREACH_DOMAIN_PTR(assembler_->postAssemble, ctx.get());
    assembler_->fix(collectInactiveInteriorParts(), {}, ctx.get(), {});

    // Prescribed-velocity boundaries: the value is already in the field, so
    // the correction there is zero. The pressure step skips the same nodes
    // (eulerEulerModel::correctPhaseVelocity), so the three stay consistent.
    for (const auto& domain : domainVector_)
    {
        const auto dirichletParts =
            model_->velocityDirichletParts(domain, phaseIndex_);
        if (!dirichletParts.empty())
        {
            assembler_->fix(dirichletParts, {}, ctx.get(), {});
        }
    }

    if (!model_->controlsRef()
             .solverRef()
             .solverControl_.expertParameters_.disableMomentumPredictor_)
    {
        if (ctx->getGraph()->isGraphMember())
        {
            linearSystem::solve();
        }
        messager::barrier();

        for (const auto& domain : domainVector_)
        {
            correctField_<linearSystem::BLOCKSIZE, SPATIAL_DIM>(
                domain.get(),
                ctx->coeffs().get(),
                stk::topology::NODE_RANK,
                model_->URef(phaseIndex_).stkFieldRef());
            model_->synchronizeVelocity(domain, phaseIndex_);
        }
        model_->updateVelocityScale(phaseIndex_);
    }
}

void phasicNavierStokesEquation::preTimeStep()
{
}

void phasicNavierStokesEquation::setResidualScales_()
{
    auto& scales = linearSystem::getContext()->getResidualScales();
    const scalar inverseScale =
        1.0 / (model_->URef(phaseIndex_).scale() + SMALL);
    std::fill(scales.begin(), scales.end(), inverseScale);
}

void phasicNavierStokesEquation::printScales()
{
    if (messager::master())
    {
        std::cout << model_->URef(phaseIndex_).name()
                  << " scale: " << std::scientific
                  << model_->URef(phaseIndex_).scale() << '\n';
    }
}


} // namespace accel
