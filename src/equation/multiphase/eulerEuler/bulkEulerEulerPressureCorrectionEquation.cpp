// File       : bulkEulerEulerPressureCorrectionEquation.cpp
// Description: Common pressure correction assembled from all phase fluxes

#include "bulkEulerEulerPressureCorrectionEquation.h"

namespace accel
{

bulkEulerEulerPressureCorrectionEquation::
    bulkEulerEulerPressureCorrectionEquation(realm* realm,
                                              eulerEulerModel* model)
    : equation("Pressure Correction", true),
      linearSystem(realm->simulationRef()), model_(model),
      assembler_(
          std::make_unique<bulkEulerEulerPressureCorrectionAssembler>(model))
{
    setEquationName({"P-Mass"});
    model_->pRef().setURF(
        model_->controlsRef()
            .solverRef()
            .solverControl_.basicSettings_.convergenceControl_
            .relaxationParameters_.pressureRelaxationFactor_);
    subIters_ = model_->controlsRef()
                    .solverRef()
                    .solverControl_.advancedOptions_.equationControls_
                    .subIterations_.pressureCorrection_;
}

void bulkEulerEulerPressureCorrectionEquation::checkDomain(
    const std::shared_ptr<domain> domain)
{
    assert(domain->type() == domainType::fluid);
}

bool bulkEulerEulerPressureCorrectionEquation::isConverged() const
{
    return linearSystem::isConverged();
}

void bulkEulerEulerPressureCorrectionEquation::setup()
{
    FOREACH_DOMAIN(model_->setupPressure);
    FOREACH_DOMAIN(model_->setupPressureCorrection);
    assembler_->setup(&model_->pRef(), null, domainVector_, nullptr);
    linearSystem::setupSolver(name(),
                              model_->meshRef(),
                              domainZones_(),
                              fallbackName());
    isCreated_ = true;
}

void bulkEulerEulerPressureCorrectionEquation::initialize()
{
    FOREACH_DOMAIN(model_->initializePressure);
    FOREACH_DOMAIN(model_->updatePressureGradientField);
    model_->updatePressureScale();
    isInitialized_ = true;
}

void bulkEulerEulerPressureCorrectionEquation::postInitialize()
{
}

void bulkEulerEulerPressureCorrectionEquation::preSolve()
{
    FOREACH_DOMAIN(model_->updatePressurePrevIterField);
    FOREACH_DOMAIN(model_->updatePressure);
}

stk::mesh::PartVector
bulkEulerEulerPressureCorrectionEquation::collectSpecifiedPressureParts_() const
{
    stk::mesh::PartVector parts;
    for (const auto& domain : domainVector_)
    {
        const zone* zonePtr = domain->zonePtr();
        for (label iBoundary = 0; iBoundary < zonePtr->nBoundaries();
             ++iBoundary)
        {
            const auto& boundaryRef = zonePtr->boundaryRef(iBoundary);
            const auto type = boundaryRef.type();
            if (type != boundaryPhysicalType::outlet &&
                type != boundaryPhysicalType::opening)
            {
                continue;
            }
            // Only when the pressure is genuinely prescribed. A mass-flow or
            // velocity-driven outlet leaves the level free and must not be
            // pinned here.
            const auto bcType =
                model_->pRef()
                    .boundaryConditionRef(domain->index(), iBoundary)
                    .type();
            if (bcType != boundaryConditionType::staticPressure &&
                bcType != boundaryConditionType::totalPressure)
            {
                continue;
            }
            for (auto* part : boundaryRef.parts())
            {
                parts.push_back(part);
            }
        }
    }
    return parts;
}

void bulkEulerEulerPressureCorrectionEquation::storePressureCorrection_(
    const domain* domain,
    ::linearSolver::coefficients<linearSystem::BLOCKSIZE>* coefficients)
{
    auto& correctionField = model_->pCorrRef();
    auto& field = correctionField.stkFieldRef();
    const Vector& correction = coefficients->getXVector();
    const auto& mesh = model_->meshRef();
    const auto& metaData = mesh.metaDataRef();
    const auto& bulkData = mesh.bulkDataRef();
    const auto selection =
        metaData.locally_owned_part() &
        stk::mesh::selectUnion(domain->zonePtr()->interiorParts());
    const auto& buckets =
        bulkData.get_buckets(stk::topology::NODE_RANK, selection);
    for (const stk::mesh::Bucket* bucketPtr : buckets)
    {
        const auto& bucket = *bucketPtr;
        scalar* values = stk::mesh::field_data(field, bucket);
        for (stk::mesh::Bucket::size_type nodeIndex = 0;
             nodeIndex < bucket.size();
             ++nodeIndex)
        {
            const int64_t row = coefficients->getGraph()->localToRow(
                bulkData.local_id(bucket[nodeIndex]));
            values[nodeIndex] = row < 0 ? 0.0 : correction[row];
        }
    }
    correctionField.synchronizeGhostedEntities(domain->index());
    correctionField.updateGradientField(domain->index());
}

void bulkEulerEulerPressureCorrectionEquation::solve()
{
    auto ctx = linearSystem::getContext();
    ctx->zeroSystemStorage();
    FOREACH_DOMAIN_PTR(assembler_->preAssemble, ctx.get());
    FOREACH_DOMAIN_PTR(assembler_->assemble, ctx.get());
    FOREACH_DOMAIN_PTR(assembler_->postAssemble, ctx.get());
    FOREACH_DOMAIN_PTR_IF(assembler_->adjustMatrixForPressureReference,
                          domain->pressureLevelRequired(),
                          ctx.get());
    assembler_->fix(collectInactiveInteriorParts(), {}, ctx.get(), {});

    // Pin p' = 0 where the pressure is prescribed. domain::pressureLevelRequired
    // is false as soon as an inlet/outlet/opening exists, on the assumption
    // that such a boundary anchors the level -- so if we do not impose the
    // Dirichlet here, nothing does, and the level drifts monotonically.
    {
        const auto pressureParts = collectSpecifiedPressureParts_();
        if (!pressureParts.empty())
        {
            assembler_->fix(pressureParts, {}, ctx.get(), {});
        }
    }

    if (ctx->getGraph()->isGraphMember())
    {
        linearSystem::solve();
    }
    messager::barrier();

    for (const auto& domain : domainVector_)
    {
        equation::correctField_<linearSystem::BLOCKSIZE>(
            domain.get(),
            ctx->coeffs().get(),
            stk::topology::NODE_RANK,
            model_->pRef().stkFieldRef(),
            model_->pRef().urf());
        model_->pRef().synchronizeGhostedEntities(domain->index());
        storePressureCorrection_(domain.get(), ctx->coeffs().get());
    }
    model_->pRef().updateScale();
}

void bulkEulerEulerPressureCorrectionEquation::postSolve()
{
}

void bulkEulerEulerPressureCorrectionEquation::preTimeStep()
{
    FOREACH_DOMAIN(model_->updatePressurePrevTimeField);
}

void bulkEulerEulerPressureCorrectionEquation::setResidualScales_()
{
    linearSystem::getContext()->getResidualScales() = {
        1.0 / (model_->pRef().scale() + SMALL)};
}

void bulkEulerEulerPressureCorrectionEquation::printScales()
{
    if (messager::master())
    {
        std::cout << model_->pRef().name() << " scale: " << std::scientific
                  << model_->pRef().scale() << '\n';
    }
}

} // namespace accel
