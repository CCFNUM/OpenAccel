// File       : phasicContinuityEquation.cpp
// Description: Conservative Euler-Euler phase-fraction transport

#include "phasicContinuityEquation.h"

namespace accel
{

phasicContinuityEquation::phasicContinuityEquation(
    realm* realm,
    eulerEulerModel* model,
    label phaseIndex)
    : equation("Volume Fraction - " +
                   realm->simulationRef().materialName(phaseIndex),
               true),
      linearSystem(realm->simulationRef()), model_(model),
      phaseIndex_(phaseIndex),
      assembler_(
          std::make_unique<phasicContinuityAssembler>(model, phaseIndex))
{
    setEquationName(
        {"alpha." + realm->simulationRef().materialName(phaseIndex)});
    subIters_ = model_->controlsRef()
                    .solverRef()
                    .solverControl_.advancedOptions_.equationControls_
                    .subIterations_.volumeFraction_;
}

void phasicContinuityEquation::checkDomain(
    const std::shared_ptr<domain> domain)
{
    assert(domain->type() == domainType::fluid);
    assert(model_->isPhaseActive(domain, phaseIndex_));
    assert(!model_->isPrimaryPhase(domain, phaseIndex_));
}

bool phasicContinuityEquation::isConverged() const
{
    return linearSystem::isConverged();
}

void phasicContinuityEquation::setup()
{
    assembler_->setup(
        &model_->alphaRef(phaseIndex_), advection, domainVector_);
    linearSystem::setupSolver(name(),
                              model_->meshRef(),
                              domainZones_(),
                              fallbackName());
    isCreated_ = true;
}

void phasicContinuityEquation::initialize()
{
    FOREACH_DOMAIN(model_->updateVolumeFractionGradientField, phaseIndex_);
    FOREACH_DOMAIN(model_->updateVolumeFractionBlendingFactorField,
                   phaseIndex_);
    model_->alphaRef(phaseIndex_).updateScale();
}

void phasicContinuityEquation::postInitialize()
{
    isInitialized_ = true;
}

void phasicContinuityEquation::preSolve()
{
    FOREACH_DOMAIN(model_->updateVolumeFractionPrevIterField, phaseIndex_);
    FOREACH_DOMAIN(model_->updateVolumeFraction, phaseIndex_);
}

void phasicContinuityEquation::solve()
{
    auto ctx = linearSystem::getContext();
    ctx->zeroSystemStorage();
    FOREACH_DOMAIN_PTR(assembler_->preAssemble, ctx.get());
    FOREACH_DOMAIN_PTR(assembler_->assemble, ctx.get());
    FOREACH_DOMAIN_PTR(assembler_->postAssemble, ctx.get());
    assembler_->fix(collectInactiveInteriorParts(), {}, ctx.get(), {});

    if (ctx->getGraph()->isGraphMember())
    {
        linearSystem::solve();
    }
    messager::barrier();

    const scalar dt = model_->controlsRef().getTimestep();

    for (const auto& domain : domainVector_)
    {
        const scalar alphaMin =
            domain->multiphase_.residualVolumeFraction_;
        const scalar alphaMax = 1.0 - alphaMin;
        correctField_<linearSystem::BLOCKSIZE, 1, 0, 1>(
            domain.get(),
            ctx->coeffs().get(),
            stk::topology::NODE_RANK,
            model_->alphaRef(phaseIndex_).stkFieldRef(),
            1.0,
            0.0,
            1.0);

        auto& meshRef = model_->meshRef();
        auto& metaData = meshRef.metaDataRef();
        auto& bulkData = meshRef.bulkDataRef();

        auto& alpha = model_->alphaRef(phaseIndex_);
        auto& alphaField = alpha.stkFieldRef();
        const auto& alphaOldField = alpha.prevTimeRef().stkFieldRef();
        const auto& rhoField = model_->rhoRef(phaseIndex_).stkFieldRef();
        const auto& fluxDivergenceField =
            model_->mDotRef(phaseIndex_).divRef().stkFieldRef();
        const auto* dualVolumeField = metaData.get_field<scalar>(
            stk::topology::NODE_RANK, mesh::dual_nodal_volume_ID);
        STK_ThrowRequire(dualVolumeField != nullptr);

        const auto selection =
            metaData.locally_owned_part() &
            stk::mesh::selectUnion(domain->zonePtr()->interiorParts());
        const auto& buckets =
            bulkData.get_buckets(stk::topology::NODE_RANK, selection);

        alpha.synchronizeGhostedEntities(domain->index());
        model_->updateVolumeFraction(domain, phaseIndex_);
        // The pressure correction has already produced a conservative
        // intrinsic flux.  Alpha transport must keep it and refresh only the
        // phase weighting; reconstructing from U loses the p' correction.
        model_->updatePhaseMassFluxWeighting(domain, phaseIndex_);
        model_->updatePhaseFluxDivergence(domain, phaseIndex_);

        if (domain->multiphase_.eulerEulerFCT_.fluxCorrectedTransport_)
        {
            // Locally bounded, exactly conservative correction -- see
            // eulerEulerModel::correctVolumeFractionFCT. Deliberately not
            // followed by the projection loop below: that loop's non-local,
            // shared-factor redistribution is exactly the mechanism FCT
            // replaces, so running both would undo the fix.
            model_->correctVolumeFractionFCT(domain, phaseIndex_);
            model_->updateVolumeFraction(domain, phaseIndex_);
            model_->updatePhaseMassFluxWeighting(domain, phaseIndex_);
            model_->updatePhaseFluxDivergence(domain, phaseIndex_);
            continue;
        }

        // Retain the transported phase's old-time maximum. This is a global
        // discrete maximum principle, analogous to using the old solution as
        // the admissible range for a bounded flux correction. It contains no
        // case- or phase-specific value.
        scalar oldMaximum = alphaMin;
        for (const stk::mesh::Bucket* bucketPtr : buckets)
        {
            const auto& bucket = *bucketPtr;
            const scalar* oldAlphaValues =
                stk::mesh::field_data(alphaOldField, bucket);
            for (stk::mesh::Bucket::size_type node = 0;
                 node < bucket.size();
                 ++node)
            {
                oldMaximum = std::max(oldMaximum, oldAlphaValues[node]);
            }
        }
        messager::maxReduce(oldMaximum);
        const scalar projectionMax =
            std::max(alphaMin, std::min(alphaMax, oldMaximum));

        scalar oldMass = 0.0;
        scalar currentMass = 0.0;
        scalar boundaryRate = 0.0;
        scalar minimumMass = 0.0;
        scalar maximumMass = 0.0;
        for (const stk::mesh::Bucket* bucketPtr : buckets)
        {
            const auto& bucket = *bucketPtr;
            scalar* alphaValues =
                stk::mesh::field_data(alphaField, bucket);
            const scalar* oldAlphaValues =
                stk::mesh::field_data(alphaOldField, bucket);
            const scalar* rhoValues =
                stk::mesh::field_data(rhoField, bucket);
            const scalar* divergenceValues =
                stk::mesh::field_data(fluxDivergenceField, bucket);
            const scalar* volumes =
                stk::mesh::field_data(*dualVolumeField, bucket);
            for (stk::mesh::Bucket::size_type node = 0;
                 node < bucket.size();
                 ++node)
            {
                alphaValues[node] =
                    std::max(alphaMin,
                             std::min(projectionMax, alphaValues[node]));
                const scalar rhoVolume = rhoValues[node] * volumes[node];
                oldMass += rhoVolume * oldAlphaValues[node];
                currentMass += rhoVolume * alphaValues[node];
                boundaryRate += divergenceValues[node];
                minimumMass += rhoVolume * alphaMin;
                maximumMass += rhoVolume * projectionMax;
            }
        }
        messager::sumReduce(oldMass);
        messager::sumReduce(currentMass);
        messager::sumReduce(boundaryRate);
        messager::sumReduce(minimumMass);
        messager::sumReduce(maximumMass);

        // The generic scalar boundedness correction is required for stable
        // segregated phase velocities, but for div(U_k) != 0 it does not
        // conserve phase volume.  Project the bounded solution back onto the
        // exact global phase inventory implied by old-time alpha and the
        // discrete boundary flux. Positive corrections restore only the local
        // deficit relative to old-time alpha, weighted linearly by the local
        // old phase inventory. This preserves the transported profile instead
        // of creating phase in unrelated low-alpha regions. Negative
        // corrections are proportional to the removable phase inventory.
        // Residual-alpha and old-time maximum bounds are preserved throughout.
        const scalar targetMass = std::max(
            minimumMass,
            std::min(maximumMass, oldMass - dt * boundaryRate));

        for (label projection = 0; projection < 8; ++projection)
        {
            const scalar remaining = targetMass - currentMass;
            if (std::abs(remaining) <=
                1.0e-12 * (std::abs(targetMass) + 1.0))
            {
                break;
            }

            scalar activeWeight = 0.0;
            for (const stk::mesh::Bucket* bucketPtr : buckets)
            {
                const auto& bucket = *bucketPtr;
                const scalar* alphaValues =
                    stk::mesh::field_data(alphaField, bucket);
                const scalar* oldAlphaValues =
                    stk::mesh::field_data(alphaOldField, bucket);
                const scalar* rhoValues =
                    stk::mesh::field_data(rhoField, bucket);
                const scalar* volumes =
                    stk::mesh::field_data(*dualVolumeField, bucket);
                for (stk::mesh::Bucket::size_type node = 0;
                     node < bucket.size();
                     ++node)
                {
                    const scalar profileWeight =
                        remaining > 0.0
                            ? (alphaValues[node] < projectionMax
                                   ? std::max(oldAlphaValues[node] -
                                                  alphaValues[node],
                                              0.0) *
                                         std::max(oldAlphaValues[node],
                                                  alphaMin)
                                   : 0.0)
                            : std::max(alphaValues[node] - alphaMin, 0.0);
                    activeWeight +=
                        rhoValues[node] * volumes[node] * profileWeight;
                }
            }
            messager::sumReduce(activeWeight);
            if (activeWeight <= SMALL)
            {
                break;
            }

            const scalar factor = remaining / activeWeight;
            currentMass = 0.0;
            for (const stk::mesh::Bucket* bucketPtr : buckets)
            {
                const auto& bucket = *bucketPtr;
                scalar* alphaValues =
                    stk::mesh::field_data(alphaField, bucket);
                const scalar* oldAlphaValues =
                    stk::mesh::field_data(alphaOldField, bucket);
                const scalar* rhoValues =
                    stk::mesh::field_data(rhoField, bucket);
                const scalar* volumes =
                    stk::mesh::field_data(*dualVolumeField, bucket);
                for (stk::mesh::Bucket::size_type node = 0;
                     node < bucket.size();
                     ++node)
                {
                    const scalar profileWeight =
                        remaining > 0.0
                            ? (alphaValues[node] < projectionMax
                                   ? std::max(oldAlphaValues[node] -
                                                  alphaValues[node],
                                              0.0) *
                                         std::max(oldAlphaValues[node],
                                                  alphaMin)
                                   : 0.0)
                            : std::max(alphaValues[node] - alphaMin, 0.0);
                    alphaValues[node] = std::max(
                        alphaMin,
                        std::min(projectionMax,
                                 alphaValues[node] + factor * profileWeight));
                    currentMass += rhoValues[node] * volumes[node] *
                                   alphaValues[node];
                }
            }
            messager::sumReduce(currentMass);

            alpha.synchronizeGhostedEntities(domain->index());
        }
        model_->updateVolumeFraction(domain, phaseIndex_);
    }

    FOREACH_DOMAIN(model_->updateVolumeFractionGradientField, phaseIndex_);
    FOREACH_DOMAIN(model_->updateVolumeFractionBlendingFactorField,
                   phaseIndex_);
    model_->alphaRef(phaseIndex_).updateScale();

    // Do not feed the algebraic alpha-transport residual into momentum as a
    // mass-transfer source.  Boundedness/projection can alter alpha after
    // its linear solve, so ddt(alpha*rho)+div(mDot) is not a physical source
    // unless both are updated by one fully consistent transport operator.
}

void phasicContinuityEquation::preTimeStep()
{
}

void phasicContinuityEquation::setResidualScales_()
{
    linearSystem::getContext()->getResidualScales() = {
        1.0 / (model_->alphaRef(phaseIndex_).scale() + SMALL)};
}

void phasicContinuityEquation::printScales()
{
    if (messager::master())
    {
        std::cout << model_->alphaRef(phaseIndex_).name()
                  << " scale: " << std::scientific
                  << model_->alphaRef(phaseIndex_).scale() << '\n';
    }
}


} // namespace accel
