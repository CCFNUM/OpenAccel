// File       : eulerEulerModel.cpp
// Description: Inhomogeneous Euler-Euler multiphase model infrastructure

#include "eulerEulerModel.h"

#include "Common.h"
#include "simulation.h"
#include "zoneTransformation.h"

namespace accel
{

eulerEulerModel::eulerEulerModel(realm* realm) : multiphaseModel(realm)
{
    for (const auto& domain : realm->simulationRef().domainVector())
    {
        if (domain->multiphase_.option_ != multiphaseModelOption::eulerEuler)
        {
            continue;
        }

        for (const auto& material : domain->materialVector())
        {
            const label globalIndex =
                realm->simulationRef().materialIndex(material.name_);
            const auto alreadyRegistered =
                std::find_if(phases_.begin(),
                             phases_.end(),
                             [globalIndex](const phase& current)
            { return current.index_ == globalIndex; });

            if (alreadyRegistered == phases_.end())
            {
                phases_.emplace_back(globalIndex, material.name_);
            }
        }
    }

    for (const auto& currentPhase : phases_)
    {
        const label phaseIndex = currentPhase.index_;
        URef(phaseIndex);
        alphaRef(phaseIndex);
        rhoRef(phaseIndex);
        muRef(phaseIndex);
        mDotRef(phaseIndex);
        intrinsicMDotRef(phaseIndex);
        phaseMassCoefficientRef(phaseIndex);
        interphaseMomentumSourceRef(phaseIndex);
        dragDiagonalRef(phaseIndex);
        // must be declared here, before the STK meta data is committed
        phaseBodyForceRef(phaseIndex);
        duRef(phaseIndex);
        duTildeRef(phaseIndex);

        if (meshRef().anyZoneFrameRotating() ||
            meshRef().anyZoneMeshTransforming() ||
            meshRef().anyZoneMeshDeforming())
        {
            UrRef(phaseIndex);
        }
    }

    phaseSmoothedAlpha_.reserve(phases_.size());
    for (const auto& currentPhase : phases_)
    {
        phaseSmoothedAlpha_.push_back(std::make_unique<nodeField<1>>(
            realm->meshPtr(),
            "euler_euler_smoothed_volume_fraction." + currentPhase.name_,
            1,
            false));
    }
}

nodeField<1>& eulerEulerModel::phaseSmoothedAlphaRef(label phaseIndex)
{
    for (label phase = 0; phase < nPhases(); ++phase)
    {
        if (this->phaseIndex(phase) == phaseIndex)
        {
            return *phaseSmoothedAlpha_[phase];
        }
    }
    STK_ThrowRequireMsg(false, "Unknown Euler-Euler phase index");
    return *phaseSmoothedAlpha_.front();
}

void eulerEulerModel::updatePhaseSmoothedVolumeFraction(
    const std::shared_ptr<domain> domain)
{
    const auto& metaData = meshRef().metaDataRef();
    const auto& bulkData = meshRef().bulkDataRef();
    const auto& partVec = domain->zonePtr()->interiorParts();
    auto* rawField = FOriginalSTKFieldPtr_;
    auto* smoothedField = FSTKFieldPtr_;
    assert(rawField && smoothedField);

    const bool redistribute =
        controlsRef()
            .solverRef()
            .solverControl_.expertParameters_.bodyForceRedistribution_;

    const auto selection =
        metaData.universal_part() & stk::mesh::selectUnion(partVec);
    const auto& buckets =
        bulkData.get_buckets(stk::topology::NODE_RANK, selection);

    for (const label phaseIndex : activePhaseIndices(domain))
    {
        const auto& alphaField = alphaRef(phaseIndex).stkFieldRef();
        auto& targetRef = phaseSmoothedAlphaRef(phaseIndex);
        auto& targetField = targetRef.stkFieldRef();

        if (!redistribute)
        {
            // Body-force redistribution off: the drag keeps the pointwise
            // alpha, exactly as the body force does.
            for (const stk::mesh::Bucket* bucketPtr : buckets)
            {
                const auto& bucket = *bucketPtr;
                const scalar* alpha =
                    stk::mesh::field_data(alphaField, bucket);
                scalar* target = stk::mesh::field_data(targetField, bucket);
                for (stk::mesh::Bucket::size_type node = 0;
                     node < bucket.size();
                     ++node)
                {
                    target[node] = alpha[node];
                }
            }
            targetRef.synchronizeGhostedEntities(domain->index());
            continue;
        }

        // Stage alpha through the shared body-force scratch, one phase at a
        // time, exactly as updatePhaseBodyForces stages its force.
        for (const stk::mesh::Bucket* bucketPtr : buckets)
        {
            const auto& bucket = *bucketPtr;
            scalar* seed = stk::mesh::field_data(*rawField, bucket);
            const scalar* alpha = stk::mesh::field_data(alphaField, bucket);
            for (stk::mesh::Bucket::size_type node = 0; node < bucket.size();
                 ++node)
            {
                for (label component = 0; component < SPATIAL_DIM;
                     ++component)
                {
                    seed[SPATIAL_DIM * node + component] = alpha[node];
                }
            }
        }

        flowModel::redistributeBodyForces(domain);

        for (const stk::mesh::Bucket* bucketPtr : buckets)
        {
            const auto& bucket = *bucketPtr;
            const scalar* raw = stk::mesh::field_data(*rawField, bucket);
            const scalar* smoothed =
                stk::mesh::field_data(*smoothedField, bucket);
            scalar* target = stk::mesh::field_data(targetField, bucket);
            for (stk::mesh::Bucket::size_type node = 0; node < bucket.size();
                 ++node)
            {
                // Every seeded component carried alpha, so any component of
                // the result is the smoothed alpha. Same limiter as the body
                // force, so the two stay in step across a free surface.
                const scalar rawValue = raw[SPATIAL_DIM * node];
                const scalar limit = redistributionLimit_ * std::abs(rawValue);
                const scalar correction =
                    smoothed[SPATIAL_DIM * node] - rawValue;
                target[node] =
                    rawValue + std::max(-limit, std::min(limit, correction));
            }
        }
        targetRef.synchronizeGhostedEntities(domain->index());
    }
}

std::vector<label> eulerEulerModel::activePhaseIndices(
    const std::shared_ptr<domain> domain) const
{
    return activePhaseIndices(domain.get());
}

std::vector<label> eulerEulerModel::activePhaseIndices(
    const domain* domain) const
{
    std::vector<label> indices;
    if (domain->multiphase_.option_ != multiphaseModelOption::eulerEuler)
    {
        return indices;
    }

    indices.reserve(domain->nMaterials());
    for (label localIndex = 0; localIndex < domain->nMaterials(); ++localIndex)
    {
        indices.push_back(domain->localToGlobalMaterialIndex(localIndex));
    }
    return indices;
}

bool eulerEulerModel::isPhaseActive(const std::shared_ptr<domain> domain,
                                    label globalPhaseIndex) const
{
    const auto indices = activePhaseIndices(domain);
    return std::find(indices.begin(), indices.end(), globalPhaseIndex) !=
           indices.end();
}

bool eulerEulerModel::isPrimaryPhase(const std::shared_ptr<domain> domain,
                                     label globalPhaseIndex) const
{
    return domain->multiphase_.option_ == multiphaseModelOption::eulerEuler &&
           domain->multiphase_.primaryPhaseGlobalIndex_ == globalPhaseIndex;
}

const fluidPairModel* eulerEulerModel::findFluidPairModel(
    const std::shared_ptr<domain> domain,
    label globalPhaseIndexA,
    label globalPhaseIndexB) const
{
    for (const auto& pair : domain->fluidPairModels_)
    {
        if ((pair.globalIndexA_ == globalPhaseIndexA &&
             pair.globalIndexB_ == globalPhaseIndexB) ||
            (pair.globalIndexA_ == globalPhaseIndexB &&
             pair.globalIndexB_ == globalPhaseIndexA))
        {
            return &pair;
        }
    }
    return nullptr;
}

void eulerEulerModel::setupPhaseFields(const std::shared_ptr<domain> domain)
{
    for (const label phaseIndex : activePhaseIndices(domain))
    {
        setupVelocity(domain, phaseIndex);
        setupVolumeFraction(domain, phaseIndex);
        setupDensity(domain, phaseIndex);
        setupDynamicViscosity(domain, phaseIndex);
        const label localPhaseIndex =
            domain->globalToLocalMaterialIndex(phaseIndex);
        if (domain->isMaterialCompressible(localPhaseIndex))
        {
            psiRef(phaseIndex);
            setupCompressibility(domain, phaseIndex);
        }
        setupMassFlowRate(domain, phaseIndex);
        setupMomentumCorrection(domain, phaseIndex);
        setupEulerEulerAuxiliaryFields_(domain, phaseIndex);
    }
}

void eulerEulerModel::initializePhaseFields(
    const std::shared_ptr<domain> domain)
{
    for (const label phaseIndex : activePhaseIndices(domain))
    {
        fieldBroker::initializeVelocity(domain, phaseIndex);
        initializeVolumeFraction(domain, phaseIndex);
        fieldBroker::initializeDensity(domain, phaseIndex);
        fieldBroker::initializeDynamicViscosity(domain, phaseIndex);
        const label localPhaseIndex =
            domain->globalToLocalMaterialIndex(phaseIndex);
        if (domain->isMaterialCompressible(localPhaseIndex))
        {
            fieldBroker::initializeCompressibility(domain, phaseIndex);
        }
    }
    applyVolumeFractionClosure(domain);

    for (const label phaseIndex : activePhaseIndices(domain))
    {
        phaseMassCoefficientRef(phaseIndex).initialize(domain->index(), true);
        ops::zero(interphaseMomentumSourceRef(phaseIndex).stkFieldPtr(),
                  domain->zonePtr()->interiorParts());
        ops::zero(dragDiagonalRef(phaseIndex).stkFieldPtr(),
                  domain->zonePtr()->interiorParts());
        updatePhaseMassCoefficient(domain, phaseIndex);
    }
    updateInterphaseMomentumSources(domain);
    for (const label phaseIndex : activePhaseIndices(domain))
    {
        updatePhaseMassFlux(domain, phaseIndex, false);
        updatePhaseFluxDivergence(domain, phaseIndex);
    }
}

void eulerEulerModel::applyVolumeFractionClosure(
    const std::shared_ptr<domain> domain)
{
    const auto phaseIndices = activePhaseIndices(domain);
    const label primaryIndex = domain->multiphase_.primaryPhaseGlobalIndex_;
    if (phaseIndices.empty() || primaryIndex < 0)
    {
        errorMsg("Euler-Euler volume-fraction closure has no primary phase in "
                 "domain " +
                 domain->name());
    }

    auto& bulkData = meshRef().bulkDataRef();
    auto& metaData = meshRef().metaDataRef();
    const auto& partVec = domain->zonePtr()->interiorParts();
    const stk::mesh::Selector selection =
        metaData.universal_part() & stk::mesh::selectUnion(partVec);
    const auto& buckets =
        bulkData.get_buckets(stk::topology::NODE_RANK, selection);

    std::vector<label> secondaryIndices;
    for (const label phaseIndex : phaseIndices)
    {
        if (phaseIndex != primaryIndex)
        {
            secondaryIndices.push_back(phaseIndex);
        }
    }

    STKScalarField* primaryField = alphaRef(primaryIndex).stkFieldPtr();
    for (const stk::mesh::Bucket* bucketPtr : buckets)
    {
        const auto& bucket = *bucketPtr;
        scalar* primary = stk::mesh::field_data(*primaryField, bucket);
        std::vector<scalar*> secondary;
        secondary.reserve(secondaryIndices.size());
        for (const label phaseIndex : secondaryIndices)
        {
            secondary.push_back(
                stk::mesh::field_data(*alphaRef(phaseIndex).stkFieldPtr(),
                                      bucket));
        }

        for (stk::mesh::Bucket::size_type node = 0; node < bucket.size();
             ++node)
        {
            // Bound every secondary phase into [alphaMin, 1] first.
            const scalar alphaMin =
                domain->multiphase_.residualVolumeFraction_;
            scalar sum = 0.0;
            for (scalar* alpha : secondary)
            {
                alpha[node] = std::max(alphaMin, std::min(1.0, alpha[node]));
                sum += alpha[node];
            }

            // The primary phase takes up the remainder, so the secondaries may
            // occupy at most 1 - alphaMin for it to keep its own residual
            // fraction. Rescaling preserves the ratios between secondaries.
            const scalar sumMax = 1.0 - alphaMin;
            if (sum > sumMax)
            {
                const scalar rescale = sumMax / sum;
                for (scalar* alpha : secondary)
                {
                    alpha[node] *= rescale;
                }
                sum = sumMax;
            }
            primary[node] = 1.0 - sum;
        }
    }

    for (const label phaseIndex : phaseIndices)
    {
        alphaRef(phaseIndex).synchronizeGhostedEntities(domain->index());
    }
}

void eulerEulerModel::updatePhaseBodyForces(
    const std::shared_ptr<domain> domain)
{
    if (domain->buoyancy_.option_ == buoyancyOption::nonBuoyant)
    {
        return;
    }

    const auto& metaData = meshRef().metaDataRef();
    const auto& bulkData = meshRef().bulkDataRef();
    const auto& partVec = domain->zonePtr()->interiorParts();
    // flowModel::redistributeBodyForces zeroes `body_forces` and gathers the
    // raw nodal force from `body_forces_original`, scattering the redistributed
    // result back into `body_forces`. So the seed must go into FOriginal; the
    // answer is read back out of F.
    auto* rawForceField = FOriginalSTKFieldPtr_;
    auto* bodyForceField = FSTKFieldPtr_;
    assert(rawForceField && bodyForceField);

    const auto& gravity = domain->buoyancy_.gravity_;
    const scalar referenceDensity = domain->buoyancy_.referenceDensity_;

    const auto selection =
        metaData.universal_part() & stk::mesh::selectUnion(partVec);
    const auto& buckets =
        bulkData.get_buckets(stk::topology::NODE_RANK, selection);

    for (const label phaseIndex : activePhaseIndices(domain))
    {
        // Seed the shared scratch field with this phase's raw body force ...
        const auto& alphaField = alphaRef(phaseIndex).stkFieldRef();
        const auto& rhoField = rhoRef(phaseIndex).stkFieldRef();
        for (const stk::mesh::Bucket* bucketPtr : buckets)
        {
            const auto& bucket = *bucketPtr;
            scalar* force = stk::mesh::field_data(*rawForceField, bucket);
            const scalar* alpha = stk::mesh::field_data(alphaField, bucket);
            const scalar* rho = stk::mesh::field_data(rhoField, bucket);
            for (stk::mesh::Bucket::size_type node = 0; node < bucket.size();
                 ++node)
            {
                // Seed the full reduced-pressure body force. Keeping alpha_k
                // inside the force before redistribution gives every phase
                // the same compact treatment used by the pressure gradient.
                //
                // The redistribution is limited below, because across a free
                // surface alpha jumps two orders of magnitude and the
                // redistribution would otherwise haul force from the gas side
                // onto the liquid side.
                const scalar weight =
                    alpha[node] * (rho[node] - referenceDensity);
                for (label component = 0; component < SPATIAL_DIM; ++component)
                {
                    force[SPATIAL_DIM * node + component] =
                        weight * gravity[component];
                }
            }
        }

        // ... redistribute it with the shared (unmodified) flowModel kernel ...
        flowModel::redistributeBodyForces(domain);

        // ... and store it, limiting how far the redistribution may move the
        // force away from the raw value at each node.
        //
        // Where alpha is smooth the correction is a fraction of a percent and
        // the limiter never binds, so the hydrostatic balance is preserved
        // exactly as before. Across the interface the correction is an order
        // of magnitude larger than the local force -- that is the case that
        // produced a 508 N/m^3 body force at a node whose raw value was 31 --
        // and there it is clipped.
        auto& phaseForceField = phaseBodyForceRef(phaseIndex).stkFieldRef();
        for (const stk::mesh::Bucket* bucketPtr : buckets)
        {
            const auto& bucket = *bucketPtr;
            const scalar* force =
                stk::mesh::field_data(*bodyForceField, bucket);
            const scalar* raw = stk::mesh::field_data(*rawForceField, bucket);
            scalar* phaseForce =
                stk::mesh::field_data(phaseForceField, bucket);
            for (stk::mesh::Bucket::size_type node = 0; node < bucket.size();
                 ++node)
            {
                for (label component = 0; component < SPATIAL_DIM; ++component)
                {
                    const label off = SPATIAL_DIM * node + component;
                    const scalar rawValue = raw[off];
                    const scalar limit =
                        redistributionLimit_ * std::abs(rawValue);
                    const scalar correction = force[off] - rawValue;
                    phaseForce[off] =
                        rawValue +
                        std::max(-limit, std::min(limit, correction));
                }
            }
        }
        phaseBodyForceRef(phaseIndex)
            .synchronizeGhostedEntities(domain->index());
    }
}

void eulerEulerModel::setupEulerEulerAuxiliaryFields_(
    const std::shared_ptr<domain> domain,
    label phaseIndex)
{
    auto& phaseBodyForce = phaseBodyForceRef(phaseIndex);
    if (phaseBodyForce.isZoneUnset(domain->index()))
    {
        phaseBodyForce.setZone(domain->index());
    }
    auto& intrinsicFlux = intrinsicMDotRef(phaseIndex);
    if (intrinsicFlux.isZoneUnset(domain->index()))
    {
        intrinsicFlux.setZone(domain->index());
        intrinsicFlux.divRef().setZone(domain->index());

        for (interface* interf : domain->interfacesRef())
        {
            if (interf->isFluidSolidType())
            {
                continue;
            }

            if (interf->isInternal())
            {
                intrinsicFlux.registerSideFieldsForInterfaceSide(
                    interf->index(), true);
                intrinsicFlux.registerSideFieldsForInterfaceSide(
                    interf->index(), false);
            }
            else
            {
                intrinsicFlux.registerSideFieldsForInterfaceSide(
                    interf->index(), interf->isMasterZone(domain->index()));
            }
        }

        for (label iBoundary = 0;
             iBoundary < domain->zonePtr()->nBoundaries();
             ++iBoundary)
        {
            const auto type =
                domain->zonePtr()->boundaryRef(iBoundary).type();
            if (type == boundaryPhysicalType::inlet ||
                type == boundaryPhysicalType::outlet ||
                type == boundaryPhysicalType::opening)
            {
                intrinsicFlux.registerSideField(domain->index(), iBoundary);
            }
        }
    }

    if (phaseMassCoefficientRef(phaseIndex).isZoneUnset(domain->index()))
    {
        phaseMassCoefficientRef(phaseIndex).setZone(domain->index());
    }
    if (interphaseMomentumSourceRef(phaseIndex)
            .isZoneUnset(domain->index()))
    {
        interphaseMomentumSourceRef(phaseIndex).setZone(domain->index());
    }
    if (dragDiagonalRef(phaseIndex).isZoneUnset(domain->index()))
    {
        dragDiagonalRef(phaseIndex).setZone(domain->index());
    }
}

void eulerEulerModel::updatePhaseMassCoefficient(
    const std::shared_ptr<domain> domain,
    label phaseIndex)
{
    const auto& mesh = meshRef();
    const auto& metaData = mesh.metaDataRef();
    const auto& bulkData = mesh.bulkDataRef();
    const auto selection =
        metaData.universal_part() &
        stk::mesh::selectUnion(domain->zonePtr()->interiorParts());

    auto& coefficientField = phaseMassCoefficientRef(phaseIndex).stkFieldRef();
    const auto& alphaField = alphaRef(phaseIndex).stkFieldRef();
    const auto& rhoField = rhoRef(phaseIndex).stkFieldRef();

    const auto& buckets =
        bulkData.get_buckets(stk::topology::NODE_RANK, selection);
    for (const stk::mesh::Bucket* bucketPtr : buckets)
    {
        const auto& bucket = *bucketPtr;
        scalar* coefficient =
            stk::mesh::field_data(coefficientField, bucket);
        const scalar* alpha = stk::mesh::field_data(alphaField, bucket);
        const scalar* rho = stk::mesh::field_data(rhoField, bucket);
        for (stk::mesh::Bucket::size_type node = 0; node < bucket.size();
             ++node)
        {
            coefficient[node] = alpha[node] * rho[node];
        }
    }

    phaseMassCoefficientRef(phaseIndex)
        .synchronizeGhostedEntities(domain->index());
}

void eulerEulerModel::updatePhaseMassCoefficientPrevTime(
    const std::shared_ptr<domain> domain,
    label phaseIndex)
{
    phaseMassCoefficientRef(phaseIndex)
        .updatePrevTimeField(domain->index());
}

void eulerEulerModel::updateInterphaseMomentumSources(
    const std::shared_ptr<domain> domain)
{
    const auto phaseIndices = activePhaseIndices(domain);
    const auto& parts = domain->zonePtr()->interiorParts();
    for (const label phaseIndex : phaseIndices)
    {
        ops::zero(interphaseMomentumSourceRef(phaseIndex).stkFieldPtr(), parts);
        ops::zero(dragDiagonalRef(phaseIndex).stkFieldPtr(), parts);
    }

    const auto& mesh = meshRef();
    const auto& metaData = mesh.metaDataRef();
    const auto& bulkData = mesh.bulkDataRef();
    const auto selection =
        metaData.universal_part() & stk::mesh::selectUnion(parts);
    const auto& buckets =
        bulkData.get_buckets(stk::topology::NODE_RANK, selection);

    for (const auto& pair : domain->fluidPairModels_)
    {
        if (pair.drag_.option_ == dragModelOption::none)
        {
            continue;
        }

        const label phaseA = pair.globalIndexA_;
        const label phaseB = pair.globalIndexB_;
        auto& sourceAField =
            interphaseMomentumSourceRef(phaseA).stkFieldRef();
        auto& sourceBField =
            interphaseMomentumSourceRef(phaseB).stkFieldRef();
        auto& diagonalAField = dragDiagonalRef(phaseA).stkFieldRef();
        auto& diagonalBField = dragDiagonalRef(phaseB).stkFieldRef();
        const auto& velocityAField = URef(phaseA).stkFieldRef();
        const auto& velocityBField = URef(phaseB).stkFieldRef();
        // Smoothed alpha, matching the body force -- see
        // updatePhaseSmoothedVolumeFraction for why the drag must use the
        // same volume fraction the buoyancy does.
        const auto& alphaAField = phaseSmoothedAlphaRef(phaseA).stkFieldRef();
        const auto& alphaBField = phaseSmoothedAlphaRef(phaseB).stkFieldRef();
        const auto& alphaGradientAField =
            alphaRef(phaseA).gradRef().stkFieldRef();
        const auto& alphaGradientBField =
            alphaRef(phaseB).gradRef().stkFieldRef();
        const auto& rhoAField = rhoRef(phaseA).stkFieldRef();
        const auto& rhoBField = rhoRef(phaseB).stkFieldRef();
        const auto& muAField = muRef(phaseA).stkFieldRef();
        const auto& muBField = muRef(phaseB).stkFieldRef();
        const auto& volumeField = *metaData.get_field<scalar>(
            stk::topology::NODE_RANK, mesh::dual_nodal_volume_ID);

        // The nominal roles from the input. With blending enabled they hold
        // only where the named continuous phase really is continuous; the
        // inverted arrangement takes over elsewhere.
        const label dispersedPhase =
            pair.drag_.dispersedPhaseGlobalIndex_;
        const label continuousPhase =
            dispersedPhase == phaseA ? phaseB : phaseA;

        // `constant` drag carries no dispersed phase, so the roles are
        // meaningless there and the property fields are never read. Resolve
        // them to phaseA/phaseB in that case so the lookups stay in range.
        const bool rolesDefined = dispersedPhase >= 0;
        const label dispersedLookup = rolesDefined ? dispersedPhase : phaseA;
        const label continuousLookup = rolesDefined ? continuousPhase : phaseB;

        // Both phases' properties are needed: either can be the continuous one
        // depending on the local volume fraction.
        const auto& rhoDispersedField = rhoRef(dispersedLookup).stkFieldRef();
        const auto& muDispersedField = muRef(dispersedLookup).stkFieldRef();
        const auto& rhoContinuousField = rhoRef(continuousLookup).stkFieldRef();
        const auto& muContinuousField = muRef(continuousLookup).stkFieldRef();

        const bool blendRoles =
            rolesDefined &&
            pair.drag_.blending_ == dragBlendingOption::linear;
        const scalar aPartly = pair.drag_.minPartlyContinuous_;
        const scalar aFully = pair.drag_.minFullyContinuous_;
        const scalar dNominal = pair.drag_.diameter_;
        const scalar dInverted = pair.drag_.invertedDiameter_ > 0.0
                                     ? pair.drag_.invertedDiameter_
                                     : pair.drag_.diameter_;
        // The two linear dispersed weights leave a positive gap precisely
        // when the blending criteria permit a segregated regime.
        const bool segregatedEnabled =
            blendRoles && aFully + aPartly > 1.0 + SMALL;
        const scalar segregatedM = pair.drag_.segregatedM_;
        const scalar segregatedN = pair.drag_.segregatedN_;
        const scalar residualAlpha =
            domain->multiphase_.residualVolumeFraction_;

        for (const stk::mesh::Bucket* bucketPtr : buckets)
        {
            const auto& bucket = *bucketPtr;
            scalar* sourceA =
                stk::mesh::field_data(sourceAField, bucket);
            scalar* sourceB =
                stk::mesh::field_data(sourceBField, bucket);
            scalar* diagonalA =
                stk::mesh::field_data(diagonalAField, bucket);
            scalar* diagonalB =
                stk::mesh::field_data(diagonalBField, bucket);
            const scalar* velocityA =
                stk::mesh::field_data(velocityAField, bucket);
            const scalar* velocityB =
                stk::mesh::field_data(velocityBField, bucket);
            const scalar* alphaA =
                stk::mesh::field_data(alphaAField, bucket);
            const scalar* alphaB =
                stk::mesh::field_data(alphaBField, bucket);
            const scalar* alphaGradientA =
                stk::mesh::field_data(alphaGradientAField, bucket);
            const scalar* alphaGradientB =
                stk::mesh::field_data(alphaGradientBField, bucket);
            const scalar* rhoA = stk::mesh::field_data(rhoAField, bucket);
            const scalar* rhoB = stk::mesh::field_data(rhoBField, bucket);
            const scalar* muA = stk::mesh::field_data(muAField, bucket);
            const scalar* muB = stk::mesh::field_data(muBField, bucket);
            const scalar* volume =
                stk::mesh::field_data(volumeField, bucket);
            const scalar* rhoDispersed =
                stk::mesh::field_data(rhoDispersedField, bucket);
            const scalar* muDispersed =
                stk::mesh::field_data(muDispersedField, bucket);
            const scalar* rhoContinuous =
                stk::mesh::field_data(rhoContinuousField, bucket);
            const scalar* muContinuous =
                stk::mesh::field_data(muContinuousField, bucket);

            for (stk::mesh::Bucket::size_type node = 0;
                 node < bucket.size();
                 ++node)
            {
                scalar relativeSpeedSquared = 0.0;
                for (label component = 0; component < SPATIAL_DIM;
                     ++component)
                {
                    const scalar relativeVelocity =
                        velocityB[SPATIAL_DIM * node + component] -
                        velocityA[SPATIAL_DIM * node + component];
                    relativeSpeedSquared +=
                        relativeVelocity * relativeVelocity;
                }
                const scalar relativeSpeed =
                    std::sqrt(relativeSpeedSquared);

                scalar coefficient = pair.drag_.coefficient_;
                if (pair.drag_.option_ ==
                    dragModelOption::schillerNaumann)
                {
                    const scalar alphaContinuous =
                        continuousLookup == phaseA ? alphaA[node]
                                                   : alphaB[node];
                    const scalar alphaDispersed =
                        dispersedLookup == phaseA ? alphaA[node]
                                                  : alphaB[node];

                    // Schiller-Naumann for one choice of continuous phase.
                    auto schillerNaumannK =
                        [&](scalar rhoC,
                            scalar muC,
                            scalar diameter,
                            scalar alphaD)
                    {
                        const scalar reynolds = rhoC * relativeSpeed *
                                                diameter / (muC + SMALL);
                        const scalar Cd =
                            reynolds <= 1000.0
                                ? 24.0 / (reynolds + SMALL) *
                                      (1.0 +
                                       0.15 * std::pow(reynolds, 0.687))
                                : 0.44;
                        // OpenFOAM dispersedDragModel::K() multiplies Ki by
                        // the dispersed phase fraction only.  The continuous
                        // fraction is not an additional multiplier.  Keeping
                        // it here weakens drag around the interface and makes
                        // the phasic velocity extrema oscillatory.
                        return 0.75 * Cd * rhoC *
                               std::max(alphaD, residualAlpha) *
                               relativeSpeed / diameter;
                    };

                    // Weights for the two dispersed configurations. When the
                    // optional segregated model is active, their unoccupied
                    // fraction is assigned to that model exactly as in
                    // OpenFOAM's BlendedInterfacialModel.
                    scalar nominalWeight = 1.0;
                    scalar invertedWeight = 0.0;
                    if (blendRoles)
                    {
                        nominalWeight =
                            (alphaContinuous - aPartly) / (aFully - aPartly);
                        nominalWeight =
                            std::min(1.0, std::max(0.0, nominalWeight));
                        if (segregatedEnabled)
                        {
                            invertedWeight =
                                (alphaDispersed - aPartly) /
                                (aFully - aPartly);
                            invertedWeight =
                                std::min(1.0,
                                         std::max(0.0, invertedWeight));
                        }
                        else
                        {
                            invertedWeight = 1.0 - nominalWeight;
                        }
                    }

                    coefficient =
                        nominalWeight *
                        schillerNaumannK(rhoContinuous[node],
                                         muContinuous[node],
                                         dNominal,
                                         alphaDispersed);
                    if (invertedWeight > 0.0)
                    {
                        coefficient +=
                            invertedWeight *
                            schillerNaumannK(rhoDispersed[node],
                                             muDispersed[node],
                                             dInverted,
                                             alphaContinuous);
                    }

                    const scalar segregatedWeight = std::max(
                        0.0, 1.0 - nominalWeight - invertedWeight);
                    if (segregatedWeight > 0.0)
                    {
                        scalar gradientMagnitudeA2 = 0.0;
                        scalar gradientMagnitudeB2 = 0.0;
                        for (label component = 0;
                             component < SPATIAL_DIM;
                             ++component)
                        {
                            const label offset =
                                SPATIAL_DIM * node + component;
                            gradientMagnitudeA2 +=
                                alphaGradientA[offset] *
                                alphaGradientA[offset];
                            gradientMagnitudeB2 +=
                                alphaGradientB[offset] *
                                alphaGradientB[offset];
                        }
                        const scalar characteristicLength = std::pow(
                            std::max(volume[node], SMALL),
                            1.0 / static_cast<scalar>(SPATIAL_DIM));
                        const scalar gradientMagnitude = std::max(
                            (rhoB[node] * std::sqrt(gradientMagnitudeA2) +
                             rhoA[node] * std::sqrt(gradientMagnitudeB2)) /
                                (rhoA[node] + rhoB[node] + SMALL),
                            residualAlpha / (2.0 * characteristicLength));
                        const scalar limitedAlphaA = std::max(
                            alphaA[node], residualAlpha);
                        const scalar limitedAlphaB = std::max(
                            alphaB[node], residualAlpha);
                        const scalar alphaWeightedViscosity =
                            limitedAlphaA * muA[node] * limitedAlphaB *
                            muB[node] /
                            (limitedAlphaB * muA[node] +
                             limitedAlphaA * muB[node] + SMALL);
                        const scalar mixtureDensity =
                            alphaA[node] * rhoA[node] +
                            alphaB[node] * rhoB[node];
                        const scalar segregatedCoefficient =
                            segregatedM * mixtureDensity * relativeSpeed *
                                gradientMagnitude /
                                (limitedAlphaA * limitedAlphaB) +
                            segregatedN * alphaWeightedViscosity *
                                gradientMagnitude * gradientMagnitude;
                        coefficient +=
                            segregatedWeight * segregatedCoefficient;
                    }
                }

                coefficient = std::max(0.0, coefficient);
                diagonalA[node] += coefficient;
                diagonalB[node] += coefficient;
                for (label component = 0; component < SPATIAL_DIM;
                     ++component)
                {
                    const scalar force =
                        coefficient *
                        (velocityB[SPATIAL_DIM * node + component] -
                         velocityA[SPATIAL_DIM * node + component]);
                    sourceA[SPATIAL_DIM * node + component] += force;
                    sourceB[SPATIAL_DIM * node + component] -= force;
                }
            }
        }
    }

    for (const label phaseIndex : phaseIndices)
    {
        interphaseMomentumSourceRef(phaseIndex)
            .synchronizeGhostedEntities(domain->index());
        dragDiagonalRef(phaseIndex)
            .synchronizeGhostedEntities(domain->index());
    }
}

void eulerEulerModel::fillPhaseEffectiveViscosity(
    const domain* domain,
    label phaseIndex,
    STKScalarField& field) const
{
    const auto& mesh = meshRef();
    const auto& metaData = mesh.metaDataRef();
    const auto& bulkData = mesh.bulkDataRef();
    const auto selection =
        metaData.universal_part() &
        stk::mesh::selectUnion(domain->zonePtr()->interiorParts());
    const auto& alphaField = alphaRef(phaseIndex).stkFieldRef();
    const auto& muField = muRef(phaseIndex).stkFieldRef();

    const auto& buckets =
        bulkData.get_buckets(stk::topology::NODE_RANK, selection);
    for (const stk::mesh::Bucket* bucketPtr : buckets)
    {
        const auto& bucket = *bucketPtr;
        scalar* effectiveMu = stk::mesh::field_data(field, bucket);
        const scalar* alpha = stk::mesh::field_data(alphaField, bucket);
        const scalar* mu = stk::mesh::field_data(muField, bucket);
        for (stk::mesh::Bucket::size_type node = 0; node < bucket.size();
             ++node)
        {
            effectiveMu[node] = alpha[node] * mu[node];
        }
    }
}

void eulerEulerModel::correctPhaseVelocity(
    const std::shared_ptr<domain> domain,
    label phaseIndex)
{
    const bool consistent = controlsRef()
                                .solverRef()
                                .solverControl_.expertParameters_.consistent_;
    const auto& correctionCoefficient =
        consistent ? duTildeRef(phaseIndex).stkFieldRef()
                   : duRef(phaseIndex).stkFieldRef();
    const auto& pressureCorrectionGradient =
        pCorrRef().gradRef().stkFieldRef();
    auto& velocity = URef(phaseIndex).stkFieldRef();

    const auto& metaData = meshRef().metaDataRef();
    const auto& bulkData = meshRef().bulkDataRef();
    // A prescribed velocity has zero correction. Correcting those nodes
    // silently undoes the boundary condition and, because the boundary mass
    // flux is built from the prescribed value while the nodal velocity is not,
    // leaves the two inconsistent.
    const auto dirichletParts = velocityDirichletParts(domain, phaseIndex);
    auto selection = metaData.universal_part() &
                     stk::mesh::selectUnion(domain->zonePtr()->interiorParts());
    if (!dirichletParts.empty())
    {
        selection &= !stk::mesh::selectUnion(dirichletParts);
    }
    const auto& buckets =
        bulkData.get_buckets(stk::topology::NODE_RANK, selection);
    for (const stk::mesh::Bucket* bucketPtr : buckets)
    {
        const auto& bucket = *bucketPtr;
        scalar* phaseVelocity = stk::mesh::field_data(velocity, bucket);
        const scalar* coefficient =
            stk::mesh::field_data(correctionCoefficient, bucket);
        const scalar* gradient =
            stk::mesh::field_data(pressureCorrectionGradient, bucket);
        for (stk::mesh::Bucket::size_type node = 0; node < bucket.size();
             ++node)
        {
            for (label component = 0; component < SPATIAL_DIM; ++component)
            {
                const label offset = SPATIAL_DIM * node + component;
                phaseVelocity[offset] -= coefficient[offset] * gradient[offset];
            }
        }
    }
    URef(phaseIndex).synchronizeGhostedEntities(domain->index());
}

void eulerEulerModel::updatePhaseMassFlux(
    const std::shared_ptr<domain> domain,
    label phaseIndex,
    bool includeRhieChow)
{
    updatePhaseMassFluxInterior_(domain, phaseIndex, includeRhieChow);

    for (const interface* interf : domain->interfacesRef())
    {
        if (interf->isFluidSolidType())
        {
            continue;
        }
        if (interf->isInternal())
        {
            updatePhaseMassFluxSideParts_(
                domain, phaseIndex, interf->masterInfoRef().currentPartVec_);
            updatePhaseMassFluxSideParts_(
                domain, phaseIndex, interf->slaveInfoRef().currentPartVec_);
        }
        else
        {
            updatePhaseMassFluxSideParts_(
                domain,
                phaseIndex,
                interf->interfaceSideInfoPtr(domain->index())
                    ->currentPartVec_);
        }
    }

    for (label iBoundary = 0;
         iBoundary < domain->zonePtr()->nBoundaries();
         ++iBoundary)
    {
        updatePhaseMassFluxSideParts_(
            domain,
            phaseIndex,
            domain->zonePtr()->boundaryRef(iBoundary).parts());
    }
    updatePhaseFlowReversalFlags_(domain, phaseIndex);
}

void eulerEulerModel::updatePhaseMassFluxInterior_(
    const std::shared_ptr<domain> domain,
    label phaseIndex,
    bool includeRhieChow)
{
    const auto& mesh = meshRef();
    const auto& metaData = mesh.metaDataRef();
    const auto& bulkData = mesh.bulkDataRef();

    auto& totalFluxField = mDotRef(phaseIndex).stkFieldRef();
    auto& intrinsicFluxField = intrinsicMDotRef(phaseIndex).stkFieldRef();
    const auto& velocityField = URef(phaseIndex).stkFieldRef();
    const auto& alphaField = alphaRef(phaseIndex).stkFieldRef();
    const auto& rhoField = rhoRef(phaseIndex).stkFieldRef();
    const auto& pressureField = pRef().stkFieldRef();
    const auto& pressureGradientField = pRef().gradRef().stkFieldRef();
    const auto& duField = duRef(phaseIndex).stkFieldRef();
    const auto& coordinatesField = *metaData.get_field<scalar>(
        stk::topology::NODE_RANK, mesh::coordinates_ID);

    const auto selection =
        metaData.universal_part() &
        stk::mesh::selectUnion(domain->zonePtr()->interiorParts());
    const auto& buckets =
        bulkData.get_buckets(stk::topology::ELEMENT_RANK, selection);

    std::vector<scalar> shapeFunction;
    std::vector<scalar> coordinates;
    std::vector<scalar> velocity;
    std::vector<scalar> pressureGradient;
    std::vector<scalar> du;
    std::vector<scalar> alpha;
    std::vector<scalar> rho;
    std::vector<scalar> pressure;
    std::vector<scalar> areaVector;
    std::vector<scalar> dndx;
    std::vector<scalar> deriv;
    std::vector<scalar> detJ;

    for (const stk::mesh::Bucket* bucketPtr : buckets)
    {
        const auto& bucket = *bucketPtr;
        MasterElement* meSCS =
            MasterElementRepo::get_surface_master_element(bucket.topology());
        const label nodesPerElement = meSCS->nodesPerElement_;
        const label numScsIp = meSCS->numIntPoints_;

        shapeFunction.resize(numScsIp * nodesPerElement);
        coordinates.resize(nodesPerElement * SPATIAL_DIM);
        velocity.resize(nodesPerElement * SPATIAL_DIM);
        pressureGradient.resize(nodesPerElement * SPATIAL_DIM);
        du.resize(nodesPerElement * SPATIAL_DIM);
        alpha.resize(nodesPerElement);
        rho.resize(nodesPerElement);
        pressure.resize(nodesPerElement);
        areaVector.resize(numScsIp * SPATIAL_DIM);
        dndx.resize(numScsIp * nodesPerElement * SPATIAL_DIM);
        deriv.resize(numScsIp * nodesPerElement * SPATIAL_DIM);
        detJ.resize(numScsIp);

        if (URef(phaseIndex).isShifted())
        {
            meSCS->shifted_shape_fcn(shapeFunction.data());
        }
        else
        {
            meSCS->shape_fcn(shapeFunction.data());
        }

        for (stk::mesh::Bucket::size_type elementIndex = 0;
             elementIndex < bucket.size();
             ++elementIndex)
        {
            const auto* nodes = bucket.begin_nodes(elementIndex);
            for (label nodeIndex = 0; nodeIndex < nodesPerElement;
                 ++nodeIndex)
            {
                const auto node = nodes[nodeIndex];
                const scalar* nodeCoordinates =
                    stk::mesh::field_data(coordinatesField, node);
                const scalar* nodeVelocity =
                    stk::mesh::field_data(velocityField, node);
                const scalar* nodePressureGradient =
                    stk::mesh::field_data(pressureGradientField, node);
                const scalar* nodeDu = stk::mesh::field_data(duField, node);
                alpha[nodeIndex] =
                    *stk::mesh::field_data(alphaField, node);
                rho[nodeIndex] = *stk::mesh::field_data(rhoField, node);
                pressure[nodeIndex] =
                    *stk::mesh::field_data(pressureField, node);
                for (label component = 0; component < SPATIAL_DIM;
                     ++component)
                {
                    const label offset =
                        SPATIAL_DIM * nodeIndex + component;
                    coordinates[offset] = nodeCoordinates[component];
                    velocity[offset] = nodeVelocity[component];
                    pressureGradient[offset] =
                        nodePressureGradient[component];
                    du[offset] = nodeDu[component];
                }
            }

            scalar scsError = 0.0;
            meSCS->determinant(
                1, coordinates.data(), areaVector.data(), &scsError);
            meSCS->grad_op(1,
                           coordinates.data(),
                           dndx.data(),
                           deriv.data(),
                           detJ.data(),
                           &scsError);

            scalar* totalFlux = stk::mesh::field_data(
                totalFluxField, bucket, elementIndex);
            scalar* intrinsicFlux = stk::mesh::field_data(
                intrinsicFluxField, bucket, elementIndex);
            for (label ip = 0; ip < numScsIp; ++ip)
            {
                scalar alphaIp = 0.0;
                scalar rhoIp = 0.0;
                scalar normalVelocity = 0.0;
                for (label component = 0; component < SPATIAL_DIM;
                     ++component)
                {
                    scalar velocityIp = 0.0;
                    scalar gradientIp = 0.0;
                    scalar compactGradient = 0.0;
                    scalar duIp = 0.0;
                    for (label nodeIndex = 0;
                         nodeIndex < nodesPerElement;
                         ++nodeIndex)
                    {
                        const scalar shape =
                            shapeFunction[ip * nodesPerElement + nodeIndex];
                        const label offset =
                            SPATIAL_DIM * nodeIndex + component;
                        velocityIp += shape * velocity[offset];
                        gradientIp += shape * pressureGradient[offset];
                        duIp += shape * du[offset];
                        compactGradient +=
                            dndx[(ip * nodesPerElement + nodeIndex) *
                                      SPATIAL_DIM +
                                  component] *
                            pressure[nodeIndex];
                    }
                    const scalar area =
                        areaVector[ip * SPATIAL_DIM + component];
                    normalVelocity += velocityIp * area;
                    if (includeRhieChow)
                    {
                        normalVelocity -=
                            duIp * (compactGradient - gradientIp) * area;
                    }
                }

                for (label nodeIndex = 0;
                     nodeIndex < nodesPerElement;
                     ++nodeIndex)
                {
                    const scalar shape =
                        shapeFunction[ip * nodesPerElement + nodeIndex];
                    alphaIp += shape * alpha[nodeIndex];
                    rhoIp += shape * rho[nodeIndex];
                }
                intrinsicFlux[ip] = rhoIp * normalVelocity;
                totalFlux[ip] = alphaIp * intrinsicFlux[ip];
            }
        }
    }
}

void eulerEulerModel::updatePhaseMassFluxSideParts_(
    const std::shared_ptr<domain> domain,
    label phaseIndex,
    const stk::mesh::PartVector& parts)
{
    if (parts.empty())
    {
        return;
    }

    auto& totalFlux = mDotRef(phaseIndex);
    auto& intrinsicFlux = intrinsicMDotRef(phaseIndex);
    if (!totalFlux.sideFieldPtr() || !intrinsicFlux.sideFieldPtr())
    {
        return;
    }

    const auto& mesh = meshRef();
    const auto& metaData = mesh.metaDataRef();
    const auto& bulkData = mesh.bulkDataRef();
    auto& totalSideField = totalFlux.sideFieldRef().stkFieldRef();
    auto& intrinsicSideField = intrinsicFlux.sideFieldRef().stkFieldRef();
    const auto& velocityField = URef(phaseIndex).stkFieldRef();
    const auto& alphaField = alphaRef(phaseIndex).stkFieldRef();
    const auto& rhoField = rhoRef(phaseIndex).stkFieldRef();
    const auto* sideVelocityField =
        URef(phaseIndex).sideFieldPtr()
            ? URef(phaseIndex).sideFieldRef().stkFieldPtr()
            : nullptr;
    const auto* sideAlphaField =
        alphaRef(phaseIndex).sideFieldPtr()
            ? alphaRef(phaseIndex).sideFieldRef().stkFieldPtr()
            : nullptr;
    const auto& areaField = *metaData.get_field<scalar>(
        metaData.side_rank(), mesh::exposed_area_vector_ID);

    const auto selection =
        metaData.universal_part() & stk::mesh::selectUnion(parts);
    const auto& buckets =
        bulkData.get_buckets(metaData.side_rank(), selection);
    std::vector<scalar> shapeFunction;
    std::vector<scalar> nodalVelocity;
    std::vector<scalar> nodalAlpha;
    std::vector<scalar> nodalRho;

    for (const stk::mesh::Bucket* bucketPtr : buckets)
    {
        const auto& bucket = *bucketPtr;
        MasterElement* meFC =
            MasterElementRepo::get_surface_master_element(bucket.topology());
        const label nodesPerSide = bucket.topology().num_nodes();
        const label numScsIp = meFC->numIntPoints_;
        shapeFunction.resize(numScsIp * nodesPerSide);
        nodalVelocity.resize(nodesPerSide * SPATIAL_DIM);
        nodalAlpha.resize(nodesPerSide);
        nodalRho.resize(nodesPerSide);
        if (URef(phaseIndex).isShifted())
        {
            meFC->shifted_shape_fcn(shapeFunction.data());
        }
        else
        {
            meFC->shape_fcn(shapeFunction.data());
        }

        for (stk::mesh::Bucket::size_type sideIndex = 0;
             sideIndex < bucket.size();
             ++sideIndex)
        {
            const auto side = bucket[sideIndex];
            scalar* total = stk::mesh::field_data(totalSideField, side);
            scalar* intrinsic =
                stk::mesh::field_data(intrinsicSideField, side);
            if (!total || !intrinsic)
            {
                continue;
            }

            const auto* nodes = bucket.begin_nodes(sideIndex);
            for (label nodeIndex = 0; nodeIndex < nodesPerSide; ++nodeIndex)
            {
                const scalar* nodeVelocity =
                    stk::mesh::field_data(velocityField, nodes[nodeIndex]);
                nodalAlpha[nodeIndex] =
                    *stk::mesh::field_data(alphaField, nodes[nodeIndex]);
                nodalRho[nodeIndex] =
                    *stk::mesh::field_data(rhoField, nodes[nodeIndex]);
                for (label component = 0; component < SPATIAL_DIM;
                     ++component)
                {
                    nodalVelocity[SPATIAL_DIM * nodeIndex + component] =
                        nodeVelocity[component];
                }
            }

            const scalar* specifiedVelocity =
                sideVelocityField
                    ? stk::mesh::field_data(*sideVelocityField, side)
                    : nullptr;
            const scalar* specifiedAlpha =
                sideAlphaField ? stk::mesh::field_data(*sideAlphaField, side)
                               : nullptr;
            const scalar* area = stk::mesh::field_data(areaField, side);
            for (label ip = 0; ip < numScsIp; ++ip)
            {
                scalar rhoIp = 0.0;
                scalar nodalAlphaIp = 0.0;
                std::array<scalar, SPATIAL_DIM> velocityIp{};
                for (label nodeIndex = 0; nodeIndex < nodesPerSide;
                     ++nodeIndex)
                {
                    const scalar shape =
                        shapeFunction[ip * nodesPerSide + nodeIndex];
                    rhoIp += shape * nodalRho[nodeIndex];
                    nodalAlphaIp += shape * nodalAlpha[nodeIndex];
                    if (!specifiedVelocity)
                    {
                        for (label component = 0;
                             component < SPATIAL_DIM;
                             ++component)
                        {
                            velocityIp[component] +=
                                shape *
                                nodalVelocity[SPATIAL_DIM * nodeIndex +
                                              component];
                        }
                    }
                }
                if (specifiedVelocity)
                {
                    for (label component = 0; component < SPATIAL_DIM;
                         ++component)
                    {
                        velocityIp[component] =
                            specifiedVelocity[ip * SPATIAL_DIM + component];
                    }
                }

                scalar normalVelocity = 0.0;
                for (label component = 0; component < SPATIAL_DIM;
                     ++component)
                {
                    normalVelocity +=
                        velocityIp[component] *
                        area[ip * SPATIAL_DIM + component];
                }
                // A side alpha value on an inlet/opening is an entrainment
                // value, not a permanent boundary value. Outflow must carry
                // the interior phase fraction; using the prescribed value on
                // outflow can differ from the nodal alpha by six orders of
                // magnitude near a pure-phase boundary.
                const scalar alphaIp =
                    specifiedAlpha && normalVelocity < 0.0
                        ? specifiedAlpha[ip]
                        : nodalAlphaIp;
                intrinsic[ip] = rhoIp * normalVelocity;
                total[ip] = alphaIp * intrinsic[ip];
            }
        }
    }
}

void eulerEulerModel::updatePhaseFlowReversalFlags_(
    const std::shared_ptr<domain> domain,
    label phaseIndex)
{
    auto* reversalField = URef(phaseIndex).reversalFlagPtr();
    auto* fluxSideField = mDotRef(phaseIndex).sideFieldPtr();
    if (!reversalField || !fluxSideField)
    {
        return;
    }

    const auto& metaData = meshRef().metaDataRef();
    const auto& bulkData = meshRef().bulkDataRef();
    auto& reversalSTKField = reversalField->stkFieldRef();
    auto& fluxSTKField = fluxSideField->stkFieldRef();

    for (label iBoundary = 0;
         iBoundary < domain->zonePtr()->nBoundaries();
         ++iBoundary)
    {
        const auto& boundary = domain->zonePtr()->boundaryRef(iBoundary);
        const auto type = boundary.type();
        if (type != boundaryPhysicalType::inlet &&
            type != boundaryPhysicalType::outlet &&
            type != boundaryPhysicalType::opening)
        {
            continue;
        }

        const auto selection =
            metaData.universal_part() &
            stk::mesh::selectUnion(boundary.parts());
        const auto& buckets =
            bulkData.get_buckets(metaData.side_rank(), selection);
        for (const stk::mesh::Bucket* bucketPtr : buckets)
        {
            const auto& bucket = *bucketPtr;
            MasterElement* meFC =
                MasterElementRepo::get_surface_master_element(
                    bucket.topology());
            const label numScsIp = meFC->numIntPoints_;
            for (stk::mesh::Bucket::size_type sideIndex = 0;
                 sideIndex < bucket.size();
                 ++sideIndex)
            {
                const auto side = bucket[sideIndex];
                scalar* flux = stk::mesh::field_data(fluxSTKField, side);
                label* reversed =
                    stk::mesh::field_data(reversalSTKField, side);
                if (!flux || !reversed)
                {
                    continue;
                }

                if (type == boundaryPhysicalType::outlet)
                {
                    scalar sideFlux = 0.0;
                    for (label ip = 0; ip < numScsIp; ++ip)
                    {
                        sideFlux += flux[ip];
                        flux[ip] = std::max(flux[ip], 0.0);
                    }
                    const label flag = sideFlux > 0.0 ? 0 : 1;
                    for (label ip = 0; ip < numScsIp; ++ip)
                    {
                        reversed[ip] = flag;
                        if (flag == 1)
                        {
                            flux[ip] = 0.0;
                        }
                    }
                }
                else if (type == boundaryPhysicalType::opening)
                {
                    for (label ip = 0; ip < numScsIp; ++ip)
                    {
                        reversed[ip] = flux[ip] < 0.0 ? 1 : 0;
                    }
                }
                else
                {
                    // Positive flux is reversed flow at an inlet, whose
                    // outward face normal gives normal inflow a negative sign.
                    for (label ip = 0; ip < numScsIp; ++ip)
                    {
                        reversed[ip] = flux[ip] > 0.0 ? 1 : 0;
                    }
                }
            }
        }
    }
}

void eulerEulerModel::updatePhaseFluxDivergence(
    const std::shared_ptr<domain> domain,
    label phaseIndex)
{
    auto accumulate = [&](massFlowRate& flux)
    {
        auto& divergence = flux.divRef();
        divergence.setToValue({0.0}, domain->zonePtr()->interiorParts());

        const auto& mesh = meshRef();
        const auto& metaData = mesh.metaDataRef();
        const auto& bulkData = mesh.bulkDataRef();
        auto& divergenceField = divergence.stkFieldRef();
        const auto& fluxField = flux.stkFieldRef();
        const auto elementSelection =
            metaData.universal_part() &
            stk::mesh::selectUnion(domain->zonePtr()->interiorParts());
        const auto& elementBuckets = bulkData.get_buckets(
            stk::topology::ELEMENT_RANK, elementSelection);

        for (const stk::mesh::Bucket* bucketPtr : elementBuckets)
        {
            const auto& bucket = *bucketPtr;
            MasterElement* meSCS = MasterElementRepo::get_surface_master_element(
                bucket.topology());
            const label numScsIp = meSCS->numIntPoints_;
            const label* adjacentNodes = meSCS->adjacentNodes();
            for (stk::mesh::Bucket::size_type elementIndex = 0;
                 elementIndex < bucket.size();
                 ++elementIndex)
            {
                const scalar* elementFlux =
                    stk::mesh::field_data(fluxField, bucket, elementIndex);
                const auto* nodes = bucket.begin_nodes(elementIndex);
                for (label ip = 0; ip < numScsIp; ++ip)
                {
                    scalar* left = stk::mesh::field_data(
                        divergenceField, nodes[adjacentNodes[2 * ip]]);
                    scalar* right = stk::mesh::field_data(
                        divergenceField, nodes[adjacentNodes[2 * ip + 1]]);
                    *left += elementFlux[ip];
                    *right -= elementFlux[ip];
                }
            }
        }

        if (flux.sideFieldPtr())
        {
            const auto& sideFluxField = flux.sideFieldRef().stkFieldRef();
            stk::mesh::PartVector sideParts;
            for (const interface* interf : domain->interfacesRef())
            {
                if (interf->isFluidSolidType())
                {
                    continue;
                }
                if (interf->isInternal())
                {
                    const auto& masterParts =
                        interf->masterInfoRef().currentPartVec_;
                    sideParts.insert(sideParts.end(),
                                     masterParts.begin(),
                                     masterParts.end());
                    const auto& slaveParts =
                        interf->slaveInfoRef().currentPartVec_;
                    sideParts.insert(sideParts.end(),
                                     slaveParts.begin(),
                                     slaveParts.end());
                }
                else
                {
                    const auto& interfaceParts =
                        interf->interfaceSideInfoPtr(domain->index())
                            ->currentPartVec_;
                    sideParts.insert(sideParts.end(),
                                     interfaceParts.begin(),
                                     interfaceParts.end());
                }
            }
            for (label iBoundary = 0;
                 iBoundary < domain->zonePtr()->nBoundaries();
                 ++iBoundary)
            {
                const auto& boundaryParts =
                    domain->zonePtr()->boundaryRef(iBoundary).parts();
                sideParts.insert(sideParts.end(),
                                 boundaryParts.begin(),
                                 boundaryParts.end());
            }

            if (!sideParts.empty())
            {
                const auto sideSelection =
                    metaData.universal_part() &
                    stk::mesh::selectUnion(sideParts);
                const auto& sideBuckets = bulkData.get_buckets(
                    metaData.side_rank(), sideSelection);
                for (const stk::mesh::Bucket* bucketPtr : sideBuckets)
                {
                    const auto& bucket = *bucketPtr;
                    MasterElement* meFC =
                        MasterElementRepo::get_surface_master_element(
                            bucket.topology());
                    const label numScsIp = meFC->numIntPoints_;
                    const label* ipNodeMap = meFC->ipNodeMap();
                    for (stk::mesh::Bucket::size_type sideIndex = 0;
                         sideIndex < bucket.size();
                         ++sideIndex)
                    {
                        const auto side = bucket[sideIndex];
                        const scalar* sideFlux =
                            stk::mesh::field_data(sideFluxField, side);
                        if (!sideFlux)
                        {
                            continue;
                        }
                        const auto* nodes = bucket.begin_nodes(sideIndex);
                        for (label ip = 0; ip < numScsIp; ++ip)
                        {
                            scalar* nodeDivergence = stk::mesh::field_data(
                                divergenceField, nodes[ipNodeMap[ip]]);
                            *nodeDivergence += sideFlux[ip];
                        }
                    }
                }
            }
        }

        for (const interface* interf : domain->interfacesRef())
        {
            if (interf->isConformalTreatment() &&
                interf->isMasterZone(domain->index()))
            {
                divergence.transfer(
                    interf->index(), dataTransferType::move, true);
            }
        }
        divergence.synchronizeGhostedEntities(domain->index());
    };

    accumulate(intrinsicMDotRef(phaseIndex));
    accumulate(mDotRef(phaseIndex));
}


void eulerEulerModel::updatePressure(const std::shared_ptr<domain> domain)
{
    // raw update: pressure::update handles staticPressure boundaries itself
    // and is a no-op for totalPressure ("model-related", filled below)
    fieldBroker::updatePressure(domain);

    updatePhasicPressureSideFields_(domain);
}

void eulerEulerModel::updatePhasicPressureSideFields_(
    const std::shared_ptr<domain> domain)
{
    const zone* zonePtr = meshRef().zonePtr(domain->index());
    for (label iBoundary = 0; iBoundary < zonePtr->nBoundaries(); ++iBoundary)
    {
        const boundary* boundaryPtr = zonePtr->boundaryPtr(iBoundary);
        if (boundaryPtr->type() != boundaryPhysicalType::opening)
        {
            continue;
        }
        const auto& bcType =
            pRef().boundaryConditionRef(domain->index(), iBoundary).type();
        if (bcType != boundaryConditionType::totalPressure)
        {
            continue;
        }
        updatePhasicPressureSideFieldOpening_(domain, boundaryPtr);
    }
}

void eulerEulerModel::updatePhasicPressureSideFieldOpening_(
    const std::shared_ptr<domain> domain,
    const boundary* boundary)
{
    auto& bc = pRef().boundaryConditionRef(domain->index(), boundary->index());
    auto& data = bc.data<1>("value");

    if (data.type() == inputDataType::null)
    {
        return;
    }
    STK_ThrowRequireMsg(
        data.type() == inputDataType::constant ||
            data.type() == inputDataType::timeTable,
        "Euler-Euler opening `"
            << boundary->name()
            << "`: only constant or time-table opening pressure is supported");

    const scalar openValue =
        data.type() == inputDataType::constant
            ? *data.value()
            : data.interpolate(this->meshRef().controlsRef().time)[0];

    const auto phaseIndices = activePhaseIndices(domain);
    if (phaseIndices.empty())
    {
        return;
    }

    // A total-pressure condition is defined with one velocity, flux, and
    // density. OpenFOAM's multiphaseEuler cases make that phase explicit
    // (for the bubble column: U.air, phi.air, rho.air). Retain the historical
    // mixture reconstruction when no reference phase was requested.
    label pressurePhaseIndex = -1;
    if (bc.isInputDataAdded("reference_phase"))
    {
        const std::string phaseName = bc.rawStringValue("reference_phase");
        pressurePhaseIndex = realmRef().simulationRef().materialIndex(phaseName);
        STK_ThrowRequireMsg(
            std::find(phaseIndices.begin(),
                      phaseIndices.end(),
                      pressurePhaseIndex) != phaseIndices.end(),
            "Euler-Euler opening `"
                << boundary->name() << "`: reference phase `" << phaseName
                << "` is not active in domain `" << domain->name() << "`");
    }

    // Opening direction. Every phase registers the same
    // `normal_to_boundary_condition` field at an opening, so the first phase
    // that has one is representative.
    const STKScalarField* sideFlowDirectionSTKFieldPtr = nullptr;
    for (const label phaseIndex : phaseIndices)
    {
        if (URef(phaseIndex).hasSideFlowDirectionFields())
        {
            sideFlowDirectionSTKFieldPtr =
                URef(phaseIndex).sideFlowDirectionFieldRef().stkFieldPtr();
            break;
        }
    }

    auto& sidePSTKFieldRef = pRef().sideFieldRef().stkFieldRef();

    auto& mesh = this->meshRef();
    const stk::mesh::BulkData& bulkData = mesh.bulkDataRef();
    const stk::mesh::MetaData& metaData = mesh.metaDataRef();

    const auto& exposedAreaVecSTKFieldRef = *metaData.get_field<scalar>(
        metaData.side_rank(), this->getExposedAreaVectorID_(domain));

    const bool isShifted = pRef().isShifted();

    std::vector<stk::topology> parentTopo;
    std::vector<scalar> ws_rhoPressure; // nodesPerSide
    std::vector<scalar> ws_UPressure;   // nodesPerSide*SPATIAL_DIM
    std::vector<scalar> ws_shape;       // numScsBip*nodesPerSide
    std::vector<scalar> pressureSideFlux; // numScsBip

    stk::mesh::Selector selAllSides =
        metaData.universal_part() & stk::mesh::selectUnion(boundary->parts());
    stk::mesh::BucketVector const& sideBuckets =
        bulkData.get_buckets(metaData.side_rank(), selAllSides);

    for (const stk::mesh::Bucket* sideBucketPtr : sideBuckets)
    {
        const stk::mesh::Bucket& sideBucket = *sideBucketPtr;

        MasterElement* meFC = MasterElementRepo::get_surface_master_element(
            sideBucket.topology());
        const label nodesPerSide = meFC->nodesPerElement_;
        const label numScsBip = meFC->numIntPoints_;

        ws_rhoPressure.resize(nodesPerSide);
        ws_UPressure.resize(nodesPerSide * SPATIAL_DIM);
        ws_shape.resize(numScsBip * nodesPerSide);
        pressureSideFlux.resize(numScsBip);

        if (isShifted)
        {
            meFC->shifted_shape_fcn(ws_shape.data());
        }
        else
        {
            meFC->shape_fcn(ws_shape.data());
        }

        for (stk::mesh::Bucket::size_type iSide = 0; iSide < sideBucket.size();
             ++iSide)
        {
            stk::mesh::Entity side = sideBucket[iSide];

            scalar* pbc = stk::mesh::field_data(sidePSTKFieldRef, side);
            if (!pbc)
            {
                continue;
            }
            const scalar* areaVec =
                stk::mesh::field_data(exposedAreaVecSTKFieldRef, side);
            const scalar* dir =
                sideFlowDirectionSTKFieldPtr
                    ? stk::mesh::field_data(*sideFlowDirectionSTKFieldPtr, side)
                    : nullptr;

            // The reference-phase flux determines inlet/outlet switching.
            // Without an explicit phase, use the mixture mass flux as before.
            std::fill(pressureSideFlux.begin(), pressureSideFlux.end(), 0.0);
            for (const label phaseIndex : phaseIndices)
            {
                if (pressurePhaseIndex >= 0 &&
                    phaseIndex != pressurePhaseIndex)
                {
                    continue;
                }
                if (!mDotRef(phaseIndex).sideFieldPtr())
                {
                    continue;
                }
                const scalar* phaseFlux = stk::mesh::field_data(
                    mDotRef(phaseIndex).sideFieldRef().stkFieldRef(), side);
                if (!phaseFlux)
                {
                    continue;
                }
                for (label ip = 0; ip < numScsBip; ++ip)
                {
                    pressureSideFlux[ip] += phaseFlux[ip];
                }
            }

            // Nodal density and velocity used by the kinetic-pressure term.
            stk::mesh::Entity const* sideNodeRels = bulkData.begin_nodes(side);
            for (label ni = 0; ni < nodesPerSide; ++ni)
            {
                stk::mesh::Entity node = sideNodeRels[ni];

                if (pressurePhaseIndex >= 0)
                {
                    ws_rhoPressure[ni] = *stk::mesh::field_data(
                        rhoRef(pressurePhaseIndex).stkFieldRef(), node);
                    const scalar* U = stk::mesh::field_data(
                        URef(pressurePhaseIndex).stkFieldRef(), node);
                    for (label j = 0; j < SPATIAL_DIM; ++j)
                    {
                        ws_UPressure[ni * SPATIAL_DIM + j] = U[j];
                    }
                    continue;
                }

                scalar rhoMix = 0.0;
                scalar momentum[SPATIAL_DIM] = {};
                for (const label phaseIndex : phaseIndices)
                {
                    const scalar alpha = *stk::mesh::field_data(
                        alphaRef(phaseIndex).stkFieldRef(), node);
                    const scalar rho = *stk::mesh::field_data(
                        rhoRef(phaseIndex).stkFieldRef(), node);
                    const scalar* U = stk::mesh::field_data(
                        URef(phaseIndex).stkFieldRef(), node);
                    const scalar alphaRho = alpha * rho;
                    rhoMix += alphaRho;
                    for (label j = 0; j < SPATIAL_DIM; ++j)
                    {
                        momentum[j] += alphaRho * U[j];
                    }
                }
                ws_rhoPressure[ni] = rhoMix;
                for (label j = 0; j < SPATIAL_DIM; ++j)
                {
                    ws_UPressure[ni * SPATIAL_DIM + j] =
                        momentum[j] / (rhoMix + SMALL);
                }
            }

            for (label ip = 0; ip < numScsBip; ++ip)
            {
                if (pressureSideFlux[ip] > 0.0)
                {
                    // outflow: p = p_spec
                    pbc[ip] = openValue;
                    continue;
                }

                // inflow: p = p0 - 1/2 rho_m U^2
                const label faceOffSet = ip * SPATIAL_DIM;
                const label offSetSF_face = ip * nodesPerSide;

                scalar asq = 0.0;
                for (label j = 0; j < SPATIAL_DIM; ++j)
                {
                    const scalar axj = areaVec[faceOffSet + j];
                    asq += axj * axj;
                }
                const scalar amag = std::sqrt(asq);

                scalar rhoBip = 0.0;
                scalar uBip[SPATIAL_DIM] = {};
                for (label ic = 0; ic < nodesPerSide; ++ic)
                {
                    const scalar r = ws_shape[offSetSF_face + ic];
                    rhoBip += r * ws_rhoPressure[ic];
                    for (label j = 0; j < SPATIAL_DIM; ++j)
                    {
                        uBip[j] +=
                            r * ws_UPressure[ic * SPATIAL_DIM + j];
                    }
                }

                scalar num = 0.0;
                scalar den = 0.0;
                scalar d2 = 0.0;
                for (label i = 0; i < SPATIAL_DIM; ++i)
                {
                    const scalar nxi = areaVec[faceOffSet + i] / amag;
                    // Opening directions point into the domain; without a
                    // direction field fall back to the inward face normal.
                    const scalar d = dir ? dir[ip * SPATIAL_DIM + i] : -nxi;
                    num += uBip[i] * nxi;
                    den += d * nxi;
                    d2 += d * d;
                }

                const scalar Iu_ipI = num / (den + SMALL);
                pbc[ip] = openValue - 0.5 * rhoBip * Iu_ipI * Iu_ipI * d2;
            }
        }
    }

    // interpolate to the node-side field and, if requested, to the nodes
    pRef().nodeSideFieldRef().interpolate(
        pRef().sideFieldRef(), domain->index(), boundary->index());

    if (pRef().correctedBoundaryNodeValues())
    {
        pRef().correctBoundaryNodes(domain->index(), boundary->index());
    }
}


stk::mesh::PartVector eulerEulerModel::velocityDirichletParts(
    const std::shared_ptr<domain> domain,
    label phaseIndex) const
{
    stk::mesh::PartVector parts;
    const zone* zonePtr = domain->zonePtr();
    for (label iBoundary = 0; iBoundary < zonePtr->nBoundaries(); ++iBoundary)
    {
        const auto& boundaryRef = zonePtr->boundaryRef(iBoundary);
        const auto bcType = URef(phaseIndex)
                                .boundaryConditionRef(domain->index(),
                                                      iBoundary)
                                .type();
        const auto physical = boundaryRef.type();
        const bool prescribed =
            (physical == boundaryPhysicalType::inlet &&
             bcType == boundaryConditionType::specifiedValue) ||
            (physical == boundaryPhysicalType::wall &&
             (bcType == boundaryConditionType::noSlip ||
              bcType == boundaryConditionType::specifiedValue));
        if (!prescribed)
        {
            continue;
        }
        for (auto* part : boundaryRef.parts())
        {
            parts.push_back(part);
        }
    }
    return parts;
}

void eulerEulerModel::applyPhaseVelocityDirichlet(
    const std::shared_ptr<domain> domain,
    label phaseIndex)
{
    auto& U = URef(phaseIndex);
    if (U.nodeSideFieldPtr() == nullptr)
    {
        return;
    }
    const zone* zonePtr = domain->zonePtr();
    bool wrote = false;
    for (label iBoundary = 0; iBoundary < zonePtr->nBoundaries(); ++iBoundary)
    {
        const auto& boundaryRef = zonePtr->boundaryRef(iBoundary);
        const auto bcType =
            U.boundaryConditionRef(domain->index(), iBoundary).type();
        const auto physical = boundaryRef.type();
        const bool prescribed =
            (physical == boundaryPhysicalType::inlet &&
             bcType == boundaryConditionType::specifiedValue) ||
            (physical == boundaryPhysicalType::wall &&
             (bcType == boundaryConditionType::noSlip ||
              bcType == boundaryConditionType::specifiedValue));
        if (!prescribed ||
            !U.nodeSideFieldRef().definedOnBoundary(domain->index(), iBoundary))
        {
            continue;
        }
        U.correctBoundaryNodes(domain->index(), iBoundary);
        wrote = true;
    }
    if (wrote)
    {
        U.synchronizeGhostedEntities(domain->index());
    }
}


} // namespace accel
