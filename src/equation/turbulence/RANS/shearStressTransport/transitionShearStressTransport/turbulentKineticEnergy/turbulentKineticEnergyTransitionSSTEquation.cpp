// File       : turbulentKineticEnergyTransitionSSTEquation.cpp
// Created    : Mon Jan 14 2025
// Author     : Adam Fares
// Description:
// Copyright 2025 CCFNUM HSLU T&A. All Rights Reserved.

#include "turbulentKineticEnergyTransitionSSTEquation.h"
#include "realm.h"
#include "scaling.h"

namespace accel
{

turbulentKineticEnergyTransitionSSTEquation::
    turbulentKineticEnergyTransitionSSTEquation(
        realm* realm,
        transitionShearStressTransportModel* model)
    : turbulentKineticEnergyEquation(realm, model), model_(model),
      assembler_(std::make_unique<Assembler>(model))
{
}

void turbulentKineticEnergyTransitionSSTEquation::setup()
{
    turbulentKineticEnergyEquation::setup();

    using Bucket = typename Assembler::Bucket;
    using BucketVector = typename Assembler::BucketVector;

    assembler_->setup(&model_->kRef(),
                      advectionDiffusion,
                      domainVector_,
                      // anonymous function to compute Gamma for tke equation:
                      [this](const domain* domain, STKScalarField& Gamma)
    {
        const auto& mesh = model_->meshRef();
        const stk::mesh::MetaData& metaData = mesh.metaDataRef();
        const stk::mesh::BulkData& bulkData = mesh.bulkDataRef();

        // define selector for domain
        const stk::mesh::PartVector& partVec =
            domain->zonePtr()->interiorParts();
        stk::mesh::Selector selAllNodes =
            metaData.universal_part() & stk::mesh::selectUnion(partVec);

        const BucketVector& nodeBuckets =
            bulkData.get_buckets(stk::topology::NODE_RANK, selAllNodes);

        using stk::mesh::field_data;
        const STKScalarField& muSTKFieldRef = model_->muRef().stkFieldRef();
        const STKScalarField& mutSTKFieldRef = model_->mutRef().stkFieldRef();

        const STKScalarField& fOneBlendingSTKFieldRef =
            model_->F1Ref().stkFieldRef();
        scalar sigmaOne = model_->sigmaKOne();
        scalar sigmaTwo = model_->sigmaKTwo();

        for (size_t ib = 0; ib < nodeBuckets.size(); ib++)
        {
            const Bucket& nodeBucket = *nodeBuckets[ib];
            const Bucket::size_type nNodesPerBucket = nodeBucket.size();

            // field chunks in bucket
            scalar* Gammab = field_data(Gamma, nodeBucket);
            const scalar* mu = field_data(muSTKFieldRef, nodeBucket);
            const scalar* mut = field_data(mutSTKFieldRef, nodeBucket);
            const scalar* fOneBlend =
                field_data(fOneBlendingSTKFieldRef, nodeBucket);

            for (Bucket::size_type i = 0; i < nNodesPerBucket; i++)
            {
                const scalar blendedConstant =
                    fOneBlend[i] * sigmaOne + (1.0 - fOneBlend[i]) * sigmaTwo;
                Gammab[i] = mu[i] + mut[i] * blendedConstant;
            }
        }
    });
}

void turbulentKineticEnergyTransitionSSTEquation::initialize()
{
    turbulentKineticEnergyEquation::initialize();
}

void turbulentKineticEnergyTransitionSSTEquation::postInitialize()
{
}

void turbulentKineticEnergyTransitionSSTEquation::solve()
{
    auto ctx = linearSystem::getContext();
    ctx->zeroSystemStorage();

    // assembly
    linearSystem::simulationRef().getProfiler().push("linear_system_assembly");

    // assemble
    FOREACH_DOMAIN_PTR(assembler_->preAssemble, ctx.get());
    FOREACH_DOMAIN_PTR(assembler_->assemble, ctx.get());
    FOREACH_DOMAIN_PTR(assembler_->postAssemble, ctx.get());

    // fix system in domains where the model is not active
    assembler_->fix(this->collectInactiveInteriorParts(), {}, ctx.get(), {});

    linearSystem::simulationRef().getProfiler().pop();

    // overrides of a plain simulation act on the assembled system
    this->applyOverrides();

    // solve linear system
    if (ctx->getGraph()->isGraphMember())
    {
        linearSystem::solve();
    }
    messager::barrier();

    // correction
    // clip values in source field below `lower_clip_value` to
    // `lower_clip_value`
    static constexpr int CLIP = 1; // true
    const scalar relax_value = 1.0;
    const scalar lower_clip_value = accel::SMALL;

    auto& k = model_->kRef();
    for (const auto& domain : domainVector_)
    {
        this->template correctField_<linearSystem::BLOCKSIZE, 1, 0, CLIP>(
            domain.get(),
            ctx->coeffs().get(),
            stk::topology::NODE_RANK,
            k.stkFieldRef(),
            relax_value,
            lower_clip_value);

        // synchronize
        k.synchronizeGhostedEntities(domain->index());
    }

    // post correction

    // gradient must only be updated after solving the whole k-omega model
    // ...

    // update scale
    k.updateScale();
}

void turbulentKineticEnergyTransitionSSTEquation::applyOverrides()
{
    masking* masks = this->maskingPtr();
    if (!masks)
    {
        return;
    }

    auto ctx = linearSystem::getContext();

    auto& phi = model_->kRef().stkFieldRef();
    const auto& U = model_->URef().stkFieldRef();
    const auto& rho = model_->rhoRef().stkFieldRef();
    const auto& mu = model_->muRef().stkFieldRef();

    for (label r = 0; r < masks->regionCount(); ++r)
    {
        maskedRegion& region = masks->regionRef(r);

        // no turbulence inside a region
        const std::vector<scalar> zero(region.coveredNodes().size(), 0.0);
        assembler_->constrain(region.coveredNodes(), zero, phi, ctx.get());

        if (region.wallTreatment() != maskWallTreatment::wallFunction)
        {
            continue;
        }

        // equilibrium log-layer state of the first ring outside it
        std::vector<scalar> ring(region.ringNodes().size(), 0.0);
        for (size_t n = 0; n < region.ringNodes().size(); ++n)
        {
            const scalar uTau = region.frictionVelocity(
                n, U, rho, mu, masks->kappa(), masks->B());
            const scalar distance = std::max(region.ringDistance(n), SMALL);
            const scalar sqrtCmu = std::sqrt(masks->Cmu());
            (void)distance;
            (void)sqrtCmu;

            ring[n] = uTau * uTau / sqrtCmu;
        }

        assembler_->constrain(region.ringNodes(), ring, phi, ctx.get());
    }

    // k and omega carry no information inside a region, so their ratio must not
    // reach the neighbouring fluid nodes through the effective viscosity
    masks->setAtCovered(model_->mutRef().stkFieldRef(), 0.0);
    masks->copyAtCovered(model_->muEffRef().stkFieldRef(),
                         model_->muRef().stkFieldRef());
}

} /* namespace accel */
