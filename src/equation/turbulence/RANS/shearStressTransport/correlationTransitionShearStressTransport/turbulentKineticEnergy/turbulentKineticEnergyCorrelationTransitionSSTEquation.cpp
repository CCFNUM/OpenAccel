// File       : turbulentKineticEnergyCorrelationTransitionSSTEquation.cpp
// Created    : Sun Dec 29 2024
// Author     : Adam Fares
// Description: TKE equation for correlation-based transition SST (Menter 2015)
// Copyright (c) 2024 CCFNUM, Lucerne University of Applied Sciences and Arts.
// SPDX-License-Identifier: BSD-3-Clause

#include "turbulentKineticEnergyCorrelationTransitionSSTEquation.h"
#include "realm.h"

namespace accel
{

turbulentKineticEnergyCorrelationTransitionSSTEquation::
    turbulentKineticEnergyCorrelationTransitionSSTEquation(
        realm* realm,
        correlationTransitionShearStressTransportModel* model)
    : turbulentKineticEnergyEquation(realm, model), model_(model),
      assembler_(std::make_unique<Assembler>(model))
{
}

void turbulentKineticEnergyCorrelationTransitionSSTEquation::setup()
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

void turbulentKineticEnergyCorrelationTransitionSSTEquation::initialize()
{
    turbulentKineticEnergyEquation::initialize();
}

void turbulentKineticEnergyCorrelationTransitionSSTEquation::postInitialize()
{
}

void turbulentKineticEnergyCorrelationTransitionSSTEquation::solve()
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
            ctx->getXVector(),
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

} /* namespace accel */
