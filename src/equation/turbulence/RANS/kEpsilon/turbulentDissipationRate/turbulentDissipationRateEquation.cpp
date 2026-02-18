// File : turbulentDissipationRateEquation.cpp
// Created : Thu Feb 22 2025 13:38:51 (+0100)
// Author : ACHRAF NAGIHI
// Description:
// Copyright (c) 2025 CCFNUM, Lucerne University of Applied Sciences and Arts.
// SPDX-License-Identifier: BSD-3-Clause

#include "turbulentDissipationRateEquation.h"
#include "kEpsilonModel.h"
#include "realm.h"

namespace accel
{

turbulentDissipationRateEquation::turbulentDissipationRateEquation(
    realm* realm,
    kEpsilonModel* model)
    : equation("Turbulent Dissipation Rate"),
      linearSystem(realm->simulationRef()), model_(model),
      assembler_(std::make_unique<Assembler>(model))
{
    this->setEquationName({"TDR"});

    // set relaxation for epsilon in kEpsilon
    model_->epsilonRef().setURF(
        model_->meshRef()
            .controlsRef()
            .solverRef()
            .solverControl_.basicSettings_.convergenceControl_
            .relaxationParameters_.turbulenceRelaxationFactor_);
}

void turbulentDissipationRateEquation::checkDomain(
    const std::shared_ptr<domain> domain)
{
    assert(domain->type() == domainType::fluid);
}

bool turbulentDissipationRateEquation::isConverged() const
{
    return linearSystem::isConverged();
}

void turbulentDissipationRateEquation::setup()
{
    // setup of fields (initial values and boundary conditions)
    FOREACH_DOMAIN(model_->setupTurbulentDissipationRate);

    using Bucket = typename Assembler::Bucket;
    using BucketVector = typename Assembler::BucketVector;

    // linear solver
    linearSystem::setupSolver(this->name(), model_->meshRef());

    assembler_->setup(&model_->epsilonRef(),
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

        scalar sigmaEpsilon = model_->sigmaEpsilon();

        for (size_t ib = 0; ib < nodeBuckets.size(); ib++)
        {
            const Bucket& nodeBucket = *nodeBuckets[ib];
            const Bucket::size_type nNodesPerBucket = nodeBucket.size();

            // field chunks in bucket
            scalar* Gammab = field_data(Gamma, nodeBucket);
            const scalar* mu = field_data(muSTKFieldRef, nodeBucket);
            const scalar* mut = field_data(mutSTKFieldRef, nodeBucket);
            for (Bucket::size_type i = 0; i < nNodesPerBucket; i++)
            {
                Gammab[i] = mu[i] + mut[i] / sigmaEpsilon;
            }
        }
    });

    equation::isCreated_ = true;
}

void turbulentDissipationRateEquation::initialize()
{
    // initialization of fields (initial values and boundary conditions)
    FOREACH_DOMAIN(model_->initializeTurbulentDissipationRate);

    // post initialization

    // 1) update gradient
    FOREACH_DOMAIN(model_->updateTurbulentDissipationRateGradientField);

    // 2) update high-res fields
    FOREACH_DOMAIN(model_->updateTurbulentDissipationRateBlendingFactorField);

    // 3) update scale
    model_->epsilonRef().updateScale();

    equation::isInitialized_ = true;
}

void turbulentDissipationRateEquation::postInitialize()
{
}

void turbulentDissipationRateEquation::preSolve()
{
    // raw update
    FOREACH_DOMAIN(model_->updateTurbulentDissipationRatePrevIterField);
    FOREACH_DOMAIN(model_->updateTurbulentDissipationRate);

    // update epsilon values on walls
    FOREACH_DOMAIN(model_->updateEpsilonAtWalls);
}

void turbulentDissipationRateEquation::solve()
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

    auto& epsilon = model_->epsilonRef();
    for (const auto& domain : domainVector_)
    {
        this->template correctField_<linearSystem::BLOCKSIZE, 1, 0, CLIP>(
            domain.get(),
            ctx->getXVector(),
            stk::topology::NODE_RANK,
            epsilon.stkFieldRef(),
            relax_value,
            lower_clip_value);

        // synchronize
        epsilon.synchronizeGhostedEntities(domain->index());
    }

    // post correction

    // gradient must only be updated after solving the whole k-epsilon model
    // ...

    // update scale
    epsilon.updateScale();
}

void turbulentDissipationRateEquation::preTimeStep()
{
    FOREACH_DOMAIN(model_->updateTurbulentDissipationRatePrevTimeField);
}

void turbulentDissipationRateEquation::setResidualScales_()
{
    auto& residual_scales = linearSystem::getContext()->getResidualScales();
    residual_scales = {1.0 / (model_->epsilonRef().scale() + ::accel::SMALL)};
}

} /* namespace accel */
