// File : transitionOnsetReynoldsNumberTransitionSSTEquation.cpp
// Created : Tue Jan 14 2025
// Author : Adam Fares
// Description:
// Copyright (c) 2025 CCFNUM, Lucerne University of Applied Sciences and Arts.
// SPDX-License-Identifier: BSD-3-Clause

#include "transitionOnsetReynoldsNumberTransitionSSTEquation.h"
#include "realm.h"

namespace accel
{

transitionOnsetReynoldsNumberTransitionSSTEquation::
    transitionOnsetReynoldsNumberTransitionSSTEquation(
        realm* realm,
        transitionShearStressTransportModel* model)
    : equation("Transition Onset Reynolds Number", true),
      linearSystem(realm->simulationRef()), model_(model),
      assembler_(std::make_unique<Assembler>(model))
{
    this->setEquationName({"TOR"});

    // set relaxation for ReTheta
    model_->ReThetaRef().setURF(
        model_->meshRef()
            .controlsRef()
            .solverRef()
            .solverControl_.basicSettings_.convergenceControl_
            .relaxationParameters_.turbulenceRelaxationFactor_);
}

void transitionOnsetReynoldsNumberTransitionSSTEquation::checkDomain(
    const std::shared_ptr<domain> domain)
{
    assert(domain->type() == domainType::fluid);
}

bool transitionOnsetReynoldsNumberTransitionSSTEquation::isConverged() const
{
    return linearSystem::isConverged();
}

void transitionOnsetReynoldsNumberTransitionSSTEquation::setup()
{
    // setup of fields (initial values and boundary conditions)
    FOREACH_DOMAIN(model_->setupTransitionOnsetReynoldsNumber);

    using Bucket = typename Assembler::Bucket;
    using BucketVector = typename Assembler::BucketVector;

    assembler_->setup(&model_->ReThetaRef(),
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

        scalar sigmaThetat = model_->sigmaThetat();

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
                Gammab[i] = (mu[i] + mut[i]) * sigmaThetat;
            }
        }
    });

    // linear solver
    linearSystem::setupSolver(this->name(), model_->meshRef());

    equation::isCreated_ = true;
}

void transitionOnsetReynoldsNumberTransitionSSTEquation::initialize()
{
    // initialization of fields (initial values and boundary conditions)
    FOREACH_DOMAIN(model_->initializeTransitionOnsetReynoldsNumber);

    // post initialization

    // 1) update gradient
    FOREACH_DOMAIN(model_->updateTransitionOnsetReynoldsNumberGradientField);

    // 2) update high-res fields
    FOREACH_DOMAIN(
        model_->updateTransitionOnsetReynoldsNumberBlendingFactorField);

    // 3) update scale
    model_->ReThetaRef().updateScale();

    equation::isInitialized_ = true;
}

void transitionOnsetReynoldsNumberTransitionSSTEquation::postInitialize()
{
}

void transitionOnsetReynoldsNumberTransitionSSTEquation::preSolve()
{
    // raw update
    FOREACH_DOMAIN(model_->updateTransitionOnsetReynoldsNumberPrevIterField);
    FOREACH_DOMAIN(model_->updateTransitionOnsetReynoldsNumber);
}

void transitionOnsetReynoldsNumberTransitionSSTEquation::solve()
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

    auto& ReTheta = model_->ReThetaRef();
    for (const auto& domain : domainVector_)
    {
        this->template correctField_<linearSystem::BLOCKSIZE, 1, 0, CLIP>(
            domain.get(),
            ctx->getXVector(),
            stk::topology::NODE_RANK,
            ReTheta.stkFieldRef(),
            relax_value);

        // synchronise
        ReTheta.synchronizeGhostedEntities(domain->index());
    }

    // Update
    // ReTheta.updateBlendingFactorField();
    ReTheta.updateScale();
}

void transitionOnsetReynoldsNumberTransitionSSTEquation::preTimeStep()
{
    FOREACH_DOMAIN(model_->updateTransitionOnsetReynoldsNumberPrevTimeField);
}

void transitionOnsetReynoldsNumberTransitionSSTEquation::setResidualScales_()
{
    auto& residual_scales = linearSystem::getContext()->getResidualScales();
    residual_scales = {1.0 / (model_->ReThetaRef().scale() + ::accel::SMALL)};
}

void transitionOnsetReynoldsNumberTransitionSSTEquation::printScales()
{
    if (messager::master())
    {
        std::cout << model_->ReThetaRef().name() << " scales:" << std::endl;
        std::cout << "\tmin: " << std::scientific << std::setprecision(8)
                  << model_->ReThetaRef().minScale() << std::endl;
        std::cout << "\tmax: " << std::scientific << std::setprecision(8)
                  << model_->ReThetaRef().maxScale() << std::endl;
        std::cout << "\tscale: " << std::scientific << std::setprecision(8)
                  << model_->ReThetaRef().scale() << std::endl
                  << std::endl;
    }
}

} /* namespace accel */
