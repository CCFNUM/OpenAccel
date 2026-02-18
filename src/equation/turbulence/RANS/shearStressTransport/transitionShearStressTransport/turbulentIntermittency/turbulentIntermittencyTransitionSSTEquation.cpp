// File : turbulentIntermittencyTransitionSSTEquation.cpp
// Created : Tue Jan 14 2025
// Author : Adam Fares
// Description:
// Copyright (c) 2025 CCFNUM, Lucerne University of Applied Sciences and Arts.
// SPDX-License-Identifier: BSD-3-Clause

#include "turbulentIntermittencyTransitionSSTEquation.h"
#include "realm.h"

namespace accel
{

turbulentIntermittencyTransitionSSTEquation::
    turbulentIntermittencyTransitionSSTEquation(
        realm* realm,
        transitionShearStressTransportModel* model)
    : equation("Turbulent Intermittency", true),
      linearSystem(realm->simulationRef()), model_(model),
      assembler_(std::make_unique<Assembler>(model))
{
    this->setEquationName({"TI"});

    // set relaxation for gamma
    model_->gammaRef().setURF(
        model_->meshRef()
            .controlsRef()
            .solverRef()
            .solverControl_.basicSettings_.convergenceControl_
            .relaxationParameters_.turbulenceRelaxationFactor_);
}

void turbulentIntermittencyTransitionSSTEquation::checkDomain(
    const std::shared_ptr<domain> domain)
{
    assert(domain->type() == domainType::fluid);
}

bool turbulentIntermittencyTransitionSSTEquation::isConverged() const
{
    return linearSystem::isConverged();
}

void turbulentIntermittencyTransitionSSTEquation::setup()
{
    // setup of fields (initial values and boundary conditions)
    FOREACH_DOMAIN(model_->setupTurbulentIntermittency);

    using Bucket = typename Assembler::Bucket;
    using BucketVector = typename Assembler::BucketVector;

    assembler_->setup(&model_->gammaRef(),
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

        scalar sigmaF = model_->sigmaF();

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
                Gammab[i] = mu[i] + mut[i] / sigmaF;
            }
        }
    });

    // linear solver
    linearSystem::setupSolver(this->name(), model_->meshRef());

    equation::isCreated_ = true;
}

void turbulentIntermittencyTransitionSSTEquation::initialize()
{
    // initialization of fields (initial values and boundary conditions)
    FOREACH_DOMAIN(model_->initializeTurbulentIntermittency);

    // post initialization

    // 1) update gradient
    FOREACH_DOMAIN(model_->updateTurbulentIntermittencyGradientField);

    // 2) update high-res fields
    FOREACH_DOMAIN(model_->updateTurbulentIntermittencyBlendingFactorField);

    // 3) update scale
    model_->gammaRef().updateScale();

    equation::isInitialized_ = true;
}

void turbulentIntermittencyTransitionSSTEquation::postInitialize()
{
}

void turbulentIntermittencyTransitionSSTEquation::preSolve()
{
    // raw update
    FOREACH_DOMAIN(model_->updateTurbulentIntermittencyPrevIterField);
    FOREACH_DOMAIN(model_->updateTurbulentIntermittency);
}

void turbulentIntermittencyTransitionSSTEquation::solve()
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

    auto& gamma = model_->gammaRef();
    for (const auto& domain : domainVector_)
    {
        this->template correctField_<linearSystem::BLOCKSIZE, 1, 0, CLIP>(
            domain.get(),
            ctx->getXVector(),
            stk::topology::NODE_RANK,
            gamma.stkFieldRef(),
            relax_value,
            lower_clip_value);

        // synchronise
        gamma.synchronizeGhostedEntities(domain->index());
    }

    // Update
    // WARNING: gamma.updateGradientField() is deferred to SST driver equation
    // gamma.updateBlendingFactorField();
    gamma.updateScale();
}

void turbulentIntermittencyTransitionSSTEquation::preTimeStep()
{
    FOREACH_DOMAIN(model_->updateTurbulentIntermittencyPrevTimeField);
}

void turbulentIntermittencyTransitionSSTEquation::setResidualScales_()
{
    auto& residual_scales = linearSystem::getContext()->getResidualScales();
    residual_scales = {1.0 / (model_->gammaRef().scale() + ::accel::SMALL)};
}

void turbulentIntermittencyTransitionSSTEquation::applyDependencyUpdates_(
    const domain* domain,
    const stk::mesh::EntityRank entityRank,
    STKScalarField& /*stk_dst*/)
{
    assert(entityRank == stk::mesh::EntityRank::NODE_RANK);

    // if this is a fluid-solid interface, then it is necessary to transfer the
    // gamma field to the solid, as this will be used in the hybrid heat
    // transfer appraoch
#ifdef HAS_INTERFACE
    for (const interface* interf : domain->interfacesRef())
    {
        if (interf->isFluidSolidType())
        {
            model_->gammaRef().transfer(interf->index(),
                                        dataTransferType::copy,
                                        !interf->isMasterZone(domain->index()));
        }
    }
#endif /* HAS_INTERFACE */
}

void turbulentIntermittencyTransitionSSTEquation::printScales()
{
    if (messager::master())
    {
        std::cout << model_->gammaRef().name() << " scales:" << std::endl;
        std::cout << "\tmin: " << std::scientific << std::setprecision(8)
                  << model_->gammaRef().minScale() << std::endl;
        std::cout << "\tmax: " << std::scientific << std::setprecision(8)
                  << model_->gammaRef().maxScale() << std::endl;
        std::cout << "\tscale: " << std::scientific << std::setprecision(8)
                  << model_->gammaRef().scale() << std::endl
                  << std::endl;
    }
}

} /* namespace accel */
