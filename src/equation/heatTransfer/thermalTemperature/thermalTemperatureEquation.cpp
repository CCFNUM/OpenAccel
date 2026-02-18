// File : thermalTemperatureEquation.cpp
// Created : Thu Feb 22 2024 08:42:10 (+0100)
// Author : Fabian Wermelinger
// Description: Thermal energy equation implementation details
// Copyright (c) 2024 CCFNUM, Lucerne University of Applied Sciences and Arts.
// SPDX-License-Identifier: BSD-3-Clause

#ifdef WITH_THERMAL_TEMPERATURE

#include "thermalTemperatureEquation.h"

namespace accel
{

thermalTemperatureEquation::thermalTemperatureEquation(realm* realm)
    : equation("Thermal Energy"), heatTransferModel(realm),
      linearSystem(realm->simulationRef()),
      assembler_(std::make_unique<Assembler>(this))
{
    this->setEquationName({"T"});

    // set relaxation factor for temperature field. Default is 1
    this->TRef().setURF(controlsRef()
                            .solverRef()
                            .solverControl_.basicSettings_.convergenceControl_
                            .relaxationParameters_.energyRelaxationFactor_);
}

bool thermalTemperatureEquation::isConverged() const
{
    return linearSystem::isConverged();
}

void thermalTemperatureEquation::setup()
{
    // setup of fields on defined domains: rho field
    // might have been already initialized over other domains
    FOREACH_DOMAIN(setupDensity);
    FOREACH_DOMAIN(setupVelocity);
    FOREACH_DOMAIN(setupMassFlowRate);
    FOREACH_DOMAIN(setupSpecificHeatCapacity);
    FOREACH_DOMAIN(setupThermalConductivity);
    FOREACH_DOMAIN(setupThermalExpansivity);
    FOREACH_DOMAIN(setupTemperature);
    FOREACH_DOMAIN_IF(setupCompressibility,
                      domain->type() == domainType::fluid &&
                          domain->isMaterialCompressible());
    FOREACH_DOMAIN(setupHeatFlowRate);

    using Bucket = typename Assembler::Bucket;
    using BucketVector = typename Assembler::BucketVector;

    // setup assembler
    assembler_->setup(&TRef(),
                      advectionDiffusion,
                      domainVector_,
                      // anonymous function to compute Gamma for temperature
                      // equation (assuming cp and lambda are constants):
                      [this](const domain* domain, STKScalarField& Gamma)
    {
        const auto& mesh = this->meshRef();
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
        const STKScalarField& lambdaEff = lambdaEffRef().stkFieldRef();
        for (size_t ib = 0; ib < nodeBuckets.size(); ib++)
        {
            const Bucket& nodeBucket = *nodeBuckets[ib];
            const Bucket::size_type nNodesPerBucket = nodeBucket.size();

            // field chunks in bucket
            scalar* Gammab = field_data(Gamma, nodeBucket);
            const scalar* lambdaEffb = field_data(lambdaEff, nodeBucket);

            for (Bucket::size_type i = 0; i < nNodesPerBucket; i++)
            {
                Gammab[i] = lambdaEffb[i];
            }
        }
    });

    // setup linear solver
    // FIXME: Consider passing mesh argument or
    // connectivity arrays passed to initialize() directly is more flexible
    // rather than this->meshRef() which is set through simulation object
    // obtained via realm in fieldBroker
    linearSystem::setupSolver(this->name(), fieldBroker::meshRef());

    equation::isCreated_ = true;
}

void thermalTemperatureEquation::initialize()
{
    // raw initialization of temperature
    FOREACH_DOMAIN(initializeTemperature);

    // 1) update temperature gradient
    FOREACH_DOMAIN(updateTemperatureGradientField);

    // 2) update high-res fields
    FOREACH_DOMAIN(updateTemperatureBlendingFactorField);

    // update velocity scale
    TRef().updateScale();

    // raw initialization of velocity if not yet: this is a situation when a
    // flow model is not enabled. The check is done in
    // nodeField::initialize(label iZone)
    FOREACH_DOMAIN(initializeVelocity);
}

void thermalTemperatureEquation::postInitialize()
{
    // property initialization: density will be initialized here for solid
    // domains, because density in fluid domains has already been initialized
    FOREACH_DOMAIN(initializeDensity);
    FOREACH_DOMAIN(initializeSpecificHeatCapacity);
    FOREACH_DOMAIN(initializeThermalConductivity);
    FOREACH_DOMAIN(initializeThermalExpansivity);
    FOREACH_DOMAIN_IF(initializeCompressibility,
                      domain->type() == domainType::fluid &&
                          domain->isMaterialCompressible());

    cpRef().updateScale();
    lambdaRef().updateScale();

    equation::isInitialized_ = true;
}

void thermalTemperatureEquation::preSolve()
{
    // temperature updates: temperature field must be updated with a raw update
    // because we are solving for it and we only intend to update its side
    // fields
    FOREACH_DOMAIN(updateTemperaturePrevIterField);
    FOREACH_DOMAIN(updateTemperature);

    // update density only for solid domains: for fluid domains, it will have
    // already been updated by the flow model
    FOREACH_DOMAIN_IF(updateDensity, domain->type() == domainType::solid);

    // other property updates
    FOREACH_DOMAIN(updateSpecificHeatCapacity);
    FOREACH_DOMAIN(updateSpecificHeatCapacityGradientField);
    FOREACH_DOMAIN(updateThermalConductivity);
    FOREACH_DOMAIN(updateEffectiveThermalConductivity); // laminar/solid
    FOREACH_DOMAIN(updateThermalExpansivity);
}

void thermalTemperatureEquation::solve()
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
    // clip values in source field to `lower_clip_value` and `upper_clip_value`
    scalar relaxCorrection = 1.0;
    static constexpr int CLIP = 1; // true
    const scalar lower_clip_value = 0;
    for (const auto& domain : domainVector_)
    {
        scalar effectiveRelaxationFactor = relaxCorrection;

        // find the ramp value (explicit relaxation) for compressible domains
        if (domain->isMaterialCompressible())
        {
            scalar iter = controlsRef().iter;
            label rampIter = 20;
            if (iter <= rampIter)
            {
                effectiveRelaxationFactor *=
                    std::max(scalar(iter) / scalar(rampIter), 0.1);
            }
        }

        this->template correctField_<linearSystem::BLOCKSIZE, 1, 0, CLIP>(
            domain.get(),
            ctx->getXVector(),
            stk::topology::NODE_RANK,
            TRef().stkFieldRef(),
            effectiveRelaxationFactor,
            lower_clip_value);

        // synchronize
        TRef().synchronizeGhostedEntities(domain->index());
    }

    // post correction

    // 1) update gradient
    FOREACH_DOMAIN(updateTemperatureGradientField);

    // 2) update high-res fields
    FOREACH_DOMAIN(updateTemperatureBlendingFactorField);

    // 3) update scale
    TRef().updateScale();
    cpRef().updateScale();
    lambdaRef().updateScale();

#ifdef HAS_INTERFACE
    // 4) Post-processed quantities
    FOREACH_DOMAIN(updateInterfaceHeatImbalance_);
#endif /* HAS_INTERFACE */
}

void thermalTemperatureEquation::preTimeStep()
{
    FOREACH_DOMAIN(updateTemperaturePrevTimeField);
    FOREACH_DOMAIN_IF(updateDensityPrevTimeField,
                      domain->type() == domainType::solid);
    FOREACH_DOMAIN(updateSpecificHeatCapacityPrevTimeField);
    FOREACH_DOMAIN(updateCompressibilityPrevTimeField);
}

void thermalTemperatureEquation::postSolve()
{
    this->reportHeatData_();

    FOREACH_DOMAIN_IF(updateCompressibility,
                      domain->type() == domainType::fluid &&
                          domain->isMaterialCompressible());
    FOREACH_DOMAIN(updateSpecificHeatCapacity);

    // calculate post-processed quantities
    FOREACH_DOMAIN_IF(updateTotalTemperatureField_,
                      domain->type() == domainType::fluid);
}

void thermalTemperatureEquation::setResidualScales_()
{
    auto& residual_scales = linearSystem::getContext()->getResidualScales();
    residual_scales = {1.0 / (TRef().scale() + ::accel::SMALL)};
}

void thermalTemperatureEquation::printScales()
{
    if (messager::master())
    {
        std::cout << this->TRef().name() << " scales:" << std::endl;
        std::cout << "\tmin: " << std::scientific << std::setprecision(8)
                  << TRef().minScale() << std::endl;
        std::cout << "\tmax: " << std::scientific << std::setprecision(8)
                  << TRef().maxScale() << std::endl;
        std::cout << "\tscale: " << std::scientific << std::setprecision(8)
                  << TRef().scale() << std::endl
                  << std::endl;

        std::cout << lambdaRef().name() << " scales:" << std::endl;
        std::cout << "\tmin: " << std::scientific << std::setprecision(8)
                  << lambdaRef().minScale() << std::endl;
        std::cout << "\tmax: " << std::scientific << std::setprecision(8)
                  << lambdaRef().maxScale() << std::endl;
        std::cout << "\tscale: " << std::scientific << std::setprecision(8)
                  << lambdaRef().scale() << std::endl
                  << std::endl;

        std::cout << cpRef().name() << " scales:" << std::endl;
        std::cout << "\tmin: " << std::scientific << std::setprecision(8)
                  << cpRef().minScale() << std::endl;
        std::cout << "\tmax: " << std::scientific << std::setprecision(8)
                  << cpRef().maxScale() << std::endl;
        std::cout << "\tscale: " << std::scientific << std::setprecision(8)
                  << cpRef().scale() << std::endl
                  << std::endl;
    }
}

} /* namespace accel */

#endif
