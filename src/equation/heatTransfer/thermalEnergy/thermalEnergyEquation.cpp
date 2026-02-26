// File       : thermalEnergyEquation.cpp
// Created    : Mon Dec 01 2025 08:42:10 (+0100)
// Author     : Mhamad Mahdi Alloush
// Description: Thermal energy equation implementation details
// Copyright (c) 2025 CCFNUM, Lucerne University of Applied Sciences and Arts.
// SPDX-License-Identifier: BSD-3-Clause

#ifndef WITH_THERMAL_TEMPERATURE

#include "thermalEnergyEquation.h"

namespace accel
{

thermalEnergyEquation::thermalEnergyEquation(realm* realm)
    : equation("Thermal Energy"), heatTransferModel(realm),
      linearSystem(realm->simulationRef()),
      assembler_(std::make_unique<Assembler>(this))
{
    this->setEquationName({"h"});

    // set relaxation factor for enthalpy field. Default is 1
    this->hRef().setURF(controlsRef()
                            .solverRef()
                            .solverControl_.basicSettings_.convergenceControl_
                            .relaxationParameters_.energyRelaxationFactor_);
}

bool thermalEnergyEquation::isConverged() const
{
    return linearSystem::isConverged();
}

void thermalEnergyEquation::setup()
{
    // setup of fields on defined domains: rho field
    // might have been already initialized over other domains
    FOREACH_DOMAIN(setupDensity);
    FOREACH_DOMAIN_IF(setupVelocity, domain->type() == domainType::fluid);
    FOREACH_DOMAIN_IF(setupMassFlowRate, domain->type() == domainType::fluid);
    FOREACH_DOMAIN(setupSpecificHeatCapacity);
    FOREACH_DOMAIN(setupThermalConductivity);
    FOREACH_DOMAIN_IF(setupThermalExpansivity,
                      domain->type() == domainType::fluid &&
                          !domain->isMaterialCompressible() &&
                          domain->buoyancy_.option_ == buoyancyOption::buoyant);
    FOREACH_DOMAIN(setupTemperature);      // derived from enthalpy
    FOREACH_DOMAIN(setupSpecificEnthalpy); // primary variable
    FOREACH_DOMAIN_IF(setupSpecificTotalEnthalpy,
                      domain->type() == domainType::fluid);
    FOREACH_DOMAIN_IF(setupCompressibility,
                      domain->type() == domainType::fluid &&
                          domain->isMaterialCompressible());
    FOREACH_DOMAIN(setupHeatFlowRate);

    using Bucket = typename Assembler::Bucket;
    using BucketVector = typename Assembler::BucketVector;

    // setup assembler
    assembler_->setup(&hRef(),
                      advectionDiffusion,
                      domainVector_,
                      // anonymous function to compute Gamma for enthalpy
                      // equation: Gamma = lambda_eff / cp (thermal diffusivity)
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
        const STKScalarField& cp = cpRef().stkFieldRef();
        for (size_t ib = 0; ib < nodeBuckets.size(); ib++)
        {
            const Bucket& nodeBucket = *nodeBuckets[ib];
            const Bucket::size_type nNodesPerBucket = nodeBucket.size();

            // field chunks in bucket
            scalar* Gammab = field_data(Gamma, nodeBucket);
            const scalar* lambdaEffb = field_data(lambdaEff, nodeBucket);
            const scalar* cpb = field_data(cp, nodeBucket);

            for (Bucket::size_type i = 0; i < nNodesPerBucket; i++)
            {
                Gammab[i] = lambdaEffb[i] / cpb[i];
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

void thermalEnergyEquation::initialize()
{
    // raw initialization of temperature (user provides T in initial conditions)
    FOREACH_DOMAIN(initializeTemperature);

    // specific heat capacity may be initialized right after temperature in case
    // it varies with temperature, but also to initialize specific enthalpy
    FOREACH_DOMAIN(initializeSpecificHeatCapacity);

    // compute enthalpy from temperature
    FOREACH_DOMAIN(initializeSpecificEnthalpy);

    // 1) update enthalpy gradient
    FOREACH_DOMAIN(updateSpecificEnthalpyGradientField);

    // 2) update high-res fields for enthalpy
    FOREACH_DOMAIN_IF(updateSpecificEnthalpyBlendingFactorField,
                      domain->type() == domainType::fluid);

    // update enthalpy scale
    hRef().updateScale();

    // raw initialization of velocity if not yet: this is a situation when a
    // flow model is not enabled. The check is done inside
    FOREACH_DOMAIN_IF(initializeVelocity, domain->type() == domainType::fluid);

    // compute specific total enthalpy from specific enthalpy and velocity.
    // Placed here to make sure all the independent variables are initialized
    FOREACH_DOMAIN_IF(initializeSpecificTotalEnthalpy,
                      domain->type() == domainType::fluid);
}

void thermalEnergyEquation::postInitialize()
{
    // property initialization: density will be initialized here for solid
    // domains, because density in fluid domains has already been initialized
    FOREACH_DOMAIN(initializeDensity);
    FOREACH_DOMAIN(initializeThermalConductivity);
    FOREACH_DOMAIN_IF(initializeThermalExpansivity,
                      domain->type() == domainType::fluid &&
                          !domain->isMaterialCompressible() &&
                          domain->buoyancy_.option_ == buoyancyOption::buoyant);
    FOREACH_DOMAIN_IF(initializeCompressibility,
                      domain->type() == domainType::fluid &&
                          domain->isMaterialCompressible());

    cpRef().updateScale();
    lambdaRef().updateScale();

    equation::isInitialized_ = true;
}

void thermalEnergyEquation::preSolve()
{
    // enthalpy updates: enthalpy field must be updated with a raw update
    // because we are solving for it and we only intend to update its side
    // fields
    FOREACH_DOMAIN(updateSpecificEnthalpyPrevIterField);
    FOREACH_DOMAIN(updateSpecificEnthalpy);
    FOREACH_DOMAIN_IF(updateSpecificTotalEnthalpy,
                      domain->type() == domainType::fluid);

    // extract temperature from enthalpy to update properties
    FOREACH_DOMAIN(updateTemperature);

    // update density only for solid domains: for fluid domains, it will have
    // already been updated by the flow model
    FOREACH_DOMAIN_IF(updateDensity, domain->type() == domainType::solid);

    // other property updates
    FOREACH_DOMAIN(updateSpecificHeatCapacity);
    FOREACH_DOMAIN(updateSpecificHeatCapacityGradientField);
    FOREACH_DOMAIN(updateThermalConductivity);
    FOREACH_DOMAIN(updateEffectiveThermalConductivity); // laminar/solid
    FOREACH_DOMAIN_IF(updateThermalExpansivity,
                      domain->type() == domainType::fluid &&
                          !domain->isMaterialCompressible() &&
                          domain->buoyancy_.option_ == buoyancyOption::buoyant);
}

void thermalEnergyEquation::solve()
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
        if (domain->type() == domainType::fluid &&
            domain->isMaterialCompressible())
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
            hRef().stkFieldRef(),
            effectiveRelaxationFactor,
            lower_clip_value);

        // synchronize
        hRef().synchronizeGhostedEntities(domain->index());
    }

    // post correction

    // 1) extract temperature from solved enthalpy
    FOREACH_DOMAIN(updateTemperature);
    FOREACH_DOMAIN_IF(updateSpecificTotalEnthalpy,
                      domain->type() == domainType::fluid);

    // 2) update enthalpy gradient
    FOREACH_DOMAIN(updateSpecificEnthalpyGradientField);

    // 3) update high-res fields
    FOREACH_DOMAIN_IF(updateSpecificEnthalpyBlendingFactorField,
                      domain->type() == domainType::fluid);

    // 4) update scale
    hRef().updateScale();
    cpRef().updateScale();
    lambdaRef().updateScale();

#ifdef HAS_INTERFACE
    // 5) Post-processed quantities
    FOREACH_DOMAIN(updateInterfaceHeatImbalance_);
#endif /* HAS_INTERFACE */
}

void thermalEnergyEquation::preTimeStep()
{
    FOREACH_DOMAIN(updateSpecificEnthalpyPrevTimeField);
    FOREACH_DOMAIN_IF(updateSpecificTotalEnthalpyPrevTimeField,
                      domain->type() == domainType::fluid);
    FOREACH_DOMAIN_IF(updateDensityPrevTimeField,
                      domain->type() == domainType::solid);
    FOREACH_DOMAIN(updateSpecificHeatCapacityPrevTimeField);
    FOREACH_DOMAIN_IF(updateCompressibilityPrevTimeField,
                      domain->type() == domainType::fluid &&
                          domain->isMaterialCompressible());
}

void thermalEnergyEquation::postSolve()
{
    this->reportHeatData_();

    FOREACH_DOMAIN_IF(updateCompressibility,
                      domain->type() == domainType::fluid &&
                          domain->isMaterialCompressible());
    FOREACH_DOMAIN(updateSpecificHeatCapacity);

    // calculate post-processed quantities
    FOREACH_DOMAIN(updateTotalTemperatureField_);
}

void thermalEnergyEquation::setResidualScales_()
{
    auto& residual_scales = linearSystem::getContext()->getResidualScales();
    residual_scales = {1.0 / (hRef().scale() + ::accel::SMALL)};
}

void thermalEnergyEquation::printScales()
{
    if (messager::master())
    {
        std::cout << this->hRef().name() << " scales:" << std::endl;
        std::cout << "\tmin: " << std::scientific << std::setprecision(8)
                  << hRef().minScale() << std::endl;
        std::cout << "\tmax: " << std::scientific << std::setprecision(8)
                  << hRef().maxScale() << std::endl;
        std::cout << "\tscale: " << std::scientific << std::setprecision(8)
                  << hRef().scale() << std::endl
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
