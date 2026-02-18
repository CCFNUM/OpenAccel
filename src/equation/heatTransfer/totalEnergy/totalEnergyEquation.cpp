// File : totalEnergyEquation.cpp
// Created : Thu Mar 27 2025 10:42:10 (+0100)
// Author : Fabian Wermelinger
// Description:
// Copyright (c) 2025 CCFNUM, Lucerne University of Applied Sciences and Arts.
// SPDX-License-Identifier: BSD-3-Clause

#include "totalEnergyEquation.h"

namespace accel
{

totalEnergyEquation::totalEnergyEquation(realm* realm)
    : equation("Total Energy"), heatTransferModel(realm),
      linearSystem(realm->simulationRef()),
      assembler_(std::make_unique<Assembler>(this))
{
    this->setEquationName({"h0"});

    // set relaxation factor for wall scale field. Default is 1
    this->h0Ref().setURF(controlsRef()
                             .solverRef()
                             .solverControl_.basicSettings_.convergenceControl_
                             .relaxationParameters_.energyRelaxationFactor_);
}

void totalEnergyEquation::checkDomain(const std::shared_ptr<domain> domain)
{
    // allowed for fluids and solids
}

bool totalEnergyEquation::isConverged() const
{
    return linearSystem::isConverged();
}

void totalEnergyEquation::setup()
{
    // setup of fields on defined domains: rho field
    // might have been already initialized over other domains
    FOREACH_DOMAIN(setupDensity);
    FOREACH_DOMAIN_IF(setupVelocity, domain->type() == domainType::fluid);
    FOREACH_DOMAIN_IF(setupMassFlowRate, domain->type() == domainType::fluid);
    FOREACH_DOMAIN(setupSpecificHeatCapacity);
    FOREACH_DOMAIN(setupThermalConductivity);
    FOREACH_DOMAIN(setupThermalExpansivity);
    FOREACH_DOMAIN(setupTemperature);
    FOREACH_DOMAIN(setupSpecificEnthalpy);
    FOREACH_DOMAIN(setupSpecificTotalEnthalpy);
    FOREACH_DOMAIN_IF(setupCompressibility,
                      domain->type() == domainType::fluid &&
                          domain->isMaterialCompressible());
    FOREACH_DOMAIN(setupHeatFlowRate);

    using Bucket = typename Assembler::Bucket;
    using BucketVector = typename Assembler::BucketVector;

    // setup assembler
    assembler_->setup(&h0Ref(),
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
    // FIXME: [2024-03-13] Consider passing mesh argument or
    // connectivity arrays passed to initialize() directly is more flexible
    // rather than this->meshRef() which is set through simulation object
    // obtained via realm in fieldBroker
    linearSystem::setupSolver(this->name(), fieldBroker::meshRef());

    equation::isCreated_ = true;
}

void totalEnergyEquation::initialize()
{
    // raw initialization of temperature
    FOREACH_DOMAIN(initializeTemperature);

    // update temperature gradient
    FOREACH_DOMAIN(updateTemperatureGradientField);

    // update velocity scale
    TRef().updateScale();
}

void totalEnergyEquation::postInitialize()
{
    // property initialization: density will be initialized here for solid
    // domains only (implicitly), because density in fluid domains has already
    // been initialized
    FOREACH_DOMAIN(initializeDensity);
    FOREACH_DOMAIN(initializeSpecificHeatCapacity);
    FOREACH_DOMAIN(initializeThermalConductivity);
    FOREACH_DOMAIN(initializeThermalExpansivity);
    FOREACH_DOMAIN_IF(initializeCompressibility,
                      domain->type() == domainType::fluid &&
                          domain->isMaterialCompressible());

    // initialization of necessary fields
    FOREACH_DOMAIN(initializeSpecificEnthalpy);
    FOREACH_DOMAIN(initializeSpecificTotalEnthalpy);

    // update gradients of necessary fields
    FOREACH_DOMAIN(updateSpecificEnthalpyGradientField);
    FOREACH_DOMAIN(updateSpecificTotalEnthalpyGradientField);

    // update high-res field of specific total enthalpy
    FOREACH_DOMAIN_IF(updateSpecificTotalEnthalpyBlendingFactorField,
                      domain->type() == domainType::fluid);

    // update scales
    hRef().updateScale();
    h0Ref().updateScale();
    cpRef().updateScale();
    lambdaRef().updateScale();

    equation::isInitialized_ = true;
}

void totalEnergyEquation::preSolve()
{
    // update density only for solid domains: for fluid domains, it will have
    // already been updated by the flow model
    FOREACH_DOMAIN_IF(updateDensity, domain->type() == domainType::solid);

    // other property updates
    FOREACH_DOMAIN(updateSpecificHeatCapacity);
    FOREACH_DOMAIN(updateThermalConductivity);
    FOREACH_DOMAIN(updateEffectiveThermalConductivity); // laminar/solid
    FOREACH_DOMAIN(updateThermalExpansivity);
    FOREACH_DOMAIN(updateSpecificEnthalpyPrevIterField);
    FOREACH_DOMAIN(updateSpecificTotalEnthalpyPrevIterField);
}

void totalEnergyEquation::solve()
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
    scalar relaxCorrection = 1.0;
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

        this->template correctField_<linearSystem::BLOCKSIZE, 1>(
            domain.get(),
            ctx->getXVector(),
            stk::topology::NODE_RANK,
            h0Ref().stkFieldRef(),
            effectiveRelaxationFactor);

        // synchronize
        h0Ref().synchronizeGhostedEntities(domain->index());
    }

    // post correction

    // 1) update gradient
    FOREACH_DOMAIN(updateSpecificTotalEnthalpyGradientField);

    // 2) update high-res fields
    FOREACH_DOMAIN_IF(updateSpecificTotalEnthalpyBlendingFactorField,
                      domain->type() == domainType::fluid);

    // 3) update scale
    hRef().updateScale();
    h0Ref().updateScale();
    TRef().updateScale();
    cpRef().updateScale();
    lambdaRef().updateScale();
}

void totalEnergyEquation::postSolve()
{
    this->reportHeatData_();

    FOREACH_DOMAIN(updateSpecificEnthalpy);
    FOREACH_DOMAIN(updateTemperature);
    FOREACH_DOMAIN_IF(updateCompressibility,
                      domain->type() == domainType::fluid &&
                          domain->isMaterialCompressible());
    FOREACH_DOMAIN(updateSpecificHeatCapacity);
    FOREACH_DOMAIN(updateSpecificTotalEnthalpy); // side-fields

    FOREACH_DOMAIN(updateTemperatureGradientField);
    FOREACH_DOMAIN(updateSpecificEnthalpyGradientField);

    // calculate post-processed quantities
    FOREACH_DOMAIN(updateTotalTemperatureField_);
}

void totalEnergyEquation::preTimeStep()
{
    FOREACH_DOMAIN(updateSpecificEnthalpyPrevTimeField);
    FOREACH_DOMAIN(updateSpecificTotalEnthalpyPrevTimeField);
    FOREACH_DOMAIN(updateTemperaturePrevTimeField);
    FOREACH_DOMAIN(updateDensityPrevTimeField);
    FOREACH_DOMAIN(updateSpecificHeatCapacityPrevTimeField);
    FOREACH_DOMAIN_IF(updateCompressibilityPrevTimeField,
                      domain->type() == domainType::fluid &&
                          domain->isMaterialCompressible());
}

void totalEnergyEquation::setResidualScales_()
{
    auto& residual_scales = linearSystem::getContext()->getResidualScales();
    residual_scales = {1.0 / (h0Ref().scale() + ::accel::SMALL)};
}

void totalEnergyEquation::printScales()
{
    if (messager::master())
    {
        std::cout << h0Ref().name() << " scales:" << std::endl;
        std::cout << "\tmin: " << std::scientific << std::setprecision(8)
                  << h0Ref().minScale() << std::endl;
        std::cout << "\tmax: " << std::scientific << std::setprecision(8)
                  << h0Ref().maxScale() << std::endl;
        std::cout << "\tscale: " << std::scientific << std::setprecision(8)
                  << h0Ref().scale() << std::endl
                  << std::endl;

        std::cout << TRef().name() << " scales:" << std::endl;
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
