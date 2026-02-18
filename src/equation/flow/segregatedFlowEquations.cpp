// File : segregatedFlowEquations.cpp
// Created : Fri Mar 15 2024 15:06:38 (+0100)
// Author : Fabian Wermelinger
// Description: Segregated Navier-Stokes equations implementation details
// Copyright (c) 2024 CCFNUM, Lucerne University of Applied Sciences and Arts.
// SPDX-License-Identifier: BSD-3-Clause

#include "segregatedFlowEquations.h"
#include "realm.h"

namespace accel
{

segregatedFlowEquations::segregatedFlowEquations(realm* realm)
    : equation("Segregated Flow"), flowModel(realm)
{
    U_eq_ = std::make_unique<navierStokesEquation>(realm, this);
    pCorr_eq_ = std::make_unique<pressureCorrectionEquation>(realm, this);

    // set relaxation factor for mass flux field to 0.75 for steady-state cases,
    // if and only if, not specified by user
    if (controlsRef()
            .solverRef()
            .solverControl_.basicSettings_.convergenceControl_
            .relaxationParameters_.relaxMass_ < 0.9999)
    {
        this->mDotRef().setURF(
            this->controlsRef()
                .solverRef()
                .solverControl_.basicSettings_.convergenceControl_
                .relaxationParameters_.relaxMass_);
    }

    // set sub-iterations
    subIters_ = this->controlsRef()
                    .solverRef()
                    .solverControl_.advancedOptions_.equationControls_
                    .subIterations_.segregatedFlow_;
}

void segregatedFlowEquations::addDomain(std::shared_ptr<domain> domain)
{
    equation::addDomain(domain);

    U_eq_->addDomain(domain);
    pCorr_eq_->addDomain(domain);
}

bool segregatedFlowEquations::isConverged() const
{
    return U_eq_->isConverged() && pCorr_eq_->isConverged();
}

void segregatedFlowEquations::setup()
{
    U_eq_->setup();
    pCorr_eq_->setup();

    equation::isCreated_ = U_eq_->isCreated() && pCorr_eq_->isCreated();
}

void segregatedFlowEquations::initialize()
{
    U_eq_->initialize();
    pCorr_eq_->initialize();
}

void segregatedFlowEquations::postInitialize()
{
    U_eq_->postInitialize();
    pCorr_eq_->postInitialize();

    equation::isInitialized_ =
        U_eq_->isInitialized() && pCorr_eq_->isInitialized();
}

void segregatedFlowEquations::solve()
{
    // predictor step: solve momentum
    {
        U_eq_->preSolve();
        U_eq_->solve();
        U_eq_->postSolve();
    }

    // corrrector step: solve pressure correction
    {
        if (messager::master() && pCorr_eq_->subIters() > 1)
        {
            std::cout << std::endl
                      << "Sub-iterating " + pCorr_eq_->name() << ":\n";
        }

        for (label subIter = 1; subIter <= pCorr_eq_->subIters(); subIter++)
        {
            // only print sub-iter if it is active (> 1)
            if (messager::master() && pCorr_eq_->subIters() > 1)
            {
                std::cout << std::endl << " sub-iter: " << subIter << "\n";
            }

            pCorr_eq_->preSolve();
            pCorr_eq_->solve();

            // correct velocity field step 1
            FOREACH_DOMAIN_RAW({
                stk::mesh::BulkData& bulkData = meshRef().bulkDataRef();
                stk::mesh::MetaData& metaData = meshRef().metaDataRef();

                // SIMPLE-Consistent
                const bool consistent =
                    controlsRef()
                        .solverRef()
                        .solverControl_.expertParameters_.consistent_;

                const auto* duSTKFieldPtr =
                    consistent
                        ? metaData.get_field<scalar>(stk::topology::NODE_RANK,
                                                     flowModel::duTilde_ID)
                        : metaData.get_field<scalar>(stk::topology::NODE_RANK,
                                                     flowModel::du_ID);

                const auto& gradPSTKFieldRef = pRef().gradRef().stkFieldRef();
                auto& USTKFieldRef = URef().stkFieldRef();

                // get relaxation factor
                scalar urf = pRef().urf();

                // find the ramp value (explicit relaxation) for compressible
                // domains
                if (domain->isMaterialCompressible())
                {
                    scalar iter = controlsRef().iter;
                    label rampIter = 20;
                    if (iter <= rampIter)
                    {
                        urf *= std::max(scalar(iter) / scalar(rampIter), 0.1);
                    }
                }

                stk::mesh::Selector selAllNodes =
                    this->meshRef().metaDataRef().universal_part() &
                    stk::mesh::selectUnion(domain->zonePtr()->interiorParts());
                stk::mesh::BucketVector const& nodeBuckets =
                    this->meshRef().bulkDataRef().get_buckets(
                        stk::topology::NODE_RANK, selAllNodes);
                for (stk::mesh::BucketVector::const_iterator ib =
                         nodeBuckets.begin();
                     ib != nodeBuckets.end();
                     ++ib)
                {
                    stk::mesh::Bucket& nodeBucket = **ib;
                    const stk::mesh::Bucket::size_type nNodesPerBucket =
                        nodeBucket.size();

                    scalar* Ub =
                        stk::mesh::field_data(USTKFieldRef, nodeBucket);
                    const scalar* dub =
                        stk::mesh::field_data(*duSTKFieldPtr, nodeBucket);
                    const scalar* dpdxb =
                        stk::mesh::field_data(gradPSTKFieldRef, nodeBucket);

                    for (stk::mesh::Bucket::size_type iNode = 0;
                         iNode < nNodesPerBucket;
                         ++iNode)
                    {
                        for (label i = 0; i < SPATIAL_DIM; i++)
                        {
                            Ub[SPATIAL_DIM * iNode + i] +=
                                1.0 / urf * dub[SPATIAL_DIM * iNode + i] *
                                dpdxb[SPATIAL_DIM * iNode + i];
                        }
                    }
                }
            });

            // update pressure gradient
            FOREACH_DOMAIN(updatePressureGradientField);

            // correct velocity field step 2
            FOREACH_DOMAIN_RAW({
                stk::mesh::BulkData& bulkData = meshRef().bulkDataRef();
                stk::mesh::MetaData& metaData = meshRef().metaDataRef();

                // SIMPLE-Consistent
                const bool consistent =
                    controlsRef()
                        .solverRef()
                        .solverControl_.expertParameters_.consistent_;

                const auto* duSTKFieldPtr =
                    consistent
                        ? metaData.get_field<scalar>(stk::topology::NODE_RANK,
                                                     flowModel::duTilde_ID)
                        : metaData.get_field<scalar>(stk::topology::NODE_RANK,
                                                     flowModel::du_ID);

                const auto& gradPSTKFieldRef = pRef().gradRef().stkFieldRef();
                auto& USTKFieldRef = URef().stkFieldRef();

                // get relaxation factor
                scalar urf = pRef().urf();

                // find the ramp value (explicit relaxation) for compressible
                // domains
                if (domain->isMaterialCompressible())
                {
                    scalar iter = controlsRef().iter;
                    label rampIter = 20;
                    if (iter <= rampIter)
                    {
                        urf *= std::max(scalar(iter) / scalar(rampIter), 0.1);
                    }
                }

                stk::mesh::Selector selAllNodes =
                    this->meshRef().metaDataRef().universal_part() &
                    stk::mesh::selectUnion(domain->zonePtr()->interiorParts());
                stk::mesh::BucketVector const& nodeBuckets =
                    this->meshRef().bulkDataRef().get_buckets(
                        stk::topology::NODE_RANK, selAllNodes);
                for (stk::mesh::BucketVector::const_iterator ib =
                         nodeBuckets.begin();
                     ib != nodeBuckets.end();
                     ++ib)
                {
                    stk::mesh::Bucket& nodeBucket = **ib;
                    const stk::mesh::Bucket::size_type nNodesPerBucket =
                        nodeBucket.size();

                    scalar* Ub =
                        stk::mesh::field_data(USTKFieldRef, nodeBucket);
                    const scalar* dub =
                        stk::mesh::field_data(*duSTKFieldPtr, nodeBucket);
                    const scalar* dpdxb =
                        stk::mesh::field_data(gradPSTKFieldRef, nodeBucket);

                    for (stk::mesh::Bucket::size_type iNode = 0;
                         iNode < nNodesPerBucket;
                         ++iNode)
                    {
                        for (label i = 0; i < SPATIAL_DIM; i++)
                        {
                            Ub[SPATIAL_DIM * iNode + i] -=
                                1.0 / urf * dub[SPATIAL_DIM * iNode + i] *
                                dpdxb[SPATIAL_DIM * iNode + i];
                        }
                    }
                }
            });

            // update velocity scale: raw
            this->URef().updateScale();
            this->updatePressureScale();

            // correct density (only if compressible)
            FOREACH_DOMAIN_IF(updateDensity, domain->isMaterialCompressible());

            // update density-related fields
            FOREACH_DOMAIN_IF(updateDensityGradientField,
                              domain->isMaterialCompressible());
            FOREACH_DOMAIN_IF(updateDensityBlendingFactorField,
                              domain->isMaterialCompressible());

            // update velocity gradient
            FOREACH_DOMAIN(updateVelocityGradientField);

            // update velocity high-res fields
            FOREACH_DOMAIN(updateVelocityBlendingFactorField);

            // correct mass flux and do updates:
            // 1) mass divergence
            // 2) flow reversal flag side field
            FOREACH_DOMAIN(updateSideMassFlowRateFraction);
            FOREACH_DOMAIN(transformMassFlowRateToAbsolute);
            FOREACH_DOMAIN(updateMassFlowRate);
            FOREACH_DOMAIN(transformMassFlowRateToRelative);
            FOREACH_DOMAIN(updateFlowReversalFlag);
            FOREACH_DOMAIN(updateMassDivergenceField);

            pCorr_eq_->postSolve();

            // if converged .. break sub-iter loop
            if (pCorr_eq_->isConverged())
            {
                break;
            }
        }
    }
}

void segregatedFlowEquations::postSolve()
{
    // calculate post-processed quantities

    FOREACH_DOMAIN(updateCourantNumberField_);
    FOREACH_DOMAIN(updateMachNumberField_);
    FOREACH_DOMAIN(updateTotalPressureField_);
    FOREACH_DOMAIN(updateRelativeVelocityField_);
    FOREACH_DOMAIN(updateWallShearStress); // laminar
    FOREACH_DOMAIN(updateMassImbalance_);

    this->reportFlowData_();
}

void segregatedFlowEquations::preTimeStep()
{
    resetCourantNumber();

    U_eq_->preTimeStep();
    pCorr_eq_->preTimeStep();
}

void segregatedFlowEquations::printScales()
{
    U_eq_->printScales();
    pCorr_eq_->printScales();
}

} /* namespace accel */
