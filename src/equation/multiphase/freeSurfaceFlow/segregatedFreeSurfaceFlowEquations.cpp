// File : segregatedFreeSurfaceFlowEquations.cpp
// Created : Sun Jan 26 2025 22:53:38 (+0100)
// Author : Mhamad Mahdi Alloush
// Description:
// Copyright (c) 2025 CCFNUM, Lucerne University of Applied Sciences and Arts.
// SPDX-License-Identifier: BSD-3-Clause

#include "segregatedFreeSurfaceFlowEquations.h"

namespace accel
{

segregatedFreeSurfaceFlowEquations::segregatedFreeSurfaceFlowEquations(
    realm* realm)
    : equation("Segregated Free Surface Flow"), freeSurfaceFlowModel(realm)
{
    U_eq_ = std::make_unique<bulkNavierStokesEquation>(realm, this);
    pCorr_eq_ = std::make_unique<bulkPressureCorrectionEquation>(realm, this);

    for (label iPhase = 0; iPhase < nPhases(); iPhase++)
    {
        if (!this->phaseRef(iPhase).primaryPhase_)
        {
            std::unique_ptr alpha_eq = std::make_unique<volumeFractionEquation>(
                realm, this, phaseIndex(iPhase));
            alpha_eq_.push_back(std::move(alpha_eq));
        }
        else
        {
            alpha_eq_.push_back(nullptr);
        }
    }

    // set sub-iterations
    subIters_ = this->controlsRef()
                    .solverRef()
                    .solverControl_.advancedOptions_.equationControls_
                    .subIterations_.segregatedFlow_;
}

void segregatedFreeSurfaceFlowEquations::addDomain(
    std::shared_ptr<domain> domain)
{
    equation::addDomain(domain);

    U_eq_->addDomain(domain);
    pCorr_eq_->addDomain(domain);

    for (label iPhase = 0; iPhase < nPhases(); iPhase++)
    {
        if (!this->phaseRef(iPhase).primaryPhase_)
        {
            assert(alpha_eq_[iPhase]);
            alpha_eq_[iPhase]->addDomain(domain);
        }
    }
}

bool segregatedFreeSurfaceFlowEquations::isConverged() const
{
    bool converged = U_eq_->isConverged() && pCorr_eq_->isConverged();
    for (label iPhase = 0; iPhase < nPhases(); iPhase++)
    {
        if (!this->phaseRef(iPhase).primaryPhase_)
        {
            assert(alpha_eq_[iPhase]);
            converged = converged && alpha_eq_[iPhase]->isConverged();
        }
    }

    return converged;
}

void segregatedFreeSurfaceFlowEquations::setup()
{
    U_eq_->setup();
    pCorr_eq_->setup();
    for (label iPhase = 0; iPhase < nPhases(); iPhase++)
    {
        if (!this->phaseRef(iPhase).primaryPhase_)
        {
            assert(alpha_eq_[iPhase]);
            alpha_eq_[iPhase]->setup();
        }
        else
        {
            FOREACH_DOMAIN(setupVolumeFraction, phaseIndex(iPhase));
        }
    }

    equation::isCreated_ = U_eq_->isCreated() && pCorr_eq_->isCreated();
    for (label iPhase = 0; iPhase < nPhases(); iPhase++)
    {
        if (!this->phaseRef(iPhase).primaryPhase_)
        {
            equation::isCreated_ =
                equation::isCreated_ && alpha_eq_[iPhase]->isCreated();
        }
    }
}

void segregatedFreeSurfaceFlowEquations::initialize()
{
    // initialize velocity and pressure
    U_eq_->initialize();
    pCorr_eq_->initialize();

    // initialize volume fractions
    for (label iPhase = 0; iPhase < nPhases(); iPhase++)
    {
        if (!this->phaseRef(iPhase).primaryPhase_)
        {
            assert(alpha_eq_[iPhase]);
            alpha_eq_[iPhase]->initialize();
        }
        else
        {
            // primary phase: raw initialization and update gradient. No HR
            // fields are updated
            FOREACH_DOMAIN(initializeVolumeFraction, phaseIndex(iPhase));

            // 1) update gradient
            FOREACH_DOMAIN(updateVolumeFractionGradientField,
                           phaseIndex(iPhase));

            // 2) update scale
            this->alphaRef(phaseIndex(iPhase)).updateScale();
        }
    }
}

void segregatedFreeSurfaceFlowEquations::postInitialize()
{
    // initialize phase and bulk properties, and initialize phase mass flux
    U_eq_->postInitialize();
    pCorr_eq_->postInitialize();

    equation::isInitialized_ =
        U_eq_->isInitialized() && pCorr_eq_->isInitialized();
    for (label iPhase = 0; iPhase < nPhases(); iPhase++)
    {
        if (!this->phaseRef(iPhase).primaryPhase_)
        {
            assert(alpha_eq_[iPhase]);
            alpha_eq_[iPhase]->postInitialize();

            equation::isInitialized_ =
                equation::isInitialized_ && alpha_eq_[iPhase]->isInitialized();
        }
    }
}

void segregatedFreeSurfaceFlowEquations::solve()
{
    // predictor step: solve bulk momentum
    {
        U_eq_->preSolve();
        U_eq_->solve();
        U_eq_->postSolve();
    }

    // corrector step: solve bulk pressure correction
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
                stk::mesh::BulkData& bulkData = this->meshRef().bulkDataRef();
                stk::mesh::MetaData& metaData = this->meshRef().metaDataRef();

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

                const auto& duSTKFieldRef = *metaData.get_field<scalar>(
                    stk::topology::NODE_RANK, freeSurfaceFlowModel::du_ID);

                const auto& gradPSTKFieldRef =
                    this->pRef().gradRef().stkFieldRef();
                auto& USTKFieldRef = this->URef().stkFieldRef();

                // get relaxation factor
                scalar urf = this->pRef().urf();

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
                        stk::mesh::field_data(duSTKFieldRef, nodeBucket);
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

            // Update pressure gradient
            FOREACH_DOMAIN(updatePressureGradientField);

            // correct velocity field step 2
            FOREACH_DOMAIN_RAW({
                stk::mesh::BulkData& bulkData = this->meshRef().bulkDataRef();
                stk::mesh::MetaData& metaData = this->meshRef().metaDataRef();

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

                const auto& gradPSTKFieldRef =
                    this->pRef().gradRef().stkFieldRef();
                auto& USTKFieldRef = this->URef().stkFieldRef();

                // get relaxation factor
                scalar urf = this->pRef().urf();

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

            // correct phase density (only if compressible)
            for (label iPhase = 0; iPhase < nPhases(); iPhase++)
            {
                FOREACH_DOMAIN_IF(
                    fieldBroker::updateDensity,
                    domain->isMaterialCompressible(
                        domain->globalToLocalMaterialIndex(phaseIndex(iPhase))),
                    phaseIndex(iPhase));
            }

            // update velocity gradient
            FOREACH_DOMAIN(updateVelocityGradientField);

            // update high-res fields
            FOREACH_DOMAIN(updateVelocityBlendingFactorField);

            // correct mass flux and do updates:
            // 1) mass divergence
            // 2) flow reversal flag side field
            for (label iPhase = 0; iPhase < nPhases(); iPhase++)
            {
                FOREACH_DOMAIN(updateSideMassFlowRateFraction,
                               phaseIndex(iPhase));
                FOREACH_DOMAIN(transformMassFlowRateToAbsolute,
                               phaseIndex(iPhase));
                FOREACH_DOMAIN(updateMassFlowRate, phaseIndex(iPhase));
                FOREACH_DOMAIN(transformMassFlowRateToRelative,
                               phaseIndex(iPhase));
                FOREACH_DOMAIN(updateFlowReversalFlag, phaseIndex(iPhase));
                FOREACH_DOMAIN(updateMassDivergenceField, phaseIndex(iPhase));
            }

            pCorr_eq_->postSolve();

            // if converged .. break sub-iter loop
            if (pCorr_eq_->isConverged())
            {
                break;
            }
        }
    }

    // solve volume fraction
    {
        for (label iPhase = 0; iPhase < nPhases(); iPhase++)
        {
            if (!this->phaseRef(iPhase).primaryPhase_)
            {
                assert(alpha_eq_[iPhase]);
                alpha_eq_[iPhase]->preSolve();
            }
            else
            {
                // update primary phase (i.e. boundary conditions)
                FOREACH_DOMAIN(updateVolumeFraction, phaseIndex(iPhase));
            }
        }

        for (label iPhase = 0; iPhase < nPhases(); iPhase++)
        {
            if (!this->phaseRef(iPhase).primaryPhase_)
            {
                assert(alpha_eq_[iPhase]);
                alpha_eq_[iPhase]->solve();
            }
        }

        // deduce the volume fraction of the last phase from volume conservation
        FOREACH_DOMAIN(applyVolumeConservation);

        // update gradient and beta field for primary phases
        for (label iPhase = 0; iPhase < nPhases(); iPhase++)
        {
            if (this->phaseRef(iPhase).primaryPhase_)
            {
                // 1) update alpha gradient
                FOREACH_DOMAIN(updateVolumeFractionGradientField,
                               phaseIndex(iPhase));

                // 2) update alpha scale
                this->alphaRef(phaseIndex(iPhase)).updateScale();
            }
        }

        for (label iPhase = 0; iPhase < nPhases(); iPhase++)
        {
            if (!this->phaseRef(iPhase).primaryPhase_)
            {
                assert(alpha_eq_[iPhase]);
                alpha_eq_[iPhase]->postSolve();
            }
        }

        // bulk properties update
        FOREACH_DOMAIN(updateDensity);
        FOREACH_DOMAIN(updateDynamicViscosity);

        // update bulk mass flux: flow reversal field is updated above and in
        // case of any zero mass flux at outlet this will be implicitly be
        // calculated upon the following update
        FOREACH_DOMAIN(updateMassFlowRate);

        // update div for the mixture mass flux and others
        FOREACH_DOMAIN(flowModel::updateMassDivergenceField);
    }
}

void segregatedFreeSurfaceFlowEquations::postSolve()
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

void segregatedFreeSurfaceFlowEquations::preTimeStep()
{
    resetCourantNumber();

    U_eq_->preTimeStep();
    pCorr_eq_->preTimeStep();

    for (label iPhase = 0; iPhase < nPhases(); iPhase++)
    {
        if (!this->phaseRef(iPhase).primaryPhase_)
        {
            assert(alpha_eq_[iPhase]);
            alpha_eq_[iPhase]->preTimeStep();
        }
        else
        {
            FOREACH_DOMAIN(updateVolumeFractionPrevTimeField,
                           phaseIndex(iPhase));
        }
    }
}

void segregatedFreeSurfaceFlowEquations::printScales()
{
    U_eq_->printScales();
    pCorr_eq_->printScales();

    for (label iPhase = 0; iPhase < nPhases(); iPhase++)
    {
        if (!this->phaseRef(iPhase).primaryPhase_)
        {
            assert(alpha_eq_[iPhase]);
            alpha_eq_[iPhase]->printScales();
        }
    }
}

} /* namespace accel */
