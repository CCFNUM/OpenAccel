// File : bulkNavierStokesEquation.cpp
// Created : Thu Jan 30 2025 10:02:56 (+0100)
// Author : Mhamad Mahdi Alloush
// Description:
// Copyright (c) 2025 CCFNUM, Lucerne University of Applied Sciences and Arts.
// SPDX-License-Identifier: BSD-3-Clause

#include "bulkNavierStokesEquation.h"
#include "realm.h"

namespace accel
{

bulkNavierStokesEquation::bulkNavierStokesEquation(realm* realm,
                                                   freeSurfaceFlowModel* model)
    : navierStokesEquation(realm, model), model_(model)
{
}

void bulkNavierStokesEquation::setup()
{
    // setup bulk fields
    navierStokesEquation::setup();

    // setup phase fields
    for (label iPhase = 0; iPhase < model_->nPhases(); iPhase++)
    {
        label phaseIndex = model_->phaseIndex(iPhase);

        FOREACH_DOMAIN(model_->fieldBroker::setupMassFlowRate, phaseIndex);
        FOREACH_DOMAIN(model_->fieldBroker::setupDensity, phaseIndex);
        FOREACH_DOMAIN(model_->fieldBroker::setupDynamicViscosity, phaseIndex);

        if (model_->rhoRef(phaseIndex).blendingFactorPtr() == nullptr)
        {
            model_->rhoRef(phaseIndex)
                .setupBlendingFactorField(
                    /*enable only a dummy field*/ true);
        }
    }
}

void bulkNavierStokesEquation::postInitialize()
{
    // initialize phase fields
    for (label iPhase = 0; iPhase < model_->nPhases(); iPhase++)
    {
        label phaseIndex = model_->phaseIndex(iPhase);

        // initialize phase density and viscosity
        FOREACH_DOMAIN(model_->fieldBroker::initializeDensity, phaseIndex);
        FOREACH_DOMAIN_IF(model_->fieldBroker::updateDensityGradientField,
                          phaseIndex,
                          domain->isMaterialCompressible(
                              domain->globalToLocalMaterialIndex(phaseIndex)));
        FOREACH_DOMAIN_IF(model_->fieldBroker::updateDensityBlendingFactorField,
                          phaseIndex,
                          domain->isMaterialCompressible(
                              domain->globalToLocalMaterialIndex(phaseIndex)));
        FOREACH_DOMAIN(model_->fieldBroker::initializeDynamicViscosity,
                       phaseIndex);

        // update scale of phase density and viscosity
        model_->rhoRef(phaseIndex).updateScale();
        model_->muRef(phaseIndex).updateScale();

        // initialize phase mass flux
        FOREACH_DOMAIN(model_->fieldBroker::initializeMassFlowRate, phaseIndex);

        FOREACH_DOMAIN(model_->transformMassFlowRateToRelative, phaseIndex);

        // in case of transient, old density must be updated before div update
        if (model_->controlsRef().isTransient())
        {
            FOREACH_DOMAIN(model_->updateDensityPrevTimeField, phaseIndex);
        }

        // update div for mass of the current phase
        FOREACH_DOMAIN(model_->updateMassDivergenceField, phaseIndex);
    }

    // initialize bulk density and viscosity
    FOREACH_DOMAIN(model_->initializeDensity);
    FOREACH_DOMAIN(model_->initializeDynamicViscosity);

    // update scale of bulk density and viscosity
    model_->rhoRef().updateScale();
    model_->muRef().updateScale();

    // update bulk mass flux: we do here in order to make use of all
    // transformations done just before
    FOREACH_DOMAIN(model_->updateMassFlowRate);

    // in case of transient, old density must be updated before div update
    if (model_->controlsRef().isTransient())
    {
        FOREACH_DOMAIN(model_->updateDensityPrevTimeField);
    }

    // update div for the bulk mass flux
    FOREACH_DOMAIN(model_->flowModel::updateMassDivergenceField);
}

void bulkNavierStokesEquation::preSolve()
{
    // update phase properties: this is in case any material has its property
    // changing due to some physics (ideal gas law in a compressible material)
    for (label iPhase = 0; iPhase < model_->nPhases(); iPhase++)
    {
        label phaseIndex = model_->phaseIndex(iPhase);

        FOREACH_DOMAIN(model_->fieldBroker::updateDensity, phaseIndex);
        FOREACH_DOMAIN(model_->fieldBroker::updateDynamicViscosity, phaseIndex);

        // update density-related fields
        FOREACH_DOMAIN_IF(model_->fieldBroker::updateDensityGradientField,
                          phaseIndex,
                          domain->isMaterialCompressible(
                              domain->globalToLocalMaterialIndex(phaseIndex)));
        FOREACH_DOMAIN_IF(model_->fieldBroker::updateDensityBlendingFactorField,
                          phaseIndex,
                          domain->isMaterialCompressible(
                              domain->globalToLocalMaterialIndex(phaseIndex)));
    }

    // update velocity and bulk properties
    navierStokesEquation::preSolve();
}

void bulkNavierStokesEquation::preTimeStep()
{
    for (label iPhase = 0; iPhase < model_->nPhases(); iPhase++)
    {
        label phaseIndex = model_->phaseIndex(iPhase);

        FOREACH_DOMAIN(model_->updateDensityPrevTimeField, phaseIndex);
    }

    navierStokesEquation::preTimeStep();
}

} /* namespace accel */
