// File       : fieldBrokerPhasic.cpp
// Created    : Mon Apr 21 2025 16:09:24 (+0100)
// Author     : Mhamad Mahdi Alloush
// Description:
// Copyright (c) 2025 CCFNUM, Lucerne University of Applied Sciences and Arts.
// SPDX-License-Identifier: BSD-3-Clause

// code
#include "dataHandler.h"
#include "domain.h"
#include "fieldBroker.h"
#include "initialConditions.h"
#include "macros.h"
#include "realm.h"
#include "simulation.h"
#include "types.h"

namespace accel
{

void fieldBroker::setupVolumeFraction(const std::shared_ptr<domain> domain,
                                      label iPhase)
{
    // set zone
    if (alphaRef(iPhase).isZoneUnset(domain->index()))
    {
        alphaRef(iPhase).setZone(domain->index());

        // set initial conditions from input
        initialCondition::
            setupFluidSpecificFieldInitializationOverDomainFromInput(
                alphaRef(iPhase),
                realm::alpha_ID,
                realmPtr_->simulationRef().materialName(iPhase),
                domain);

        // boundary conditions for this domain and phase
        setupBoundaryConditions_(
            domain,
            iPhase,
            // anonymous function to set boundary conditions for this model
            [this](const ::accel::domain* domain,
                   const label iBoundary,
                   const label iPhase,
                   const boundaryPhysicalType bc_type,
                   const YAML::Node fluidValues)
        {
            auto& bc = alphaRef(iPhase).boundaryConditionRef(domain->index(),
                                                             iBoundary);

#ifndef NDEBUG
            if (messager::master())
            {
                std::cout << "Setting boundary conditions:\n";
                std::cout << "\tdomain name:   " << domain->name() << "\n";
                std::cout << "\tdomain index:  " << domain->index() << "\n";
                std::cout << "\tpatch index:   " << iBoundary << "\n";
                std::cout << "\tphase index:  " << iPhase << "\n";
                std::cout << "\tBoundary type: " << ::accel::toString(bc_type)
                          << "\n";
                std::cout << "\tYAML values:\n" << fluidValues << "\n\n";
            }
#endif /* NDEBUG */

            switch (bc_type)
            {
                case boundaryPhysicalType::inlet:
                    {
                        const auto& fluidValuesForMaterialNode =
                            fluidValues[realmPtr_->simulationRef().materialName(
                                iPhase)];
                        if (fluidValuesForMaterialNode["volume_fraction"])
                        {
                            const auto& volumeFractionNode =
                                fluidValuesForMaterialNode["volume_fraction"];
                            std::string option =
                                volumeFractionNode["option"]
                                    .template as<std::string>();
                            if (option == "value")
                            {
                                bc.setType(
                                    boundaryConditionType::specifiedValue);
                                bc.query<1>(volumeFractionNode,
                                            "value",
                                            "volume_fraction");
                                alphaRef(iPhase).registerSideFields(
                                    domain->index(), iBoundary);
                            }
                            else
                            {
                                errorMsg("Unrecognized volume fraction option "
                                         "for boundary");
                            }
                        }
                        else
                        {
                            errorMsg("phase specification at inlet must be "
                                     "specified");
                        }
                    }
                    break;

                case boundaryPhysicalType::opening:
                    {
                        const auto& fluidValuesForMaterialNode =
                            fluidValues[realmPtr_->simulationRef().materialName(
                                iPhase)];
                        if (fluidValuesForMaterialNode["volume_fraction"])
                        {
                            const auto& volumeFractionNode =
                                fluidValuesForMaterialNode["volume_fraction"];
                            std::string option =
                                volumeFractionNode["option"]
                                    .template as<std::string>();

                            if (option == "value")
                            {
                                // setting specified value here is only for the
                                // purpose of updating the side fields,
                                // otherwise the boundary condition of volume
                                // fraction is not specifid in fact, but
                                // zero-gradient while specified only in case of
                                // a reversed flow
                                bc.setType(
                                    boundaryConditionType::specifiedValue);
                                bc.query<1>(volumeFractionNode,
                                            "value",
                                            "volume_fraction");
                                alphaRef(iPhase).registerSideFields(
                                    domain->index(), iBoundary);
                            }
                            else
                            {
                                errorMsg("Unrecognized volume fraction option "
                                         "for boundary");
                            }
                        }
                        else
                        {
                            errorMsg("phase specification at inlet must be "
                                     "specified");
                        }
                    }
                    break;

                case boundaryPhysicalType::wall:
                case boundaryPhysicalType::symmetry:
                    bc.setType(boundaryConditionType::zeroGradient);
                    break;

                default:
                    break;
            }
        });
    }
}

void fieldBroker::setupDensity(const std::shared_ptr<domain> domain,
                               label iPhase)
{
    if (realmPtr_->rhoVector_[iPhase])
    {
        label localPhaseIndexInDomain =
            domain->globalToLocalMaterialIndex(iPhase);

        SET_DENSITY_PROPERTY(this->rhoRef(iPhase),
                             domain->getYAMLMaterial(localPhaseIndexInDomain),
                             localPhaseIndexInDomain);
    }
}

void fieldBroker::setupDynamicViscosity(const std::shared_ptr<domain> domain,
                                        label iPhase)
{
    if (realmPtr_->muVector_[iPhase])
    {
        label localPhaseIndexInDomain =
            domain->globalToLocalMaterialIndex(iPhase);

        SET_DYNAMIC_VISCOSITY_PROPERTY(
            this->muRef(iPhase),
            domain->getYAMLMaterial(localPhaseIndexInDomain),
            localPhaseIndexInDomain);
    }
}

void fieldBroker::setupSpecificHeatCapacity(
    const std::shared_ptr<domain> domain,
    label iPhase)
{
    if (realmPtr_->cpVector_[iPhase])
    {
        label localPhaseIndexInDomain =
            domain->globalToLocalMaterialIndex(iPhase);

        SET_SPECIFIC_HEAT_CAPACITY_PROPERTY(
            this->cpRef(iPhase),
            domain->getYAMLMaterial(localPhaseIndexInDomain),
            localPhaseIndexInDomain);
    }
}

void fieldBroker::setupThermalConductivity(const std::shared_ptr<domain> domain,
                                           label iPhase)
{
    if (realmPtr_->lambdaVector_[iPhase])
    {
        label localPhaseIndexInDomain =
            domain->globalToLocalMaterialIndex(iPhase);

        SET_THERMAL_CONDUCTIVITY_PROPERTY(
            this->lambdaRef(iPhase),
            domain->getYAMLMaterial(localPhaseIndexInDomain),
            localPhaseIndexInDomain);
    }
}

void fieldBroker::setupThermalExpansivity(const std::shared_ptr<domain> domain,
                                          label iPhase)
{
    if (realmPtr_->betaVector_[iPhase])
    {
        label localPhaseIndexInDomain =
            domain->globalToLocalMaterialIndex(iPhase);

        SET_THERMAL_EXPANSIVITY_PROPERTY(
            this->betaRef(iPhase),
            domain->getYAMLMaterial(localPhaseIndexInDomain),
            localPhaseIndexInDomain);
    }
}

void fieldBroker::setupCompressibility(const std::shared_ptr<domain> domain,
                                       label iPhase)
{
    if (realmPtr_->psiVector_[iPhase])
    {
        psiRef(iPhase).setZone(domain->index());
    }
}

void fieldBroker::setupMassFlowRate(const std::shared_ptr<domain> domain,
                                    label iPhase)
{
    if (this->mDotRef(iPhase).isZoneUnset(domain->index()))
    {
        this->mDotRef(iPhase).setZone(domain->index());
        this->mDotRef(iPhase).divRef().setZone(domain->index());

#ifdef HAS_INTERFACE
        // register mass flux side field for interfaces in fluid domain
        for (interface* interf : domain->interfacesRef())
        {
            if (interf->isInternal())
            {
                this->mDotRef(iPhase).registerSideFieldsForInterfaceSide(
                    interf->index(), true);
                this->mDotRef(iPhase).registerSideFieldsForInterfaceSide(
                    interf->index(), false);
            }
            else
            {
                if (interf->isFluidSolidType())
                {
                    // do nothing
                }
                else
                {
                    this->mDotRef(iPhase).registerSideFieldsForInterfaceSide(
                        interf->index(), interf->isMasterZone(domain->index()));
                }
            }
        }
#endif /* HAS_INTERFACE */

        // Register mass flux side field for patches in fluid domain
        for (label iBoundary = 0; iBoundary < domain->zonePtr()->nBoundaries();
             iBoundary++)
        {
            const auto& boundaryRef = domain->zonePtr()->boundaryRef(iBoundary);

            boundaryPhysicalType type = boundaryRef.type();

            switch (type)
            {
                case boundaryPhysicalType::inlet:
                case boundaryPhysicalType::outlet:
                case boundaryPhysicalType::opening:
                    this->mDotRef(iPhase).registerSideField(domain->index(),
                                                            iBoundary);
                    break;

                default:
                    break;
            }
        }
    }
}

void fieldBroker::initializeVolumeFraction(const std::shared_ptr<domain> domain,
                                           label iPhase)
{
    if (realmPtr_->alphaVector_[iPhase])
    {
        alphaRef(iPhase).initialize(domain->index());
    }
}

void fieldBroker::initializeDensity(const std::shared_ptr<domain> domain,
                                    label iPhase)
{
    if (realmPtr_->rhoVector_[iPhase])
    {
        rhoRef(iPhase).initialize(domain->index());
    }
}

void fieldBroker::initializeDynamicViscosity(
    const std::shared_ptr<domain> domain,
    label iPhase)
{
    if (realmPtr_->muVector_[iPhase])
    {
        muRef(iPhase).initialize(domain->index());
    }
}

void fieldBroker::initializeSpecificHeatCapacity(
    const std::shared_ptr<domain> domain,
    label iPhase)
{
    if (realmPtr_->cpVector_[iPhase])
    {
        cpRef(iPhase).initialize(domain->index());
    }
}

void fieldBroker::initializeThermalConductivity(
    const std::shared_ptr<domain> domain,
    label iPhase)
{
    if (realmPtr_->lambdaVector_[iPhase])
    {
        lambdaRef(iPhase).initialize(domain->index());
    }
}

void fieldBroker::initializeThermalExpansivity(
    const std::shared_ptr<domain> domain,
    label iPhase)
{
    if (realmPtr_->betaVector_[iPhase])
    {
        betaRef(iPhase).initialize(domain->index());
    }
}

void fieldBroker::initializeCompressibility(
    const std::shared_ptr<domain> domain,
    label iPhase)
{
    if (realmPtr_->psiVector_[iPhase])
    {
        psiRef(iPhase).initialize(domain->index());
    }
}

void fieldBroker::initializeMassFlowRate(const std::shared_ptr<domain> domain,
                                         label iPhase)
{
    // interior
    initializeMassFlowRateInterior_(domain, iPhase);

#ifdef HAS_INTERFACE
    // Interfaces
    for (const interface* interf : domain->zonePtr()->interfacesRef())
    {
        if (interf->isInternal())
        {
            initializeMassFlowRateInterfaceSideField_(
                domain, interf->masterInfoPtr(), iPhase);
            initializeMassFlowRateInterfaceSideField_(
                domain, interf->slaveInfoPtr(), iPhase);
        }
        else if (!interf->isFluidSolidType())
        {
            // get interface side that is sitting in this domain
            const auto* interfaceSideInfoPtr =
                interf->interfaceSideInfoPtr(domain->index());

            initializeMassFlowRateInterfaceSideField_(
                domain, interfaceSideInfoPtr, iPhase);
        }
    }
#endif /* HAS_INTERFACE */

    // Boundary
    for (label iBoundary = 0; iBoundary < domain->zonePtr()->nBoundaries();
         iBoundary++)
    {
        initializeMassFlowRateBoundaryField_(
            domain, domain->zonePtr()->boundaryPtr(iBoundary), iPhase);
    }
}

void fieldBroker::resetVolumeFraction(const std::shared_ptr<domain> domain,
                                      label iPhase)
{
    if (realmPtr_->alphaVector_[iPhase])
    {
        alphaRef(iPhase).initialize(domain->index(), true);
    }
}

void fieldBroker::resetDensity(const std::shared_ptr<domain> domain,
                               label iPhase)
{
    if (realmPtr_->rhoVector_[iPhase])
    {
        rhoRef(iPhase).initialize(domain->index(), true);
    }
}

void fieldBroker::resetDynamicViscosity(const std::shared_ptr<domain> domain,
                                        label iPhase)
{
    if (realmPtr_->muVector_[iPhase])
    {
        muRef(iPhase).initialize(domain->index(), true);
    }
}

void fieldBroker::resetSpecificHeatCapacity(
    const std::shared_ptr<domain> domain,
    label iPhase)
{
    if (realmPtr_->cpVector_[iPhase])
    {
        cpRef(iPhase).initialize(domain->index(), true);
    }
}

void fieldBroker::resetThermalConductivity(const std::shared_ptr<domain> domain,
                                           label iPhase)
{
    if (realmPtr_->lambdaVector_[iPhase])
    {
        lambdaRef(iPhase).initialize(domain->index(), true);
    }
}

void fieldBroker::resetThermalExpansivity(const std::shared_ptr<domain> domain,
                                          label iPhase)
{
    if (realmPtr_->betaVector_[iPhase])
    {
        betaRef(iPhase).initialize(domain->index(), true);
    }
}

void fieldBroker::resetCompressibility(const std::shared_ptr<domain> domain,
                                       label iPhase)
{
    if (realmPtr_->psiVector_[iPhase])
    {
        psiRef(iPhase).initialize(domain->index(), true);
    }
}

void fieldBroker::updateVolumeFraction(const std::shared_ptr<domain> domain,
                                       label iPhase)
{
    if (realmPtr_->alphaVector_[iPhase])
    {
        alphaRef(iPhase).update(domain->index());
    }
}

void fieldBroker::updateDensity(const std::shared_ptr<domain> domain,
                                label iPhase)
{
    if (realmPtr_->rhoVector_[iPhase])
    {
        rhoRef(iPhase).update(domain->index());
    }
}

void fieldBroker::updateDynamicViscosity(const std::shared_ptr<domain> domain,
                                         label iPhase)
{
    if (realmPtr_->muVector_[iPhase])
    {
        muRef(iPhase).update(domain->index());
    }
}

void fieldBroker::updateSpecificHeatCapacity(
    const std::shared_ptr<domain> domain,
    label iPhase)
{
    if (realmPtr_->cpVector_[iPhase])
    {
        cpRef(iPhase).update(domain->index());
    }
}

void fieldBroker::updateThermalConductivity(
    const std::shared_ptr<domain> domain,
    label iPhase)
{
    if (realmPtr_->lambdaVector_[iPhase])
    {
        lambdaRef(iPhase).update(domain->index());
    }
}

void fieldBroker::updateThermalExpansivity(const std::shared_ptr<domain> domain,
                                           label iPhase)
{
    if (realmPtr_->betaVector_[iPhase])
    {
        betaRef(iPhase).update(domain->index());
    }
}

void fieldBroker::updateCompressibility(const std::shared_ptr<domain> domain,
                                        label iPhase)
{
    if (realmPtr_->psiVector_[iPhase])
    {
        psiRef(iPhase).update(domain->index());
    }
}

void fieldBroker::updateMassFlowRate(const std::shared_ptr<domain> domain,
                                     label iPhase)
{
    // interior
    updateMassFlowRateInterior_(domain, iPhase);

#ifdef HAS_INTERFACE
    // Interfaces
    for (const interface* interf : domain->interfacesRef())
    {
        if (interf->isInternal())
        {
            updateMassFlowRateInterfaceSideField_(
                domain, interf->masterInfoPtr(), iPhase);
            updateMassFlowRateInterfaceSideField_(
                domain, interf->slaveInfoPtr(), iPhase);
        }
        else if (!interf->isFluidSolidType())
        {
            // get interface side that is sitting in this domain
            const auto* interfaceSideInfoPtr =
                interf->interfaceSideInfoPtr(domain->index());

            updateMassFlowRateInterfaceSideField_(
                domain, interfaceSideInfoPtr, iPhase);
        }
    }
#endif /* HAS_INTERFACE */

    // Boundary
    for (label iBoundary = 0; iBoundary < domain->zonePtr()->nBoundaries();
         iBoundary++)
    {
        updateMassFlowRateBoundaryField_(
            domain, domain->zonePtr()->boundaryPtr(iBoundary), iPhase);
    }
}

void fieldBroker::updateVolumeFractionGradientField(
    const std::shared_ptr<domain> domain,
    label iPhase)
{
    if (realmPtr_->alphaVector_[iPhase])
    {
        alphaRef(iPhase).updateGradientField(domain->index());
    }
}

void fieldBroker::updateDensityGradientField(
    const std::shared_ptr<domain> domain,
    label iPhase)
{
    if (realmPtr_->rhoVector_[iPhase])
    {
        rhoRef(iPhase).updateGradientField(domain->index());
    }
}

void fieldBroker::updateDynamicViscosityGradientField(
    const std::shared_ptr<domain> domain,
    label iPhase)
{
    if (realmPtr_->muVector_[iPhase])
    {
        muRef(iPhase).updateGradientField(domain->index());
    }
}

void fieldBroker::updateSpecificHeatCapacityGradientField(
    const std::shared_ptr<domain> domain,
    label iPhase)
{
    if (realmPtr_->cpVector_[iPhase])
    {
        cpRef(iPhase).updateGradientField(domain->index());
    }
}

void fieldBroker::updateThermalConductivityGradientField(
    const std::shared_ptr<domain> domain,
    label iPhase)
{
    if (realmPtr_->lambdaVector_[iPhase])
    {
        lambdaRef(iPhase).updateGradientField(domain->index());
    }
}

void fieldBroker::updateVolumeFractionBlendingFactorField(
    const std::shared_ptr<domain> domain,
    label iPhase)
{
    if (realmPtr_->alphaVector_[iPhase])
    {
        alphaRef(iPhase).updateBlendingFactorField(domain->index());
    }
}

void fieldBroker::updateDensityBlendingFactorField(
    const std::shared_ptr<domain> domain,
    label iPhase)
{
    if (realmPtr_->rhoVector_[iPhase])
    {
        rhoRef(iPhase).updateBlendingFactorField(domain->index());
    }
}

void fieldBroker::updateDynamicViscosityBlendingFactorField(
    const std::shared_ptr<domain> domain,
    label iPhase)
{
    if (realmPtr_->muVector_[iPhase])
    {
        muRef(iPhase).updateBlendingFactorField(domain->index());
    }
}

void fieldBroker::updateSpecificHeatCapacityBlendingFactorField(
    const std::shared_ptr<domain> domain,
    label iPhase)
{
    if (realmPtr_->cpVector_[iPhase])
    {
        cpRef(iPhase).updateBlendingFactorField(domain->index());
    }
}

void fieldBroker::updateThermalConductivityBlendingFactorField(
    const std::shared_ptr<domain> domain,
    label iPhase)
{
    if (realmPtr_->lambdaVector_[iPhase])
    {
        lambdaRef(iPhase).updateBlendingFactorField(domain->index());
    }
}

void fieldBroker::updateVolumeFractionPrevIterField(
    const std::shared_ptr<domain> domain,
    label iPhase)
{
    if (realmPtr_->alphaVector_[iPhase])
    {
        alphaRef(iPhase).updatePrevIterField(domain->index());
    }
}

void fieldBroker::updateDensityPrevIterField(
    const std::shared_ptr<domain> domain,
    label iPhase)
{
    if (realmPtr_->rhoVector_[iPhase])
    {
        rhoRef(iPhase).updatePrevIterField(domain->index());
    }
}

void fieldBroker::updateDynamicViscosityPrevIterField(
    const std::shared_ptr<domain> domain,
    label iPhase)
{
    if (realmPtr_->muVector_[iPhase])
    {
        muRef(iPhase).updatePrevIterField(domain->index());
    }
}

void fieldBroker::updateSpecificHeatCapacityPrevIterField(
    const std::shared_ptr<domain> domain,
    label iPhase)
{
    if (realmPtr_->cpVector_[iPhase])
    {
        cpRef(iPhase).updatePrevIterField(domain->index());
    }
}

void fieldBroker::updateThermalConductivityPrevIterField(
    const std::shared_ptr<domain> domain,
    label iPhase)
{
    if (realmPtr_->lambdaVector_[iPhase])
    {
        lambdaRef(iPhase).updatePrevIterField(domain->index());
    }
}

void fieldBroker::updateVolumeFractionPrevTimeField(
    const std::shared_ptr<domain> domain,
    label iPhase)
{
    if (realmPtr_->alphaVector_[iPhase])
    {
        alphaRef(iPhase).updatePrevTimeField(domain->index());
    }
}

void fieldBroker::updateDensityPrevTimeField(
    const std::shared_ptr<domain> domain,
    label iPhase)
{
    if (realmPtr_->rhoVector_[iPhase])
    {
        rhoRef(iPhase).updatePrevTimeField(domain->index());
    }
}

void fieldBroker::updateDynamicViscosityPrevTimeField(
    const std::shared_ptr<domain> domain,
    label iPhase)
{
    if (realmPtr_->muVector_[iPhase])
    {
        muRef(iPhase).updatePrevTimeField(domain->index());
    }
}

void fieldBroker::updateSpecificHeatCapacityPrevTimeField(
    const std::shared_ptr<domain> domain,
    label iPhase)
{
    if (realmPtr_->cpVector_[iPhase])
    {
        cpRef(iPhase).updatePrevTimeField(domain->index());
    }
}

void fieldBroker::updateThermalConductivityPrevTimeField(
    const std::shared_ptr<domain> domain,
    label iPhase)
{
    if (realmPtr_->lambdaVector_[iPhase])
    {
        lambdaRef(iPhase).updatePrevTimeField(domain->index());
    }
}

// Public access

density& fieldBroker::rhoRef(label iPhase)
{
    RETURN_REF_ALLOC_MATERIAL(
        realmPtr_->rhoVector_,
        iPhase,
        density,
        realmPtr_,
        realm::rho_ID + std::string(".") +
            realmPtr_->simulationRef().materialName(iPhase),
        n_states,
        high_res,
        true);
}

const density& fieldBroker::rhoRef(label iPhase) const
{
    RETURN_REF_ALLOC_MATERIAL(
        realmPtr_->rhoVector_,
        iPhase,
        density,
        realmPtr_,
        realm::rho_ID + std::string(".") +
            realmPtr_->simulationRef().materialName(iPhase),
        n_states,
        high_res,
        true);
}

massFlowRate& fieldBroker::mDotRef(label iPhase)
{
    RETURN_REF_ALLOC_MATERIAL(
        realmPtr_->mDotVector_,
        iPhase,
        massFlowRate,
        realmPtr_,
        realm::mDot_ID + std::string(".") +
            realmPtr_->simulationRef().materialName(iPhase),
        n_states);
}

const massFlowRate& fieldBroker::mDotRef(label iPhase) const
{
    RETURN_REF_ALLOC_MATERIAL(
        realmPtr_->mDotVector_,
        iPhase,
        massFlowRate,
        realmPtr_,
        realm::mDot_ID + std::string(".") +
            realmPtr_->simulationRef().materialName(iPhase),
        n_states);
}

// Protected access

volumeFraction& fieldBroker::alphaRef(label iPhase)
{
    RETURN_REF_ALLOC_MATERIAL(
        realmPtr_->alphaVector_,
        iPhase,
        volumeFraction,
        realmPtr_,
        realm::alpha_ID + std::string(".") +
            realmPtr_->simulationRef().materialName(iPhase),
        n_states,
        high_res);
}

const volumeFraction& fieldBroker::alphaRef(label iPhase) const
{
    RETURN_REF_ALLOC_MATERIAL(
        realmPtr_->alphaVector_,
        iPhase,
        volumeFraction,
        realmPtr_,
        realm::alpha_ID + std::string(".") +
            realmPtr_->simulationRef().materialName(iPhase),
        n_states,
        high_res);
}

simpleScalarField& fieldBroker::alphaSmoothRef(label iPhase)
{
    RETURN_REF_ALLOC_AUX_MATERIAL(
        realmPtr_->alphaSmoothVector_,
        iPhase,
        simpleScalarField,
        realmPtr_,
        realm::alpha_smooth_ID + std::string(".") +
            realmPtr_->simulationRef().materialName(iPhase),
        stk::topology::NODE_RANK);
}

const simpleScalarField& fieldBroker::alphaSmoothRef(label iPhase) const
{
    RETURN_REF_ALLOC_AUX_MATERIAL(
        realmPtr_->alphaSmoothVector_,
        iPhase,
        simpleScalarField,
        realmPtr_,
        realm::alpha_smooth_ID + std::string(".") +
            realmPtr_->simulationRef().materialName(iPhase),
        stk::topology::NODE_RANK);
}

simpleScalarField& fieldBroker::rhsSmoothRef(label iPhase)
{
    RETURN_REF_ALLOC_AUX_MATERIAL(
        realmPtr_->rhsSmoothVector_,
        iPhase,
        simpleScalarField,
        realmPtr_,
        realm::rhs_smooth_ID + std::string(".") +
            realmPtr_->simulationRef().materialName(iPhase),
        stk::topology::NODE_RANK);
}

const simpleScalarField& fieldBroker::rhsSmoothRef(label iPhase) const
{
    RETURN_REF_ALLOC_AUX_MATERIAL(
        realmPtr_->rhsSmoothVector_,
        iPhase,
        simpleScalarField,
        realmPtr_,
        realm::rhs_smooth_ID + std::string(".") +
            realmPtr_->simulationRef().materialName(iPhase),
        stk::topology::NODE_RANK);
}

dynamicViscosity& fieldBroker::muRef(label iPhase)
{
    RETURN_REF_ALLOC_MATERIAL(
        realmPtr_->muVector_,
        iPhase,
        dynamicViscosity,
        realmPtr_,
        realm::mu_ID + std::string(".") +
            realmPtr_->simulationRef().materialName(iPhase),
        n_states,
        false,
        false);
}

const dynamicViscosity& fieldBroker::muRef(label iPhase) const
{
    RETURN_REF_ALLOC_MATERIAL(
        realmPtr_->muVector_,
        iPhase,
        dynamicViscosity,
        realmPtr_,
        realm::mu_ID + std::string(".") +
            realmPtr_->simulationRef().materialName(iPhase),
        n_states,
        false,
        false);
}

specificHeatCapacity& fieldBroker::cpRef(label iPhase)
{
    RETURN_REF_ALLOC_MATERIAL(
        realmPtr_->cpVector_,
        iPhase,
        specificHeatCapacity,
        realmPtr_,
        realm::cp_ID + std::string(".") +
            realmPtr_->simulationRef().materialName(iPhase),
        n_states,
        false,
        false);
}

const specificHeatCapacity& fieldBroker::cpRef(label iPhase) const
{
    RETURN_REF_ALLOC_MATERIAL(
        realmPtr_->cpVector_,
        iPhase,
        specificHeatCapacity,
        realmPtr_,
        realm::cp_ID + std::string(".") +
            realmPtr_->simulationRef().materialName(iPhase),
        n_states,
        false,
        false);
}

thermalConductivity& fieldBroker::lambdaRef(label iPhase)
{
    RETURN_REF_ALLOC_MATERIAL(
        realmPtr_->lambdaVector_,
        iPhase,
        thermalConductivity,
        realmPtr_,
        realm::lambda_ID + std::string(".") +
            realmPtr_->simulationRef().materialName(iPhase),
        n_states,
        false,
        false);
}

const thermalConductivity& fieldBroker::lambdaRef(label iPhase) const
{
    RETURN_REF_ALLOC_MATERIAL(
        realmPtr_->lambdaVector_,
        iPhase,
        thermalConductivity,
        realmPtr_,
        realm::lambda_ID + std::string(".") +
            realmPtr_->simulationRef().materialName(iPhase),
        n_states,
        false,
        false);
}

thermalExpansivity& fieldBroker::betaRef(label iPhase)
{
    RETURN_REF_ALLOC_MATERIAL(
        realmPtr_->betaVector_,
        iPhase,
        thermalExpansivity,
        realmPtr_,
        realm::beta_ID + std::string(".") +
            realmPtr_->simulationRef().materialName(iPhase),
        n_states,
        false,
        false);
}

const thermalExpansivity& fieldBroker::betaRef(label iPhase) const
{
    RETURN_REF_ALLOC_MATERIAL(
        realmPtr_->betaVector_,
        iPhase,
        thermalExpansivity,
        realmPtr_,
        realm::beta_ID + std::string(".") +
            realmPtr_->simulationRef().materialName(iPhase),
        n_states,
        false,
        false);
}

const compressibility& fieldBroker::psiRef(label iPhase) const
{
    RETURN_REF_ALLOC_MATERIAL(
        realmPtr_->psiVector_,
        iPhase,
        compressibility,
        realmPtr_,
        realm::psi_ID + std::string(".") +
            realmPtr_->simulationRef().materialName(iPhase),
        n_states,
        false,
        false);
}

compressibility& fieldBroker::psiRef(label iPhase)
{
    RETURN_REF_ALLOC_MATERIAL(
        realmPtr_->psiVector_,
        iPhase,
        compressibility,
        realmPtr_,
        realm::psi_ID + std::string(".") +
            realmPtr_->simulationRef().materialName(iPhase),
        n_states,
        false,
        false);
}

simpleVectorField& fieldBroker::nHatRef(label iPhase)
{
    RETURN_REF_ALLOC_AUX_MATERIAL(
        realmPtr_->nHatVector_,
        iPhase,
        simpleVectorField,
        realmPtr_,
        realm::nHat_ID + std::string(".") +
            realmPtr_->simulationRef().materialName(iPhase),
        stk::topology::NODE_RANK);
};

const simpleVectorField& fieldBroker::nHatRef(label iPhase) const
{
    RETURN_REF_ALLOC_AUX_MATERIAL(
        realmPtr_->nHatVector_,
        iPhase,
        simpleVectorField,
        realmPtr_,
        realm::nHat_ID + std::string(".") +
            realmPtr_->simulationRef().materialName(iPhase),
        stk::topology::NODE_RANK);
};

void fieldBroker::setupBoundaryConditions_(
    const std::shared_ptr<domain> domain,
    const label iPhase,
    std::function<void(const accel::domain* domain,
                       const label iBoundary,
                       const label iPhase,
                       const boundaryPhysicalType bc_type,
                       const YAML::Node fluidValues)> setupBC)
{
    const YAML::Node bc_list = domain->getYAMLBoundaryConditions();
    for (label iBoundary = 0; iBoundary < bc_list.size(); iBoundary++)
    {
        const YAML::Node patch = bc_list[iBoundary];
        if (!patch["type"])
        {
            errorMsg("fieldBroker: boundary patch " +
                     std::to_string(iBoundary) + " for domain `" +
                     domain->name() + "` does not define a `type`");
        }
        const boundaryPhysicalType boundary_type =
            convertBoundaryPhysicalTypeFromString(
                patch["type"].template as<std::string>());
        YAML::Node bc_values;
        if (patch["fluid_values"])
        {
            bc_values = patch["fluid_values"];
        }

        setupBC(domain.get(), iBoundary, iPhase, boundary_type, bc_values);
    }
}

} /* namespace accel */
