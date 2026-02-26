// File       : fieldBroker.h
// Created    : Tue Feb 20 2024 12:55:24 (+0100)
// Author     : Fabian Wermelinger
// Description: Field broker is the interface between a realm and the fields
// owned by the realm
// Copyright (c) 2024 CCFNUM, Lucerne University of Applied Sciences and Arts.
// SPDX-License-Identifier: BSD-3-Clause

#ifndef FIELDBROKER_H
#define FIELDBROKER_H

// code
#include "controls.h"
#include "domain.h"
#include "realm.h"

namespace accel
{

class fieldBroker
{
public:
    // Constructors

    fieldBroker() = delete;

    fieldBroker(realm* realm);

    virtual ~fieldBroker() = default;

    // Access

    controls& controlsRef()
    {
        return realmPtr_->meshRef().controlsRef();
    }

    const controls& controlsRef() const
    {
        return realmPtr_->meshRef().controlsRef();
    }

    mesh& meshRef()
    {
        return realmPtr_->meshRef();
    }

    const mesh& meshRef() const
    {
        return realmPtr_->meshRef();
    }

    realm& realmRef()
    {
        return *realmPtr_;
    }

    const realm& realmRef() const
    {
        return *realmPtr_;
    }

    // public API for transport fields
    velocity& URef();

    const velocity& URef() const;

    simpleVectorField& UrRef();

    const simpleVectorField& UrRef() const;

    density& rhoRef();

    const density& rhoRef() const;

    massFlowRate& mDotRef();

    const massFlowRate& mDotRef() const;

    nodeVectorField& UmRef();

    const nodeVectorField& UmRef() const;

    nodeScalarField& divUmRef();

    const nodeScalarField& divUmRef() const;

    // public API for phasic transport fields (iPhase is global to the
    // simulation: i.e. material index)

    density& rhoRef(label iPhase);

    const density& rhoRef(label iPhase) const;

    massFlowRate& mDotRef(label iPhase);

    const massFlowRate& mDotRef(label iPhase) const;

    // public API for fluid transport
    // TODO: Provide convenience methods for transport
    // field initialization and facilitation of code re-usage
    //
    // Having further separation for fluidBroker (and solidBroker) would make
    // more clear that transport is a phenomenon that occurs in fluids only.
    // Can be decided upon later -- this implementation is consistent with
    // discussion 2024-05-24

    // public API for field setup

    virtual void setupVelocity(const std::shared_ptr<domain> domain);

    virtual void setupPressure(const std::shared_ptr<domain> domain);

    virtual void setupTemperature(const std::shared_ptr<domain> domain);

    virtual void setupSpecificEnthalpy(const std::shared_ptr<domain> domain);

    virtual void
    setupSpecificTotalEnthalpy(const std::shared_ptr<domain> domain);

    virtual void setupWallScale(const std::shared_ptr<domain> domain);

    virtual void
    setupTurbulentKineticEnergy(const std::shared_ptr<domain> domain);

    virtual void
    setupTurbulentEddyFrequency(const std::shared_ptr<domain> domain);

    virtual void
    setupTurbulentDissipationRate(const std::shared_ptr<domain> domain);

    virtual void
    setupTransitionOnsetReynoldsNumber(const std::shared_ptr<domain> domain);

    virtual void
    setupTurbulentIntermittency(const std::shared_ptr<domain> domain);

    virtual void setupDensity(const std::shared_ptr<domain> domain);

    virtual void setupDynamicViscosity(const std::shared_ptr<domain> domain);

    virtual void
    setupSpecificHeatCapacity(const std::shared_ptr<domain> domain);

    virtual void setupThermalConductivity(const std::shared_ptr<domain> domain);

    virtual void setupThermalExpansivity(const std::shared_ptr<domain> domain);

    virtual void setupCompressibility(const std::shared_ptr<domain> domain);

    virtual void setupDisplacement(const std::shared_ptr<domain> domain);

    virtual void setupTotalDisplacement(const std::shared_ptr<domain> domain);

    virtual void setupMeshVelocity(const std::shared_ptr<domain> domain);

    virtual void setupDivMeshVelocity(const std::shared_ptr<domain> domain);

    virtual void setupMassFlowRate(const std::shared_ptr<domain> domain);

    virtual void setupHeatFlowRate(const std::shared_ptr<domain> domain);

    virtual void setupYoungModulus(const std::shared_ptr<domain> domain);

    virtual void setupPoissonRatio(const std::shared_ptr<domain> domain);

    // public API for phasic field setup (iPhase is global to the simulation:
    // i.e. material index)

    virtual void setupVolumeFraction(const std::shared_ptr<domain> domain,
                                     label iPhase);

    virtual void setupDensity(const std::shared_ptr<domain> domain,
                              label iPhase);

    virtual void setupDynamicViscosity(const std::shared_ptr<domain> domain,
                                       label iPhase);

    virtual void setupSpecificHeatCapacity(const std::shared_ptr<domain> domain,
                                           label iPhase);

    virtual void setupThermalConductivity(const std::shared_ptr<domain> domain,
                                          label iPhase);

    virtual void setupThermalExpansivity(const std::shared_ptr<domain> domain,
                                         label iPhase);

    virtual void setupCompressibility(const std::shared_ptr<domain> domain,
                                      label iPhase);

    virtual void setupMassFlowRate(const std::shared_ptr<domain> domain,
                                   label iPhase);

    // public API for field raw initialization

    virtual void initializeVelocity(const std::shared_ptr<domain> domain);

    virtual void initializePressure(const std::shared_ptr<domain> domain);

    virtual void initializeTemperature(const std::shared_ptr<domain> domain);

    virtual void
    initializeSpecificEnthalpy(const std::shared_ptr<domain> domain);

    virtual void
    initializeSpecificTotalEnthalpy(const std::shared_ptr<domain> domain);

    virtual void initializeWallScale(const std::shared_ptr<domain> domain);

    virtual void
    initializeTurbulentKineticEnergy(const std::shared_ptr<domain> domain);

    virtual void
    initializeTurbulentEddyFrequency(const std::shared_ptr<domain> domain);

    virtual void
    initializeTurbulentDissipationRate(const std::shared_ptr<domain> domain);

    virtual void initializeTransitionOnsetReynoldsNumber(
        const std::shared_ptr<domain> domain);

    virtual void
    initializeTurbulentIntermittency(const std::shared_ptr<domain> domain);

    virtual void initializeYoungModulus(const std::shared_ptr<domain> domain);

    virtual void initializePoissonRatio(const std::shared_ptr<domain> domain);

    virtual void initializeDensity(const std::shared_ptr<domain> domain);

    virtual void
    initializeDynamicViscosity(const std::shared_ptr<domain> domain);

    virtual void
    initializeSpecificHeatCapacity(const std::shared_ptr<domain> domain);

    virtual void
    initializeThermalConductivity(const std::shared_ptr<domain> domain);

    virtual void
    initializeThermalExpansivity(const std::shared_ptr<domain> domain);

    virtual void
    initializeCompressibility(const std::shared_ptr<domain> domain);

    virtual void initializeDisplacement(const std::shared_ptr<domain> domain);

    virtual void
    initializeTotalDisplacement(const std::shared_ptr<domain> domain);

    virtual void initializeMeshVelocity(const std::shared_ptr<domain> domain);

    virtual void
    initializeDivMeshVelocity(const std::shared_ptr<domain> domain);

    virtual void initializeMassFlowRate(const std::shared_ptr<domain> domain);

    // public API for phasic field raw initialization (iPhase is global to the
    // simulation: i.e. material index)

    virtual void initializeVolumeFraction(const std::shared_ptr<domain> domain,
                                          label iPhase);

    virtual void initializeDensity(const std::shared_ptr<domain> domain,
                                   label iPhase);

    virtual void
    initializeDynamicViscosity(const std::shared_ptr<domain> domain,
                               label iPhase);

    virtual void
    initializeSpecificHeatCapacity(const std::shared_ptr<domain> domain,
                                   label iPhase);

    virtual void
    initializeThermalConductivity(const std::shared_ptr<domain> domain,
                                  label iPhase);

    virtual void
    initializeThermalExpansivity(const std::shared_ptr<domain> domain,
                                 label iPhase);

    virtual void initializeCompressibility(const std::shared_ptr<domain> domain,
                                           label iPhase);

    virtual void initializeMassFlowRate(const std::shared_ptr<domain> domain,
                                        label iPhase);

    // public API for field raw resetting

    virtual void resetVelocity(const std::shared_ptr<domain> domain);

    virtual void resetPressure(const std::shared_ptr<domain> domain);

    virtual void resetTemperature(const std::shared_ptr<domain> domain);

    virtual void resetSpecificEnthalpy(const std::shared_ptr<domain> domain);

    virtual void
    resetSpecificTotalEnthalpy(const std::shared_ptr<domain> domain);

    virtual void resetWallScale(const std::shared_ptr<domain> domain);

    virtual void
    resetTurbulentKineticEnergy(const std::shared_ptr<domain> domain);

    virtual void
    resetTurbulentEddyFrequency(const std::shared_ptr<domain> domain);

    virtual void
    resetTurbulentDissipationRate(const std::shared_ptr<domain> domain);

    virtual void
    resetTransitionOnsetReynoldsNumber(const std::shared_ptr<domain> domain);

    virtual void
    resetTurbulentIntermittency(const std::shared_ptr<domain> domain);

    virtual void resetDensity(const std::shared_ptr<domain> domain);

    virtual void resetDynamicViscosity(const std::shared_ptr<domain> domain);

    virtual void
    resetSpecificHeatCapacity(const std::shared_ptr<domain> domain);

    virtual void resetThermalConductivity(const std::shared_ptr<domain> domain);

    virtual void resetThermalExpansivity(const std::shared_ptr<domain> domain);

    virtual void resetCompressibility(const std::shared_ptr<domain> domain);

    virtual void resetDisplacement(const std::shared_ptr<domain> domain);

    virtual void resetTotalDisplacement(const std::shared_ptr<domain> domain);

    virtual void resetMeshVelocity(const std::shared_ptr<domain> domain);

    virtual void resetDivMeshVelocity(const std::shared_ptr<domain> domain);

    // public API for phasic field raw initialization (iPhase is global to the
    // simulation: i.e. material index)

    virtual void resetVolumeFraction(const std::shared_ptr<domain> domain,
                                     label iPhase);

    virtual void resetDensity(const std::shared_ptr<domain> domain,
                              label iPhase);

    virtual void resetDynamicViscosity(const std::shared_ptr<domain> domain,
                                       label iPhase);

    virtual void resetSpecificHeatCapacity(const std::shared_ptr<domain> domain,
                                           label iPhase);

    virtual void resetThermalConductivity(const std::shared_ptr<domain> domain,
                                          label iPhase);

    virtual void resetThermalExpansivity(const std::shared_ptr<domain> domain,
                                         label iPhase);

    virtual void resetCompressibility(const std::shared_ptr<domain> domain,
                                      label iPhase);

    // public API for field raw update

    virtual void updateVelocity(const std::shared_ptr<domain> domain);

    virtual void updatePressure(const std::shared_ptr<domain> domain);

    virtual void updateTemperature(const std::shared_ptr<domain> domain);

    virtual void updateSpecificEnthalpy(const std::shared_ptr<domain> domain);

    virtual void
    updateSpecificTotalEnthalpy(const std::shared_ptr<domain> domain);

    virtual void updateWallScale(const std::shared_ptr<domain> domain);

    virtual void
    updateTurbulentKineticEnergy(const std::shared_ptr<domain> domain);

    virtual void
    updateTurbulentEddyFrequency(const std::shared_ptr<domain> domain);

    virtual void
    updateTurbulentDissipationRate(const std::shared_ptr<domain> domain);

    virtual void
    updateTransitionOnsetReynoldsNumber(const std::shared_ptr<domain> domain);

    virtual void
    updateTurbulentIntermittency(const std::shared_ptr<domain> domain);

    virtual void updateYoungModulus(const std::shared_ptr<domain> domain);

    virtual void updatePoissonRatio(const std::shared_ptr<domain> domain);

    virtual void updateDensity(const std::shared_ptr<domain> domain);

    virtual void updateDynamicViscosity(const std::shared_ptr<domain> domain);

    virtual void
    updateSpecificHeatCapacity(const std::shared_ptr<domain> domain);

    virtual void
    updateThermalConductivity(const std::shared_ptr<domain> domain);

    virtual void updateThermalExpansivity(const std::shared_ptr<domain> domain);

    virtual void updateCompressibility(const std::shared_ptr<domain> domain);

    virtual void updateDisplacement(const std::shared_ptr<domain> domain);

    virtual void updateTotalDisplacement(const std::shared_ptr<domain> domain);

    virtual void updateMeshVelocity(const std::shared_ptr<domain> domain);

    virtual void updateDivMeshVelocity(const std::shared_ptr<domain> domain);

    virtual void updateMassFlowRate(const std::shared_ptr<domain> domain);

    // public API for phasic field raw initialization (iPhase is global to the
    // simulation: i.e. material index)

    virtual void updateVolumeFraction(const std::shared_ptr<domain> domain,
                                      label iPhase);

    virtual void updateDensity(const std::shared_ptr<domain> domain,
                               label iPhase);

    virtual void updateDynamicViscosity(const std::shared_ptr<domain> domain,
                                        label iPhase);

    virtual void
    updateSpecificHeatCapacity(const std::shared_ptr<domain> domain,
                               label iPhase);

    virtual void updateThermalConductivity(const std::shared_ptr<domain> domain,
                                           label iPhase);

    virtual void updateThermalExpansivity(const std::shared_ptr<domain> domain,
                                          label iPhase);

    virtual void updateCompressibility(const std::shared_ptr<domain> domain,
                                       label iPhase);

    virtual void updateMassFlowRate(const std::shared_ptr<domain> domain,
                                    label iPhase);

    // gradient

    virtual void
    updateVelocityGradientField(const std::shared_ptr<domain> domain);

    virtual void
    updatePressureGradientField(const std::shared_ptr<domain> domain);

    virtual void
    updateTemperatureGradientField(const std::shared_ptr<domain> domain);

    virtual void
    updateSpecificEnthalpyGradientField(const std::shared_ptr<domain> domain);

    virtual void updateSpecificTotalEnthalpyGradientField(
        const std::shared_ptr<domain> domain);

    virtual void
    updateWallScaleGradientField(const std::shared_ptr<domain> domain);

    virtual void updateTurbulentKineticEnergyGradientField(
        const std::shared_ptr<domain> domain);

    virtual void updateTurbulentEddyFrequencyGradientField(
        const std::shared_ptr<domain> domain);

    virtual void updateTurbulentDissipationRateGradientField(
        const std::shared_ptr<domain> domain);

    virtual void updateTransitionOnsetReynoldsNumberGradientField(
        const std::shared_ptr<domain> domain);

    virtual void updateTurbulentIntermittencyGradientField(
        const std::shared_ptr<domain> domain);

    virtual void
    updateDensityGradientField(const std::shared_ptr<domain> domain);

    virtual void
    updateDynamicViscosityGradientField(const std::shared_ptr<domain> domain);

    virtual void updateSpecificHeatCapacityGradientField(
        const std::shared_ptr<domain> domain);

    virtual void updateThermalConductivityGradientField(
        const std::shared_ptr<domain> domain);

    virtual void
    updateDisplacementGradientField(const std::shared_ptr<domain> domain);

    virtual void
    updateTotalDisplacementGradientField(const std::shared_ptr<domain> domain);

    virtual void
    updateMeshVelocityGradientField(const std::shared_ptr<domain> domain);

    // gradient for phasic fields (iPhase is global to the simulation: i.e.
    // material index)

    virtual void
    updateVolumeFractionGradientField(const std::shared_ptr<domain> domain,
                                      label iPhase);

    virtual void
    updateDensityGradientField(const std::shared_ptr<domain> domain,
                               label iPhase);

    virtual void
    updateDynamicViscosityGradientField(const std::shared_ptr<domain> domain,
                                        label iPhase);

    virtual void updateSpecificHeatCapacityGradientField(
        const std::shared_ptr<domain> domain,
        label iPhase);

    virtual void
    updateThermalConductivityGradientField(const std::shared_ptr<domain> domain,
                                           label iPhase);

    // beta

    virtual void
    updateVelocityBlendingFactorField(const std::shared_ptr<domain> domain);

    virtual void
    updatePressureBlendingFactorField(const std::shared_ptr<domain> domain);

    virtual void
    updateTemperatureBlendingFactorField(const std::shared_ptr<domain> domain);

    virtual void updateSpecificEnthalpyBlendingFactorField(
        const std::shared_ptr<domain> domain);

    virtual void updateSpecificTotalEnthalpyBlendingFactorField(
        const std::shared_ptr<domain> domain);

    virtual void updateTurbulentKineticEnergyBlendingFactorField(
        const std::shared_ptr<domain> domain);

    virtual void updateTurbulentEddyFrequencyBlendingFactorField(
        const std::shared_ptr<domain> domain);

    virtual void updateTransitionOnsetReynoldsNumberBlendingFactorField(
        const std::shared_ptr<domain> domain);

    virtual void updateTurbulentIntermittencyBlendingFactorField(
        const std::shared_ptr<domain> domain);

    virtual void updateTurbulentDissipationRateBlendingFactorField(
        const std::shared_ptr<domain> domain);

    virtual void
    updateDensityBlendingFactorField(const std::shared_ptr<domain> domain);

    virtual void updateDynamicViscosityBlendingFactorField(
        const std::shared_ptr<domain> domain);

    virtual void updateSpecificHeatCapacityBlendingFactorField(
        const std::shared_ptr<domain> domain);

    virtual void updateThermalConductivityBlendingFactorField(
        const std::shared_ptr<domain> domain);

    // beta for phasic fields (iPhase is global to the simulation: i.e. material
    // index)

    virtual void updateVolumeFractionBlendingFactorField(
        const std::shared_ptr<domain> domain,
        label iPhase);

    virtual void
    updateDensityBlendingFactorField(const std::shared_ptr<domain> domain,
                                     label iPhase);

    virtual void updateDynamicViscosityBlendingFactorField(
        const std::shared_ptr<domain> domain,
        label iPhase);

    virtual void updateSpecificHeatCapacityBlendingFactorField(
        const std::shared_ptr<domain> domain,
        label iPhase);

    virtual void updateThermalConductivityBlendingFactorField(
        const std::shared_ptr<domain> domain,
        label iPhase);

    // prev-iter

    virtual void
    updateVelocityPrevIterField(const std::shared_ptr<domain> domain);

    virtual void
    updatePressurePrevIterField(const std::shared_ptr<domain> domain);

    virtual void
    updateTemperaturePrevIterField(const std::shared_ptr<domain> domain);

    virtual void
    updateSpecificEnthalpyPrevIterField(const std::shared_ptr<domain> domain);

    virtual void updateSpecificTotalEnthalpyPrevIterField(
        const std::shared_ptr<domain> domain);

    virtual void
    updateWallScalePrevIterField(const std::shared_ptr<domain> domain);

    virtual void updateTurbulentKineticEnergyPrevIterField(
        const std::shared_ptr<domain> domain);

    virtual void updateTurbulentEddyFrequencyPrevIterField(
        const std::shared_ptr<domain> domain);

    virtual void updateTurbulentDissipationRatePrevIterField(
        const std::shared_ptr<domain> domain);

    virtual void updateTransitionOnsetReynoldsNumberPrevIterField(
        const std::shared_ptr<domain> domain);

    virtual void updateTurbulentIntermittencyPrevIterField(
        const std::shared_ptr<domain> domain);

    virtual void
    updateDensityPrevIterField(const std::shared_ptr<domain> domain);

    virtual void
    updateDynamicViscosityPrevIterField(const std::shared_ptr<domain> domain);

    virtual void updateSpecificHeatCapacityPrevIterField(
        const std::shared_ptr<domain> domain);

    virtual void updateThermalConductivityPrevIterField(
        const std::shared_ptr<domain> domain);

    virtual void
    updateDisplacementPrevIterField(const std::shared_ptr<domain> domain);

    virtual void
    updateTotalDisplacementPrevIterField(const std::shared_ptr<domain> domain);

    // prev-iter for phasic fields (iPhase is global to the simulation: i.e.
    // material index)

    virtual void
    updateVolumeFractionPrevIterField(const std::shared_ptr<domain> domain,
                                      label iPhase);

    virtual void
    updateDensityPrevIterField(const std::shared_ptr<domain> domain,
                               label iPhase);

    virtual void
    updateDynamicViscosityPrevIterField(const std::shared_ptr<domain> domain,
                                        label iPhase);

    virtual void updateSpecificHeatCapacityPrevIterField(
        const std::shared_ptr<domain> domain,
        label iPhase);

    virtual void
    updateThermalConductivityPrevIterField(const std::shared_ptr<domain> domain,
                                           label iPhase);

    // prev-time

    virtual void
    updateVelocityPrevTimeField(const std::shared_ptr<domain> domain);

    virtual void
    updatePressurePrevTimeField(const std::shared_ptr<domain> domain);

    virtual void
    updateTemperaturePrevTimeField(const std::shared_ptr<domain> domain);

    virtual void
    updateSpecificEnthalpyPrevTimeField(const std::shared_ptr<domain> domain);

    virtual void updateSpecificTotalEnthalpyPrevTimeField(
        const std::shared_ptr<domain> domain);

    virtual void
    updateWallScalePrevTimeField(const std::shared_ptr<domain> domain);

    virtual void updateTurbulentKineticEnergyPrevTimeField(
        const std::shared_ptr<domain> domain);

    virtual void updateTurbulentEddyFrequencyPrevTimeField(
        const std::shared_ptr<domain> domain);

    virtual void updateTurbulentDissipationRatePrevTimeField(
        const std::shared_ptr<domain> domain);

    virtual void updateTransitionOnsetReynoldsNumberPrevTimeField(
        const std::shared_ptr<domain> domain);

    virtual void updateTurbulentIntermittencyPrevTimeField(
        const std::shared_ptr<domain> domain);

    virtual void
    updateDensityPrevTimeField(const std::shared_ptr<domain> domain);

    virtual void
    updateDynamicViscosityPrevTimeField(const std::shared_ptr<domain> domain);

    virtual void updateSpecificHeatCapacityPrevTimeField(
        const std::shared_ptr<domain> domain);

    virtual void updateThermalConductivityPrevTimeField(
        const std::shared_ptr<domain> domain);

    virtual void
    updateCompressibilityPrevTimeField(const std::shared_ptr<domain> domain);

    virtual void
    updateDisplacementPrevTimeField(const std::shared_ptr<domain> domain);

    virtual void
    updateTotalDisplacementPrevTimeField(const std::shared_ptr<domain> domain);

    // prev-time for phasic fields (iPhase is global to the simulation: i.e.
    // material index)

    virtual void
    updateVolumeFractionPrevTimeField(const std::shared_ptr<domain> domain,
                                      label iPhase);

    virtual void
    updateDensityPrevTimeField(const std::shared_ptr<domain> domain,
                               label iPhase);

    virtual void
    updateDynamicViscosityPrevTimeField(const std::shared_ptr<domain> domain,
                                        label iPhase);

    virtual void updateSpecificHeatCapacityPrevTimeField(
        const std::shared_ptr<domain> domain,
        label iPhase);

    virtual void
    updateThermalConductivityPrevTimeField(const std::shared_ptr<domain> domain,
                                           label iPhase);

protected:
    // Protected auxiliary initializations

    virtual void
    initializeMassFlowRateInterior_(const std::shared_ptr<domain> domain)
    {
        errorMsg("Must not reach here");
    }

#ifdef HAS_INTERFACE
    virtual void initializeMassFlowRateInterfaceSideField_(
        const std::shared_ptr<domain> domain,
        const interfaceSideInfo* interfaceSideInfoPtr)
    {
        errorMsg("Must not reach here");
    }
#endif /* HAS_INTERFACE */

    virtual void
    initializeMassFlowRateBoundaryField_(const std::shared_ptr<domain> domain,
                                         const boundary* boundary)
    {
        errorMsg("Must not reach here");
    }

    virtual void
    updateMassFlowRateInterior_(const std::shared_ptr<domain> domain)
    {
        errorMsg("Must not reach here");
    }

#ifdef HAS_INTERFACE
    virtual void updateMassFlowRateInterfaceSideField_(
        const std::shared_ptr<domain> domain,
        const interfaceSideInfo* interfaceSideInfoPtr)
    {
        errorMsg("Must not reach here");
    }
#endif /* HAS_INTERFACE */

    virtual void
    updateMassFlowRateBoundaryField_(const std::shared_ptr<domain> domain,
                                     const boundary* boundary)
    {
        errorMsg("Must not reach here");
    }

    // Protected auxiliary initializations for phase fields

    virtual void
    initializeMassFlowRateInterior_(const std::shared_ptr<domain> domain,
                                    label iPhase)
    {
        errorMsg("Must not reach here");
    }

#ifdef HAS_INTERFACE
    virtual void initializeMassFlowRateInterfaceSideField_(
        const std::shared_ptr<domain> domain,
        const interfaceSideInfo* interfaceSideInfoPtr,
        label iPhase)
    {
        errorMsg("Must not reach here");
    }
#endif /* HAS_INTERFACE */

    virtual void
    initializeMassFlowRateBoundaryField_(const std::shared_ptr<domain> domain,
                                         const boundary* boundary,
                                         label iPhase)
    {
        errorMsg("Must not reach here");
    }

    virtual void
    updateMassFlowRateInterior_(const std::shared_ptr<domain> domain,
                                label iPhase)
    {
        errorMsg("Must not reach here");
    }

#ifdef HAS_INTERFACE
    virtual void updateMassFlowRateInterfaceSideField_(
        const std::shared_ptr<domain> domain,
        const interfaceSideInfo* interfaceSideInfoPtr,
        label iPhase)
    {
        errorMsg("Must not reach here");
    }
#endif /* HAS_INTERFACE */

    virtual void
    updateMassFlowRateBoundaryField_(const std::shared_ptr<domain> domain,
                                     const boundary* boundary,
                                     label iPhase)
    {
        errorMsg("Must not reach here");
    }

    // protected API for setting boundary conditions (mainly used in
    // initialize_* methods)
    void setupBoundaryConditions_(
        const std::shared_ptr<domain> domain,
        std::function<void(const ::accel::domain* domain,
                           const label iBoundary,
                           const boundaryPhysicalType bc_type,
                           const YAML::Node bc_values)> set_bc);

    // devoted to fields that are fluid-specific
    void setupBoundaryConditions_(
        const std::shared_ptr<domain> domain,
        const label iPhase,
        std::function<void(const ::accel::domain* domain,
                           const label iBoundary,
                           const label iPhase,
                           const boundaryPhysicalType bc_type,
                           const YAML::Node bc_values)> set_bc);

protected:
    realm* realmPtr_;

    // protected field API for field broker, to be used selectively in models
    // that require these fields to solve the equations (these must not have an
    // `_` suffix in order to forward them to public scope with `using`
    // statements in derived classes)
    pressure& pRef();

    const pressure& pRef() const;

    simpleScalarField& p0Ref();

    const simpleScalarField& p0Ref() const;

    simpleScalarField& MaRef();

    const simpleScalarField& MaRef() const;

    simpleScalarField& CoRef();

    const simpleScalarField& CoRef() const;

    temperature& TRef();

    const temperature& TRef() const;

    simpleScalarField& T0Ref();

    const simpleScalarField& T0Ref() const;

    specificEnthalpy& hRef();

    const specificEnthalpy& hRef() const;

    specificTotalEnthalpy& h0Ref();

    const specificTotalEnthalpy& h0Ref() const;

    nodeScalarField& yScaleRef();

    const nodeScalarField& yScaleRef() const;

    turbulentKineticEnergy& kRef();

    const turbulentKineticEnergy& kRef() const;

    turbulentEddyFrequency& omegaRef();

    const turbulentEddyFrequency& omegaRef() const;

    turbulentDissipationRate& epsilonRef();

    const turbulentDissipationRate& epsilonRef() const;

    transitionOnsetReynoldsNumber& ReThetaRef();

    const transitionOnsetReynoldsNumber& ReThetaRef() const;

    turbulentIntermittency& gammaRef();

    const turbulentIntermittency& gammaRef() const;

    heatFlowRate& qDotRef();

    const heatFlowRate& qDotRef() const;

    youngModulus& ERef();

    const youngModulus& ERef() const;

    poissonRatio& nuRef();

    const poissonRatio& nuRef() const;

    dynamicViscosity& muRef();

    const dynamicViscosity& muRef() const;

    specificHeatCapacity& cpRef();

    const specificHeatCapacity& cpRef() const;

    thermalConductivity& lambdaRef();

    const thermalConductivity& lambdaRef() const;

    thermalExpansivity& betaRef();

    const thermalExpansivity& betaRef() const;

    compressibility& psiRef();

    const compressibility& psiRef() const;

    turbulentViscosity& mutRef();

    const turbulentViscosity& mutRef() const;

    turbulentThermalConductivity& lambdatRef();

    const turbulentThermalConductivity& lambdatRef() const;

    nodeScalarField& muEffRef();

    const nodeScalarField& muEffRef() const;

    nodeScalarField& lambdaEffRef();

    const nodeScalarField& lambdaEffRef() const;

    simpleScalarField& uTauRef();

    const simpleScalarField& uTauRef() const;

    simpleScalarField& yMinRef();

    const simpleScalarField& yMinRef() const;

    simpleScalarField& F1Ref();

    const simpleScalarField& F1Ref() const;

    simpleScalarField& PkRef();

    const simpleScalarField& PkRef() const;

    simpleScalarField& yPlusRef();

    const simpleScalarField& yPlusRef() const;

    simpleScalarField& uPlusRef();

    const simpleScalarField& uPlusRef() const;

    simpleScalarField& TPlusRef();

    const simpleScalarField& TPlusRef() const;

    simpleScalarField& yStarRef();

    const simpleScalarField& yStarRef() const;

    simpleScalarField& uStarRef();

    const simpleScalarField& uStarRef() const;

    simpleScalarField& duPlusdyPlusRef();

    const simpleScalarField& duPlusdyPlusRef() const;

    simpleScalarField& uWallCoeffsRef();

    const simpleScalarField& uWallCoeffsRef() const;

    simpleScalarField& TWallCoeffsRef();

    const simpleScalarField& TWallCoeffsRef() const;

    simpleVectorField& wallShearStressRef();

    const simpleVectorField& wallShearStressRef() const;

    nodeVectorField& DtRef();

    const nodeVectorField& DtRef() const;

    displacement& DRef();

    const displacement& DRef() const;

    simpleTensorField& stressRef();

    const simpleTensorField& stressRef() const;

    simpleTensorField& strainRef();

    const simpleTensorField& strainRef() const;

    // Protected access for phasic fields (iPhase is global to the simulation)

    volumeFraction& alphaRef(label iPhase);

    const volumeFraction& alphaRef(label iPhase) const;

    simpleScalarField& alphaSmoothRef(label iPhase);

    const simpleScalarField& alphaSmoothRef(label iPhase) const;

    simpleScalarField& rhsSmoothRef(label iPhase);

    const simpleScalarField& rhsSmoothRef(label iPhase) const;

    dynamicViscosity& muRef(label iPhase);

    const dynamicViscosity& muRef(label iPhase) const;

    specificHeatCapacity& cpRef(label iPhase);

    const specificHeatCapacity& cpRef(label iPhase) const;

    thermalConductivity& lambdaRef(label iPhase);

    const thermalConductivity& lambdaRef(label iPhase) const;

    thermalExpansivity& betaRef(label iPhase);

    const thermalExpansivity& betaRef(label iPhase) const;

    compressibility& psiRef(label iPhase);

    const compressibility& psiRef(label iPhase) const;

    simpleVectorField& nHatRef(label iPhase);

    const simpleVectorField& nHatRef(label iPhase) const;

protected:
    // helper methods for auxiliary fields
    template <typename T>
    stk::mesh::Field<T>* newSTKField_(const std::string& name,
                                      const unsigned ncomponents,
                                      stk::topology::rank_t topology) const
    {
        stk::mesh::Field<T>* f =
            &const_cast<stk::mesh::MetaData&>(meshRef().metaDataRef())
                 .declare_field<T>(topology, name);
        for (auto zonePtr : realmPtr_->meshRef().zoneVector())
        {
            stk::mesh::put_field_on_mesh(
                *f,
                stk::mesh::selectUnion(zonePtr->interiorParts()),
                ncomponents,
                nullptr);
        }
        return f;
    }

    template <typename T>
    stk::mesh::Field<T>* newSTKWallField_(const std::string& name,
                                          const unsigned ncomponents) const
    {
        stk::mesh::Field<T>* f = nullptr;

        for (label iZone = 0; iZone < realmPtr_->meshRef().zoneVector().size();
             iZone++)
        {
            for (label iBoundary = 0;
                 iBoundary < realmPtr_->meshRef().zonePtr(iZone)->nBoundaries();
                 iBoundary++)
            {
                const auto& boundaryRef =
                    realmPtr_->meshRef().zonePtr(iZone)->boundaryRef(iBoundary);
                boundaryPhysicalType type = boundaryRef.type();
                switch (type)
                {
                    case boundaryPhysicalType::wall:
                        {
                            f = &const_cast<stk::mesh::MetaData&>(
                                     meshRef().metaDataRef())
                                     .declare_field<T>(
                                         meshRef().metaDataRef().side_rank(),
                                         name);
                            for (auto* part : boundaryRef.parts())
                            {
                                for (const stk::mesh::Part* subPart :
                                     part->subsets())
                                {
                                    MasterElement* meFC = MasterElementRepo::
                                        get_surface_master_element(
                                            subPart->topology());
                                    const label numScsIp = meFC->numIntPoints_;
                                    stk::mesh::put_field_on_mesh(*f,
                                                                 *subPart,
                                                                 ncomponents *
                                                                     numScsIp,
                                                                 nullptr);
                                }
                            }
                        }
                        break;

                    default:
                        break;
                }
            }
        }

#ifdef HAS_INTERFACE
        for (label iInterface = 0; iInterface < meshRef().nInterfaces();
             iInterface++)
        {
            const auto& interf = meshRef().interfaceRef(iInterface);

            // In case of a very rare situation where no walls exist in the
            // simulation, f needs to be declared since not yet declared
            if (f == nullptr)
            {
                f = &const_cast<stk::mesh::MetaData&>(meshRef().metaDataRef())
                         .declare_field<T>(meshRef().metaDataRef().side_rank(),
                                           name);
            }

            if (interf.isFluidSolidType())
            {
                for (auto* part : interf.masterInfoRef().currentPartVec_)
                {
                    for (const stk::mesh::Part* subPart : part->subsets())
                    {
                        MasterElement* meFC =
                            MasterElementRepo::get_surface_master_element(
                                subPart->topology());
                        const label numScsIp = meFC->numIntPoints_;
                        stk::mesh::put_field_on_mesh(
                            *f, *subPart, ncomponents * numScsIp, nullptr);
                    }
                }

                for (auto* part : interf.slaveInfoRef().currentPartVec_)
                {
                    for (const stk::mesh::Part* subPart : part->subsets())
                    {
                        MasterElement* meFC =
                            MasterElementRepo::get_surface_master_element(
                                subPart->topology());
                        const label numScsIp = meFC->numIntPoints_;
                        stk::mesh::put_field_on_mesh(
                            *f, *subPart, ncomponents * numScsIp, nullptr);
                    }
                }
            }
        }
#endif /* HAS_INTERFACE */

        // For very rare situations where no walls exist, f still needs to be
        // instantiated for consistency
        if (f == nullptr)
        {
            f = &const_cast<stk::mesh::MetaData&>(meshRef().metaDataRef())
                     .declare_field<T>(meshRef().metaDataRef().side_rank(),
                                       name);
        }
        return f;
    }
};

#define RETURN_REF_ALLOC_AUX(ptr, type, realm, name, topology)                 \
    do                                                                         \
    {                                                                          \
        if ((ptr) == nullptr)                                                  \
        {                                                                      \
            auto* stk_field =                                                  \
                this->template newSTKField_<typename type::DataType>(          \
                    name, type::NComponents, topology);                        \
            assert(stk_field != nullptr);                                      \
            (ptr) = std::make_unique<type>(&(realm)->meshRef(), stk_field);    \
        }                                                                      \
        return *(ptr);                                                         \
    } while (0)

#define RETURN_REF_ALLOC_AUX_WALL(ptr, type, realm, name)                      \
    do                                                                         \
    {                                                                          \
        if ((ptr) == nullptr)                                                  \
        {                                                                      \
            auto* stk_field =                                                  \
                this->template newSTKWallField_<typename type::DataType>(      \
                    name, type::NComponents);                                  \
            assert(stk_field != nullptr);                                      \
            (ptr) = std::make_unique<type>(&(realm)->meshRef(), stk_field);    \
        }                                                                      \
        return *(ptr);                                                         \
    } while (0)

#define RETURN_REF_ALLOC(ptr, type, ...)                                       \
    do                                                                         \
    {                                                                          \
        if ((ptr) == nullptr)                                                  \
        {                                                                      \
            const label n_states =                                             \
                meshRef().controlsRef().getNumberOfStates();                   \
            const bool high_res = meshRef().controlsRef().isHighResolution();  \
            (ptr) = std::make_unique<type>(__VA_ARGS__);                       \
        }                                                                      \
        return *(ptr);                                                         \
    } while (0)

#define RETURN_REF_ALLOC_TURB(ptr, type, ...)                                  \
    do                                                                         \
    {                                                                          \
        if ((ptr) == nullptr)                                                  \
        {                                                                      \
            const label n_states =                                             \
                meshRef().controlsRef().getNumberOfStates();                   \
            const bool high_res =                                              \
                meshRef().controlsRef().isHighResolutionTurbulenceNumerics();  \
            (ptr) = std::make_unique<type>(__VA_ARGS__);                       \
        }                                                                      \
        return *(ptr);                                                         \
    } while (0)

#define SET_DENSITY_PROPERTY(rhoField, materialBlock, ...)                     \
    do                                                                         \
    {                                                                          \
        if ((rhoField).isZoneUnset((domain)->index()))                         \
        {                                                                      \
            (rhoField).setZone((domain)->index());                             \
            if (domain->materialRef(__VA_ARGS__)                               \
                    .thermodynamicProperties_.equationOfState_.option_ ==      \
                equationOfStateOption::value)                                  \
            {                                                                  \
                initialCondition::                                             \
                    setupFieldInitializationOverDomainFromConfig(              \
                        (rhoField),                                            \
                        (domain)->index(),                                     \
                        (materialBlock)["thermodynamic_properties"]            \
                                       ["equation_of_state"]);                 \
            }                                                                  \
        }                                                                      \
    } while (0)

#define SET_DYNAMIC_VISCOSITY_PROPERTY(muField, materialBlock, ...)            \
    do                                                                         \
    {                                                                          \
        if ((muField).isZoneUnset((domain)->index()))                          \
        {                                                                      \
            (muField).setZone((domain)->index());                              \
            if (domain->materialRef(__VA_ARGS__)                               \
                    .transportProperties_.dynamicViscosity_.option_ ==         \
                dynamicViscosityOption::value)                                 \
            {                                                                  \
                initialCondition::                                             \
                    setupFieldInitializationOverDomainFromConfig(              \
                        (muField),                                             \
                        (domain)->index(),                                     \
                        (materialBlock)["transport_properties"]                \
                                       ["dynamic_viscosity"]);                 \
            }                                                                  \
        }                                                                      \
    } while (0)

#define SET_SPECIFIC_HEAT_CAPACITY_PROPERTY(cpField, materialBlock, ...)       \
    do                                                                         \
    {                                                                          \
        if ((cpField).isZoneUnset((domain)->index()))                          \
        {                                                                      \
            (cpField).setZone((domain)->index());                              \
            if (domain->materialRef(__VA_ARGS__)                               \
                    .thermodynamicProperties_.specificHeatCapacity_.option_ == \
                specificHeatCapacityOption::value)                             \
            {                                                                  \
                initialCondition::                                             \
                    setupFieldInitializationOverDomainFromConfig(              \
                        (cpField),                                             \
                        (domain)->index(),                                     \
                        (materialBlock)["thermodynamic_properties"]            \
                                       ["specific_heat_capacity"]);            \
            }                                                                  \
        }                                                                      \
    } while (0)

#define SET_THERMAL_CONDUCTIVITY_PROPERTY(lambdaField, materialBlock, ...)     \
    do                                                                         \
    {                                                                          \
        if ((lambdaField).isZoneUnset((domain)->index()))                      \
        {                                                                      \
            (lambdaField).setZone((domain)->index());                          \
            if (domain->materialRef(__VA_ARGS__)                               \
                    .transportProperties_.thermalConductivity_.option_ ==      \
                thermalConductivityOption::value)                              \
            {                                                                  \
                initialCondition::                                             \
                    setupFieldInitializationOverDomainFromConfig(              \
                        (lambdaField),                                         \
                        (domain)->index(),                                     \
                        (materialBlock)["transport_properties"]                \
                                       ["thermal_conductivity"]);              \
            }                                                                  \
        }                                                                      \
    } while (0)

#define SET_THERMAL_EXPANSIVITY_PROPERTY(betaField, materialBlock, ...)        \
    do                                                                         \
    {                                                                          \
        if ((betaField).isZoneUnset((domain)->index()))                        \
        {                                                                      \
            (betaField).setZone((domain)->index());                            \
            if (domain->materialRef(__VA_ARGS__)                               \
                    .buoyancyProperties_.thermalExpansivity_.option_ ==        \
                thermalExpansivityOption::value)                               \
            {                                                                  \
                initialCondition::                                             \
                    setupFieldInitializationOverDomainFromConfig(              \
                        (betaField),                                           \
                        (domain)->index(),                                     \
                        (materialBlock)["buoyancy_properties"]                 \
                                       ["thermal_expansivity"]);               \
            }                                                                  \
        }                                                                      \
    } while (0)

#define SET_ELECTRICAL_CONDUCTIVITY_PROPERTY(sigmaField, materialBlock, ...)   \
    do                                                                         \
    {                                                                          \
        if ((sigmaField).isZoneUnset((domain)->index()))                       \
        {                                                                      \
            (sigmaField).setZone((domain)->index());                           \
            if (domain->materialRef(__VA_ARGS__)                               \
                    .electromagneticProperties_.electricalConductivity_        \
                    .option_ == electricalConductivityOption::value)           \
            {                                                                  \
                initialCondition::                                             \
                    setupFieldInitializationOverDomainFromConfig(              \
                        (sigmaField),                                          \
                        (domain)->index(),                                     \
                        (materialBlock)["electromagnetic_properties"]          \
                                       ["electrical_conductivity"]);           \
            }                                                                  \
        }                                                                      \
    } while (0)

#define SET_YOUNG_MODULUS_PROPERTY(EField, materialBlock, ...)                 \
    do                                                                         \
    {                                                                          \
        if ((EField).isZoneUnset((domain)->index()))                           \
        {                                                                      \
            (EField).setZone((domain)->index());                               \
            if (domain->materialRef(__VA_ARGS__)                               \
                    .mechanicalProperties_.youngModulus_.option_ ==            \
                youngModulusOption::value)                                     \
            {                                                                  \
                initialCondition::                                             \
                    setupFieldInitializationOverDomainFromConfig(              \
                        (EField),                                              \
                        (domain)->index(),                                     \
                        (materialBlock)["mechanical_properties"]               \
                                       ["young_modulus"]);                     \
            }                                                                  \
        }                                                                      \
    } while (0)

#define SET_POISSON_RATIO_PROPERTY(nuField, materialBlock, ...)                \
    do                                                                         \
    {                                                                          \
        if ((nuField).isZoneUnset((domain)->index()))                          \
        {                                                                      \
            (nuField).setZone((domain)->index());                              \
            if (domain->materialRef(__VA_ARGS__)                               \
                    .mechanicalProperties_.poissonRatio_.option_ ==            \
                poissonRatioOption::value)                                     \
            {                                                                  \
                initialCondition::                                             \
                    setupFieldInitializationOverDomainFromConfig(              \
                        (nuField),                                             \
                        (domain)->index(),                                     \
                        (materialBlock)["mechanical_properties"]               \
                                       ["poisson_ratio"]);                     \
            }                                                                  \
        }                                                                      \
    } while (0)

#define RETURN_REF_ALLOC_MATERIAL(ptrVec, idx, type, ...)                      \
    do                                                                         \
    {                                                                          \
        if (ptrVec.empty())                                                    \
        {                                                                      \
            ptrVec.resize(realmPtr_->simulationRef().nRegisteredMaterials());  \
        }                                                                      \
        assert(idx < ptrVec.size());                                           \
        if ((ptrVec[idx]) == nullptr)                                          \
        {                                                                      \
            const label n_states =                                             \
                meshRef().controlsRef().getNumberOfStates();                   \
            const bool high_res = meshRef().controlsRef().isHighResolution();  \
            (ptrVec[idx]) = std::make_unique<type>(__VA_ARGS__);               \
        }                                                                      \
        return *(ptrVec[idx]);                                                 \
    } while (0)

#define RETURN_REF_ALLOC_AUX_MATERIAL(                                         \
    ptrVec, idx, type, realm, name, topology)                                  \
    do                                                                         \
    {                                                                          \
        if (ptrVec.empty())                                                    \
        {                                                                      \
            ptrVec.resize(realmPtr_->simulationRef().nRegisteredMaterials());  \
        }                                                                      \
        assert(idx < ptrVec.size());                                           \
        if ((ptrVec[idx]) == nullptr)                                          \
        {                                                                      \
            auto* stk_field =                                                  \
                this->template newSTKField_<typename type::DataType>(          \
                    name, type::NComponents, topology);                        \
            assert(stk_field != nullptr);                                      \
            (ptrVec[idx]) =                                                    \
                std::make_unique<type>(&(realm)->meshRef(), stk_field);        \
        }                                                                      \
        return *(ptrVec[idx]);                                                 \
    } while (0)

} /* namespace accel */

#endif // FIELDBROKER_H
