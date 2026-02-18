// File : freeSurfaceFlowModel.h
// Created : Sun Jan 26 2025 22:06:38 (+0100)
// Author : Mhamad Mahdi Alloush
// Description: VOF-based free surface flow model with FCT and surface tension
// Copyright (c) 2025 CCFNUM, Lucerne University of Applied Sciences and Arts.
// SPDX-License-Identifier: BSD-3-Clause

#ifndef FREESURFACEFLOWMODEL_H
#define FREESURFACEFLOWMODEL_H

// code
#include "multiphaseModel.h"

namespace accel
{

class freeSurfaceFlowModel : public multiphaseModel
{
private:
    // parameters for smoothing alpha

    scalar dxMin_ = 1.0e16; // initial value

    scalar Fo_ = 0.25;

    // cMULES iterative lambda finding parameters

    // tol = 0 means all nIter iterations are performed (no early exit)
    scalar lambdaTolerance_ = 0.0;

    label maxLambdaIterations_ = 3;

    void computeFL_(const std::shared_ptr<domain> domain, label iPhase);

    void computeFH_(const std::shared_ptr<domain> domain, label iPhase);

    void computeA_(const std::shared_ptr<domain> domain, label iPhase);

    void computeLambda_(const std::shared_ptr<domain> domain, label iPhase);

    void computeLambdaNode_(const std::shared_ptr<domain> domain, label iPhase);

    scalar computeLambdaNodeWithBounds_(const std::shared_ptr<domain> domain,
                                        label iPhase,
                                        label iter);

    void computeLambdaIP_(const std::shared_ptr<domain> domain, label iPhase);

    void computeP_(const std::shared_ptr<domain> domain, label iPhase);

    void computeQ_(const std::shared_ptr<domain> domain, label iPhase);

    void updateA_(const std::shared_ptr<domain> domain, label iPhase);

    void computeSumA_(const std::shared_ptr<domain> domain, label iPhase);

    void updateAlpha_(const std::shared_ptr<domain> domain,
                      label iPhase,
                      label iter);

    void
    updateFL_(const std::shared_ptr<domain> domain, label iPhase, label iter);

public:
    // Constructors

    freeSurfaceFlowModel(realm* realm);

    // Enable public use

    using fieldBroker::nHatRef;

    // Body forces override (adds CSF surface tension after base class forces)
    void computeBodyForces(const std::shared_ptr<domain> domain) override;

    // body force redistribution with harmonic density-weighted interpolation
    void redistributeBodyForces(const std::shared_ptr<domain> domain) override;

    // methods (iPhase is global to the simulation)

    void transformMassFlowRateToRelative(const std::shared_ptr<domain> domain,
                                         label iPhase);

    void transformMassFlowRateToAbsolute(const std::shared_ptr<domain> domain,
                                         label iPhase);

    void updateFlowReversalFlag(const std::shared_ptr<domain> domain,
                                label iPhase);

    void updateMassDivergenceField(const std::shared_ptr<domain> domain,
                                   label iPhase);

    void updateSideMassFlowRateFraction(const std::shared_ptr<domain> domain,
                                        label iPhase);

    void updateInterfaceNormal(const std::shared_ptr<domain> domain,
                               label iPhase);

    void applyVOFSmoothing(const std::shared_ptr<domain> domain, label iPhase);

    void applyVolumeConservation(const std::shared_ptr<domain> domain);

    // Flux Corrected Transport FCT

    void setupFCTFields(const std::shared_ptr<domain> domain, label iPhase);

    void correctFCT(const std::shared_ptr<domain> domain, label iPhase);

    // setup (iPhase is global to the simulation)

    virtual void setupVolumeFraction(const std::shared_ptr<domain> domain,
                                     label iPhase) override;

    virtual void setupDensity(const std::shared_ptr<domain> domain) override;

    virtual void
    setupDynamicViscosity(const std::shared_ptr<domain> domain) override;

    virtual void
    setupSpecificHeatCapacity(const std::shared_ptr<domain> domain) override;

    virtual void
    setupThermalConductivity(const std::shared_ptr<domain> domain) override;

    // initialize (iPhase is global to the simulation)

    virtual void
    initializeDensity(const std::shared_ptr<domain> domain) override;

    virtual void initializeDensity(const std::shared_ptr<domain> domain,
                                   label iPhase) override;

    virtual void
    initializeDynamicViscosity(const std::shared_ptr<domain> domain) override;

    virtual void
    initializeDynamicViscosity(const std::shared_ptr<domain> domain,
                               label iPhase) override;

    virtual void
    initializeMassFlowRate(const std::shared_ptr<domain> domain) override;

    // update (iPhase is global to the simulation)

    virtual void updateDensity(const std::shared_ptr<domain> domain) override;

    virtual void updateDensity(const std::shared_ptr<domain> domain,
                               label iPhase) override;

    virtual void
    updateDynamicViscosity(const std::shared_ptr<domain> domain) override;

    virtual void updateDynamicViscosity(const std::shared_ptr<domain> domain,
                                        label iPhase) override;

    // Access

    scalar gamma(const domain* domain) const
    {
        if (domain->multiphase_.freeSurfaceModel_.interfaceCompressionLevel_ ==
            0)
        {
            return 0.0;
        }
        else if (domain->multiphase_.freeSurfaceModel_
                     .interfaceCompressionLevel_ == 1)
        {
            return 0.5;
        }
        else if (domain->multiphase_.freeSurfaceModel_
                     .interfaceCompressionLevel_ == 2)
        {
            return 1.0;
        }

        return 0.5;
    }

protected:
    // Surface tension: per-pair curvature fields
    void computeCurvature_(const std::shared_ptr<domain> domain,
                           label iPhase,
                           STKScalarField* kappaFieldPtr);
    std::map<std::string, STKScalarField*> kappaSTKFieldPtrs_;

    void initializeMassFlowRateInterior_(const std::shared_ptr<domain> domain,
                                         label iPhase) override;

    void
    initializeMassFlowRateBoundaryField_(const std::shared_ptr<domain> domain,
                                         const boundary* boundary,
                                         label iPhase) override;

    void
    updateMassFlowRateInterior_(const std::shared_ptr<domain> domain) override;

    void updateMassFlowRateBoundaryField_(const std::shared_ptr<domain> domain,
                                          const boundary* boundary) override;

    void updateMassFlowRateInterior_(const std::shared_ptr<domain> domain,
                                     label iPhase) override;

    void updateMassFlowRateBoundaryField_(const std::shared_ptr<domain> domain,
                                          const boundary* boundary,
                                          label iPhase) override;

    // Protected methods (iPhase is global to the simulation)

    void computeSmoothRHS_(const std::shared_ptr<domain> domain, label iPhase);

    void assembleSmoothingTerm_(const std::shared_ptr<domain> domain,
                                label iPhase);

    // FCT Fields
    STKScalarField* FLSTKFieldPtr_ = nullptr;
    STKScalarField* sideFLSTKFieldPtr_ = nullptr;
    STKScalarField* FHSTKFieldPtr_ = nullptr;
    STKScalarField* sideFHSTKFieldPtr_ = nullptr;
    STKScalarField* ASTKFieldPtr_ = nullptr;
    STKScalarField* sideASTKFieldPtr_ = nullptr;
    STKScalarField* lambdaSTKFieldPtr_ = nullptr;
    STKScalarField* sideLambdaSTKFieldPtr_ = nullptr;
    STKScalarField* QplusSTKFieldPtr_ = nullptr;
    STKScalarField* QminusSTKFieldPtr_ = nullptr;
    STKScalarField* PplusSTKFieldPtr_ = nullptr;
    STKScalarField* PminusSTKFieldPtr_ = nullptr;
    STKScalarField* sumAPlusSTKFieldPtr_ = nullptr;
    STKScalarField* sumAMinusSTKFieldPtr_ = nullptr;
    STKScalarField* cMULESLimiterPlusSTKFieldPtr_ = nullptr;
    STKScalarField* cMULESLimiterMinusSTKFieldPtr_ = nullptr;
    STKScalarField* cMULESLimiterPlusPrevSTKFieldPtr_ = nullptr;
    STKScalarField* cMULESLimiterMinusPrevSTKFieldPtr_ = nullptr;
};

} /* namespace accel */

#endif // FREESURFACEFLOWMODEL_H
