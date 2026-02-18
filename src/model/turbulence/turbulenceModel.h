// File : turbulenceModel.h
// Created : Fri Mar 15 2024 15:06:38 (+0100)
// Author : Mhamad Mahdi Alloush
// Description: Base turbulence model with wall function parameters
// Copyright (c) 2024 CCFNUM, Lucerne University of Applied Sciences and Arts.
// SPDX-License-Identifier: BSD-3-Clause

#ifndef TURBULENCEMODEL_H
#define TURBULENCEMODEL_H

// code
#include "model.h"

namespace accel
{

class turbulenceModel : public model
{
private:
    scalar Cmu_ = 0.09;

    // Wall function-related parameters

    scalar elog_ = 9.8;

    scalar yplusCrit_ = 11.63;

    scalar kappa_ = 0.41;

    scalar B_ = 5.2;

public:
    // Constructors

    turbulenceModel(realm* realm);

    // Initialize

    void initializeTurbulentKineticEnergy(
        const std::shared_ptr<domain> domain) override;

    void initializeTurbulentEddyFrequency(
        const std::shared_ptr<domain> domain) override;

    void initializeTurbulentDissipationRate(
        const std::shared_ptr<domain> domain) override;

    // Update

    void
    updateTurbulentKineticEnergy(const std::shared_ptr<domain> domain) override;

    void
    updateTurbulentEddyFrequency(const std::shared_ptr<domain> domain) override;

    void updateTurbulentDissipationRate(
        const std::shared_ptr<domain> domain) override;

    // Enable public use

    using fieldBroker::cpRef;
    using fieldBroker::duPlusdyPlusRef;
    using fieldBroker::epsilonRef;
    using fieldBroker::F1Ref;
    using fieldBroker::gammaRef;
    using fieldBroker::kRef;
    using fieldBroker::lambdaEffRef;
    using fieldBroker::lambdaRef;
    using fieldBroker::lambdatRef;
    using fieldBroker::muEffRef;
    using fieldBroker::muRef;
    using fieldBroker::mutRef;
    using fieldBroker::omegaRef;
    using fieldBroker::PkRef;
    using fieldBroker::ReThetaRef;
    using fieldBroker::TPlusRef;
    using fieldBroker::TRef;
    using fieldBroker::TWallCoeffsRef;
    using fieldBroker::uPlusRef;
    using fieldBroker::uStarRef;
    using fieldBroker::uTauRef;
    using fieldBroker::uWallCoeffsRef;
    using fieldBroker::wallShearStressRef;
    using fieldBroker::yMinRef;
    using fieldBroker::yPlusRef;
    using fieldBroker::yScaleRef;
    using fieldBroker::yStarRef;

    // Access

    scalar Cmu() const
    {
        return Cmu_;
    }

    scalar elog() const
    {
        return elog_;
    }

    scalar yplusCrit() const
    {
        return yplusCrit_;
    }

    scalar kappa() const
    {
        return kappa_;
    }

    scalar B() const
    {
        return B_;
    }

protected:
    // Wall-functions

    void updateYPlus(const std::shared_ptr<domain> domain);

    void updateUPlus(const std::shared_ptr<domain> domain);

    void updateYStar(const std::shared_ptr<domain> domain);

    void updateUStar(const std::shared_ptr<domain> domain);

    void updateUTau(const std::shared_ptr<domain> domain);

    void updateDuPlusDyPlus(const std::shared_ptr<domain> domain);

    void updateUWallCoeffs(const std::shared_ptr<domain> domain);

    void updateWallShearStress(const std::shared_ptr<domain> domain);

    void updateTPlus(const std::shared_ptr<domain> domain);

    void updateTWallCoeffs(const std::shared_ptr<domain> domain);

    // Protected operations

    virtual void
    updateTurbulentDynamicViscosity(const std::shared_ptr<domain> domain) = 0;

    virtual void
    updateEffectiveDynamicViscosity(const std::shared_ptr<domain> domain);

    virtual void
    updateTurbulentThermalConductivity(const std::shared_ptr<domain> domain);

    virtual void
    updateEffectiveThermalConductivity(const std::shared_ptr<domain> domain);

    virtual void clipMinDistToWall(const std::shared_ptr<domain> domain);

    // Other protected methods

    const stk::mesh::PartVector
    collectNoSlipWallParts_(const std::shared_ptr<domain> domain);

private:
    // Private Wall-Function methods
    void updateYStarScalable_(const std::shared_ptr<domain> domain);

    void updateUStarScalable_(const std::shared_ptr<domain> domain);

    void updateYPlusScalable_(const std::shared_ptr<domain> domain);

    void updateUTauScalable_(const std::shared_ptr<domain> domain);

    void updateYStarAutomatic_(const std::shared_ptr<domain> domain);

    void updateUStarAutomatic_(const std::shared_ptr<domain> domain);

    void updateYPlusAutomatic_(const std::shared_ptr<domain> domain);

    void updateUTauAutomatic_(const std::shared_ptr<domain> domain);

    void updateTPlusScalable_(const std::shared_ptr<domain> domain);

    void updateTPlusAutomatic_(const std::shared_ptr<domain> domain);

    // Side Updates

    void updateTurbulentKineticEnergySideFields_(
        const std::shared_ptr<domain> domain);

    void updateTurbulentEddyFrequencySideFields_(
        const std::shared_ptr<domain> domain);

    void updateTurbulentDissipationRateSideFields_(
        const std::shared_ptr<domain> domain);

    void
    updateTurbulentKineticEnergyBoundarySideFieldInletIntensityAndLengthScale_(
        const std::shared_ptr<domain> domain,
        const boundary* boundary);

    void
    updateTurbulentEddyFrequencyBoundarySideFieldInletIntensityAndLengthScale_(
        const std::shared_ptr<domain> domain,
        const boundary* boundary);

    void
    updateTurbulentDissipationRateBoundarySideFieldInletIntensityAndLengthScale_(
        const std::shared_ptr<domain> domain,
        const boundary* boundary);

    void
    updateTurbulentKineticEnergyBoundarySideFieldInletIntensityAndEddyViscosityRatio_(
        const std::shared_ptr<domain> domain,
        const boundary* boundary);

    void
    updateTurbulentEddyFrequencyBoundarySideFieldInletIntensityAndEddyViscosityRatio_(
        const std::shared_ptr<domain> domain,
        const boundary* boundary);

    void
    updateTurbulentDissipationRateBoundarySideFieldInletIntensityAndEddyViscosityRatio_(
        const std::shared_ptr<domain> domain,
        const boundary* boundary);
};

} /* namespace accel */

#endif // TURBULENCEMODEL_H
