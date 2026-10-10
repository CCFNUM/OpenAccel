// File       : shearStressTransportModel.h
// Created    : Fri Mar 15 2024 15:06:38 (+0100)
// Author     : Mhamad Mahdi Alloush
// Description: Menter SST k-omega turbulence model coefficients and blending
// Copyright 2024 CCFNUM HSLU T&A. All Rights Reserved.

#ifndef SHEARSTRESSTRANSPORTMODEL_H
#define SHEARSTRESSTRANSPORTMODEL_H

#include "RANSModel.h"
#include "nodeTensorField.h"

namespace accel
{

class shearStressTransportModel : public RANSModel
{
private:
    // Model parameters: default SST

    scalar sigmaKOne_ = 0.85;

    scalar sigmaKTwo_ = 1.0;

    scalar sigmaWOne_ = 0.5;

    scalar sigmaWTwo_ = 0.856;

    scalar betaOne_ = 0.075;

    scalar betaTwo_ = 0.0828;

    scalar gammaOne_ = 0.5532;

    scalar gammaTwo_ = 0.4403;

    scalar aOne_ = 0.31;

    // Curvature correction (Smirnov-Menter 2009) constants, SST values
    scalar cR1_ = 1.0;

    scalar cR2_ = 2.0;

    scalar cR3_ = 1.0;

    // strain-rate tensor S_ij with gradient dS_ij/dx_k, for DS_ij/Dt; only
    // allocated when curvature correction is enabled
    std::unique_ptr<nodeTensorField> SijPtr_ = nullptr;

public:
    shearStressTransportModel(realm* realm);

    // Wall-functions

    void
    updateTurbulentEddyFrequencyAtWalls(const std::shared_ptr<domain> domain);

    // Access

    scalar sigmaKOne() const
    {
        return sigmaKOne_;
    }

    scalar sigmaKTwo() const
    {
        return sigmaKTwo_;
    }

    scalar sigmaWOne() const
    {
        return sigmaWOne_;
    }

    scalar sigmaWTwo() const
    {
        return sigmaWTwo_;
    }

    scalar betaOne() const
    {
        return betaOne_;
    }

    scalar betaTwo() const
    {
        return betaTwo_;
    }

    scalar gammaOne() const
    {
        return gammaOne_;
    }

    scalar gammaTwo() const
    {
        return gammaTwo_;
    }

    scalar aOne() const
    {
        return aOne_;
    }

protected:
    void clipValues(const std::shared_ptr<domain> domain);

    virtual void updateFOneBlending(const std::shared_ptr<domain> domain);

    virtual void
    updateTurbulentProduction(const std::shared_ptr<domain> domain);

    // Curvature correction (Smirnov-Menter 2009)

    void setupCurvatureCorrection(const std::shared_ptr<domain> domain);

    void updateCurvatureCorrection(const std::shared_ptr<domain> domain);

    // Other

    void updateTurbulentDynamicViscosity(
        const std::shared_ptr<domain> domain) override;
};

} /* namespace accel */

#endif // SHEARSTRESSTRANSPORTMODEL_H
