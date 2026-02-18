// File : correlationTransitionShearStressTransportModel.h
// Created : Sun Dec 29 2024
// Author : Adam Fares
// Description: Correlation-based transition SST model (Menter 2015)
// Solves only gamma transport equation with local correlations
// Copyright (c) 2024 CCFNUM, Lucerne University of Applied Sciences and Arts.
// SPDX-License-Identifier: BSD-3-Clause

#ifndef CORRELATIONTRANSITIONSHEARSTRESSTRANSPORTMODEL_H
#define CORRELATIONTRANSITIONSHEARSTRESSTRANSPORTMODEL_H

#include "shearStressTransportModel.h"

namespace accel
{

class realm;

/**
 * @brief Correlation-based transition SST model based on Menter et al. (2015)
 *
 * Reference:
 * Menter, F.R., Smirnov, P.E., Liu, T. et al. A One-Equation Local
 * Correlation-Based Transition Model. Flow Turbulence Combust 95, 583-619
 * (2015). https://doi.org/10.1007/s10494-015-9622-4
 *
 * This model solves only one additional transport equation for gamma
 * (intermittency) and uses local correlations to determine transition onset.
 * Unlike the full Langtry-Menter (2009) model, it does not solve for ReTheta.
 */
class correlationTransitionShearStressTransportModel
    : public shearStressTransportModel
{
private:
    // Model parameters: Menter 2015 correlation-based transition model

    // Gamma equation source term constants (Menter 2015)
    scalar flength_ = 100.0; // Production length scale (absorbs ca1)
    scalar caTwo_ = 0.06;    // Destruction coefficient
    scalar ceTwo_ = 50.0;    // Destruction coefficient

    // Correlation constants for Re_theta_c (critical Re_theta)
    scalar Ctu1_ = 100.0;
    scalar Ctu2_ = 1000.0;
    scalar Ctu3_ = 1.0;

    // Pressure gradient function constants (FPG)
    scalar CPG1_ = 14.68;
    scalar CPG2_ = -7.34;
    scalar CPG3_ = 0.0;
    scalar CPG1_lim_ = 1.5;
    scalar CPG2_lim_ = 3.0;

    // k-equation modification constants
    scalar Ck_BLT_ = 1.0;     // Boundary layer transition coefficient
    scalar CSEP_ = 1.0;       // Separation coefficient
    scalar Retclim_ = 1100.0; // Limiting Re_theta

    // Freestream turbulence intensity (if specified, -1 means use local Tu)
    scalar fsti_ = -1.0;

    // Gamma equation diffusion coefficient (sigma_gamma)
    scalar sigmaGamma_ = 1.0;

public:
    correlationTransitionShearStressTransportModel(realm* realm);

    // Access functions for model constants

    scalar flength() const
    {
        return flength_;
    }

    scalar caTwo() const
    {
        return caTwo_;
    }

    scalar ceTwo() const
    {
        return ceTwo_;
    }

    scalar Ctu1() const
    {
        return Ctu1_;
    }

    scalar Ctu2() const
    {
        return Ctu2_;
    }

    scalar Ctu3() const
    {
        return Ctu3_;
    }

    scalar CPG1() const
    {
        return CPG1_;
    }

    scalar CPG2() const
    {
        return CPG2_;
    }

    scalar CPG3() const
    {
        return CPG3_;
    }

    scalar CPG1_lim() const
    {
        return CPG1_lim_;
    }

    scalar CPG2_lim() const
    {
        return CPG2_lim_;
    }

    scalar Ck_BLT() const
    {
        return Ck_BLT_;
    }

    scalar CSEP() const
    {
        return CSEP_;
    }

    scalar Retclim() const
    {
        return Retclim_;
    }

    scalar fsti() const
    {
        return fsti_;
    }

    void setFsti(scalar fsti)
    {
        fsti_ = fsti;
    }

    scalar sigmaGamma() const
    {
        return sigmaGamma_;
    }

    /**
     * @brief Compute the pressure gradient function FPG(lambda_0L)
     * @param lambda0L Pressure gradient parameter
     * @return FPG value
     */
    static scalar FPG(scalar lambda0L,
                      scalar CPG1,
                      scalar CPG2,
                      scalar CPG3,
                      scalar CPG1_lim,
                      scalar CPG2_lim);

protected:
    void updateFOneBlending(const std::shared_ptr<domain> domain) override;

    void
    updateTurbulentProduction(const std::shared_ptr<domain> domain) override;
};

} /* namespace accel */

#endif // CORRELATIONTRANSITIONSHEARSTRESSTRANSPORTMODEL_H
