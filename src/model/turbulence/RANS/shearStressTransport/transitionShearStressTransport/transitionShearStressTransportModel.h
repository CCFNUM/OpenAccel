// File : transitionShearStressTransportModel.h
// Created : Mon Jan 14 2025
// Author : Adam Fares
// Description: Langtry-Menter transition SST model coefficients and
// correlations
// Copyright (c) 2025 CCFNUM, Lucerne University of Applied
// Sciences and Arts.
// SPDX-License-Identifier: BSD-3-Clause

#ifndef TRANSITIONSHEARSTRESSTRANSPORTMODEL_H
#define TRANSITIONSHEARSTRESSTRANSPORTMODEL_H

#include "shearStressTransportModel.h"

namespace accel
{

class realm;

class transitionShearStressTransportModel : public shearStressTransportModel
{
private:
    // Model parameters: default Transition SST (Langtry Menter 2003)

    scalar ca1_ = 2.0;

    scalar ca2_ = 0.06;

    scalar ce1_ = 1.0;

    scalar ce2_ = 50.0;

    scalar cThetat_ = 0.03;

    scalar s1_ = 2.0;

    scalar sigmaF_ = 1.0;

    scalar sigmaThetat_ = 2.0;

public:
    transitionShearStressTransportModel(realm* realm);

    // initialize

    virtual void initializeTransitionOnsetReynoldsNumber(
        const std::shared_ptr<domain> domain) override;

    // update

    virtual void updateTransitionOnsetReynoldsNumber(
        const std::shared_ptr<domain> domain) override;

    // Access

    scalar ca1() const
    {
        return ca1_;
    }

    scalar ca2() const
    {
        return ca2_;
    }

    scalar ce1() const
    {
        return ce1_;
    }

    scalar ce2() const
    {
        return ce2_;
    }

    scalar cThetat() const
    {
        return cThetat_;
    }

    scalar s1() const
    {
        return s1_;
    }

    scalar sigmaF() const
    {
        return sigmaF_;
    }

    scalar sigmaThetat() const
    {
        return sigmaThetat_;
    }

protected:
    void updateFOneBlending(const std::shared_ptr<domain> domain) override;

    void
    updateTurbulentProduction(const std::shared_ptr<domain> domain) override;

private:
    void updateTransitionOnsetReynoldsNumberSideFields_(
        const std::shared_ptr<domain> domain);
};

} /* namespace accel */

#endif // TRANSITIONSHEARSTRESSTRANSPORTMODEL_H
