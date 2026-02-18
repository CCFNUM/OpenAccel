// File : kEpsilonModel.h
// Created : Thu Feb 22 2025 13:38:51 (+0100)
// Author : Achraf Nagihi
// Description: Standard k-epsilon turbulence model coefficients and wall
// functions
// Copyright (c) 2025 CCFNUM, Lucerne University of Applied Sciences
// and Arts.
// SPDX-License-Identifier: BSD-3-Clause

#ifndef KEPSILONMODEL_H
#define KEPSILONMODEL_H

#include "RANSModel.h"

namespace accel
{

class kEpsilonModel : public RANSModel
{
private:
    scalar sigmaK_ = 1;

    scalar sigmaEpsilon_ = 1.3;

    scalar CEpsilonOne_ = 1.45;

    scalar CEpsilonTwo_ = 1.9;

public:
    kEpsilonModel(realm* realm);

    // Wall-functions

    void updateEpsilonAtWalls(const std::shared_ptr<domain> domain);

    // Access

    scalar sigmaK() const
    {
        return sigmaK_;
    }

    scalar sigmaEpsilon() const
    {
        return sigmaEpsilon_;
    }

    scalar CEpsilonOne() const
    {
        return CEpsilonOne_;
    }

    scalar CEpsilonTwo() const
    {
        return CEpsilonTwo_;
    }

protected:
    // protected methods

    void clipValues(const std::shared_ptr<domain> domain);

    virtual void
    updateTurbulentProduction(const std::shared_ptr<domain> domain);

    virtual void updateTurbulentDynamicViscosity(
        const std::shared_ptr<domain> domain) override;

    virtual void
    clipMinDistToWall(const std::shared_ptr<domain> domain) override;
};

} /* namespace accel */

#endif // KEPSILONMODEL_H
