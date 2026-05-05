// File       : segregatedKEpsilonEquations.h
// Created    : Thu Feb 22 2025 13:38:51 (+0100)
// Author     : Achraf Nagihi
// Description: Segregated solution of the k-epsilon turbulence equations
// Copyright 2025 CCFNUM HSLU T&A. All Rights Reserved.

#ifndef SEGREGATEDKEPSILONEQUATIONS_H
#define SEGREGATEDKEPSILONEQUATIONS_H

// code
#include "kEpsilonModel.h"
#include "turbulentDissipationRateEquation.h"
#include "turbulentKineticEnergyKEpsilonEquation.h"

namespace accel
{

class segregatedKEpsilonEquations : public equation, public kEpsilonModel
{
public:
    static constexpr equationID ID = equationID::segregatedKEpsilon;

    segregatedKEpsilonEquations(realm* realm);

    void addDomain(std::shared_ptr<domain> domain) override;

    bool isConverged() const override;

    void setup() override;

    void initialize() override;

    void postInitialize() override;

    void preSolve() override;

    void solve() override;

    void postSolve() override;

    void preTimeStep() override;

    equationID getID() override
    {
        return ID;
    }

private:
    // instances for segregated equations

    std::unique_ptr<turbulentKineticEnergyKEpsilonEquation> tke_eq_;

    std::unique_ptr<turbulentDissipationRateEquation> tdr_eq_;
};

} /* namespace accel */

#endif // SEGREGATEDKEPSILONEQUATIONS_H
