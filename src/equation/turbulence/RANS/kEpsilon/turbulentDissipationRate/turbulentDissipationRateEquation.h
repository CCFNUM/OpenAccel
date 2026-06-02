// File       : turbulentDissipationRateEquation.h
// Created    : Thu Feb 22 2025 13:38:51 (+0100)
// Author     : ACHRAF NAGIHI
// Description: Turbulent dissipation rate equation for the k-epsilon model
// Copyright 2025 CCFNUM HSLU T&A. All Rights Reserved.

#ifndef TURBULENTDISSIPATIONRATEEQUATION_H
#define TURBULENTDISSIPATIONRATEEQUATION_H

#include "equation.h"
#include "kEpsilonModel.h"
#include "linearSystem.h"
#include "macros.h"
#include "turbulentDissipationRateAssembler.h"
#include "types.h"

namespace accel
{

class turbulentDissipationRateEquation : public equation, public linearSystem<1>
{
private:
    kEpsilonModel* model_;

    using Assembler = turbulentDissipationRateAssembler;

public:
    static constexpr equationID ID = equationID::turbulentDissipationRate;

    turbulentDissipationRateEquation(realm* realm, kEpsilonModel* model);

    void checkDomain(const std::shared_ptr<domain> domain) override;

    bool isConverged() const override;

    void setup() override;

    void initialize() override;

    void postInitialize() override;

    void preSolve() override;

    void solve() override;

    void preTimeStep() override;

    equationID getID() override
    {
        return ID;
    }

protected:
    void setResidualScales_() override;

private:
    std::unique_ptr<Assembler> assembler_;
};

} /* namespace accel */

#endif // TURBULENTDISSIPATIONRATEEQUATION_H
