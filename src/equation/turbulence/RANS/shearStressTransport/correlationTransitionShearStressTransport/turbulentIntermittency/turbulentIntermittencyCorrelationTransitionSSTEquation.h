// File       : turbulentIntermittencyCorrelationTransitionSSTEquation.h
// Created    : Sun Dec 29 2024
// Author     : Adam Fares
// Description: Gamma equation for correlation-based transition SST (Menter
//              2015)
// Copyright 2024 CCFNUM HSLU T&A. All Rights Reserved.

#ifndef TURBULENTINTERMITTENCYCORRELATIONTRANSITIONSSTEQUATION_H
#define TURBULENTINTERMITTENCYCORRELATIONTRANSITIONSSTEQUATION_H

#include "correlationTransitionShearStressTransportModel.h"
#include "equation.h"
#include "linearSystem.h"
#include "macros.h"
#include "turbulentIntermittencyCorrelationTransitionSSTAssembler.h"
#include "types.h"

namespace accel
{

class turbulentIntermittencyCorrelationTransitionSSTEquation
    : public equation,
      public linearSystem<1>
{
private:
    correlationTransitionShearStressTransportModel* model_;

    using Assembler = turbulentIntermittencyCorrelationTransitionSSTAssembler;

public:
    static constexpr equationID ID =
        equationID::turbulentIntermittencyCorrelation;

    turbulentIntermittencyCorrelationTransitionSSTEquation(
        realm* realm,
        correlationTransitionShearStressTransportModel* model);

    void checkDomain(const std::shared_ptr<domain> domain) override;

    bool isConverged() const override;

    void setup() override;

    void initialize() override;

    void postInitialize() override;

    void preSolve() override;

    void solve() override;

    void preTimeStep() override;

    void printScales() override;

    equationID getID() override
    {
        return ID;
    }

protected:
    void setResidualScales_() override;

    void applyDependencyUpdates_(const domain* domain,
                                 const stk::mesh::EntityRank entityRank,
                                 STKScalarField& stk_dst) override;

private:
    std::unique_ptr<Assembler> assembler_;
};

} /* namespace accel */

#endif // TURBULENTINTERMITTENCYCORRELATIONTRANSITIONSSTEQUATION_H
