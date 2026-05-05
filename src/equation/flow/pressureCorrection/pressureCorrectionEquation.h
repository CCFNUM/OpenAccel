// File       : pressureCorrectionEquation.h
// Created    : Thu Mar 14 2024 12:50:04 (+0100)
// Author     : Fabian Wermelinger
// Description: Pressure correction (continuity) equation
// Copyright 2024 CCFNUM HSLU T&A. All Rights Reserved.

#ifndef PRESSURECORRECTIONEQUATION_H
#define PRESSURECORRECTIONEQUATION_H

// code
#include "equation.h"
#include "flowModel.h"
#include "linearSystem.h"
#include "macros.h"
#include "pressureCorrectionAssembler.h"

namespace accel
{

class pressureCorrectionEquation : public equation, public linearSystem<1>
{
private:
    flowModel* model_;

    using Assembler = pressureCorrectionAssembler;

public:
    static constexpr equationID ID = equationID::pressureCorrection;

    pressureCorrectionEquation(realm* realm, flowModel* model);

    void checkDomain(const std::shared_ptr<domain> domain) override;

    bool isConverged() const override;

    void setup() override;

    void initialize() override;

    void postInitialize() override;

    void preSolve() override;

    void solve() override;

    void postSolve() override;

    void preTimeStep() override;

    void printScales() override;

    equationID getID() override
    {
        return ID;
    }

protected:
    void setResidualScales_() override;

    std::unique_ptr<Assembler> assembler_;
};

} /* namespace accel */

#endif // PRESSURECORRECTIONEQUATION_H
