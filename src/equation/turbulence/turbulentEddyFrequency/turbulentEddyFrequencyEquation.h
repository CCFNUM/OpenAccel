// File       : turbulentEddyFrequencyEquation.h
// Created    : Thu Mar 14 2024 12:50:04 (+0100)
// Author     : Mhamad Mahdi Alloush
// Description: Base class for turbulent eddy frequency transport equation
// Copyright 2024 CCFNUM HSLU T&A. All Rights Reserved.

#ifndef TURBULENTEDDYFREQUENCYEQUATION_H
#define TURBULENTEDDYFREQUENCYEQUATION_H

// code
#include "equation.h"
#include "linearSystem.h"
#include "macros.h"
#include "turbulenceModel.h"
#include "turbulentEddyFrequencyAssembler.h"
#include "types.h"

namespace accel
{

class turbulentEddyFrequencyEquation : public equation, public linearSystem<1>
{
private:
    turbulenceModel* model_;

public:
    static constexpr equationID ID = equationID::turbulentEddyFrequency;

    turbulentEddyFrequencyEquation(realm* realm, turbulenceModel* model);

    void checkDomain(const std::shared_ptr<domain> domain) override;

    bool isConverged() const override;

    virtual void setup() override;

    virtual void initialize() override;

    virtual void postInitialize() override;

    virtual void preSolve() override;

    virtual void solve() override = 0;

    void preTimeStep() override;

    void printScales() override;

    equationID getID() override
    {
        return ID;
    }

protected:
    void setResidualScales_() override;
};

} /* namespace accel */

#endif // TURBULENTEDDYFREQUENCYEQUATION_H
