// File       : volumeFractionEquation.h
// Created    : Sun Jan 26 2025 22:57:38 (+0100)
// Author     : Mhamad Mahdi Alloush
// Description: Volume fraction transport equation for multiphase free-surface
// flow
// Copyright (c) 2025 CCFNUM, Lucerne University of Applied Sciences and
// Arts.
// SPDX-License-Identifier: BSD-3-Clause

#ifndef VOLUMEFRACTIONEQUATION_H
#define VOLUMEFRACTIONEQUATION_H

// code
#include "equation.h"
#include "freeSurfaceFlowModel.h"
#include "linearSystem.h"
#include "volumeFractionAssembler.h"

namespace accel
{

class volumeFractionEquation : public equation, public linearSystem<1>
{
private:
    freeSurfaceFlowModel* model_;

    using Assembler = volumeFractionAssembler;

    label phaseIndex_ = -1;

public:
    static constexpr equationID ID = equationID::volumeFraction;

    volumeFractionEquation(realm* realm,
                           freeSurfaceFlowModel* model,
                           label phaseIndex);

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

private:
    std::unique_ptr<Assembler> assembler_;
};

} /* namespace accel */

#endif // VOLUMEFRACTIONEQUATION_H
