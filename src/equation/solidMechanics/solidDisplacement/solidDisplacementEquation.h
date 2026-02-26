// File       : solidDisplacementEquation.cpp
// Created    : Thu Dec 04 2025 08:42:10 (+0100)
// Author     : Mhamad Mahdi Alloush
// Description: Solid displacement equation with Aitken acceleration for
// structural mechanics
// Copyright (c) 2025 CCFNUM, Lucerne University of Applied
// Sciences and Arts.
// SPDX-License-Identifier: BSD-3-Clause

#ifndef SOLIDDISPLACEMENTEQUATION_H
#define SOLIDDISPLACEMENTEQUATION_H

// code
#include "equation.h"
#include "linearSystem.h"
#include "macros.h"
#include "realm.h"
#include "solidDisplacementAssembler.h"
#include "solidMechanicsModel.h"
#include "types.h"

namespace accel
{

class solidDisplacementEquation : public equation,
                                  public solidMechanicsModel,
                                  public linearSystem<SPATIAL_DIM>
{
    using Assembler = solidDisplacementAssembler;

public:
    static constexpr equationID ID = equationID::solidDisplacement;

    solidDisplacementEquation(realm* realm);

    void checkDomain(const std::shared_ptr<domain> domain) override;

    bool isConverged() const override;

    scalar fsiResidualNorm() const
    {
        return fsiResidualNorm_;
    }

    bool fsiActive() const
    {
        return fsiActive_;
    }

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

    void applyDependencyUpdates_(const domain* domain,
                                 const stk::mesh::EntityRank entityRank,
                                 STKScalarField& stk_dst) override;

    stk::mesh::PartVector collectDirichletBoundaryParts_();

private:
    std::unique_ptr<Assembler> assembler_;

    // Aitken acceleration state
    bool useAitken_ = false;
    scalar aitkenOmega_ = 1.0;
    scalar aitkenOmegaInit_ = 1.0;
    scalar aitkenOmegaMin_ = 0.1;
    scalar aitkenOmegaMax_ = 1.0;
    label aitkenIter_ = 0;
    Vector aitkenResidualPrev_;

    scalar computeAitkenOmega_(const Vector& correction);

    void initializeFsiResidualFile_(label interfIdx,
                                    const std::string& interfName);

    void
    writeFsiResidualLine_(label interfIdx, scalar omega, scalar residualNorm);

    // FSI-level Aitken state (per interface, indexed by interface index)
    std::map<label, std::shared_ptr<std::ofstream>> fsiResidualStreams_;
    std::map<label, std::vector<scalar>> fsiDfluidPrev_;
    std::map<label, std::vector<scalar>> fsiResidualPrev_;
    std::map<label, scalar> fsiAitkenOmega_;
    std::map<label, scalar> fsiResidualNormMax_;
    scalar fsiResidualNorm_ = 0.0;
    bool fsiActive_ = false;
    scalar fsiOmegaInit_ = 0.4;
    scalar fsiOmegaMin_ = 0.01;
    scalar fsiOmegaMax_ = 1.0;
};

} /* namespace accel */

#endif // SOLIDDISPLACEMENTEQUATION_H
