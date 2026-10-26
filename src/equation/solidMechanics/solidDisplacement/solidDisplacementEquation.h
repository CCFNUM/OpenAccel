// File       : solidDisplacementEquation.cpp
// Created    : Thu Dec 04 2025 08:42:10 (+0100)
// Author     : Mhamad Mahdi Alloush
// Description: Solid displacement equation with Aitken acceleration for
//              structural mechanics
// Copyright 2025 CCFNUM HSLU T&A. All Rights Reserved.

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

    // true if this timestep started with a solid predictor
    bool predicted() const
    {
        return predicted_;
    }

    bool movedMeshBoundary() const override
    {
        return predicted_;
    }

protected:
    void setResidualScales_() override;

    stk::mesh::PartVector collectDirichletBoundaryParts_();

private:
    std::unique_ptr<Assembler> assembler_;

    // solid solve with the previous-timestep fluid load, as iteration 0
    void predict_();

    // not solved until all interfaces reached coupling_start_time
    bool frozen_ = false;

    bool predicted_ = false;
};

} /* namespace accel */

#endif // SOLIDDISPLACEMENTEQUATION_H
