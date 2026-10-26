// File       : physicsConvergence.h
// Created    : Thu Feb 26 2026
// Author     : Mhamad Mahdi Alloush
// Description: Physics-based convergence checks for coupled simulations
// Copyright 2026 CCFNUM HSLU T&A. All Rights Reserved.

#ifndef PHYSICSCONVERGENCE_H
#define PHYSICSCONVERGENCE_H

#include "types.h"

namespace accel
{

class equation;
class simulation;

class physicsConvergence
{
public:
    explicit physicsConvergence(simulation& sim) : sim_(sim)
    {
    }

    virtual ~physicsConvergence() = default;

    // the checks physics_convergence requests for the physics of this run
    static std::vector<std::unique_ptr<physicsConvergence>>
    create(simulation& sim);

    // start of a timestep, after the equations' preTimeStep()
    virtual void resetForTimeStep() = 0;

    // false while the coupled physics has not started
    virtual bool active() const = 0;

    // residuals of the outer iteration just solved
    virtual void update() = 0;

    // false: convergence falls back to the equation residuals
    virtual bool hasCriteria() const = 0;

    virtual bool isConverged() const = 0;

protected:
    simulation& sim_;

    const std::vector<std::unique_ptr<equation>>& equations_() const;
};

} // namespace accel

#endif // PHYSICSCONVERGENCE_H
