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

class simulation;

class physicsConvergence
{
public:
    explicit physicsConvergence(simulation& sim);

    bool enabled() const;

    void resetForTimeStep();

    void update();

    bool isConverged() const;

private:
    simulation& sim_;

    std::map<label, std::vector<scalar>> fsiInterfaceDispPrev_;
    std::map<label, scalar> fsiInterfaceResidualNormMax_;
    // running max of ||D_total|| at the interface this timestep, used for
    // the second (displacement-normalized) residual definition combined via
    // min() with the residual-normalized one above (s4f-style dual norm).
    std::map<label, scalar> fsiInterfaceMaxTotalDisplNorm_;
    std::map<label, scalar> fsiInterfaceResidualNorms_;
    scalar fsiInterfaceResidualNorm_ = 0.0;

    std::map<label, std::vector<scalar>> fsiInterfaceTractionPrev_;
    std::map<label, scalar> fsiInterfaceTractionResidualNormMax_;
    std::map<label, scalar> fsiForceResidualNorms_;
    scalar fsiForceResidualNorm_ = 0.0;

    std::map<std::string, std::map<label, std::shared_ptr<std::ofstream>>>
        residualStreams_;

    void updateFsiInterfaceResidual_(bool writeResiduals);

    void updateFsiForceResidual_(bool writeResiduals);

    void initializeResidualFile_(label interfIdx,
                                 const std::string& interfName,
                                 const std::string& criterionName);

    void writeResidualLine_(const std::string& criterionName,
                            label interfIdx,
                            scalar residualNorm);
};

} // namespace accel

#endif // PHYSICSCONVERGENCE_H
