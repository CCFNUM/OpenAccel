// File       : physicsConvergence.cpp
// Created    : Thu Feb 26 2026
// Author     : Mhamad Mahdi Alloush
// Description: Physics-based convergence checks for coupled simulations
// Copyright 2026 CCFNUM HSLU T&A. All Rights Reserved.

#include "physicsConvergence.h"
#include "domain.h"
#include "fsiConvergence.h"
#include "interface.h"
#include "messager.h"
#include "simulation.h"

namespace accel
{

std::vector<std::unique_ptr<physicsConvergence>>
physicsConvergence::create(simulation& sim)
{
    const auto& cfg = sim.controlsRef()
                          .solverRef()
                          .solverControl_.basicSettings_.convergenceCriteria_
                          .physicsConvergence_;

    bool fsi = false;
    bool forceBalanceRequested = false;
    for (const auto& domainPtr : sim.domainVector())
    {
        for (const interface* interf : domainPtr->interfacesRef())
        {
            if (interf->isFluidSolidType())
            {
                fsi = true;
                forceBalanceRequested =
                    forceBalanceRequested ||
                    domainPtr->fluidSolidCouplingRef(interf).writeForceBalance_;
            }
        }
    }

    std::vector<std::unique_ptr<physicsConvergence>> checks;

    if (!cfg.enabled_)
    {
        // FSI outputs are written by the FSI check only
        if ((forceBalanceRequested || cfg.diagnosticsEnabled_) &&
            messager::master())
        {
            warningMsg("write_force_balance / physics_convergence diagnostics "
                       "have no effect: physics_convergence is not enabled");
        }
        return checks;
    }

    if (fsi)
    {
        checks.push_back(std::make_unique<fsiConvergence>(sim));
    }

    return checks;
}

const std::vector<std::unique_ptr<equation>>&
physicsConvergence::equations_() const
{
    return sim_.equationVector_;
}

} // namespace accel
