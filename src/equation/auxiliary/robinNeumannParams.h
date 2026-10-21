// File       : robinNeumannParams.h
// Created    : Thu Oct 01 2026 20:14:34 (+0300)
// Author     : Mohammad Nasreddine
// Description: Robin-Neumann FSI helpers of the pressure-correction equations:
//              per-timestep (h_s, rho_s) and the robin_alpha.out writer
// Copyright 2026 CCFNUM HSLU T&A. All Rights Reserved.

#ifndef ROBINNEUMANNPARAMS_H
#define ROBINNEUMANNPARAMS_H

// code
#include "types.h"

#include <fstream>
#include <memory>
#include <vector>

namespace accel
{

// Forward declarations
class domain;
class robinNeumannState;
class simulation;

struct robinNeumannParams
{
    scalar hs = 0.0;
    scalar rhoSolid = 0.0;

    // false if no Robin-Neumann FSI interface exists or the solid material
    // read failed (E <= 0 or rho_s <= 0)
    bool valid = false;
};

// true if any fluid-solid interface of the given domains uses Robin-Neumann
// coupling
bool anyRobinNeumannFsi(
    const std::vector<std::shared_ptr<domain>>& domainVector);

// h_s = a_p*dt (or solid_thickness) and rho_s, uniform solid assumed;
// collective (all ranks)
robinNeumannParams computeRobinNeumannParams(
    simulation& sim,
    const std::vector<std::shared_ptr<domain>>& domainVector,
    scalar dt);

// reduce the c1 statistics and write a robin_alpha.out row; collective, file IO
// on master
void writeRobinC1Diagnostics(
    const robinNeumannState& assemblerState,
    const std::vector<std::shared_ptr<domain>>& domainVector,
    std::shared_ptr<std::ofstream>& stream,
    scalar simTime);

} // namespace accel

#endif // ROBINNEUMANNPARAMS_H
