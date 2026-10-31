// File       : robinNeumannParams.h
// Created    : Thu Oct 01 2026 20:14:34 (+0300)
// Author     : Mohammad Nasreddine
// Description: Robin-Neumann FSI helpers of the pressure-correction equations:
//              per-timestep (h_s, rho_s), the robin_alpha.out writer and the
//              equation steps shared by the single-phase and VOF equations
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
class fieldBroker;
class robinNeumannState;

struct robinNeumannParams
{
    scalar hs = 0.0;
    scalar rhoSolid = 0.0;

    // false without a Robin-Neumann interface or if E or rho_s is not > 0
    bool valid = false;
};

// true if a fluid-solid interface of these domains uses Robin-Neumann
bool anyRobinNeumannFsi(
    const std::vector<std::shared_ptr<domain>>& domainVector);

// h_s = a_p*dt (or solid_thickness) and rho_s of a uniform solid; collective
robinNeumannParams computeRobinNeumannParams(
    const fieldBroker& model,
    const std::vector<std::shared_ptr<domain>>& domainVector,
    scalar dt);

// reduce the c1 statistics and write a robin_alpha.out row; collective
void writeRobinC1Diagnostics(
    const robinNeumannState& assemblerState,
    const std::vector<std::shared_ptr<domain>>& domainVector,
    std::shared_ptr<std::ofstream>& stream,
    scalar simTime);

// Robin-Neumann steps of a pressure-correction equation (Tukovic et al. 2014)
class robinNeumannCoupling
{
public:
    robinNeumannCoupling(
        fieldBroker* model,
        robinNeumannState& state,
        const std::vector<std::shared_ptr<domain>>& domainVector);

    // Robin side fields: a_n, Q_vol (Tukovic Eq. 30) and the frozen term E
    void setup();

    // h_s and rho_s of this timestep; no-slip fallback if they are invalid
    void preTimeStep();

    // reset Q_vol and E once per FSI iteration, not per segregated_flow pass
    void startOfOuterIteration();

    // after the Q_vol accumulation of a corrector: later correctors reuse E
    void postSolve();

    // one robin_alpha.out row per timestep
    void postTimeStep();

private:
    fieldBroker* model_;

    robinNeumannState& state_;

    const std::vector<std::shared_ptr<domain>>& domainVector_;

    std::shared_ptr<std::ofstream> robinAlphaStream_;

    // timestep and FSI iteration of the last Q_vol and E reset
    label resetStep_ = -1;
    label resetIter_ = -1;
};

} // namespace accel

#endif // ROBINNEUMANNPARAMS_H
