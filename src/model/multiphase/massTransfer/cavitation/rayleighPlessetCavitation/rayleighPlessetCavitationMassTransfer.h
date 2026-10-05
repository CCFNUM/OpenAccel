// File       : rayleighPlessetCavitationMassTransfer.h
// Created    : Thu Oct 01 2026
// Author     : OpenAccel
// Description: CFX-style Rayleigh-Plesset (Zwart-type) cavitation model
// Copyright 2026 CCFNUM HSLU T&A. All Rights Reserved.

#ifndef RAYLEIGHPLESSETCAVITATIONMASSTRANSFER_H
#define RAYLEIGHPLESSETCAVITATIONMASSTRANSFER_H

#include "cavitationModel.h"

namespace accel
{

// p < p_v : mdot =  F_vap * 3 r_nuc (1-a_v) rho_v / R_nuc
//                   * sqrt(2/3 (p_v-p)/rho_l)
// p > p_v : mdot = -F_cond * 3 a_v rho_v / R_nuc * sqrt(2/3 (p-p_v)/rho_l)
class RayleighPlessetCavitationMassTransfer : public cavitationModel
{
public:
    explicit RayleighPlessetCavitationMassTransfer(
        const fluidPairModel::massTransfer& cfg)
        : cavitationModel(cfg)
    {
    }

    std::string name() const override
    {
        return "RayleighPlessetCavitationMassTransfer";
    }

    double rate(const massTransferState& state) const override;
};

} /* namespace accel */

#endif // RAYLEIGHPLESSETCAVITATIONMASSTRANSFER_H
