// File       : massTransferModel.h
// Created    : Thu Oct 01 2026
// Author     : OpenAccel
// Description: Interface of interphase mass-transfer models for a liquid/vapor
//              fluid pair. Models are pointwise kernels, independent of mesh
//              and field storage.
// Copyright 2026 CCFNUM HSLU T&A. All Rights Reserved.

#ifndef MASSTRANSFERMODEL_H
#define MASSTRANSFERMODEL_H

// std
#include <memory>
#include <span>
#include <string>

#include "domain.h"

namespace accel
{

// Local state a mass-transfer model may depend on (SI units)
struct massTransferState
{
    double pressure = 0.0;    // absolute static pressure [Pa]
    double alphaVapor = 0.0;  // vapor volume fraction [-]
    double rhoLiquid = 0.0;   // liquid density [kg/m^3]
    double rhoVapor = 0.0;    // vapor density [kg/m^3]
};

// Mass-transfer rate per unit volume mdot_lv [kg/(m^3 s)]; positive for
// liquid -> vapor, negative for vapor -> liquid.
class massTransferModel
{
public:
    virtual ~massTransferModel() = default;

    virtual std::string name() const = 0;

    // raw (un-relaxed) rate at one point
    virtual double rate(const massTransferState& state) const = 0;

    // rate under-relaxation factor omega in (0,1]
    virtual double underRelaxation() const
    {
        return 1.0;
    }

    // whether the volumetric expansion enters the pressure equation
    virtual bool includeContinuitySource() const
    {
        return false;
    }

    // Fill `mdot` for n points (raw rates).
    void compute(std::span<const double> pressure,
                 std::span<const double> alphaVapor,
                 std::span<const double> rhoLiquid,
                 std::span<const double> rhoVapor,
                 std::span<double> mdot) const;

    // mdot_new = omega*mdot_raw + (1-omega)*mdot_prev
    double relax(double mdotRaw, double mdotPrev) const
    {
        const double w = underRelaxation();
        return w * mdotRaw + (1.0 - w) * mdotPrev;
    }

    // Volumetric source of the mixture continuity equation,
    //   div(u) = mdot (1/rho_v - 1/rho_l)   [1/s]
    static double
    volumeSource(double mdot, double rhoLiquid, double rhoVapor)
    {
        if (!(rhoLiquid > 0.0) || !(rhoVapor > 0.0))
            return 0.0;
        return mdot * (1.0 / rhoVapor - 1.0 / rhoLiquid);
    }
};

// Create the model described by `cfg`; returns nullptr for `none`.
std::unique_ptr<massTransferModel>
createMassTransferModel(const fluidPairModel::massTransfer& cfg);

} /* namespace accel */

#endif // MASSTRANSFERMODEL_H
