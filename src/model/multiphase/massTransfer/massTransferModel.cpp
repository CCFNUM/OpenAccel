// File       : massTransferModel.cpp
// Created    : Thu Oct 01 2026
// Author     : OpenAccel
// Description: Mass-transfer model base helpers and factory
// Copyright 2026 CCFNUM HSLU T&A. All Rights Reserved.

#include "massTransferModel.h"
#include "cavitationModel.h"

// std
#include <stdexcept>

namespace accel
{

void massTransferModel::compute(std::span<const double> pressure,
                                std::span<const double> alphaVapor,
                                std::span<const double> rhoLiquid,
                                std::span<const double> rhoVapor,
                                std::span<double> mdot) const
{
    const std::size_t n = mdot.size();
    if (pressure.size() != n || alphaVapor.size() != n ||
        rhoLiquid.size() != n || rhoVapor.size() != n)
    {
        throw std::runtime_error(
            "massTransferModel::compute: input spans must have equal size");
    }

    for (std::size_t i = 0; i < n; ++i)
    {
        mdot[i] = rate({pressure[i], alphaVapor[i], rhoLiquid[i], rhoVapor[i]});
    }
}

std::unique_ptr<massTransferModel>
createMassTransferModel(const fluidPairModel::massTransfer& cfg)
{
    switch (cfg.option_)
    {
        case massTransferModelOption::none:
            return nullptr;

        case massTransferModelOption::cavitation:
            return createCavitationModel(cfg);
    }
    throw std::runtime_error("createMassTransferModel: unsupported model");
}

} /* namespace accel */
