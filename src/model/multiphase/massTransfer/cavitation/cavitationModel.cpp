// File       : cavitationModel.cpp
// Created    : Thu Oct 01 2026
// Author     : OpenAccel
// Description: Cavitation model base class and factory
// Copyright 2026 CCFNUM HSLU T&A. All Rights Reserved.

#include "cavitationModel.h"
#include "rayleighPlessetCavitationMassTransfer.h"

// std
#include <algorithm>
#include <cmath>
#include <stdexcept>

namespace accel
{

cavitationModel::cavitationModel(const fluidPairModel::massTransfer& cfg)
    : cfg_(cfg)
{
    if (cfg_.option_ != massTransferModelOption::cavitation)
    {
        throw std::runtime_error(
            "mass_transfer: cavitation model requires model: cavitation");
    }
    validate(cfg_, "mass_transfer");
}

void cavitationModel::validate(const fluidPairModel::massTransfer& cfg,
                               const std::string& path)
{
    auto fail = [&](const std::string& msg) {
        throw std::runtime_error(path + msg);
    };

    if (cfg.liquidPhase_.empty() || cfg.vaporPhase_.empty())
        fail(": liquid_phase and vapor_phase must be specified");
    if (cfg.liquidPhase_ == cfg.vaporPhase_)
        fail(".vapor_phase: must differ from liquid_phase");

    if (!(cfg.saturationPressure_ > 0.0) ||
        !std::isfinite(cfg.saturationPressure_))
        fail(".saturation_pressure: must be a positive (absolute) pressure "
             "[Pa]");
    if (!(cfg.nucleationSiteRadius_ > 0.0) ||
        !std::isfinite(cfg.nucleationSiteRadius_))
        fail(".nucleation_site_radius: must be positive [m]");
    if (!(cfg.nucleationSiteVolumeFraction_ >= 0.0 &&
          cfg.nucleationSiteVolumeFraction_ <= 1.0))
        fail(".nucleation_site_volume_fraction: must be in [0, 1]");
    if (!(cfg.vaporizationCoefficient_ >= 0.0) ||
        !std::isfinite(cfg.vaporizationCoefficient_))
        fail(".vaporization_coefficient: must be non-negative");
    if (!(cfg.condensationCoefficient_ >= 0.0) ||
        !std::isfinite(cfg.condensationCoefficient_))
        fail(".condensation_coefficient: must be non-negative");
    if (!(cfg.underRelaxation_ > 0.0 && cfg.underRelaxation_ <= 1.0))
        fail(".under_relaxation: must be in (0, 1]");
}

double cavitationModel::ratePressure(double pressure) const
{
    return cfg_.pressureClippingForRate_ ? std::max(pressure, 0.0) : pressure;
}

double cavitationModel::clampedAlpha(double alpha)
{
    return std::clamp(alpha, 0.0, 1.0);
}

bool cavitationModel::validDensities(const massTransferState& state)
{
    return state.rhoLiquid > 0.0 && state.rhoVapor > 0.0;
}

std::unique_ptr<cavitationModel>
createCavitationModel(const fluidPairModel::massTransfer& cfg)
{
    switch (cfg.cavitationModel_)
    {
        case cavitationModelOption::rayleighPlesset:
            return std::make_unique<RayleighPlessetCavitationMassTransfer>(
                cfg);
    }
    throw std::runtime_error("createCavitationModel: unsupported model");
}

} /* namespace accel */
