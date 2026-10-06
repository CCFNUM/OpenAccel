// File       : cavitationModel.cpp
// Created    : Thu Oct 01 2026
// Author     : OpenAccel
// Description: Validation of the cavitation parameters
// Copyright 2026 CCFNUM HSLU T&A. All Rights Reserved.

#include "cavitationModel.h"
#include "macros.h"

// std
#include <cmath>

namespace accel
{

void validateCavitationModel(const fluidPairModel::massTransfer& cfg,
                             const std::string& path)
{
    if (cfg.liquidPhase_.empty() || cfg.vaporPhase_.empty())
        errorMsg(path + ": liquid_phase and vapor_phase must be specified");
    if (cfg.liquidPhase_ == cfg.vaporPhase_)
        errorMsg(path + ".vapor_phase: must differ from liquid_phase");
    if (!(cfg.saturationPressure_ > 0.0) ||
        !std::isfinite(cfg.saturationPressure_))
        errorMsg(path + ".saturation_pressure: must be a positive absolute "
                        "pressure [Pa]");
    if (!(cfg.nucleationSiteRadius_ > 0.0) ||
        !std::isfinite(cfg.nucleationSiteRadius_))
        errorMsg(path + ".nucleation_site_radius: must be positive [m]");
    if (!(cfg.nucleationSiteVolumeFraction_ >= 0.0 &&
          cfg.nucleationSiteVolumeFraction_ <= 1.0))
        errorMsg(path + ".nucleation_site_volume_fraction: must be in [0, 1]");
    if (!(cfg.vaporizationCoefficient_ >= 0.0) ||
        !std::isfinite(cfg.vaporizationCoefficient_))
        errorMsg(path + ".vaporization_coefficient: must be non-negative");
    if (!(cfg.condensationCoefficient_ >= 0.0) ||
        !std::isfinite(cfg.condensationCoefficient_))
        errorMsg(path + ".condensation_coefficient: must be non-negative");
    if (!(cfg.underRelaxation_ > 0.0 && cfg.underRelaxation_ <= 1.0))
        errorMsg(path + ".under_relaxation: must be in (0, 1]");
    if (!(cfg.reboudCorrectionExponent_ >= 0.0) ||
        !std::isfinite(cfg.reboudCorrectionExponent_))
        errorMsg(path + ".reboud_correction_exponent: must be non-negative");
}

} /* namespace accel */
