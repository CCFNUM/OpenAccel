// File       : cavitationModel.h
// Created    : Thu Oct 01 2026
// Author     : OpenAccel
// Description: Base class of all cavitation (liquid <-> vapor) mass-transfer
//              models. Holds the parameters common to cavitation models:
//              saturation pressure, rate under-relaxation, pressure clipping
//              and the continuity-source switch. Derived classes provide the
//              rate expression.
// Copyright 2026 CCFNUM HSLU T&A. All Rights Reserved.

#ifndef CAVITATIONMODEL_H
#define CAVITATIONMODEL_H

#include "massTransferModel.h"

namespace accel
{

class cavitationModel : public massTransferModel
{
public:
    // `cfg` must describe an active cavitation model; it is validated
    explicit cavitationModel(const fluidPairModel::massTransfer& cfg);

    // Check the cavitation parameters; throws std::runtime_error with `path`
    // (YAML path, e.g. fluid_pair_models[0].mass_transfer) prefixed
    static void validate(const fluidPairModel::massTransfer& cfg,
                         const std::string& path);

    double underRelaxation() const override
    {
        return cfg_.underRelaxation_;
    }

    bool includeContinuitySource() const override
    {
        return cfg_.includeContinuitySource_;
    }

    const fluidPairModel::massTransfer& config() const
    {
        return cfg_;
    }

    double saturationPressure() const
    {
        return cfg_.saturationPressure_;
    }

protected:
    // pressure used in the rate: the true local absolute pressure unless
    // clipping for the rate is requested (then limited to p >= 0)
    double ratePressure(double pressure) const;

    // vapor volume fraction limited to [0,1]
    static double clampedAlpha(double alpha);

    // true if the densities allow evaluating a rate
    static bool validDensities(const massTransferState& state);

    const fluidPairModel::massTransfer cfg_;
};

// Create the cavitation model selected by `cfg.cavitationModel_`
std::unique_ptr<cavitationModel>
createCavitationModel(const fluidPairModel::massTransfer& cfg);

} /* namespace accel */

#endif // CAVITATIONMODEL_H
