// File       : convergenceAcceleration.h
// Created    : Thu Feb 26 2026
// Author     : Mhamad Mahdi Alloush
// Description: Generic convergence acceleration (Aitken, IQN-ILS)
// Copyright (c) 2026 CCFNUM, Lucerne University of Applied
// Sciences and Arts.
// SPDX-License-Identifier: BSD-3-Clause

#ifndef CONVERGENCEACCELERATION_H
#define CONVERGENCEACCELERATION_H

#include "types.h"

namespace accel
{

class convergenceAcceleration
{
public:
    struct Config
    {
        accelerationType type = accelerationType::none;
        scalar aitkenInitialOmega = 1.0;
        scalar aitkenOmegaMin = 0.1;
        scalar aitkenOmegaMax = 1.0;
        label iqnIlsWindow = 5;
        scalar iqnIlsRegularization = 1e-8;
    };

    explicit convergenceAcceleration(const Config& cfg);

    bool enabled() const
    {
        return type_ != accelerationType::none;
    }

    void resetForTimeStep();

    const Vector& apply(const Vector& correction,
                        scalar baseRelax,
                        std::vector<Vector>& scratch,
                        scalar& outRelaxValue);

private:
    accelerationType type_ = accelerationType::none;

    // Aitken state
    scalar aitkenOmega_ = 1.0;
    scalar aitkenOmegaInit_ = 1.0;
    scalar aitkenOmegaMin_ = 0.1;
    scalar aitkenOmegaMax_ = 1.0;
    label aitkenIter_ = 0;
    Vector aitkenResidualPrev_;

    // IQN-ILS state
    label iqnIlsWindow_ = 5;
    scalar iqnIlsRegularization_ = 1e-8;
    std::deque<Vector> iqnIlsResidualHistory_;
    std::deque<Vector> iqnIlsUpdateHistory_;

    scalar computeAitkenOmega_(const Vector& correction);

    void computeIqnIlsUpdate_(const Vector& correction,
                              scalar baseRelax,
                              Vector& update);
};

} // namespace accel

#endif // CONVERGENCEACCELERATION_H
