// File       : convergenceAcceleration.cpp
// Created    : Thu Feb 26 2026
// Author     : Mhamad Mahdi Alloush
// Description: Generic convergence acceleration (Aitken, IQN-ILS)
// Copyright (c) 2026 CCFNUM, Lucerne University of Applied
// Sciences and Arts.
// SPDX-License-Identifier: BSD-3-Clause

#include "convergenceAcceleration.h"
#include "messager.h"

#include <algorithm>
#include <cmath>

namespace accel
{

namespace
{
bool solveDenseSystem_(std::vector<scalar>& A,
                       std::vector<scalar>& b,
                       std::vector<scalar>& x,
                       size_t n)
{
    x.assign(n, 0.0);
    for (size_t k = 0; k < n; ++k)
    {
        size_t pivot = k;
        scalar maxVal = std::abs(A[k * n + k]);
        for (size_t i = k + 1; i < n; ++i)
        {
            const scalar val = std::abs(A[i * n + k]);
            if (val > maxVal)
            {
                maxVal = val;
                pivot = i;
            }
        }

        if (maxVal < SMALL)
        {
            return false;
        }

        if (pivot != k)
        {
            for (size_t j = k; j < n; ++j)
            {
                std::swap(A[k * n + j], A[pivot * n + j]);
            }
            std::swap(b[k], b[pivot]);
        }

        const scalar diag = A[k * n + k];
        for (size_t i = k + 1; i < n; ++i)
        {
            const scalar factor = A[i * n + k] / diag;
            A[i * n + k] = 0.0;
            for (size_t j = k + 1; j < n; ++j)
            {
                A[i * n + j] -= factor * A[k * n + j];
            }
            b[i] -= factor * b[k];
        }
    }

    for (size_t i = n; i-- > 0;)
    {
        scalar sum = b[i];
        for (size_t j = i + 1; j < n; ++j)
        {
            sum -= A[i * n + j] * x[j];
        }
        if (std::abs(A[i * n + i]) < SMALL)
        {
            return false;
        }
        x[i] = sum / A[i * n + i];
    }

    return true;
}
} // namespace

convergenceAcceleration::convergenceAcceleration(const Config& cfg)
    : type_(cfg.type), aitkenOmega_(cfg.aitkenInitialOmega),
      aitkenOmegaInit_(cfg.aitkenInitialOmega),
      aitkenOmegaMin_(cfg.aitkenOmegaMin), aitkenOmegaMax_(cfg.aitkenOmegaMax),
      iqnIlsWindow_(cfg.iqnIlsWindow),
      iqnIlsRegularization_(cfg.iqnIlsRegularization)
{
}

void convergenceAcceleration::resetForTimeStep()
{
    if (type_ == accelerationType::aitken)
    {
        aitkenIter_ = 0;
        aitkenOmega_ = aitkenOmegaInit_;
        aitkenResidualPrev_.clear();
    }
    if (type_ == accelerationType::iqn_ils)
    {
        iqnIlsResidualHistory_.clear();
        iqnIlsUpdateHistory_.clear();
    }
}

const Vector& convergenceAcceleration::apply(const Vector& correction,
                                             scalar baseRelax,
                                             std::vector<Vector>& scratch,
                                             scalar& outRelaxValue)
{
    outRelaxValue = baseRelax;
    if (type_ == accelerationType::none)
    {
        return correction;
    }

    if (type_ == accelerationType::aitken)
    {
        outRelaxValue = computeAitkenOmega_(correction);
        return correction;
    }

    if (type_ == accelerationType::iqn_ils)
    {
        if (scratch.empty())
        {
            scratch.emplace_back();
        }
        Vector& update = scratch[0];
        computeIqnIlsUpdate_(correction, baseRelax, update);
        outRelaxValue = 1.0;
        return update;
    }

    return correction;
}

scalar convergenceAcceleration::computeAitkenOmega_(const Vector& correction)
{
    const size_t n = correction.size();

    if (aitkenIter_ == 0)
    {
        aitkenResidualPrev_ = correction;
        aitkenIter_++;
        return aitkenOmegaInit_;
    }

    scalar normPrevSq = 0.0;
    scalar normCurrSq = 0.0;
    scalar dotPrevCurr = 0.0;

    for (size_t i = 0; i < n; ++i)
    {
        normPrevSq += aitkenResidualPrev_[i] * aitkenResidualPrev_[i];
        normCurrSq += correction[i] * correction[i];
        dotPrevCurr += aitkenResidualPrev_[i] * correction[i];
    }

    messager::sumReduce(normPrevSq);
    messager::sumReduce(normCurrSq);
    messager::sumReduce(dotPrevCurr);

    const scalar normDrSq = normCurrSq - 2.0 * dotPrevCurr + normPrevSq;
    const scalar numerator = normPrevSq - dotPrevCurr;

    if (normDrSq > SMALL)
    {
        const scalar omegaNew = aitkenOmega_ * (numerator / normDrSq);
        aitkenOmega_ = std::clamp(omegaNew, aitkenOmegaMin_, aitkenOmegaMax_);
    }

    aitkenResidualPrev_ = correction;
    aitkenIter_++;

    return aitkenOmega_;
}

void convergenceAcceleration::computeIqnIlsUpdate_(const Vector& correction,
                                                   scalar baseRelax,
                                                   Vector& update)
{
    const size_t n = correction.size();
    update.assign(n, 0.0);
    if (n == 0)
    {
        return;
    }

    if (iqnIlsWindow_ < 1 || iqnIlsResidualHistory_.size() < 2)
    {
        for (size_t i = 0; i < n; ++i)
        {
            update[i] = baseRelax * correction[i];
        }
        iqnIlsResidualHistory_.push_back(correction);
        iqnIlsUpdateHistory_.push_back(update);
        if (iqnIlsResidualHistory_.size() >
            static_cast<size_t>(iqnIlsWindow_ + 1))
        {
            iqnIlsResidualHistory_.pop_front();
            iqnIlsUpdateHistory_.pop_front();
        }
        return;
    }

    const size_t available = iqnIlsResidualHistory_.size() - 1;
    const size_t m = std::min(static_cast<size_t>(iqnIlsWindow_), available);
    if (m == 0)
    {
        for (size_t i = 0; i < n; ++i)
        {
            update[i] = baseRelax * correction[i];
        }
        iqnIlsResidualHistory_.push_back(correction);
        iqnIlsUpdateHistory_.push_back(update);
        return;
    }

    std::vector<size_t> indices(m);
    const size_t start = iqnIlsResidualHistory_.size() - m;
    for (size_t j = 0; j < m; ++j)
    {
        indices[j] = start + j;
    }

    std::vector<scalar> A(m * m, 0.0);
    std::vector<scalar> b(m, 0.0);

    std::vector<scalar> dr(m, 0.0);
    for (size_t k = 0; k < n; ++k)
    {
        for (size_t i = 0; i < m; ++i)
        {
            const size_t idx = indices[i];
            dr[i] = iqnIlsResidualHistory_[idx][k] -
                    iqnIlsResidualHistory_[idx - 1][k];
            b[i] += dr[i] * correction[k];
        }
        for (size_t i = 0; i < m; ++i)
        {
            for (size_t j = 0; j < m; ++j)
            {
                A[i * m + j] += dr[i] * dr[j];
            }
        }
    }

    messager::sumReduce(b);
    messager::sumReduce(A);

    for (size_t i = 0; i < m; ++i)
    {
        A[i * m + i] += iqnIlsRegularization_;
    }

    std::vector<scalar> coeffs;
    if (!solveDenseSystem_(A, b, coeffs, m))
    {
        for (size_t i = 0; i < n; ++i)
        {
            update[i] = baseRelax * correction[i];
        }
        iqnIlsResidualHistory_.push_back(correction);
        iqnIlsUpdateHistory_.push_back(update);
        if (iqnIlsResidualHistory_.size() >
            static_cast<size_t>(iqnIlsWindow_ + 1))
        {
            iqnIlsResidualHistory_.pop_front();
            iqnIlsUpdateHistory_.pop_front();
        }
        return;
    }

    update = correction;
    for (size_t i = 0; i < m; ++i)
    {
        const size_t idx = indices[i];
        const Vector& dx = iqnIlsUpdateHistory_[idx];
        const scalar ci = coeffs[i];
        for (size_t k = 0; k < n; ++k)
        {
            update[k] -= ci * dx[k];
        }
    }

    iqnIlsResidualHistory_.push_back(correction);
    iqnIlsUpdateHistory_.push_back(update);
    if (iqnIlsResidualHistory_.size() > static_cast<size_t>(iqnIlsWindow_ + 1))
    {
        iqnIlsResidualHistory_.pop_front();
        iqnIlsUpdateHistory_.pop_front();
    }
}

} // namespace accel
