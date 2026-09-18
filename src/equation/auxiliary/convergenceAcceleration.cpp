// File       : convergenceAcceleration.cpp
// Created    : Thu Feb 26 2026
// Author     : Mhamad Mahdi Alloush
// Description: Generic convergence acceleration (Aitken, IQN-ILS)
// Copyright 2026 CCFNUM HSLU T&A. All Rights Reserved.

#include "convergenceAcceleration.h"
#include "messager.h"
#include <iostream>

#include <algorithm>
#include <cmath>

namespace accel
{

convergenceAcceleration::convergenceAcceleration(const Config& cfg)
    : type_(cfg.type), aitkenOmega_(cfg.aitkenInitialOmega),
      aitkenOmegaInit_(cfg.aitkenInitialOmega),
      aitkenOmegaMin_(cfg.aitkenOmegaMin), aitkenOmegaMax_(cfg.aitkenOmegaMax),
      iqnIlsWindow_(cfg.iqnIlsWindow), iqnIlsFilter_(cfg.iqnIlsFilter),
      iqnIlsWindowsReused_(cfg.iqnIlsWindowsReused)
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
        // keep secant columns from recent windows, drop the stale ones
        iqnWindow_++;
        iqnHavePrev_ = false;
        while (!iqnColWindow_.empty() &&
               iqnColWindow_.back() < iqnWindow_ - iqnIlsWindowsReused_)
        {
            iqnV_.pop_back();
            iqnW_.pop_back();
            iqnColWindow_.pop_back();
        }
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

    // a ghosting update can resize the local system: restart the history,
    // decided collectively since the resize is rank-local
    scalar restartFlag = (aitkenResidualPrev_.size() != n) ? 1.0 : 0.0;
    messager::sumReduce(restartFlag);
    if (restartFlag > 0.0)
    {
        aitkenResidualPrev_ = correction;
        aitkenIter_++;
        return aitkenOmega_;
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
        if (std::isfinite(omegaNew)) // std::clamp would pass NaN through
        {
            aitkenOmega_ =
                std::clamp(omegaNew, aitkenOmegaMin_, aitkenOmegaMax_);
        }
    }

    if (messager::master())
    {
        // one line per coupling iteration: a stalled coupling (omega pinned
        // at its bound) can then be told apart from a diverging one
        std::cout << "Aitken omega: " << aitkenOmega_ << std::endl;
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

    // a resize means the dof layout changed: restart collectively
    scalar restart = (iqnHavePrev_ && iqnPrevR_.size() != n) ? 1.0 : 0.0;
    messager::sumReduce(restart);
    if (restart > 0.0)
    {
        iqnV_.clear();
        iqnW_.clear();
        iqnColWindow_.clear();
        iqnHavePrev_ = false;
    }
    // an empty rank must still join every collective below
    update.assign(n, 0.0);

    auto gdot = [](const Vector& a, const Vector& b)
    {
        scalar d = 0.0;
        for (size_t k = 0; k < a.size(); ++k)
        {
            d += a[k] * b[k];
        }
        messager::sumReduce(d);
        return d;
    };

    // grow the secant history: V = dr, W = dx~ = du_prev + dr
    if (iqnHavePrev_)
    {
        Vector v(n), w(n);
        for (size_t k = 0; k < n; ++k)
        {
            v[k] = correction[k] - iqnPrevR_[k];
            w[k] = iqnPrevUpdate_[k] + v[k];
        }
        iqnV_.push_front(std::move(v));
        iqnW_.push_front(std::move(w));
        iqnColWindow_.push_front(iqnWindow_);
        while (iqnV_.size() > static_cast<size_t>(iqnIlsWindow_))
        {
            iqnV_.pop_back();
            iqnW_.pop_back();
            iqnColWindow_.pop_back();
        }
    }

    const size_t m = iqnV_.size();
    if (m == 0)
    {
        const scalar w0 =
            (aitkenOmegaInit_ > 0.0) ? aitkenOmegaInit_ : baseRelax;
        for (size_t k = 0; k < n; ++k)
        {
            update[k] = w0 * correction[k];
        }
        iqnPrevR_ = correction;
        iqnPrevUpdate_ = update;
        iqnHavePrev_ = true;
        return;
    }

    // modified Gram-Schmidt on the V columns, newest first; a column whose
    // orthogonal remainder falls under the filter is dropped for good
    std::vector<Vector> Q;
    std::vector<std::vector<scalar>> R;
    std::vector<size_t> kept;
    size_t filtered = 0;
    for (size_t j = 0; j < iqnV_.size();)
    {
        Vector v = iqnV_[j];
        const scalar orig = std::sqrt(gdot(v, v));
        std::vector<scalar> rj(Q.size(), 0.0);
        for (size_t i = 0; i < Q.size(); ++i)
        {
            const scalar h = gdot(Q[i], v);
            rj[i] = h;
            for (size_t k = 0; k < n; ++k)
            {
                v[k] -= h * Q[i][k];
            }
        }
        const scalar nrm = std::sqrt(gdot(v, v));
        if (!(nrm > iqnIlsFilter_ * orig) || !(orig > 0.0))
        {
            iqnV_.erase(iqnV_.begin() + j);
            iqnW_.erase(iqnW_.begin() + j);
            iqnColWindow_.erase(iqnColWindow_.begin() + j);
            filtered++;
            continue;
        }
        for (size_t k = 0; k < n; ++k)
        {
            v[k] /= nrm;
        }
        rj.push_back(nrm);
        Q.push_back(std::move(v));
        R.push_back(std::move(rj));
        kept.push_back(j);
        ++j;
    }

    const size_t mk = Q.size();
    if (mk == 0)
    {
        const scalar w0 =
            (aitkenOmegaInit_ > 0.0) ? aitkenOmegaInit_ : baseRelax;
        for (size_t k = 0; k < n; ++k)
        {
            update[k] = w0 * correction[k];
        }
        iqnPrevR_ = correction;
        iqnPrevUpdate_ = update;
        iqnHavePrev_ = true;
        return;
    }

    // solve R alpha = -Q^T r by back substitution
    std::vector<scalar> beta(mk, 0.0);
    for (size_t i = 0; i < mk; ++i)
    {
        beta[i] = -gdot(Q[i], correction);
    }
    std::vector<scalar> alpha(mk, 0.0);
    for (size_t i = mk; i-- > 0;)
    {
        scalar sum = beta[i];
        for (size_t j2 = i + 1; j2 < mk; ++j2)
        {
            sum -= R[j2][i] * alpha[j2];
        }
        alpha[i] = sum / R[i][i];
    }

    // x_{k+1} = x~_k + W alpha  ->  update = r + W alpha
    update = correction;
    scalar amax = 0.0;
    for (size_t i = 0; i < mk; ++i)
    {
        const Vector& w = iqnW_[kept[i]];
        const scalar ai = alpha[i];
        amax = std::max(amax, std::abs(ai));
        for (size_t k = 0; k < n; ++k)
        {
            update[k] += ai * w[k];
        }
    }

    if (messager::master())
    {
        std::cout << "IQN-ILS: cols " << m << " kept " << mk << " filtered "
                  << filtered << " |alpha|max " << amax << std::endl;
    }

    iqnPrevR_ = correction;
    iqnPrevUpdate_ = update;
    iqnHavePrev_ = true;
}

} // namespace accel
