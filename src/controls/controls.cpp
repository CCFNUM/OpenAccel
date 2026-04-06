// File       : controls.cpp
// Created    : Fri Aug 25 2023 12:55:24 (+0100)
// Author     : Mhamad Mahdi Alloush
// Description:
// Copyright (c) 2023 CCFNUM, Lucerne University of Applied Sciences and Arts.
// SPDX-License-Identifier: BSD-3-Clause

// code
#include "controls.h"
#include "messager.h"
#include "version.h"

// std
#include <algorithm>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>

namespace accel
{

// Constructors

controls::controls(fs::path cwd)
    : workingDirectory_(cwd), profiler_(messager::comm())
{
    // required states for correct restart
    restartParameter_.set_param(
        "timeStepCount", analysisType_.timeStepCount_, false, true);
    restartParameter_.set_param("globalIter", globalIter, false, true);
    for (int i = 0; i < analysisTypeDictionary::DT_ENTRIES; i++)
    {
        const std::string dt_id("dt_");
        restartParameter_.set_param(
            dt_id + std::to_string(i), getTimestep(-i), false, true);
    }
}

// Destructor

controls::~controls()
{
}

bool controls::isReducedStencil() const
{
    return solver_.solverControl_.expertParameters_.bandwidthReduction_;
}

bool controls::isTransient() const
{
    return analysisType_.transient_;
}

bool controls::isHighResolution() const
{
    return solver_.solverControl_.basicSettings_.advectionScheme_ ==
           advectionSchemeType::highResolution;
}

bool controls::isHighResolutionTurbulenceNumerics() const
{
    return solver_.solverControl_.basicSettings_.turbulenceNumerics_ ==
           advectionSchemeType::highResolution;
}

label controls::getNumberOfStates() const
{
    label numberOfStates = 1; // 1 for steady-state, 2 for first order accurate,
    // 3 for second order accurate
    if (analysisType_.transient_)
    {
        switch (solver_.solverControl_.basicSettings_.transientScheme_)
        {
            case transientSchemeType::firstOrderBackwardEuler:
                numberOfStates = 2;
                break;

            case transientSchemeType::secondOrderBackwardEuler:
                numberOfStates = 3;
                break;

            default:
                break;
        }
    }
    assert(numberOfStates - 1 <= analysisTypeDictionary::DT_ENTRIES);
    return numberOfStates;
}

stk::util::ParameterList& controls::getRestartParam()
{
    return restartParameter_;
}

void controls::setRestartParam()
{
    restartParameter_.set_value("timeStepCount", analysisType_.timeStepCount_);
    restartParameter_.set_value("globalIter", globalIter);
    for (int i = 0; i < analysisTypeDictionary::DT_ENTRIES; i++)
    {
        const std::string dt_id("dt_");
        restartParameter_.set_value(dt_id + std::to_string(i), getTimestep(-i));
    }
}

void controls::deserializeRestartParam(
    const stk::io::StkMeshIoBroker& io_broker)
{
    if (!isTransient())
    {
        // For steady-state, globalIter is recovered from the field
        // time_restored() in nodeField::initialize(). The dt_* and
        // timeStepCount parameters are only relevant to transient runs.
        io_broker.get_global("globalIter", globalIter, false);
        return;
    }

    io_broker.get_global("timeStepCount", analysisType_.timeStepCount_);
    io_broker.get_global("globalIter", globalIter);
    for (int i = 0; i < analysisTypeDictionary::DT_ENTRIES; i++)
    {
        const std::string dt_id("dt_");
        const int idx = timestepPosition_(-i);
        io_broker.get_global(dt_id + std::to_string(i),
                             analysisType_.timestep_[idx]);
    }
}

scalar controls::getTotalTime() const
{
    return analysisType_.totalTime_;
}

void controls::setTimestep(const scalar dt)
{
    const int idx = timestepPosition_(0);
    analysisType_.timestep_[idx] = dt;
}

scalar controls::getTimestep(const int i) const
{
    const int idx = timestepPosition_(i);
    return analysisType_.timestep_[idx];
}

scalar controls::getPhysicalTimescale() const
{
    return solver_.solverControl_.basicSettings_.convergenceControl_
        .physicalTimescale_;
}

label controls::getTimeStepCount() const
{
    return analysisType_.timeStepCount_;
}

void controls::advanceAndSetTimestep()
{
    LogFileInfo info(analysisType_.timeSteps_.mode_);
    info.dt = analysisType_.initialTimestep_;

#ifndef NDEBUG
    scalar dt_before_output_adjustment = info.dt;
    bool timestep_adjusted_for_output = false;
#endif /* NDEBUG */

    if (info.mode == timestepMode::periodicInterval)
    {
        periodicIntervalTimestep_(info);
    }
    else if (info.mode == timestepMode::specifiedInterval)
    {
        specifiedIntervalTimestep_(info);
    }
    else if (info.mode == timestepMode::adaptive)
    {
        courantAdaptiveTimestep_(info);
    }

    if (solver_.outputControl_.matchFinalTime_)
    {
        // WARNING [faw 2025-01-17]: current implementation writes output based
        // on a dump frequency. Restarting simulations with input that was
        // generated with match_final_time = true results in non-uniform time
        // steps in the total series. The final time is further written to the
        // results file regardless of write frequency. Default for
        // match_final_time is `false`.
        const scalar dt_end = analysisType_.totalTime_ - this->time;
        if (dt_end < info.dt)
        {
            info.dt = dt_end;
            info.action = "match_final_time";
        }
    }

    // Adjust timestep to match output times when using time_interval
    // This applies to all timestep modes to ensure synchronized output times
    if (solver_.outputControl_.outputFrequency_.option_ ==
        outputFrequencyType::timeInterval)
    {
        const scalar interval =
            solver_.outputControl_.outputFrequency_.timeInterval_;
        if (interval > 0.0)
        {
            // Calculate next scheduled output time
            // Using floor(time/interval + 1) to always get the NEXT multiple of
            // interval
            const scalar current_time = this->time;
            const scalar next_output_time =
                std::floor(current_time / interval + 1.0) * interval;
            const scalar time_to_output = next_output_time - current_time;

            // If next timestep would overshoot the output time, reduce it
            if (time_to_output > 1.0e-12 && info.dt > time_to_output)
            {
#ifndef NDEBUG
                dt_before_output_adjustment = info.dt;
                timestep_adjusted_for_output = true;
#endif /* NDEBUG */
                info.dt = time_to_output;
                info.action = "output_sync";

#ifndef NDEBUG
                // Console output for output synchronization (debug only, master
                // only)
                if (messager::master() && info.mode == timestepMode::adaptive)
                {
                    std::cout
                        << "[Adaptive dt] Time=" << std::scientific
                        << std::setprecision(4) << this->time
                        << " | dt adjusted for output: "
                        << dt_before_output_adjustment << " -> " << info.dt
                        << " (next output at t=" << next_output_time << ")"
                        << std::endl;
                }
#endif /* NDEBUG */
            }
        }
    }

    writeTimestepLog_(info);

    ++analysisType_.timeStepCount_; // 1. advance ID
    this->setTimestep(info.dt);     // 2. set new step @ID
    resetMaxCourant();
}

int controls::timestepPosition_(const int i) const
{
    // i = 0: t_{n+1} - t_n (current timestep)
    // i = -1: t_n - t_{n-1} (previous timestep)
    // i = -2: t_{n-1} - t_{n-2}
    // i = -3: t_{n-2} - t_{n-3}
    assert(-analysisTypeDictionary::DT_ENTRIES < i && i <= 0);

    return (analysisType_.timeStepCount_ + i +
            analysisTypeDictionary::DT_ENTRIES) %
           analysisTypeDictionary::DT_ENTRIES;
}

void controls::writeTimestepLog_(const LogFileInfo& info)
{
    // Check if user enabled timestep info output on this rank
    if (!dtLogFile_)
        return;

    static bool log_file_initialized = false;

    // Initialize file with header
    auto& out = *dtLogFile_;
    if (!log_file_initialized)
    {
        auto now = std::chrono::system_clock::now();
        auto in_time_t = std::chrono::system_clock::to_time_t(now);

        out << "# Accel Time Stepping Diagnostics: "
            << std::put_time(std::localtime(&in_time_t), "%c\n");
        out << "# Version: " << accel::PROJECT_VERSION << '\n';
        out << "# Git hash: " << accel::GIT_HASH << '\n';
        out << "# Git describe: " << accel::GIT_DESCRIBE << '\n';
        out << "#\n";

        if (info.mode == timestepMode::adaptive)
        {
            out << "# Mode: Adaptive\n";
            out << "# Time\t\tTimestep\tCo_max\t\tCo_target\tAction\n";
        }
        else if (info.mode == timestepMode::periodicInterval)
        {
            out << "# Mode: Periodic interval\n";
            out << "# Time\tTimestep\tPeriod\t";
            out << "Interval_start\tInterval_length\tAction\n";
        }
        else if (info.mode == timestepMode::specifiedInterval)
        {
            out << "# Mode: Specified interval\n";
            out << "# Time\tTimestep\t";
            out << "Interval_start\tInterval_length\tAction\n";
        }
        else
        {
            out << "# Mode: Constant\n";
            out << "# Time\tTimestep\tAction\n";
        }
        log_file_initialized = true;
    }

    out << std::scientific << std::setprecision(6);
    out << this->time << "\t" << info.dt << "\t";
    if (info.mode == timestepMode::adaptive)
    {
        out << info.maxCourant << "\t"
            << analysisType_.timeSteps_.timestepAdaptation_.courantNumber_
            << "\t" << info.action;
    }
    else if (info.mode == timestepMode::periodicInterval)
    {
        out << info.period << "\t" << info.intervalStart << "\t"
            << info.intervalLength << "\t" << info.action;
    }
    else if (info.mode == timestepMode::specifiedInterval)
    {
        out << info.intervalStart << "\t" << info.intervalLength << "\t"
            << info.action;
    }
    else
    {
        out << info.action;
    }
    out << "\n";
}

void controls::periodicIntervalTimestep_(LogFileInfo& info) const
{
    const auto& start_time = analysisType_.timeSteps_.startTime_;
    const auto& interval_length = analysisType_.timeSteps_.intervalLength_;
    const auto& timestep_interval = analysisType_.timeSteps_.timestepInterval_;

    assert(start_time.size() == 1);
    assert(interval_length.size() == 1);
    assert(timestep_interval.size() == 1);

    const scalar t_start = start_time.front();
    const scalar t_length = interval_length.front();
    const scalar T = analysisType_.timeSteps_.period_;

    static int k = 0;
    scalar t_test = (this->time - t_start) - k * T;

    // periodicity
    if (t_test > T ||
        std::abs(T - t_test) < 5.0 * std::numeric_limits<scalar>::epsilon())
    {
        ++k;
        t_test -= T;
        t_test = (t_test < 0) ? -t_test : t_test; // round-off
    }

    info.period = k;
    info.intervalStart = t_start + k * T;
    info.intervalLength = T;

    if (0 <= t_test &&
        (t_length - t_test) > 5.0 * std::numeric_limits<scalar>::epsilon())
    {
        const scalar dt_interval = timestep_interval.front();
        assert(dt_interval > 0);
        info.action = "inside_interval";
        info.dt = dt_interval;
    }
    else
    {
        assert(analysisType_.initialTimestep_ > 0);
        info.action = "outside_interval";
        info.dt = analysisType_.initialTimestep_;
    }
}

void controls::specifiedIntervalTimestep_(LogFileInfo& info)
{
    auto& start_time = analysisType_.timeSteps_.startTime_;
    auto& interval_length = analysisType_.timeSteps_.intervalLength_;
    auto& timestep_interval = analysisType_.timeSteps_.timestepInterval_;

    // either a list of start times with identical interval and timestep size or
    // all three specified
    assert((start_time.size() > 0 && interval_length.size() == 1 &&
            timestep_interval.size() == 1) ||
           (start_time.size() == interval_length.size() &&
            start_time.size() == timestep_interval.size()));

    const scalar t_start = start_time.front();
    const scalar t_length = interval_length.front();
    const scalar dt_interval = timestep_interval.front();
    const scalar t_test = this->time - t_start;

    if (t_test > t_length)
    {
        // assumes that two consecutive intervals are at least
        // analysisType_.initialTimestep_ apart (or dt_interval if it is larger
        // than analysisType_.initialTimestep_)
        if (start_time.size() > 1)
        {
            start_time.pop_front();
        }
        if (interval_length.size() > 1)
        {
            interval_length.pop_front();
        }
        if (timestep_interval.size() > 1)
        {
            timestep_interval.pop_front();
        }
    }

    info.intervalStart = t_start;
    info.intervalLength = t_length;

    if (0 <= t_test &&
        (t_length - t_test) > 5.0 * std::numeric_limits<scalar>::epsilon())
    {
        assert(dt_interval > 0);
        info.action = "inside_interval";
        info.dt = dt_interval;
    }
    else
    {
        assert(analysisType_.initialTimestep_ > 0);
        info.action = "outside_interval";
        info.dt = analysisType_.initialTimestep_;
    }
}

void controls::courantAdaptiveTimestep_(LogFileInfo& info) const
{
    const scalar dt_curr = getTimestep();
    const scalar maxCo = maxCourant_;
    const auto& dtInfo = analysisType_.timeSteps_;
    const auto& adapt = dtInfo.timestepAdaptation_;

    const label freq = std::max<label>(1, dtInfo.timestepUpdateFrequency_);
    const bool updateNow = ((analysisType_.timeStepCount_ % freq) == 0);
    if (!updateNow)
    {
        info.dt = dt_curr;
        return;
    }

    if (adapt.option_ == timestepAdaptationType::rmsCourant)
    {
        errorMsg("rms_courant adaptation is not implemented");
    }

    scalar dt_new = dt_curr;

    if (maxCo > 0 && dt_curr > 0)
    {
        // Calculate desired adjustment factor
        const scalar factor = adapt.courantNumber_ / maxCo;

        // Apply asymmetric damping to prevent aggressive timestep changes
        scalar damped_factor = factor;
        if (factor < 1.0)
        {
            // Decreasing timestep - limit reduction (default: max 20%
            // decrease)
            damped_factor = std::max(factor, adapt.timestepDecreaseFactor_);
            info.action = "decreased";
        }
        else if (factor > 1.0)
        {
            // Increasing timestep - limit increase (default: max 6%
            // increase)
            damped_factor = std::min(factor, adapt.timestepIncreaseFactor_);
            info.action = "increased";
        }

        dt_new = dt_curr * damped_factor;
    }

    dt_new = std::max(adapt.minTimestep_, dt_new);
    dt_new = std::min(adapt.maxTimestep_, dt_new);

#ifndef NDEBUG
    // Console output for adaptive changes (debug only, master only)
    if (messager::master() && std::abs(dt_new - dt_curr) > 1.0e-12)
    {
        std::cout << "[Adaptive dt] Time=" << std::scientific
                  << std::setprecision(4) << this->time
                  << " | Co_max=" << std::fixed << std::setprecision(3) << maxCo
                  << " (target=" << adapt.courantNumber_ << ")"
                  << " | dt: " << std::scientific << std::setprecision(4)
                  << dt_curr << " -> " << dt_new << " (" << info.action << ")"
                  << std::endl;
    }
#endif /* NDEBUG */

    info.maxCourant = maxCo;
    info.dt = dt_new;
}

void controls::updateMaxCourant(const scalar maxCourant)
{
    maxCourant_ = std::max(maxCourant_, maxCourant);
}

scalar controls::getMaxCourant() const
{
    return maxCourant_;
}

void controls::resetMaxCourant()
{
    maxCourant_ = -1;
}

// IO paths
fs::path controls::getWorkingDirectory() const
{
    assert(fs::exists(workingDirectory_));
    return workingDirectory_;
}

fs::path controls::getPostProcessingDirectory() const
{
    fs::path dir = getWorkingDirectory() / "postProcessing";
    fs::create_directories(dir);
    return dir;
}

fs::path controls::getResidualDirectory() const
{
    fs::path dir = getPostProcessingDirectory() / "0" / "residuals";
    fs::create_directories(dir);
    return dir;
}

fs::path controls::getAdaptiveTimesteppingDirectory() const
{
    fs::path dir = getWorkingDirectory() / "adaptiveTimeStepping";
    fs::create_directories(dir);
    return dir;
}

} // namespace accel
