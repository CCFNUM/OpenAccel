// File       : Profiler.h
// Created    : Thu Nov 09 2023 16:29:40 (+0100)
// Author     : Fabian Wermelinger
// Description: Profiling agent
// Copyright (c) 2025 CCFNUM, Lucerne University of Applied Sciences and Arts.
// SPDX-License-Identifier: BSD-3-Clause

#ifndef PROFILER_H
#define PROFILER_H

#include <cassert>
#include <chrono>
#include <map>
#include <mpi.h>
#include <stack>
#include <string>
#include <vector>

/**
 * @brief Runtime profiler
 *
 * @rst
 * Used to collect runtime samples for a code section that is enclosed
 * by the ``push()`` and ``pop()`` methods.  Used for profiling.
 * @endrst
 * */
class Profiler
{
public:
    class Timer
    {
        using Clock = std::chrono::high_resolution_clock;

    public:
        /** @brief Default constructor
         *
         * Starts the timer */
        Timer() : start_(clock_.now())
        {
        }

        /**
         * @brief Restart the timer
         */
        void start()
        {
            start_ = clock_.now();
        }

        /**
         * @brief Get the currently elapsed seconds
         * @return Elapsed seconds since construction or start
         */
        double stop() const
        {
            return std::chrono::duration<double>(clock_.now() - start_).count();
        }

    private:
        Clock clock_;
        Clock::time_point start_;
    };

    class Agent
    {
    public:
        Agent() : n_samples_(0), total_time_(0), active_(false)
        {
        }

        void start()
        {
            assert(!active_);
            timer_.start();
            active_ = true;
        }

        void stop()
        {
            assert(active_);
            total_time_ += timer_.stop();
            ++n_samples_;
            active_ = false;
        }

        unsigned long int getSampleCount() const
        {
            return n_samples_;
        }

        double getTotalActiveTime() const
        {
            return total_time_;
        }

        double getAverageActiveTime() const
        {
            return total_time_ / n_samples_;
        }

        bool isActive() const
        {
            return active_;
        }

    private:
        unsigned long int n_samples_; // number of samples the agent collected
        double total_time_;           // total time the agent was active
        bool active_;                 // agent is active
        Timer timer_;                 // high-resolution clock
    };

    struct AgentSummary
    {
        std::string name;   // name of agent
        double avg_samples; // average sample count (natural number for single
                            // rank)
        double time_avg;    // total time active (average for multiple ranks)
        double time_min;    // total time active min (on some rank)
        double time_max;    // total time active max (on some rank)
        int n_relevant_ranks;
        int n_total_ranks;

        AgentSummary(const std::string& n,
                     const double as,
                     const double t,
                     const double t_min = -1,
                     const double t_max = -1,
                     const int n_relevant = 1,
                     const int n_total = 1)
            : name(n), avg_samples(as), time_avg(t), time_min(t_min),
              time_max(t_max), n_relevant_ranks(n_relevant),
              n_total_ranks(n_total)
        {
            if (time_min < 0)
            {
                time_min = time_avg;
            }
            if (time_max < 0)
            {
                time_max = time_avg;
            }
        }

        double getSampleAverageTime() const
        {
            return time_avg / avg_samples;
        }
    };

    /**
     * @brief Main constructor
     *
     * @param name Name of the profiler
     */
    Profiler(const MPI_Comm comm, const std::string& name = "Default");

    virtual ~Profiler();

    Profiler(const Profiler& c) = delete;
    Profiler& operator=(const Profiler& c) = delete;

    /**
     * @brief Activate a profiling agent
     *
     * @param name Name of the agent
     *
     * If another agent is already active the agent is stopped and resumed once
     * this active agent stops.
     */
    void push(const std::string& name)
    {
        if (stopped_agents_.size() > 0)
        {
            getAgent(stopped_agents_.top()).stop();
        }
        stopped_agents_.push(name);
        getAgent(name).start();
    }

    /** @brief Deactivate the currently active profiling agent */
    void pop()
    {
        getAgent(stopped_agents_.top()).stop();
        stopped_agents_.pop();
        if (stopped_agents_.size() > 0)
        {
            getAgent(stopped_agents_.top()).start();
        }
    }

    /**
     * @brief Get profiling agent by name
     *
     * @param name Name of the agent
     *
     * @return Reference to agent (creates new agent if it does not exist)
     */
    Agent& getAgent(const std::string& name)
    {
        auto it = agents_.find(name);
        if (it != agents_.end())
        {
            return *it->second;
        }
        Agent* new_agent = new Agent;
        agents_[name] = new_agent;
        return *new_agent;
    }

    /**
     * @brief Clear all agents
     */
    void clear()
    {
        for (auto it : agents_)
        {
            delete it.second;
        }
        agents_.clear();
    }

    /**
     * @brief Get list of agent summaries
     *
     * @param thresh_time Threshold time for agents to be considered
     *
     * @return Vector of agent summaries
     */
    std::vector<AgentSummary>
    getAgentSummaries(const double thresh_time = 1.0e-4) const;

    /**
     * @brief Print agent profiling report to standard output
     *
     * @param title Title for profiling report
     * @param thresh_time Threshold time for agents to be considered
     * @param call_overhead Overhead associated with calling push/pop pair
     */
    void printReport(const std::string title = "",
                     const double thresh_time = 1.0e-4,
                     const double call_overhead = 1.0e-6);

    /**
     * @brief Set MPI communicator
     *
     * @param comm Communicator handle
     */
    inline void setComm(const MPI_Comm comm)
    {
        comm_ = comm;
    }

private:
    MPI_Comm comm_;
    const std::string name_;
    std::stack<std::string> stopped_agents_;
    std::map<std::string, Agent*> agents_;

    static constexpr size_t NCHAR_MAX = 64;

    struct RankSummary
    {
        int rank;
        char name[NCHAR_MAX];
        double total_time;
        unsigned long int n_samples;
    };

    MPI_Datatype MPIRankSummary;

    /**
     * @brief Consolidate profile agents on multiple ranks
     *
     * @param thresh_time Threshold time for agents to be considered
     *
     * @return Vector of consolidated agent summaries
     */
    std::vector<AgentSummary>
    consolidateAgents_(const double thresh_time) const;
};

#endif /* PROFILER_H */
