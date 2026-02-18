// File       : Profiler.cpp
// Created    : Thu Nov 09 2023 16:29:40 (+0100)
// Author     : Fabian Wermelinger
// Description: Profiling agent implementation
// Copyright 2023 CCFNUM HSLU T&A. All Rights Reserved.

#include "Profiler.h"
#include <cassert>
#include <cstdio>
#include <iostream>
#include <limits>

Profiler::Profiler(const MPI_Comm comm, const std::string& name)
    : comm_(comm), name_(name)
{
    constexpr int count = 4; // number of fields in RankSummary struct
    int block_length[count] = {1, NCHAR_MAX, 1, 1};
    MPI_Aint block_disp[count] = {offsetof(RankSummary, rank),
                                  offsetof(RankSummary, name),
                                  offsetof(RankSummary, total_time),
                                  offsetof(RankSummary, n_samples)};
    MPI_Datatype block_type[count] = {
        MPI_INT, MPI_CHAR, MPI_DOUBLE, MPI_UNSIGNED_LONG};
    MPI_Datatype TempType;
    MPI_Type_create_struct(
        count, block_length, block_disp, block_type, &TempType);

    MPI_Aint lb, extent;
    MPI_Type_get_extent(TempType, &lb, &extent);
    MPI_Type_create_resized(TempType, lb, extent, &MPIRankSummary);
    MPI_Type_commit(&MPIRankSummary);
}

Profiler::~Profiler()
{
    clear();
    MPI_Type_free(&MPIRankSummary);
}

std::vector<Profiler::AgentSummary>
Profiler::getAgentSummaries(const double thresh_time) const
{
    assert(comm_ != MPI_COMM_NULL);

    int size;
    MPI_Comm_size(comm_, &size);
    if (size > 1)
    {
        return consolidateAgents_(thresh_time);
    }
    else
    {
        std::vector<AgentSummary> summary;
        summary.reserve(agents_.size());
        for (auto it : agents_)
        {
            const Agent& agent = *it.second;
            assert(!agent.isActive());
            if (agent.getTotalActiveTime() >= thresh_time)
            {
                summary.push_back(AgentSummary(it.first,
                                               agent.getSampleCount(),
                                               agent.getTotalActiveTime()));
            }
        }
        return summary;
    }
}

void Profiler::printReport(const std::string title,
                           const double thresh_time,
                           const double /* call_overhead */)
{
    assert(comm_ != MPI_COMM_NULL);

    const std::vector<AgentSummary> summary = getAgentSummaries(thresh_time);
    double total_time = 0;
    // double total_time_corrected = 0;
    for (auto& section : summary)
    {
        total_time += section.time_avg;
        // total_time_corrected += section.time_avg - section.avg_samples *
        // call_overhead;
    }

    int rank;
    MPI_Comm_rank(comm_, &rank);
    if (0 == rank)
    {
        if (summary.size() > 0)
        {
            printf("%s", title.c_str());
            printf(" [%-64s]:   perc  avg_time/sample[s]  avg_time[s]  "
                   "min_time[s]  max_time[s]\n",
                   "Name of profiled section");
        }
        for (auto& section : summary)
        {
            printf(" [%-64s]: %5.1f%%           %.3e    %.3e    %.3e    %.3e  "
                   "(%.2f samples, relevant on %d/%d ranks)\n",
                   section.name.c_str(),
                   100.0 * section.time_avg / total_time,
                   section.getSampleAverageTime(),
                   section.time_avg,
                   section.time_min,
                   section.time_max,
                   section.avg_samples,
                   section.n_relevant_ranks,
                   section.n_total_ranks);
        }
        if (summary.size() > 0)
        {
            printf("\n");
        }
    }
}

std::vector<Profiler::AgentSummary>
Profiler::consolidateAgents_(const double thresh_time) const
{
    assert(comm_ != MPI_COMM_NULL);

    int rank, size;
    MPI_Comm_rank(comm_, &rank);
    MPI_Comm_size(comm_, &size);

    // 1. Prepare send buffer
    std::vector<RankSummary> rank_summary;
    rank_summary.reserve(agents_.size());
    for (auto it : agents_)
    {
        if (it.first.size() > NCHAR_MAX)
        {
            std::cerr << "WARNING: name of profile agent `" << it.first
                      << "` exceeds " << NCHAR_MAX
                      << " characters! This may cause ambiguity in agent "
                         "statistics --> consider shortening the name\n";
        }
        const size_t lmax =
            (NCHAR_MAX < it.first.size()) ? NCHAR_MAX : it.first.size();
        const Agent& agent = *it.second;
        assert(!agent.isActive());
        RankSummary s;
        s.rank = rank;
        it.first.copy(s.name, lmax);
        s.name[lmax] = '\0'; // terminating NUL
        s.total_time = agent.getTotalActiveTime();
        s.n_samples = agent.getSampleCount();
        rank_summary.push_back(s);
    }

    // 2. Get receive counts on root rank (depending on how profiler is used
    // on different ranks, these may not be uniform on all ranks).
    std::vector<int> recv_count(size, 0);
    int rcount = static_cast<int>(rank_summary.size());
    MPI_Gather(&rcount, 1, MPI_INT, recv_count.data(), 1, MPI_INT, 0, comm_);

    // 3. Gather to root
    int total_recv = 0;
    for (int i = 0; i < size; i++)
    {
        total_recv += recv_count[i];
    }
    std::vector<int> displs(size, 0);
    for (int i = 1; i < size; i++)
    {
        displs[i] = displs[i - 1] + recv_count[i - 1];
    }
    std::vector<RankSummary> recv_buf(total_recv);
    MPI_Gatherv(rank_summary.data(),
                rank_summary.size(),
                MPIRankSummary,
                recv_buf.data(),
                recv_count.data(),
                displs.data(),
                MPIRankSummary,
                0,
                comm_);

    // 4. Post-process data
    std::map<std::string, std::vector<const RankSummary*>> all_agents;
    for (const auto& it : recv_buf)
    {
        all_agents[std::string(it.name)].push_back(&it);
    }

    std::vector<AgentSummary> summary;
    for (const auto& it : all_agents)
    {
        const auto& agents_on_ranks = it.second;
        assert(agents_on_ranks.size() != 0);

        double min = std::numeric_limits<double>::max();
        double max = std::numeric_limits<double>::min();
        double avg = 0;
        double avg_samples = 0;
        const int relevant_ranks = static_cast<int>(agents_on_ranks.size());
        for (auto agent : agents_on_ranks)
        {
            const double t_agent = agent->total_time;
            min = (t_agent < min) ? t_agent : min;
            max = (t_agent > max) ? t_agent : max;
            avg += t_agent;
            avg_samples += agent->n_samples;
        }
        avg /= relevant_ranks;
        avg_samples /= relevant_ranks;

        if (max >= thresh_time)
        {
            summary.push_back(AgentSummary(
                it.first, avg_samples, avg, min, max, relevant_ranks, size));
        }
    }
    return summary;
}
