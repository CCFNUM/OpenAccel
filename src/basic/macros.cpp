// File       : macros.cpp
// Created    : Tue Apr 30 2024 12:55:24 (+0100)
// Author     : Mhamad Mahdi Alloush
// Description:
// Copyright (c) 2024 CCFNUM, Lucerne University of Applied Sciences and Arts.
// SPDX-License-Identifier: BSD-3-Clause

// code
#include "macros.h"

namespace accel
{

void tolower(std::string& s)
{
    std::transform(s.begin(), s.end(), s.begin(), [](unsigned char c) {
        return std::tolower(c);
    });
}

void errorMsg(std::string msg)
{
    throw std::runtime_error(msg);
}

void warningMsg(std::string msg)
{
    std::cout << "Warning:" << std::endl;
    std::cout << msg << std::endl;
}

void infoMsg(std::string msg)
{
    std::cout << std::endl;
    std::cout << msg << std::endl;
}

namespace ops
{

// Ghosting

void populateGhostCommProcs(const stk::mesh::BulkData& bulk_data,
                            stk::mesh::Ghosting& ghosting,
                            std::vector<int>& ghostCommProcs)
{
    ghostCommProcs.clear();

    std::vector<stk::mesh::EntityProc> sendList;
    ghosting.send_list(sendList);

    for (const stk::mesh::EntityProc& entProc : sendList)
    {
        stk::util::insert_keep_sorted_and_unique(entProc.second,
                                                 ghostCommProcs);
    }

    std::vector<stk::mesh::EntityKey> recvList;
    ghosting.receive_list(recvList);

    for (const stk::mesh::EntityKey& key : recvList)
    {
        stk::mesh::Entity entity = bulk_data.get_entity(key);
        stk::util::insert_keep_sorted_and_unique(
            bulk_data.parallel_owner_rank(entity), ghostCommProcs);
    }
}

// Memory Diagnostics

void printMemoryDiag(const stk::mesh::BulkData& bulk)
{
    const double factor = 1024.0 * 1024 * 1024;
    size_t curr_max, curr_min, curr_avg;
    stk::get_current_memory_usage_across_processors(
        bulk.parallel(), curr_max, curr_min, curr_avg);

    if (bulk.parallel_rank() == 0)
    {
        std::cout << "Memory usage (GB): Avg. = " << 1.0 * curr_avg / factor
                  << "; Min. = " << 1.0 * curr_min / factor
                  << "; Max. = " << 1.0 * curr_max / factor << std::endl;
    }
}

void printHwmMemoryDiag(const stk::mesh::BulkData& bulk)
{
    std::ios::fmtflags cflags(std::cout.flags());
    const double factor = 1024.0 * 1024 * 1024;
    size_t curr_max, curr_min, curr_avg;
    stk::get_memory_high_water_mark_across_processors(
        bulk.parallel(), curr_max, curr_min, curr_avg);

    if (bulk.parallel_rank() == 0)
    {
        std::cout << "Memory high-water mark: Avg. = " << std::setw(6)
                  << std::fixed << std::setprecision(4) << curr_avg / factor
                  << " GB; Min. = " << std::setw(6) << std::fixed
                  << std::setprecision(4) << curr_min / factor
                  << " GB; Max. = " << std::setw(6) << std::fixed
                  << std::setprecision(4) << curr_max / factor << " GB"
                  << std::endl;
    }
    std::cout.flags(cflags);
}

} // namespace ops

} // namespace accel
