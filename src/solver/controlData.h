// File       : controlData.h
// Created    : Sun Oct 13 2024 18:21:12 (+0200)
// Author     : Fabian Wermelinger
// Description: Linear solver residual control data object
// Copyright (c) 2024 CCFNUM, Lucerne University of Applied Sciences and Arts.
// SPDX-License-Identifier: BSD-3-Clause

#ifndef CONTROLDATA_H
#define CONTROLDATA_H

#include <array>
#include <ostream>

template <size_t N>
std::ostream& operator<<(std::ostream& os, const std::array<double, N>& v)
{
    os << "[" << std::scientific << v[0];
    if (N > 1)
    {
        for (size_t i = 1; i < N; i++)
        {
            os << ", " << std::scientific << v[i];
        }
    }
    os << "]";
    return os;
}

namespace linearSolver
{

template <size_t N>
struct ControlData
{
    using Array = std::array<double, N>;
    Array scaled_initial_res = {0};
    Array scaled_final_res = {0};
    Array solver_initial_res = {0};
    Array solver_final_res = {0};
    int n_iterations = 0;

    friend std::ostream& operator<<(std::ostream& os, const ControlData& cd)
    {
        os << "Scaled initial residual: " << cd.scaled_initial_res << '\n';
        os << "Scaled final residual:   " << cd.scaled_final_res << '\n';
        os << "Solver initial residual: " << cd.solver_initial_res << '\n';
        os << "Solver final residual:   " << cd.solver_final_res << '\n';
        os << "Number of iterations:    " << cd.n_iterations << '\n';
        return os;
    }
};

} /* namespace linearSolver */

#endif /* CONTROLDATA_H */
