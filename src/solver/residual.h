// File       : residual.h
// Created    : Sun Oct 13 2024 18:48:31 (+0200)
// Author     : Fabian Wermelinger
// Description: Linear solver residual functions
// Copyright (c) 2024 CCFNUM, Lucerne University of Applied Sciences and Arts.
// SPDX-License-Identifier: BSD-3-Clause

#ifndef RESIDUAL_H
#define RESIDUAL_H

#include "CRSResidual.h"
#include "coefficients.h"
#include "controlData.h"
#include <cassert>
#include <cmath>
#include <mpi.h>

namespace linearSolver
{
namespace residual
{

template <size_t N>
using Array = typename ControlData<N>::Array;

template <size_t N>
constexpr void compute(coefficients<N>& coeff)
{
    using Matrix = typename coefficients<N>::Matrix;
    using Vector = typename coefficients<N>::Vector;

    const Matrix& A = coeff.getAMatrix();
    const Vector& b = coeff.getBVector();
    Vector& x = coeff.getXVector();
    Vector& r = coeff.getRVector();
    CRSResidual(A, b, x, r);
}

template <size_t N>
constexpr Array<N>
RMSReduceAndScale(const MPI_Comm comm, const Array<N>& rms, const double scale)
{
    Array<N> global_rms = {0};
    MPI_Allreduce(
        rms.data(), global_rms.data(), rms.size(), MPI_DOUBLE, MPI_SUM, comm);
    for (size_t k = 0; k < N; k++)
    {
        global_rms[k] = std::sqrt(global_rms[k] * scale);
    }

    return global_rms;
}

template <size_t N, typename DataType, typename Index>
constexpr Array<N> RMSError(const MPI_Comm comm,
                            const Index n_local,
                            const Index n_global,
                            const DataType* R)
{
    constexpr Index BLOCKSIZE = N;

    Array<N> rms = {0}; // root-mean-square of residual vector
    for (Index i = 0; i < n_local; i++)
    {
        for (Index k = 0; k < BLOCKSIZE; k++)
        {
            rms[k] += R[BLOCKSIZE * i + k] * R[BLOCKSIZE * i + k];
        }
    }

    return RMSReduceAndScale<BLOCKSIZE>(comm, rms, 1.0 / n_global);
}

template <size_t N>
constexpr Array<N> RMSError(const coefficients<N>& coeff)
{
    const auto& R = coeff.getRVector();
    assert(coeff.getAMatrix().nRows() == coeff.nCoefficients());
    assert(R.size() <= N * coeff.nCoefficients());
    return RMSError<N>(coeff.getCommunicator(),
                       coeff.nCoefficients(),
                       coeff.nGlobalCoefficients(),
                       R.data());
}

template <size_t N>
constexpr Array<N> RMSErrorDiagNormalized(const coefficients<N>& coeff)
{
    using Matrix = typename coefficients<N>::Matrix;
    using Vector = typename coefficients<N>::Vector;
    using Index = typename coefficients<N>::Index;
    using DataType = typename Matrix::DataType;
    constexpr Index BLOCKSIZE = N;

    const Matrix& A = coeff.getAMatrix();
    const Vector& R = coeff.getRVector();

    assert(A.nRows() == coeff.nCoefficients());
    assert(static_cast<Index>(R.size()) <= BLOCKSIZE * coeff.nCoefficients());
    const Index nrows = coeff.nCoefficients(); // rank local rows
    Array<N> rms = {0}; // root-mean-square of residual vector
    for (Index i = 0; i < nrows; i++)
    {
        const DataType* Ap_blk = A.diag(i);
        for (Index k = 0; k < BLOCKSIZE; k++)
        {
            assert(Ap_blk[BLOCKSIZE * k + k] > static_cast<DataType>(0.0));
            const DataType r = R[BLOCKSIZE * i + k] / Ap_blk[BLOCKSIZE * k + k];
            rms[k] += r * r;
        }
    }

    return RMSReduceAndScale<BLOCKSIZE>(
        coeff.getCommunicator(), rms, 1.0 / coeff.nGlobalCoefficients());
}

} /* namespace residual */
} /* namespace linearSolver */

#endif /* RESIDUAL_H */
