// File       : blockMatrixOperators.h
// Created    : Wed May 01 2024 09:39:02 (+0200)
// Author     : Fabian Wermelinger
// Description: block matrix operators
// Copyright (c) 2024 CCFNUM, Lucerne University of Applied Sciences and Arts.
// SPDX-License-Identifier: BSD-3-Clause

#ifndef BLOCKMATRIXOPERATORS_H
#define BLOCKMATRIXOPERATORS_H

#include <stdexcept>

namespace linearSolver
{
namespace BlockMatrix
{

// vector-vector (BLAS level 1)
template <int BLOCKSIZE, typename TReal>
constexpr TReal innerProduct(const TReal* x, const TReal* y)
{
    TReal dotp = static_cast<TReal>(0.0);
    for (int i = 0; i < BLOCKSIZE; i++)
    {
        dotp += x[i] * y[i];
    }
    return dotp;
}

template <int BLOCKSIZE, typename TReal>
constexpr void vectorAssign(const TReal* x, TReal* y)
{
    for (int i = 0; i < BLOCKSIZE; i++)
    {
        y[i] = x[i];
    }
}

template <int BLOCKSIZE, typename TReal>
constexpr void vectorAddInplace(const TReal* x, TReal* y)
{
    for (int i = 0; i < BLOCKSIZE; i++)
    {
        y[i] += x[i];
    }
}

template <int BLOCKSIZE, typename TReal>
constexpr void
vectorAddInplaceScaled(const TReal* x, TReal* y, const TReal scale)
{
    for (int i = 0; i < BLOCKSIZE; i++)
    {
        y[i] += scale * x[i];
    }
}

template <int BLOCKSIZE, typename TReal>
constexpr void vectorSubInplace(const TReal* x, TReal* y)
{
    for (int i = 0; i < BLOCKSIZE; i++)
    {
        y[i] -= x[i];
    }
}

template <int BLOCKSIZE, typename TReal>
constexpr void
vectorSubInplaceScaled(const TReal* x, TReal* y, const TReal scale)
{
    for (int i = 0; i < BLOCKSIZE; i++)
    {
        y[i] -= scale * x[i];
    }
}

// matrix-vector (BLAS level 2)
template <int BLOCKSIZE, typename TReal>
constexpr void matrixVectorInplace(const TReal* A, TReal* x)
{
    TReal dotp[BLOCKSIZE] = {0};
    for (int i = 0; i < BLOCKSIZE; i++)
    {
        for (int j = 0; j < BLOCKSIZE; j++)
        {
            dotp[i] += A[BLOCKSIZE * i + j] * x[j];
        }
    }
    for (int i = 0; i < BLOCKSIZE; i++)
    {
        x[i] = dotp[i];
    }
}

template <int BLOCKSIZE, typename TReal>
constexpr void
matrixVector(const TReal* A, const TReal* x, TReal* __restrict__ y)
{
    for (int i = 0; i < BLOCKSIZE; i++)
    {
        y[i] = static_cast<TReal>(0.0);
        for (int j = 0; j < BLOCKSIZE; j++)
        {
            y[i] += A[BLOCKSIZE * i + j] * x[j];
        }
    }
}

template <int BLOCKSIZE, typename TReal>
constexpr void
matrixVectorAdd(const TReal* A, const TReal* x, TReal* __restrict__ y)
{
    for (int i = 0; i < BLOCKSIZE; i++)
    {
        for (int j = 0; j < BLOCKSIZE; j++)
        {
            y[i] += A[BLOCKSIZE * i + j] * x[j];
        }
    }
}

template <int BLOCKSIZE, typename TReal>
constexpr void
matrixVectorSub(const TReal* A, const TReal* x, TReal* __restrict__ y)
{
    for (int i = 0; i < BLOCKSIZE; i++)
    {
        for (int j = 0; j < BLOCKSIZE; j++)
        {
            y[i] -= A[BLOCKSIZE * i + j] * x[j];
        }
    }
}

// matrix-matrix (BLAS level 3)
template <int BLOCKSIZE, typename TReal>
constexpr void matrixMatrixInplace(const TReal* A, const TReal* B, TReal* dst)
{
    TReal C[BLOCKSIZE * BLOCKSIZE] = {0};
    for (int i = 0; i < BLOCKSIZE; i++)
    {
        for (int k = 0; k < BLOCKSIZE; k++)
        {
            for (int j = 0; j < BLOCKSIZE; j++)
            {
                C[BLOCKSIZE * i + j] +=
                    A[BLOCKSIZE * i + k] * B[BLOCKSIZE * k + j];
            }
        }
    }
    for (int i = 0; i < BLOCKSIZE * BLOCKSIZE; i++)
    {
        dst[i] = C[i];
    }
}

template <int BLOCKSIZE, typename TReal>
constexpr void
matrixMatrix(const TReal* A, const TReal* B, TReal* __restrict__ C)
{
    for (int i = 0; i < BLOCKSIZE; i++)
    {
        for (int j = 0; j < BLOCKSIZE; j++)
        {
            C[BLOCKSIZE * i + j] = static_cast<TReal>(0.0);
        }
        for (int k = 0; k < BLOCKSIZE; k++)
        {
            for (int j = 0; j < BLOCKSIZE; j++)
            {
                C[BLOCKSIZE * i + j] +=
                    A[BLOCKSIZE * i + k] * B[BLOCKSIZE * k + j];
            }
        }
    }
}

template <int BLOCKSIZE, typename TReal>
constexpr void matrixAssign(const TReal* A, TReal* B)
{
    for (int i = 0; i < BLOCKSIZE * BLOCKSIZE; i++)
    {
        B[i] = A[i];
    }
}

template <int BLOCKSIZE, typename TReal>
constexpr void matrixAddInplace(const TReal* A, TReal* B)
{
    for (int i = 0; i < BLOCKSIZE * BLOCKSIZE; i++)
    {
        B[i] += A[i];
    }
}

template <int BLOCKSIZE, typename TReal>
constexpr void matrixSubInplace(const TReal* A, TReal* B)
{
    for (int i = 0; i < BLOCKSIZE * BLOCKSIZE; i++)
    {
        B[i] -= A[i];
    }
}

template <int BLOCKSIZE, typename TReal>
constexpr void
matrixMatrixAdd(const TReal* A, const TReal* B, TReal* __restrict__ C)
{
    for (int i = 0; i < BLOCKSIZE; i++)
    {
        for (int k = 0; k < BLOCKSIZE; k++)
        {
            for (int j = 0; j < BLOCKSIZE; j++)
            {
                C[BLOCKSIZE * i + j] +=
                    A[BLOCKSIZE * i + k] * B[BLOCKSIZE * k + j];
            }
        }
    }
}

template <int BLOCKSIZE, typename TReal>
constexpr void
matrixMatrixSub(const TReal* A, const TReal* B, TReal* __restrict__ C)
{
    for (int i = 0; i < BLOCKSIZE; i++)
    {
        for (int k = 0; k < BLOCKSIZE; k++)
        {
            for (int j = 0; j < BLOCKSIZE; j++)
            {
                C[BLOCKSIZE * i + j] -=
                    A[BLOCKSIZE * i + k] * B[BLOCKSIZE * k + j];
            }
        }
    }
}

// TODO: [2024-05-01] generic implementations for the
// following:

// determinant
template <int BLOCKSIZE, typename TReal>
inline TReal determinant(const TReal*)
{
    throw std::runtime_error(
        "BlockMatrix: determinant not implemented for generic BLOCKSIZE");
    return static_cast<TReal>(0.0);
}

// matrix inversion
template <int BLOCKSIZE, typename TReal>
inline void invertInplace(TReal*)
{
    throw std::runtime_error(
        "BlockMatrix: inplace invert not implemented for generic BLOCKSIZE");
}

template <int BLOCKSIZE, typename TReal>
inline void invert(const TReal*, TReal*)
{
    throw std::runtime_error(
        "BlockMatrix: invert not implemented for generic BLOCKSIZE");
}

} /* namespace BlockMatrix */
} /* namespace linearSolver */

// implementations/specializations
#include "blockMatrixDeterminant.hpp"
#include "blockMatrixInvert.hpp"

// TODO: [2024-05-01] AVX implementations MV and MM

#endif /* BLOCKMATRIXOPERATORS_H */
