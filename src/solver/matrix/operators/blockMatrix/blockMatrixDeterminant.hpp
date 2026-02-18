// File       : blockMatrixDeterminant.hpp
// Created    : Wed May 01 2024 10:16:27 (+0200)
// Author     : Fabian Wermelinger
// Description: block matrix determinant kernels implementation/specializations
// Copyright (c) 2024 CCFNUM, Lucerne University of Applied Sciences and Arts.
// SPDX-License-Identifier: BSD-3-Clause

// no include guard: must only be include in blockMatrixOperators.h

#include <cassert>

namespace linearSolver
{
namespace BlockMatrix
{

template <>
constexpr float determinant<1, float>(const float* src)
{
    return *src;
}

template <>
constexpr double determinant<1, double>(const double* src)
{
    return *src;
}

#define DET_2x2(T)                                                             \
    do                                                                         \
    {                                                                          \
        const T a00 = src[0];                                                  \
        const T a01 = src[1];                                                  \
        const T a10 = src[2];                                                  \
        const T a11 = src[3];                                                  \
        return a00 * a11 - a01 * a10;                                          \
    } while (0)

template <>
constexpr float determinant<2, float>(const float* src)
{
    DET_2x2(float);
}

template <>
constexpr double determinant<2, double>(const double* src)
{
    DET_2x2(double);
}

#undef DET_2x2

// clang-format off
#define DET_3x3(T)                                                             \
    do                                                                         \
    {                                                                          \
        const T a00 = src[0];                                                  \
        const T a01 = src[1];                                                  \
        const T a02 = src[2];                                                  \
        const T a10 = src[3];                                                  \
        const T a11 = src[4];                                                  \
        const T a12 = src[5];                                                  \
        const T a20 = src[6];                                                  \
        const T a21 = src[7];                                                  \
        const T a22 = src[8];                                                  \
        return a00 * a11 * a22 + a01 * a12 * a20 + a02 * a10 * a21 -           \
               a00 * a12 * a21 - a01 * a10 * a22 - a02 * a11 * a20;            \
    } while (0)

// clang-format on

template <>
constexpr float determinant<3, float>(const float* src)
{
    DET_3x3(float);
}

template <>
constexpr double determinant<3, double>(const double* src)
{
    DET_3x3(double);
}

#undef DET_3x3

// clang-format off
#define DET_4x4(T)                                                             \
    do                                                                         \
    {                                                                          \
        const T a00 = src[0];                                                  \
        const T a01 = src[1];                                                  \
        const T a02 = src[2];                                                  \
        const T a03 = src[3];                                                  \
        const T a10 = src[4];                                                  \
        const T a11 = src[5];                                                  \
        const T a12 = src[6];                                                  \
        const T a13 = src[7];                                                  \
        const T a20 = src[8];                                                  \
        const T a21 = src[9];                                                  \
        const T a22 = src[10];                                                 \
        const T a23 = src[11];                                                 \
        const T a30 = src[12];                                                 \
        const T a31 = src[13];                                                 \
        const T a32 = src[14];                                                 \
        const T a33 = src[15];                                                 \
        return a00 * a11 * a22 * a33 - a00 * a11 * a23 * a32 -                 \
               a00 * a12 * a21 * a33 + a00 * a12 * a23 * a31 +                 \
               a00 * a13 * a21 * a32 - a00 * a13 * a22 * a31 -                 \
               a01 * a10 * a22 * a33 + a01 * a10 * a23 * a32 +                 \
               a01 * a12 * a20 * a33 - a01 * a12 * a23 * a30 -                 \
               a01 * a13 * a20 * a32 + a01 * a13 * a22 * a30 +                 \
               a02 * a10 * a21 * a33 - a02 * a10 * a23 * a31 -                 \
               a02 * a11 * a20 * a33 + a02 * a11 * a23 * a30 +                 \
               a02 * a13 * a20 * a31 - a02 * a13 * a21 * a30 -                 \
               a03 * a10 * a21 * a32 + a03 * a10 * a22 * a31 +                 \
               a03 * a11 * a20 * a32 - a03 * a11 * a22 * a30 -                 \
               a03 * a12 * a20 * a31 + a03 * a12 * a21 * a30;                  \
    } while (0)

// clang-format on

template <>
constexpr float determinant<4, float>(const float* src)
{
    DET_4x4(float);
}

template <>
constexpr double determinant<4, double>(const double* src)
{
    DET_4x4(double);
}

} /* namespace BlockMatrix */
} /* namespace linearSolver */
