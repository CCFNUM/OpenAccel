// File       : blockMatrixInvert.hpp
// Created    : Wed May 01 2024 10:16:27 (+0200)
// Author     : Fabian Wermelinger
// Description: block matrix inversion kernels implementation/specializations
// Copyright (c) 2024 CCFNUM, Lucerne University of Applied Sciences and Arts.
// SPDX-License-Identifier: BSD-3-Clause

// no include guard: must only be include in blockMatrixOperators.h

#include <cassert>
#include <cmath>

namespace linearSolver
{
namespace BlockMatrix
{

// BLOCKSIZE = 1
template <>
constexpr void invertInplace<1, float>(float* src)
{
    assert(std::abs(*src) > static_cast<float>(0.0));
    *src = static_cast<float>(1.0) / *src;
}

template <>
constexpr void invertInplace<1, double>(double* src)
{
    assert(std::abs(*src) > static_cast<double>(0.0));
    *src = static_cast<double>(1.0) / *src;
}

template <>
constexpr void invert<1, float>(const float* src, float* dst)
{
    assert(std::abs(*src) > static_cast<float>(0.0));
    *dst = static_cast<float>(1.0) / *src;
}

template <>
constexpr void invert<1, double>(const double* src, double* dst)
{
    assert(std::abs(*src) > static_cast<double>(0.0));
    *dst = static_cast<double>(1.0) / *src;
}

// BLOCKSIZE = 2
#define INVERT_2x2(src, dst, T)                                                \
    do                                                                         \
    {                                                                          \
        const T a00 = src[0];                                                  \
        const T a01 = src[1];                                                  \
        const T a10 = src[2];                                                  \
        const T a11 = src[3];                                                  \
                                                                               \
        const T det = determinant<2>(src);                                     \
        assert(std::abs(det) > static_cast<T>(0.0));                           \
        const T inv_det = static_cast<T>(1.0) / det;                           \
                                                                               \
        dst[0] = inv_det * a11;                                                \
        dst[1] = -inv_det * a01;                                               \
        dst[2] = -inv_det * a10;                                               \
        dst[3] = inv_det * a00;                                                \
    } while (0)

template <>
constexpr void invertInplace<2, float>(float* src)
{
    INVERT_2x2(src, src, float); // in-place
}

template <>
constexpr void invertInplace<2, double>(double* src)
{
    INVERT_2x2(src, src, double); // in-place
}

template <>
constexpr void invert<2, float>(const float* src, float* dst)
{
    INVERT_2x2(src, dst, float); // out-of-place
}

template <>
constexpr void invert<2, double>(const double* src, double* dst)
{
    INVERT_2x2(src, dst, double); // out-of-place
}

#undef INVERT_2x2

// BLOCKSIZE = 3
#define INVERT_3x3(src, dst, T)                                                \
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
                                                                               \
        const T det = determinant<3>(src);                                     \
        assert(std::abs(det) > static_cast<T>(0.0));                           \
        const T inv_det = static_cast<T>(1.0) / det;                           \
                                                                               \
        dst[0] = inv_det * (a11 * a22 - a12 * a21);                            \
        dst[1] = -inv_det * (a01 * a22 - a02 * a21);                           \
        dst[2] = inv_det * (a01 * a12 - a02 * a11);                            \
        dst[3] = -inv_det * (a10 * a22 - a12 * a20);                           \
        dst[4] = inv_det * (a00 * a22 - a02 * a20);                            \
        dst[5] = -inv_det * (a00 * a12 - a02 * a10);                           \
        dst[6] = inv_det * (a10 * a21 - a11 * a20);                            \
        dst[7] = -inv_det * (a00 * a21 - a01 * a20);                           \
        dst[8] = inv_det * (a00 * a11 - a01 * a10);                            \
    } while (0)

template <>
constexpr void invertInplace<3, float>(float* src)
{
    INVERT_3x3(src, src, float); // in-place
}

template <>
constexpr void invertInplace<3, double>(double* src)
{
    INVERT_3x3(src, src, double); // in-place
}

template <>
constexpr void invert<3, float>(const float* src, float* dst)
{
    INVERT_3x3(src, dst, float); // out-of-place
}

template <>
constexpr void invert<3, double>(const double* src, double* dst)
{
    INVERT_3x3(src, dst, double); // out-of-place
}

#undef INVERT_3x3

// BLOCKSIZE = 4
// clang-format off
#define INVERT_4x4(src, dst, T)                                                                                   \
    do {                                                                                                          \
        const T a00 = src[0];                                                                                     \
        const T a01 = src[1];                                                                                     \
        const T a02 = src[2];                                                                                     \
        const T a03 = src[3];                                                                                     \
        const T a10 = src[4];                                                                                     \
        const T a11 = src[5];                                                                                     \
        const T a12 = src[6];                                                                                     \
        const T a13 = src[7];                                                                                     \
        const T a20 = src[8];                                                                                     \
        const T a21 = src[9];                                                                                     \
        const T a22 = src[10];                                                                                    \
        const T a23 = src[11];                                                                                    \
        const T a30 = src[12];                                                                                    \
        const T a31 = src[13];                                                                                    \
        const T a32 = src[14];                                                                                    \
        const T a33 = src[15];                                                                                    \
                                                                                                                  \
        const T det = determinant<4>(src);                                                                        \
        assert(std::abs(det) > static_cast<T>(0.0));                                                              \
        const T inv_det = static_cast<T>(1.0) / det;                                                              \
                                                                                                                  \
        dst[0]  =  inv_det * (a11*a22*a33 - a11*a23*a32 - a12*a21*a33 + a12*a23*a31 + a13*a21*a32 - a13*a22*a31); \
        dst[1]  = -inv_det * (a01*a22*a33 - a01*a23*a32 - a02*a21*a33 + a02*a23*a31 + a03*a21*a32 - a03*a22*a31); \
        dst[2]  =  inv_det * (a01*a12*a33 - a01*a13*a32 - a02*a11*a33 + a02*a13*a31 + a03*a11*a32 - a03*a12*a31); \
        dst[3]  = -inv_det * (a01*a12*a23 - a01*a13*a22 - a02*a11*a23 + a02*a13*a21 + a03*a11*a22 - a03*a12*a21); \
        dst[4]  = -inv_det * (a10*a22*a33 - a10*a23*a32 - a12*a20*a33 + a12*a23*a30 + a13*a20*a32 - a13*a22*a30); \
        dst[5]  =  inv_det * (a00*a22*a33 - a00*a23*a32 - a02*a20*a33 + a02*a23*a30 + a03*a20*a32 - a03*a22*a30); \
        dst[6]  = -inv_det * (a00*a12*a33 - a00*a13*a32 - a02*a10*a33 + a02*a13*a30 + a03*a10*a32 - a03*a12*a30); \
        dst[7]  =  inv_det * (a00*a12*a23 - a00*a13*a22 - a02*a10*a23 + a02*a13*a20 + a03*a10*a22 - a03*a12*a20); \
        dst[8]  =  inv_det * (a10*a21*a33 - a10*a23*a31 - a11*a20*a33 + a11*a23*a30 + a13*a20*a31 - a13*a21*a30); \
        dst[9]  = -inv_det * (a00*a21*a33 - a00*a23*a31 - a01*a20*a33 + a01*a23*a30 + a03*a20*a31 - a03*a21*a30); \
        dst[10] =  inv_det * (a00*a11*a33 - a00*a13*a31 - a01*a10*a33 + a01*a13*a30 + a03*a10*a31 - a03*a11*a30); \
        dst[11] = -inv_det * (a00*a11*a23 - a00*a13*a21 - a01*a10*a23 + a01*a13*a20 + a03*a10*a21 - a03*a11*a20); \
        dst[12] = -inv_det * (a10*a21*a32 - a10*a22*a31 - a11*a20*a32 + a11*a22*a30 + a12*a20*a31 - a12*a21*a30); \
        dst[13] =  inv_det * (a00*a21*a32 - a00*a22*a31 - a01*a20*a32 + a01*a22*a30 + a02*a20*a31 - a02*a21*a30); \
        dst[14] = -inv_det * (a00*a11*a32 - a00*a12*a31 - a01*a10*a32 + a01*a12*a30 + a02*a10*a31 - a02*a11*a30); \
        dst[15] =  inv_det * (a00*a11*a22 - a00*a12*a21 - a01*a10*a22 + a01*a12*a20 + a02*a10*a21 - a02*a11*a20); \
    } while (0)

// clang-format on

template <>
constexpr void invertInplace<4, float>(float* src)
{
    INVERT_4x4(src, src, float); // in-place
}

template <>
constexpr void invertInplace<4, double>(double* src)
{
    INVERT_4x4(src, src, double); // in-place
}

template <>
constexpr void invert<4, float>(const float* src, float* dst)
{
    INVERT_4x4(src, dst, float); // out-of-place
}

template <>
constexpr void invert<4, double>(const double* src, double* dst)
{
    INVERT_4x4(src, dst, double); // out-of-place
}

#undef INVERT_4x4

} /* namespace BlockMatrix */
} /* namespace linearSolver */
