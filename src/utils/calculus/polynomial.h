// File       : polynomial.h
// Created    : Wed Jun 11 2025 12:55:24 (+0100)
// Author     : Mhamad Mahdi Alloush
// Description: Templated polynomial class with Newton-Raphson solver.
// Copyright (c) 2025 CCFNUM, Lucerne University of Applied Sciences and Arts.
// SPDX-License-Identifier: BSD-3-Clause

#ifndef POLYNOMIAL_H
#define POLYNOMIAL_H

#include "types.h"

namespace accel
{
namespace utils
{

template <size_t p, typename T = scalar>
class polynomial
{
public:
    std::array<T, p + 1> coeffs_;

    // Default constructor (zero coefficients)
    polynomial() noexcept
    {
        coeffs_.fill(static_cast<T>(0));
    }

    // Constructor with coefficient array
    explicit polynomial(const std::array<T, p + 1>& coefficients) noexcept
        : coeffs_(coefficients)
    {
    }

    // Evaluate polynomial using Horner's method
    inline T evaluate(T x) const noexcept
    {
        T result = coeffs_[p];
        for (label i = static_cast<label>(p) - 1; i >= 0; --i)
        {
            result = result * x + coeffs_[i];
        }
        return result;
    }

    // Derivative using Horner's method
    inline T derivative(T x) const noexcept
    {
        if constexpr (p == 0)
        {
            return static_cast<T>(0);
        }

        T result = static_cast<T>(p) * coeffs_[p];
        for (label i = static_cast<label>(p) - 1; i >= 1; --i)
        {
            result = result * x + static_cast<T>(i) * coeffs_[i];
        }
        return result;
    }

    // Convenience operator()
    inline T operator()(T x) const noexcept
    {
        return evaluate(x);
    }

    // Newton-Raphson solver: solve P(x) = y
    T solve(T y,
            T x0 = static_cast<T>(0.0),
            T tol = static_cast<T>(1e-6),
            label maxIter = 100) const
    {
        T x = x0;

        for (label iter = 0; iter < maxIter; ++iter)
        {
            T f = evaluate(x) - y;
            T df = derivative(x);

            if (std::abs(df) < std::numeric_limits<T>::epsilon())
            {
                std::cerr << "Derivative too small, Newton-Raphson may fail.\n";
                break;
            }

            T x_new = x - f / df;

            if (std::abs(x_new - x) < tol)
                return x_new;

            x = x_new;
        }

        std::cerr << "Warning: Newton-Raphson did not converge.\n";
        return x;
    }
};

} // namespace utils
} // namespace accel

#endif // POLYNOMIAL_H
