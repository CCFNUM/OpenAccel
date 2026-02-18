// File       : linearSolver.h
// Created    : Mon Jan 22 2024 13:02:48 (+0100)
// Author     : Fabian Wermelinger
// Description: Base class for linear system solvers
// Copyright (c) 2024 CCFNUM, Lucerne University of Applied Sciences and Arts.
// SPDX-License-Identifier: BSD-3-Clause

#ifndef LINEARSOLVER_H
#define LINEARSOLVER_H

#include "CRSNodeGraph.h"
#include "linearSolverContext.h"

#include <cassert>
#include <cmath>
#include <cstdint>
#include <limits>
#include <memory>
#include <string>
#include <vector>

namespace linearSolver
{

template <size_t N>
class Base
{
public:
    using Context = ::linearSolver::Context<N>;
    using ControlData = typename Context::ControlData;
    using Array = typename Context::Array;
    using Coefficients = typename Context::Coefficients;
    using Matrix = typename Context::Matrix;
    using Vector = typename Context::Vector;
    using Index = typename Context::Index;
    using IndexVector = typename Context::IndexVector;

    static constexpr Index BLOCKSIZE = N;

    Base(const int id, const bool precond = false)
        : ctx_(nullptr), id_(id), is_preconditioner_(precond)
    {
    }

    virtual ~Base()
    {
    }

    // disallow
    Base() = delete;
    Base(const Base& c) = delete;
    Base(Base&& c) = delete;
    Base& operator=(const Base& rhs) = delete;
    Base& operator=(Base&& rhs) = delete;

    // public API
    virtual int solve() = 0;

    virtual void setup()
    {
    }

    virtual void report()
    {
        assert(ctx_);
        if (this->is_preconditioner_)
        {
            return; // only report for front-end solver instances
        }

        if (ctx_->verbose() > 1)
        {
            MemoryFootprint data, connectivity;
            const auto& coeffs = ctx_->getCoefficients();
            coeffs.getMemoryFootprint(data, connectivity);
            const bool reduced =
                coeffs.getGraph()->getLayout() & GraphLayout::Stencil__Reduced;

            // clang-format off
            auto& cout = ctx_->cout();
            cout << "Context: "
                      << static_cast<const void*>(ctx_.get()) << '\n';
            cout << "\tMatrix all rows:  " << BLOCKSIZE * coeffs.nGlobalCoefficients() << '\n';
            cout << "\tMatrix all nnz:   " << coeffs.nnzGlobal() << '\n';
            cout << "\tMatrix blocksize: " << BLOCKSIZE << '\n';
            cout << "\tMatrix fp data:   " << data.sum_byte / MemoryFootprint::BYTE_DIVIDE << " MB\n";
            cout << "\tMatrix int data:  " << connectivity.sum_byte / MemoryFootprint::BYTE_DIVIDE << " MB\n";
            cout << "\tReduced stencil:  " << std::string(reduced ? "true" : "false") << '\n';
            cout << "\tCoefficients:     " << static_cast<const void*>(ctx_->coeffs().get()) << '\n';
#ifndef NDEBUG
            // check for numerical zeros and almost zeros
            const typename Matrix::DataType zero = std::abs(0.0);
            const typename Matrix::DataType eps = std::numeric_limits<typename Matrix::DataType>::epsilon();
            unsigned long long local_zeros[2] = {0};
            unsigned long long &is_zero = local_zeros[0];
            unsigned long long &almost_zero = local_zeros[1];
            for (typename Matrix::DataType c : coeffs.valuesRef()) {
                c = std::abs(c);
                if (c == zero) {
                    ++is_zero;
                }
                else if (c <= eps) {
                    ++almost_zero;
                }
            }
            unsigned long long global_zeros[2] = {0};
            MPI_Reduce(local_zeros, global_zeros, 2, MPI_UNSIGNED_LONG_LONG, MPI_SUM, 0, ctx_->getCommunicator());
            const size_t zero_mbyte = global_zeros[0] *
                                      sizeof(typename Matrix::DataType) /
                                      MemoryFootprint::BYTE_DIVIDE;
            const size_t eps_mbyte = global_zeros[1] *
                                     sizeof(typename Matrix::DataType) /
                                     MemoryFootprint::BYTE_DIVIDE;
            cout << "\tN numerical zero: " << global_zeros[0] << " (" << zero_mbyte << " MB)\n";
            cout << "\tN machine eps.:   " << global_zeros[1] << " (" << eps_mbyte << " MB)\n";
#endif /* NDEBUG */
            cout << std::endl;
            // clang-format on
        }
    }

    virtual std::shared_ptr<Context>
    createContext(const std::string& system_name,
                  const CRSNodeGraph* graph) = 0;
    virtual std::shared_ptr<Context>
    createContext(const std::string& system_name,
                  std::shared_ptr<Coefficients> coeffs) = 0;

    int getID() const
    {
        return this->id_;
    }

    bool isPreconditioner() const
    {
        return is_preconditioner_;
    }

    std::shared_ptr<Context> getContext()
    {
        assert(ctx_);
        return ctx_;
    }

    std::shared_ptr<const Context> getContext() const
    {
        assert(ctx_);
        return std::static_pointer_cast<const Context>(ctx_);
    }

    virtual void setContext(std::shared_ptr<Context> ctx)
    {
        assert(ctx);
        ctx_ = ctx;
    }

protected:
    std::shared_ptr<Context> ctx_;

private:
    const int id_;
    const bool is_preconditioner_; // true if this instance is used for
                                   // preconditioning
};

} /* namespace linearSolver */

#endif /* LINEARSOLVER_H */
