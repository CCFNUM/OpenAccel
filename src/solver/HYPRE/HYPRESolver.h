// File       : HYPRESolver.h
// Created    : Wed Mar 27 2024 09:38:56 (+0100)
// Author     : Fabian Wermelinger
// Description: HYPRE based linear solver wrapper
// Copyright (c) 2024 CCFNUM, Lucerne University of Applied Sciences and Arts.
// SPDX-License-Identifier: BSD-3-Clause

#ifndef HYPRE_H
#define HYPRE_H

#include "ContextHYPRE.h"
#include "linearSolver.h"
#include <cassert>
#include <yaml-cpp/yaml.h>

namespace linearSolver
{

template <size_t N>
class HYPRE : public ::linearSolver::Base<N>
{
public:
    using Base = ::linearSolver::Base<N>;
    using Base::BLOCKSIZE;
    using Coefficients = typename Base::Coefficients;
    using Context = typename Base::Context;
    using Matrix = typename Context::Matrix;
    using Index = typename Context::Index;
    using IndexVector = typename Context::IndexVector;

protected:
    using Base::ctx_;

public:
    HYPRE(const YAML::Node& node, const bool precond = false)
        : Base(ID::HYPRE, precond), yaml_conf_(node), ctx_hypre_(nullptr)
    {
    }

    int solve() override;
    std::shared_ptr<Context> createContext(const std::string& system_name,
                                           const CRSNodeGraph* graph) override;
    std::shared_ptr<Context>
    createContext(const std::string& system_name,
                  std::shared_ptr<Coefficients> coeffs) override;

protected:
    const YAML::Node yaml_conf_;

private:
    ContextHYPRE<N>* ctx_hypre_;
};

template <size_t N>
std::shared_ptr<typename HYPRE<N>::Context>
HYPRE<N>::createContext(const std::string& system_name,
                        const CRSNodeGraph* graph)
{
#ifdef HAS_HYPRE
    auto ctx = std::make_shared<ContextHYPRE<N>>(
        "HYPRE", system_name, graph, &yaml_conf_);
    ctx_ = ctx;
    ctx_hypre_ = ctx.get();
#endif /* HAS_HYPRE */
    return ctx_;
}

template <size_t N>
std::shared_ptr<typename HYPRE<N>::Context>
HYPRE<N>::createContext(const std::string& system_name,
                        std::shared_ptr<Coefficients> coeffs)
{
#ifdef HAS_HYPRE
    auto ctx = std::make_shared<ContextHYPRE<N>>(
        "HYPRE", system_name, coeffs, &yaml_conf_);
    ctx_ = ctx;
    ctx_hypre_ = ctx.get();
#endif /* HAS_HYPRE */
    return ctx_;
}

template <size_t N>
int HYPRE<N>::solve()
{
#ifdef HAS_HYPRE
    assert(ctx_);
    assert(ctx_hypre_);

    if (ctx_.get() != ctx_hypre_)
    {
        ctx_hypre_ = *ctx_;
        assert(ctx_hypre_);
    }
    assert(ctx_->coeffs() == ctx_hypre_->coeffs());

    ctx_->solvePrologue(this->getID(), this->isPreconditioner());

    // Solve the system Ax = b
    ctx_->getControlData().n_iterations = ctx_hypre_->solve();

    ctx_->solveEpilogue(this->getID(), this->isPreconditioner());

    ++(*ctx_); // increment solver call count
#endif         /* HAS_HYPRE */

    return 0;
}

} /* namespace linearSolver */

#endif /* HYPRE_H */
