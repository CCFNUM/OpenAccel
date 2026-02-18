// File       : TrilinosSolver.h
// Created    : Wed Mar 05 2025 11:04:12 (+0100)
// Author     : Mhamad Mahdi Alloush
// Description: Trilinos (AztecOO) linear solver wrapper
// Copyright (c) 2025 CCFNUM, Lucerne University of Applied Sciences and Arts.
// SPDX-License-Identifier: BSD-3-Clause

#ifndef TRILINOSSOLVER_H
#define TRILINOSSOLVER_H

#include "ContextTrilinos.h"
#include "linearSolver.h"
#include <cassert>
#include <yaml-cpp/yaml.h>

namespace linearSolver
{

template <size_t N>
class Trilinos : public ::linearSolver::Base<N>
{
public:
    using Base = ::linearSolver::Base<N>;
    using Base::BLOCKSIZE;
    using Coefficients = typename Base::Coefficients;
    using Context = typename Base::Context;

protected:
    using Base::ctx_;

public:
    Trilinos(const YAML::Node& node, const bool precond = false)
        : Base(ID::Trilinos, precond), yaml_conf_(node), ctx_trilinos_(nullptr)
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
    ContextTrilinos<N>* ctx_trilinos_;
};

template <size_t N>
std::shared_ptr<typename Trilinos<N>::Context>
Trilinos<N>::createContext(const std::string& system_name,
                           const CRSNodeGraph* graph)
{
#ifdef HAS_TRILINOS
    auto ctx = std::make_shared<ContextTrilinos<N>>(
        "Trilinos", system_name, graph, &yaml_conf_);
    ctx_ = ctx;
    ctx_trilinos_ = ctx.get();
#endif /* HAS_TRILINOS */
    return ctx_;
}

template <size_t N>
std::shared_ptr<typename Trilinos<N>::Context>
Trilinos<N>::createContext(const std::string& system_name,
                           std::shared_ptr<Coefficients> coeffs)
{
#ifdef HAS_TRILINOS
    auto ctx = std::make_shared<ContextTrilinos<N>>(
        "Trilinos", system_name, coeffs, &yaml_conf_);
    ctx_ = ctx;
    ctx_trilinos_ = ctx.get();
#endif /* HAS_TRILINOS */
    return ctx_;
}

template <size_t N>
int Trilinos<N>::solve()
{
#ifdef HAS_TRILINOS
    assert(ctx_);
    assert(ctx_trilinos_);

    if (ctx_.get() != ctx_trilinos_)
    {
        ctx_trilinos_ = *ctx_;
        assert(ctx_trilinos_);
    }
    assert(ctx_->coeffs() == ctx_trilinos_->coeffs());

    ctx_->solvePrologue(this->getID(), this->isPreconditioner());
    ctx_->getControlData().n_iterations = ctx_trilinos_->solve();
    ctx_->solveEpilogue(this->getID(), this->isPreconditioner());
    ++(*ctx_);
#endif /* HAS_TRILINOS */
    return 0;
}

} /* namespace linearSolver */

#endif /* TRILINOSSOLVER_H */
