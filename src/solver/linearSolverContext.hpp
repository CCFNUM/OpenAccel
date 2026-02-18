// File       : linearSolverContext.hpp
// Created    : Fri Sep 27 2024 15:50:24 (+0200)
// Author     : Fabian Wermelinger
// Description: Linear solver context base class (implementation details)
// Copyright (c) 2024 CCFNUM, Lucerne University of Applied Sciences and Arts.
// SPDX-License-Identifier: BSD-3-Clause

// no include guard: must only be included in linearSolverContext.h

#include <iostream>

namespace linearSolver
{

template <size_t N>
Context<N>::Context(const int solver_id,
                    const std::string& family,
                    const std::string& system_name,
                    const CRSNodeGraph* graph,
                    const YAML::Node* node)
    : solver_id_(solver_id), coeffs_(std::make_shared<Coefficients>(graph)),
      profiler_(this->getCommunicator(), "Context [" + family + "]"),
      family_(family), system_name_(system_name), call_count_(0),
      call_count_relative_(0), ctx_precond_(nullptr), cout_(&::std::cout),
      active_stream_(nullptr)
{
    assert(solver_id_ < linearSolver::ID::UNDEFINED);
    setupDefaults_();
    if (node)
    {
        setOptions_(node);
    }

    setupIOStream_();
}

template <size_t N>
Context<N>::Context(const int solver_id,
                    const std::string& family,
                    const std::string& system_name,
                    std::shared_ptr<Coefficients> coeffs,
                    const YAML::Node* node)
    : solver_id_(solver_id), coeffs_(coeffs),
      profiler_(this->getCommunicator(), "Context [" + family + "]"),
      family_(family), system_name_(system_name), call_count_(0),
      call_count_relative_(0), ctx_precond_(nullptr), cout_(&::std::cout),
      active_stream_(nullptr)
{
    assert(solver_id_ < linearSolver::ID::UNDEFINED);
    setupDefaults_();
    if (node)
    {
        setOptions_(node);
    }

    setupIOStream_();
}

// protected constructor
template <size_t N>
Context<N>::Context(const std::string& family, const std::string& system_name)
    : solver_id_(linearSolver::ID::UNDEFINED),
      profiler_(MPI_COMM_NULL, "Context [" + family + "]"), family_(family),
      system_name_(system_name), ctx_precond_(nullptr), cout_(&::std::cout),
      active_stream_(nullptr), call_count_(0), call_count_relative_(0)
{
    setupDefaults_();
}

template <size_t N>
std::string Context<N>::info(std::ostream& os, const char* prefix) const
{
    std::string pfx("");
    if (prefix)
    {
        pfx = std::string(prefix);
    }

    // clang-format off
    os << pfx << "Linear solver context: " << family_ << "\n";
    os << pfx << '\t' << "Address:          " << static_cast<const void *>(this) << '\n';
    os << pfx << '\t' << "Solver name:      " << linearSolver::SolverName[solver_id_] << '\n';
    os << pfx << '\t' << "System name:      " << system_name_ << '\n';
    os << pfx << '\t' << "MPI ranks:        " << this->commSize() << '\n';
    os << pfx << '\t' << "Call count:       " << call_count_ << '\n';
    if (coeffs_) {
        os << pfx << '\t' << "Matrix all rows:  " << coeffs_->nGlobalCoefficients()<< '\n';
        os << pfx << '\t' << "Matrix all nnz:   " << coeffs_->nnzGlobal() << '\n';
        os << pfx << '\t' << "Matrix blocksize: " << BLOCKSIZE << '\n';
        MemoryFootprint data, connectivity;
        coeffs_->getMemoryFootprint(data, connectivity);
        const bool reduced = coeffs_->getGraph()->getLayout() & GraphLayout::Stencil__Reduced;
        os << pfx << '\t' << "Matrix fp data:   " << data.sum_byte / MemoryFootprint::BYTE_DIVIDE << " MB\n";
        os << pfx << '\t' << "Matrix int data:  " << connectivity.sum_byte / MemoryFootprint::BYTE_DIVIDE << " MB\n";
        os << pfx << '\t' << "Reduced stencil:  " << std::string(reduced ? "true" : "false") << '\n';
    }
    os << pfx << '\t' << "Min. iterations:  " << min_iterations_ << '\n';
    os << pfx << '\t' << "Max. iterations:  " << max_iterations_ << '\n';
    os << pfx << '\t' << "Rel. tolerance:   " << std::scientific << rtol_ << '\n';
    os << pfx << '\t' << "Abs. tolerance:   " << std::scientific << atol_ << '\n';
    os << pfx << '\t' << "Control residual: " << control_residual_ << '\n';
    // clang-format on
    return pfx;
}

template <size_t N>
void Context<N>::copyProtectedSettings(std::shared_ptr<const Context> ctx)
{
    assert(ctx);

    // override defaults
    solver_id_ = ctx->solver_id_;
    verbose_ = ctx->verbose_;
    min_iterations_ = ctx->min_iterations_;
    max_iterations_ = ctx->max_iterations_;
    rtol_ = ctx->rtol_;
    atol_ = ctx->atol_;
}

template <size_t N>
int Context<N>::setOptions_(const YAML::Node* solver_yaml)
{
    assert(solver_yaml);

    // override defaults
    const YAML::Node& s = *solver_yaml;
    if (s["verbose"])
    {
        verbose_ = s["verbose"].template as<int>();
    }
    if (s["min_iterations"])
    {
        min_iterations_ = s["min_iterations"].template as<int>();
    }
    if (s["max_iterations"])
    {
        max_iterations_ = s["max_iterations"].template as<int>();
    }
    if (s["rtol"])
    {
        rtol_ = s["rtol"].template as<double>();
    }
    if (s["atol"])
    {
        atol_ = s["atol"].template as<double>();
    }

    return 0;
}

template <size_t N>
void Context<N>::setupDefaults_()
{
    verbose_ = 0;
    min_iterations_ = 0;
    max_iterations_ = 20;
    rtol_ = 1.0e-6;
    atol_ = 1.0e-16;

    for (size_t i = 0; i < N; i++)
    {
        residual_scales_[i] = 1.0;
        control_residual_[i] = std::numeric_limits<double>::max();
    }
}

template <size_t N>
void Context<N>::setupIOStream_()
{
    active_stream_ = &dev_null_;
    if (0 == this->commRank())
    {
        active_stream_ = cout_;
    }
}

} /* namespace linearSolver */
