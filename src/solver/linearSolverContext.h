// File       : linearSolverContext.h
// Created    : Mon Mar 04 2024 11:22:21 (+0100)
// Author     : Fabian Wermelinger
// Description: Linear solver context base class
// Copyright (c) 2024 CCFNUM, Lucerne University of Applied Sciences and Arts.
// SPDX-License-Identifier: BSD-3-Clause

#ifndef LINEARSOLVERCONTEXT_H
#define LINEARSOLVERCONTEXT_H

// code
#include "CRSNodeGraph.h"
#include "Common.h"
#include "Profiler.h"
#include "coefficients.h"
#include "controlData.h"
#include "residual.h"

// std
#include <algorithm>
#include <cassert>
#include <iomanip>
#include <memory>
#include <mpi.h>
#include <ostream>
#include <stdexcept>
#include <string>
#include <vector>
#include <yaml-cpp/yaml.h>

namespace linearSolver
{

// implemented linear solver identifiers
namespace ID
{
enum : int
{
    GMRES = 0,
    FlexGMRES,
    DirectLU,
    AMG,
    PETSc,
    HYPRE,
    Trilinos,
    UNDEFINED
};
} /* namespace ID */

// must be in same order as in `ID`
static constexpr const char* SolverName[] = {"GMRES",
                                             "FelxGMRES",
                                             "DirectLU",
                                             "AMG",
                                             "PETSc",
                                             "HYPRE",
                                             "Trilinos",
                                             "UNDEFINED"};

template <size_t N>
class ContextAMGSolver;

template <size_t N>
class ContextGMRES;

template <size_t N>
class ContextPETSc;

template <size_t N>
class ContextHYPRE;

template <size_t N>
class ContextDirectLU;

template <size_t N>
class ContextTrilinos;

template <size_t N>
class Context
{
public:
    using Coefficients = coefficients<N>;
    using Matrix = typename Coefficients::Matrix;
    using Vector = typename Coefficients::Vector;
    using Index = typename Coefficients::Index;
    using IndexVector = typename Coefficients::IndexVector;
    using ControlData = ::linearSolver::ControlData<N>;
    using Array = typename ControlData::Array;

    static constexpr Index BLOCKSIZE = N;

    Context() = delete;

    Context(const int solver_id,
            const std::string& family,
            const std::string& system_name,
            const CRSNodeGraph* graph,
            const YAML::Node* node = nullptr);

    Context(const int solver_id,
            const std::string& family,
            const std::string& system_name,
            std::shared_ptr<Coefficients> coeffs,
            const YAML::Node* node = nullptr);

    virtual ~Context()
    {
    }

    Context(const Context& c) = delete;
    Context& operator=(const Context& c) = delete;
    Context(Context&& c) = delete;
    Context& operator=(Context&& c) = delete;

    friend std::ostream& operator<<(std::ostream& os, const Context& ctx)
    {
        ctx.info(os);
        return os;
    }

    // public API
    GraphLayout getLayout() const
    {
        assert(coeffs_);
        return coeffs_->getGraph()->getLayout();
    }

    MPI_Comm getCommunicator() const
    {
        assert(coeffs_);
        return coeffs_->getGraph()->getCommunicator();
    }

    int commRank() const
    {
        assert(coeffs_);
        return coeffs_->getGraph()->commRank();
    }

    int commSize() const
    {
        assert(coeffs_);
        return coeffs_->getGraph()->commSize();
    }

    size_t getCallCount() const
    {
        return call_count_;
    }

    size_t getCallCountRelative() const
    {
        return call_count_relative_;
    }

    void resetCallCountRelative()
    {
        call_count_relative_ = 0;
    }

    ControlData getControlData() const
    {
        return control_;
    }

    ControlData& getControlData()
    {
        return control_;
    }

    Array getControlResidual() const
    {
        return control_residual_;
    }

    Array& getControlResidual()
    {
        return control_residual_;
    }

    Array getResidualScales() const
    {
        return residual_scales_;
    }

    Array& getResidualScales()
    {
        return residual_scales_;
    }

    size_t operator++()
    {
        ++call_count_relative_;
        return ++call_count_;
    }

    size_t operator++(int)
    {
        const size_t old = call_count_;
        this->operator++();
        return old;
    }

    virtual void solvePrologue(const int /* solver_id */,
                               const bool preconditioner)
    {
        // initial residual
        auto& coeff = this->getCoefficients();
        auto& ctrl = this->getControlData();
        if (!preconditioner)
        {
            residual::compute(coeff);
            ctrl.scaled_initial_res = residual::RMSErrorDiagNormalized(coeff);
            ctrl.solver_initial_res = residual::RMSError(coeff);
        }
    }

    virtual void solveEpilogue(const int /* solver_id */,
                               const bool preconditioner)
    {
        // final residual
        auto& coeff = this->getCoefficients();
        auto& ctrl = this->getControlData();
        if (!preconditioner)
        {
            residual::compute(coeff);
            ctrl.scaled_final_res = residual::RMSErrorDiagNormalized(coeff);
            ctrl.solver_final_res = residual::RMSError(coeff);
        }
    }

    virtual std::string info(std::ostream& os,
                             const char* prefix = nullptr) const;

    std::ostream& cout()
    {
        assert(active_stream_);
        return *active_stream_;
    }

    std::ostream& coutAll()
    {
        return *cout_;
    }

    std::shared_ptr<std::ostream> makeFileStream(const std::string& fname) const
    {
        if (0 == this->commRank())
        {
            return std::make_shared<std::ofstream>(fname, std::ios::trunc);
        }
        else
        {
            return std::make_shared<onullstream>();
        }
    }

    static void tolower(std::string& s)
    {
        std::transform(s.begin(), s.end(), s.begin(), [](unsigned char c) {
            return std::tolower(c);
        });
    }

    // clang-format off
    // coefficient data handlers
    Coefficients &getCoefficients() { assert(coeffs_); return *coeffs_; }
    const Coefficients &getCoefficients() const { assert(coeffs_); return *coeffs_; }

    Matrix &getAMatrix() { assert(coeffs_); return coeffs_->getAMatrix(); }
    const Matrix &getAMatrix() const { assert(coeffs_); return coeffs_->getAMatrix(); }
    Vector &getXVector() { assert(coeffs_); return coeffs_->getXVector(); }
    const Vector &getXVector() const { assert(coeffs_); return coeffs_->getXVector(); }
    Vector &getBVector() { assert(coeffs_); return coeffs_->getBVector(); }
    const Vector &getBVector() const { assert(coeffs_); return coeffs_->getBVector(); }
    Vector &getRVector() { assert(coeffs_); return coeffs_->getRVector(); }
    const Vector &getRVector() const { assert(coeffs_); return coeffs_->getRVector(); }

    // clang-format on

    // coefficient management
    std::shared_ptr<Coefficients> coeffs()
    {
        assert(coeffs_);
        return coeffs_;
    }

    std::shared_ptr<const Coefficients> coeffs() const
    {
        assert(coeffs_);
        return coeffs_;
    }

    void setCoeffs(std::shared_ptr<Coefficients> coeffs)
    {
        assert(coeffs);
        this->coeffs_ = coeffs;
        this->profiler_.setComm(this->coeffs_->getCommunicator());
        this->setupIOStream_();
        if (ctx_precond_)
        {
            ctx_precond_->setCoeffs(this->coeffs_);
        }
    }

    std::shared_ptr<Coefficients>
    swapCoeffs(std::shared_ptr<Coefficients> coeffs)
    {
        std::shared_ptr<Coefficients> tmp = this->coeffs_;
        this->setCoeffs(coeffs);
        return tmp;
    }

    void zeroSystemStorage()
    {
        assert(coeffs_);
        coeffs_->resizeGraph();
        coeffs_->zeroLHS();
        coeffs_->zeroRHS();
        coeffs_->zeroSOL();
    }

    int solverID() const
    {
        return solver_id_;
    }

    int verbose() const
    {
        return verbose_;
    }

    void setVerbose(int v)
    {
        verbose_ = v;
    }

    int minIterations() const
    {
        return min_iterations_;
    }

    int maxIterations() const
    {
        return max_iterations_;
    }

    void setMinIterations(int it)
    {
        min_iterations_ = it;
    }

    void setMaxIterations(int it)
    {
        max_iterations_ = it;
    }

    double rtol() const
    {
        return rtol_;
    }

    void setRtol(double tol)
    {
        rtol_ = tol;
    }

    double atol() const
    {
        return atol_;
    }

    void setAtol(double tol)
    {
        atol_ = tol;
    }

    Profiler& getProfiler()
    {
        return profiler_;
    }

    std::string getFamily() const
    {
        return family_;
    }

    std::string getSystemName() const
    {
        return system_name_;
    }

    std::shared_ptr<Context> getPreconditioner() const
    {
        return ctx_precond_;
    }

    void setPreconditioner(std::shared_ptr<Context> ctx)
    {
        assert(ctx);
        assert(ctx->coeffs() == this->coeffs());
        ctx_precond_ = ctx;
    }

    // allowed cast operators
    operator ContextAMGSolver<N>*()
    {
        return castContextAMG_();
    }

    operator ContextGMRES<N>*()
    {
        return castContextGMRES_();
    }

    operator ContextPETSc<N>*()
    {
        return castContextPETSc_();
    }

    operator ContextHYPRE<N>*()
    {
        return castContextHYPRE_();
    }

    operator ContextDirectLU<N>*()
    {
        return castContextDirectLU_();
    }

    operator ContextTrilinos<N>*()
    {
        return castContextTrilinos_();
    }

    void copyProtectedSettings(std::shared_ptr<const Context> ctx);

protected:
    Context(const std::string& family, const std::string& system_name);

    int solver_id_;
    int verbose_;
    int min_iterations_;
    int max_iterations_;
    double rtol_;
    double atol_;

    std::shared_ptr<Coefficients> coeffs_; // system coefficients
    Profiler profiler_;

    virtual int setOptions_(const YAML::Node* solver_yaml);

    // allowed cast operators
    virtual ContextAMGSolver<N>* castContextAMG_()
    {
        throw std::runtime_error(
            "linearSolver::Context: illegal cast to ContextAMGSolver");
        return nullptr;
    }

    virtual ContextGMRES<N>* castContextGMRES_()
    {
        throw std::runtime_error(
            "linearSolver::Context: illegal cast to ContextGMRES");
        return nullptr;
    }

    virtual ContextPETSc<N>* castContextPETSc_()
    {
        throw std::runtime_error(
            "linearSolver::Context: illegal cast to ContextPETSc");
        return nullptr;
    }

    virtual ContextHYPRE<N>* castContextHYPRE_()
    {
        throw std::runtime_error(
            "linearSolver::Context: illegal cast to ContextHYPRE");
        return nullptr;
    }

    virtual ContextDirectLU<N>* castContextDirectLU_()
    {
        throw std::runtime_error(
            "linearSolver::Context: illegal cast to ContextDirectLU");
        return nullptr;
    }

    virtual ContextTrilinos<N>* castContextTrilinos_()
    {
        throw std::runtime_error(
            "linearSolver::Context: illegal cast to ContextTrilinos");
        return nullptr;
    }

private:
    const std::string family_;
    const std::string system_name_;
    size_t call_count_;
    size_t call_count_relative_;
    ControlData control_;
    Array control_residual_;
    Array residual_scales_;
    std::shared_ptr<Context> ctx_precond_; // preconditioner context

    // I/O streams
    onullstream dev_null_;
    std::ostream* cout_;
    std::ostream* active_stream_;

    void setupDefaults_();
    void setupIOStream_();
};

} /* namespace linearSolver */

#include "linearSolverContext.hpp"

#endif /* LINEARSOLVERCONTEXT_H */
