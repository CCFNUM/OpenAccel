// File       : ContextHYPRE.h
// Created    : Wed Mar 27 2024 09:44:54 (+0100)
// Author     : Fabian Wermelinger
// Description: Specialized HYPRE solver context
// Copyright (c) 2024 CCFNUM, Lucerne University of Applied Sciences and Arts.
// SPDX-License-Identifier: BSD-3-Clause

#ifndef CONTEXTHYPRE_H
#define CONTEXTHYPRE_H

#include "linearSolverContext.h"
#include "memoryLayout.h"
#include <algorithm>
#include <cassert>
#include <map>
#include <mpi.h>
#include <numeric>
#include <string>
#include <vector>

#ifdef HAS_HYPRE
#include <HYPRE.h>
#include <HYPRE_IJ_mv.h>
#include <HYPRE_krylov.h>
#include <HYPRE_parcsr_ls.h>
#endif /* HAS_HYPRE */

#define HYPRE_CALLBACK(callback, ...)                                          \
    do                                                                         \
    {                                                                          \
        if (callback)                                                          \
        {                                                                      \
            callback(__VA_ARGS__);                                             \
        }                                                                      \
    } while (0)

namespace linearSolver
{

#ifdef HAS_HYPRE
template <size_t N>
class ContextHYPRE : public Context<N>
{
    // supported HYPRE solver
    enum class SolverType
    {
        GMRES,
        FlexGMRES,
        BoomerAMG,
        MGR
    };

    // callback types
    typedef HYPRE_PtrToParSolverFcn HYPRESolverFunction;
    typedef HYPRE_Int (*HYPREGetNumIterFunction)(HYPRE_Solver, HYPRE_Int*);
    typedef HYPRE_Int (*HYPRESolverDestroy)(HYPRE_Solver);
    typedef HYPRE_Int (*HYPRESetPreconditioner)(HYPRE_Solver,
                                                HYPRE_PtrToParSolverFcn,
                                                HYPRE_PtrToParSolverFcn,
                                                HYPRE_Solver);

public:
    using Context<N>::BLOCKSIZE;
    using Coefficients = typename Context<N>::Coefficients;
    using Matrix = typename Coefficients::Matrix;
    using Vector = typename Coefficients::Vector;
    using Index = typename Coefficients::Index;
    using IndexVector = typename Coefficients::IndexVector;

    ContextHYPRE() = delete;

    ContextHYPRE(const std::string& family,
                 const std::string& system_name,
                 const CRSNodeGraph* graph,
                 const YAML::Node* node = nullptr)
        : Context<N>(linearSolver::ID::HYPRE, family, system_name, graph, node),
          options_node_(node)
    {
        setup_();
    }

    // used for preconditioning (shares coefficients)
    ContextHYPRE(const std::string& family,
                 const std::string& system_name,
                 std::shared_ptr<Coefficients> coeffs,
                 const YAML::Node* node = nullptr)
        : Context<N>(linearSolver::ID::HYPRE,
                     family,
                     system_name,
                     coeffs,
                     node),
          options_node_(node)
    {
        setup_();
    }

    ~ContextHYPRE()
    {
        HYPRE_IJMatrixDestroy(A_hypre_);
        HYPRE_IJVectorDestroy(b_hypre_);
        HYPRE_IJVectorDestroy(x_hypre_);
    }

    void solvePrologue(const int solver_id, const bool preconditioner) override
    {
        setupSolver_(options_node_);
        copyLocalToHYPRE_();
        Context<N>::solvePrologue(solver_id, preconditioner);
    }

    void solveEpilogue(const int solver_id, const bool preconditioner) override
    {
        copyHYPREToLocal_();
        destroySolver_();
        Context<N>::solveEpilogue(solver_id, preconditioner);
    }

    std::string info(std::ostream& os,
                     const char* prefix = nullptr) const override
    {
        const std::string pfx = Context<N>::info(os, prefix);
        for (const auto& option : options_)
        {
            os << pfx << '\t' << std::setw(18) << std::left << option.first
               << option.second << '\n';
        }
        return pfx;
    }

    int solve()
    {
        assert(cb_solver_setup_);
        assert(cb_solver_);
        HYPRE_CALLBACK(cb_solver_setup_, solver_hypre_, A_par_, b_par_, x_par_);
        HYPRE_CALLBACK(cb_solver_, solver_hypre_, A_par_, b_par_, x_par_);

        int iterations;
        assert(cb_solver_get_iter_);
        HYPRE_CALLBACK(cb_solver_get_iter_, solver_hypre_, &iterations);
        return iterations;
    }

protected:
    inline ContextHYPRE* castContextHYPRE_() override
    {
        return this;
    }

private:
    const YAML::Node* options_node_;
    std::map<std::string, std::string> options_;

    HYPRE_IJMatrix A_hypre_;
    HYPRE_IJVector b_hypre_;
    HYPRE_IJVector x_hypre_;
    HYPRE_Solver solver_hypre_;
    HYPRE_Solver preconditioner_hypre_;
    HYPRE_Solver mgr_coarse_hypre_;

    // helper
    HYPRE_ParCSRMatrix A_par_;
    HYPRE_ParVector b_par_;
    HYPRE_ParVector x_par_;

    // callbacks
    // solver
    HYPRESolverFunction cb_solver_;
    HYPRESolverFunction cb_solver_setup_;
    HYPRESolverDestroy cb_solver_destroy_;
    HYPREGetNumIterFunction cb_solver_get_iter_;

    // preconditioner
    HYPRESolverFunction cb_preconditioner_;
    HYPRESolverFunction cb_preconditioner_setup_;
    HYPRESolverDestroy cb_preconditioner_destroy_;
    HYPRESetPreconditioner cb_set_preconditioner_;

    // MGR only
    HYPRESolverDestroy cb_mgr_coarse_destroy_;

    // indices for global row (i) and column (j) span (for square matrices
    // i=j)
    Index i_lower_;
    Index i_upper_;
    Index j_lower_;
    Index j_upper_;
    Index n_unknowns;

    void setup_(); // at construction

    int setOptions_(const YAML::Node* solver_yaml) override
    {
        assert(solver_yaml);
        Context<N>::setOptions_(solver_yaml);
        this->setupSolver_(solver_yaml);
        return 0;
    }

    void setupSolver_(const YAML::Node* solver_yaml);
    void copyLocalToHYPRE_();
    void copyHYPREToLocal_();

    void resetCallbacks_()
    {
        cb_solver_ = nullptr;
        cb_solver_setup_ = nullptr;
        cb_solver_destroy_ = nullptr;
        cb_solver_get_iter_ = nullptr;

        cb_preconditioner_ = nullptr;
        cb_preconditioner_setup_ = nullptr;
        cb_preconditioner_destroy_ = nullptr;
        cb_set_preconditioner_ = nullptr;

        cb_mgr_coarse_destroy_ = nullptr;
    }

    void destroySolver_()
    {
        HYPRE_CALLBACK(cb_mgr_coarse_destroy_, mgr_coarse_hypre_);
        HYPRE_CALLBACK(cb_preconditioner_destroy_, preconditioner_hypre_);
        HYPRE_CALLBACK(cb_solver_destroy_, solver_hypre_);
        resetCallbacks_();
    }

    // clang-format off
    // solver/preconditioner settings
    void setSolver_(SolverType type, const YAML::Node *options = nullptr);
    void setPreconditioner_(SolverType type, const YAML::Node *options = nullptr);

    static void setGMRES_(const MPI_Comm comm, HYPRE_Solver &object, const YAML::Node *options = nullptr);
    static void setFlexGMRES_(const MPI_Comm comm, HYPRE_Solver &object, const YAML::Node *options = nullptr);
    static void setBoomerAMG_(HYPRE_Solver &object, const YAML::Node *options = nullptr);
    void setMGR_(HYPRE_Solver &object, const YAML::Node *options = nullptr);
    // clang-format on
};

template <size_t N>
void ContextHYPRE<N>::setup_()
{
    resetCallbacks_();

    const Matrix& A = this->getAMatrix();
    assert(A.getMemoryLayout() & GraphLayout::ColumnIndexOrder__Global);

    i_lower_ = A.localToGlobal(0) * BLOCKSIZE;
    i_upper_ = A.localToGlobal(A.nRows()) * BLOCKSIZE - 1;
    j_lower_ = i_lower_;
    j_upper_ = i_upper_;
    n_unknowns = i_upper_ - i_lower_ + 1;

    HYPRE_IJMatrixCreate(this->getCommunicator(),
                         i_lower_,
                         i_upper_,
                         j_lower_,
                         j_upper_,
                         &A_hypre_);
    HYPRE_IJMatrixSetObjectType(A_hypre_, HYPRE_PARCSR);

    HYPRE_IJVectorCreate(
        this->getCommunicator(), i_lower_, i_upper_, &b_hypre_);
    HYPRE_IJVectorSetObjectType(b_hypre_, HYPRE_PARCSR);
    HYPRE_IJVectorCreate(
        this->getCommunicator(), i_lower_, i_upper_, &x_hypre_);
    HYPRE_IJVectorSetObjectType(x_hypre_, HYPRE_PARCSR);

    copyLocalToHYPRE_();

    HYPRE_IJMatrixGetObject(A_hypre_, (void**)&A_par_);
    HYPRE_IJVectorGetObject(b_hypre_, (void**)&b_par_);
    HYPRE_IJVectorGetObject(x_hypre_, (void**)&x_par_);

    // populate options
    setOptions_(options_node_);
    destroySolver_();
}

template <size_t N>
void ContextHYPRE<N>::setupSolver_(const YAML::Node* node)
{
    if (!node || !(*node)["options"])
    {
        setSolver_(SolverType::GMRES); // default solver (w/o preconditioner)
    }
    else
    {
        const YAML::Node option = (*node)["options"];
        if (!option["solver"])
        {
            throw std::runtime_error(
                "ContextHYPRE: `options` node is missing required "
                "`solver` block");
        }
        if (option["solver"] && option["solver"]["type"])
        {
            const YAML::Node solver = option["solver"];
            std::string type = solver["type"].template as<std::string>();
            Context<N>::tolower(type); // set lower-case
            if (type == "gmres")
            {
                setSolver_(SolverType::GMRES, &solver);
            }
            else if (type == "flexgmres")
            {
                setSolver_(SolverType::FlexGMRES, &solver);
            }
            else if (type == "boomeramg")
            {
                setSolver_(SolverType::BoomerAMG, &solver);
            }
            else if (type == "mgr")
            {
                setSolver_(SolverType::MGR, &solver);
            }
            else
            {
                throw std::runtime_error(
                    "ContextHYPRE: solver type `" +
                    solver["type"].template as<std::string>() +
                    "` not supported");
            }
        }
        else
        {
            throw std::runtime_error(
                "ContextHYPRE: `solver` block requires field `type`");
        }

        if (option["precond"])
        {
            if (!option["precond"]["type"])
            {
                throw std::runtime_error(
                    "ContextHYPRE: `precond` block defined but is "
                    "missing required field `type`");
            }
            const YAML::Node precond = option["precond"];
            std::string type = precond["type"].template as<std::string>();
            Context<N>::tolower(type); // set lower-case
            if (type == "boomeramg")
            {
                setPreconditioner_(SolverType::BoomerAMG, &precond);
            }
            else if (type == "mgr")
            {
                setPreconditioner_(SolverType::MGR, &precond);
            }
            else
            {
                throw std::runtime_error(
                    "ContextHYPRE: preconditioner type `" +
                    precond["type"].template as<std::string>() +
                    "` not supported");
            }
        }
    }

    if (!(cb_solver_ && cb_solver_setup_ && cb_solver_get_iter_ &&
          cb_solver_destroy_))
    {
        throw std::runtime_error(
            "ContextHYPRE (internal): solver setup is inconsistent");
    }

    if (cb_set_preconditioner_ && cb_preconditioner_ &&
        cb_preconditioner_setup_ && cb_preconditioner_destroy_)
    {
        HYPRE_CALLBACK(cb_set_preconditioner_,
                       solver_hypre_,
                       cb_preconditioner_,
                       cb_preconditioner_setup_,
                       preconditioner_hypre_);
    }
}

template <size_t N>
void ContextHYPRE<N>::copyLocalToHYPRE_()
{
    const Matrix& A = this->getAMatrix();
    const Vector& b = this->getBVector();

    // TODO: [2024-10-14] assertion can be removed after
    // testing
    assert(A.getMemoryLayout() & GraphLayout::ColumnIndexOrder__Global);

    // TODO: [2024-03-27] verify HYPRE types and
    // scalar/Index are consistent
    std::vector<Index> row_nnz;
    std::vector<Index> row_idx;
    std::vector<Index> col_idx;
    std::vector<typename Matrix::DataType> values;

    // matrix
    HYPRE_IJMatrixInitialize(A_hypre_);
    for (Index i = 0; i < A.nRows(); i++)
    {
        matrixLayout::blockRowToRowMajor(i,
                                         A,
                                         A.getGraph()->rowGlobalIndices(i),
                                         row_nnz,
                                         row_idx,
                                         col_idx,
                                         values);
        HYPRE_IJMatrixSetValues(A_hypre_,
                                BLOCKSIZE,
                                row_nnz.data(),
                                row_idx.data(),
                                col_idx.data(),
                                values.data());
    }
    HYPRE_IJMatrixAssemble(A_hypre_);

    // vectors
    col_idx.resize(n_unknowns);
    std::iota(col_idx.begin(), col_idx.end(), i_lower_);
    values.resize(n_unknowns);
    std::fill(values.begin(), values.end(), 0);

    HYPRE_IJVectorInitialize(b_hypre_);
    HYPRE_IJVectorSetValues(b_hypre_, n_unknowns, col_idx.data(), b.data());
    HYPRE_IJVectorAssemble(b_hypre_);

    HYPRE_IJVectorInitialize(x_hypre_);
    HYPRE_IJVectorSetValues(
        x_hypre_, n_unknowns, col_idx.data(), values.data());
    HYPRE_IJVectorAssemble(x_hypre_);
}

template <size_t N>
void ContextHYPRE<N>::copyHYPREToLocal_()
{
    Vector& x = this->getXVector();
    std::vector<Index> col_idx(n_unknowns);
    std::iota(col_idx.begin(), col_idx.end(), i_lower_);
    HYPRE_IJVectorGetValues(x_hypre_, n_unknowns, col_idx.data(), x.data());
}

template <size_t N>
void ContextHYPRE<N>::setSolver_(SolverType type, const YAML::Node* options)
{
    if (type == SolverType::GMRES)
    {
        setGMRES_(this->getCommunicator(), solver_hypre_, options);
        cb_solver_ = HYPRE_ParCSRGMRESSolve;
        cb_solver_setup_ = HYPRE_ParCSRGMRESSetup;
        cb_solver_destroy_ = HYPRE_ParCSRGMRESDestroy;
        cb_solver_get_iter_ = HYPRE_GMRESGetNumIterations;
        cb_set_preconditioner_ = HYPRE_ParCSRGMRESSetPrecond;
        options_["HYPRE solver:"] = "GMRES";
        // enforced options (overwrite previous setting)
        HYPRE_GMRESSetMaxIter(solver_hypre_, this->max_iterations_);
        HYPRE_GMRESSetTol(solver_hypre_, this->rtol_);
        HYPRE_GMRESSetAbsoluteTol(solver_hypre_, this->atol_);
    }
    else if (type == SolverType::FlexGMRES)
    {
        setFlexGMRES_(this->getCommunicator(), solver_hypre_, options);
        cb_solver_ = HYPRE_ParCSRFlexGMRESSolve;
        cb_solver_setup_ = HYPRE_ParCSRFlexGMRESSetup;
        cb_solver_destroy_ = HYPRE_ParCSRFlexGMRESDestroy;
        cb_solver_get_iter_ = HYPRE_FlexGMRESGetNumIterations;
        cb_set_preconditioner_ = HYPRE_ParCSRFlexGMRESSetPrecond;
        options_["HYPRE solver:"] = "FlexGMRES";
        // enforced options (overwrite previous setting)
        HYPRE_FlexGMRESSetMaxIter(solver_hypre_, this->max_iterations_);
        HYPRE_FlexGMRESSetTol(solver_hypre_, this->rtol_);
        HYPRE_FlexGMRESSetAbsoluteTol(solver_hypre_, this->atol_);
    }
    else if (type == SolverType::BoomerAMG)
    {
        setBoomerAMG_(solver_hypre_, options);
        cb_solver_ = HYPRE_BoomerAMGSolve;
        cb_solver_setup_ = HYPRE_BoomerAMGSetup;
        cb_solver_destroy_ = HYPRE_BoomerAMGDestroy;
        cb_solver_get_iter_ = HYPRE_BoomerAMGGetNumIterations;
        options_["HYPRE solver:"] = "BoomerAMG";
        // enforced options (overwrite previous setting)
        HYPRE_BoomerAMGSetMaxIter(solver_hypre_, this->max_iterations_);
        HYPRE_BoomerAMGSetTol(solver_hypre_, this->rtol_);
    }
    else if (type == SolverType::MGR)
    {
        setMGR_(solver_hypre_, options);
        cb_solver_ = HYPRE_MGRSolve;
        cb_solver_setup_ = HYPRE_MGRSetup;
        cb_solver_destroy_ = HYPRE_MGRDestroy;
        cb_solver_get_iter_ = HYPRE_MGRGetNumIterations;
        options_["HYPRE solver:"] = "MGR";
        // enforced options (overwrite previous setting)
        HYPRE_MGRSetMaxIter(solver_hypre_, this->max_iterations_);
        HYPRE_MGRSetTol(solver_hypre_, this->rtol_);
    }
    else
    {
        throw std::runtime_error(
            "ContextHYPRE (internal): solver type not supported");
    }
}

template <size_t N>
void ContextHYPRE<N>::setPreconditioner_(SolverType type,
                                         const YAML::Node* options)
{
    if (type == SolverType::BoomerAMG)
    {
        setBoomerAMG_(preconditioner_hypre_, options);
        cb_preconditioner_ = HYPRE_BoomerAMGSolve;
        cb_preconditioner_setup_ = HYPRE_BoomerAMGSetup;
        cb_preconditioner_destroy_ = HYPRE_BoomerAMGDestroy;
        options_["HYPRE precond:"] = "BoomerAMG";
        // enforced options (overwrite previous setting)
        HYPRE_BoomerAMGSetMaxIter(preconditioner_hypre_, 1);
        HYPRE_BoomerAMGSetTol(preconditioner_hypre_, 0.0);
    }
    else if (type == SolverType::MGR)
    {
        setMGR_(preconditioner_hypre_, options);
        cb_preconditioner_ = HYPRE_MGRSolve;
        cb_preconditioner_setup_ = HYPRE_MGRSetup;
        cb_preconditioner_destroy_ = HYPRE_MGRDestroy;
        options_["HYPRE precond:"] = "MGR";
        // enforced options (overwrite previous setting)
        HYPRE_MGRSetMaxIter(preconditioner_hypre_, 1);
        HYPRE_MGRSetTol(preconditioner_hypre_, 0.0);
    }
    else
    {
        throw std::runtime_error(
            "ContextHYPRE (internal): preconditioner type not supported");
    }
}

template <size_t N>
void ContextHYPRE<N>::setGMRES_(const MPI_Comm comm,
                                HYPRE_Solver& object,
                                const YAML::Node* options)
{
    HYPRE_ParCSRGMRESCreate(comm, &object);
    if (options)
    {
        using Iter = YAML::const_iterator;
        for (Iter opt = options->begin(); opt != options->end(); ++opt)
        {
            std::string key = opt->first.template as<std::string>();
            Context<N>::tolower(key); // set lower-case
            // use HYPRE default if not set
            if (key == "printlevel")
            {
                HYPRE_GMRESSetPrintLevel(object,
                                         opt->second.template as<int>());
            }
            else if (key == "logging")
            {
                HYPRE_GMRESSetLogging(object, opt->second.template as<int>());
            }
            else if (key == "kdim")
            { // Krylov dimension for restarted GMRES
                HYPRE_GMRESSetKDim(object, opt->second.template as<int>());
            }
            // use OUR default if not set
        }
    }
}

template <size_t N>
void ContextHYPRE<N>::setFlexGMRES_(const MPI_Comm comm,
                                    HYPRE_Solver& object,
                                    const YAML::Node* options)
{
    HYPRE_ParCSRFlexGMRESCreate(comm, &object);
    if (options)
    {
        using Iter = YAML::const_iterator;
        for (Iter opt = options->begin(); opt != options->end(); ++opt)
        {
            std::string key = opt->first.template as<std::string>();
            Context<N>::tolower(key); // set lower-case
            // use HYPRE default if not set
            if (key == "printlevel")
            {
                HYPRE_FlexGMRESSetPrintLevel(object,
                                             opt->second.template as<int>());
            }
            else if (key == "logging")
            {
                HYPRE_FlexGMRESSetLogging(object,
                                          opt->second.template as<int>());
            }
            else if (key == "kdim")
            { // Krylov dimension for restarted GMRES
                HYPRE_FlexGMRESSetKDim(object, opt->second.template as<int>());
            }
            // use OUR default if not set
        }
    }
}

template <size_t N>
void ContextHYPRE<N>::setBoomerAMG_(HYPRE_Solver& object,
                                    const YAML::Node* options)
{
    HYPRE_BoomerAMGCreate(&object);
    if (options)
    {
        using Iter = YAML::const_iterator;
        for (Iter opt = options->begin(); opt != options->end(); ++opt)
        {
            std::string key = opt->first.template as<std::string>();
            Context<N>::tolower(key); // set lower-case
            // use HYPRE default if not set
            if (key == "printlevel")
            {
                HYPRE_BoomerAMGSetPrintLevel(object,
                                             opt->second.template as<int>());
            }
            else if (key == "logging")
            {
                HYPRE_BoomerAMGSetLogging(object,
                                          opt->second.template as<int>());
            }
            else if (key == "relaxtype")
            {
                HYPRE_BoomerAMGSetRelaxType(object,
                                            opt->second.template as<int>());
            }
            else if (key == "relaxorder")
            {
                HYPRE_BoomerAMGSetRelaxOrder(object,
                                             opt->second.template as<int>());
            }
            else if (key == "numsweeps")
            {
                HYPRE_BoomerAMGSetNumSweeps(object,
                                            opt->second.template as<int>());
            }
            else if (key == "maxiter")
            {
                HYPRE_BoomerAMGSetMaxIter(object,
                                          opt->second.template as<int>());
            }
            else if (key == "maxlevels")
            {
                HYPRE_BoomerAMGSetMaxLevels(object,
                                            opt->second.template as<int>());
            }
            else if (key == "interptype")
            {
                HYPRE_BoomerAMGSetInterpType(object,
                                             opt->second.template as<int>());
            }
            else if (key == "aggnumlevels")
            {
                HYPRE_BoomerAMGSetAggNumLevels(object,
                                               opt->second.template as<int>());
            }
            else if (key == "numpaths")
            {
                HYPRE_BoomerAMGSetNumPaths(object,
                                           opt->second.template as<int>());
            }
            else if (key == "agginterptype")
            {
                HYPRE_BoomerAMGSetAggInterpType(object,
                                                opt->second.template as<int>());
            }
            else if (key == "nodal")
            {
                HYPRE_BoomerAMGSetNodal(object, opt->second.template as<int>());
            }
            else if (key == "nodaldiag")
            {
                HYPRE_BoomerAMGSetNodalDiag(object,
                                            opt->second.template as<int>());
            }
            else if (key == "numfunctions")
            {
                HYPRE_BoomerAMGSetNumFunctions(object,
                                               opt->second.template as<int>());
            }
            else if (key == "coarsentype")
            {
                HYPRE_BoomerAMGSetCoarsenType(object,
                                              opt->second.template as<int>());
            }
            else if (key == "maxcoarsesize")
            {
                HYPRE_BoomerAMGSetMaxCoarseSize(object,
                                                opt->second.template as<int>());
            }
            else if (key == "coarsencutfactor")
            {
                HYPRE_BoomerAMGSetCoarsenCutFactor(
                    object, opt->second.template as<int>());
            }
            else if (key == "strongthreshold")
            {
                HYPRE_BoomerAMGSetStrongThreshold(
                    object,
                    opt->second.template as<typename Matrix::DataType>());
            }
            else if (key == "cycletype")
            {
                HYPRE_BoomerAMGSetCycleType(object,
                                            opt->second.template as<int>());
            }
            else if (key == "fcycle")
            {
                HYPRE_BoomerAMGSetFCycle(object,
                                         opt->second.template as<int>());
            }
            else if (key == "tol")
            {
                HYPRE_BoomerAMGSetTol(
                    object,
                    opt->second.template as<typename Matrix::DataType>());
            }
            else if (key == "truncfactor")
            {
                HYPRE_BoomerAMGSetTruncFactor(
                    object,
                    opt->second.template as<typename Matrix::DataType>());
            }
            else if (key == "aggtruncfactor")
            {
                HYPRE_BoomerAMGSetAggTruncFactor(
                    object,
                    opt->second.template as<typename Matrix::DataType>());
            }
            else if (key == "jacobitruncthreshold")
            {
                HYPRE_BoomerAMGSetJacobiTruncThreshold(
                    object,
                    opt->second.template as<typename Matrix::DataType>());
            }
            else if (key == "maxrowsum")
            {
                HYPRE_BoomerAMGSetMaxRowSum(
                    object,
                    opt->second.template as<typename Matrix::DataType>());
            }
            // use OUR default if not set
        }
    }
}

template <size_t N>
void ContextHYPRE<N>::setMGR_(HYPRE_Solver& object, const YAML::Node* options)
{
    HYPRE_MGRCreate(&object);
    if (options)
    {
        using Iter = YAML::const_iterator;

        for (Iter opt = options->begin(); opt != options->end(); ++opt)
        {
            std::string key = opt->first.template as<std::string>();
            Context<N>::tolower(key); // set lower-case
            if (key == "cpointsbyblock")
            {
                const YAML::Node cpoints_bb = opt->second;
                if (!cpoints_bb["n_levels"])
                {
                    throw std::runtime_error(
                        "ContextHYPRE: `n_levels` is required "
                        "when `CpointsByBlock` is defined");
                }
                if (!cpoints_bb["blocksize"])
                {
                    throw std::runtime_error(
                        "ContextHYPRE: `blocksize` is required "
                        "when `CpointsByBlock` is defined");
                }
                const int n_levels = cpoints_bb["n_levels"].template as<int>();
                const int blocksize =
                    cpoints_bb["blocksize"].template as<int>();
                assert(blocksize == BLOCKSIZE);

                std::vector<int> cpoints(n_levels);
                std::vector<int*> cindex(n_levels);
                std::vector<int> level_data(n_levels * blocksize, 0);
                for (int i = 0; i < n_levels; i++)
                {
                    int* idx = &level_data[i * blocksize];
                    cindex[i] = idx;
                    const std::string level_id("level_" + std::to_string(i));
                    if (!cpoints_bb[level_id])
                    {
                        throw std::runtime_error(
                            "ContextHYPRE: `" + level_id +
                            "` is required when `CpointsByBlock` is defined");
                    }
                    const std::vector<int> cpoints_level =
                        cpoints_bb[level_id].template as<std::vector<int>>();
                    assert(static_cast<int>(cpoints_level.size()) <= blocksize);
                    cpoints[i] = cpoints_level.size();
                    for (const int cpoint : cpoints_level)
                    {
                        *idx++ = cpoint;
                    }
                }
                HYPRE_MGRSetCpointsByBlock(
                    object, blocksize, n_levels, cpoints.data(), cindex.data());
            }
        }

        for (Iter opt = options->begin(); opt != options->end(); ++opt)
        {
            std::string key = opt->first.template as<std::string>();
            Context<N>::tolower(key); // set lower-case
            // use HYPRE default if not set
            if (key == "printlevel")
            {
                HYPRE_MGRSetPrintLevel(object, opt->second.template as<int>());
            }
            else if (key == "frelaxprintlevel")
            {
                HYPRE_MGRSetFrelaxPrintLevel(object,
                                             opt->second.template as<int>());
            }
            else if (key == "coarsegridprintlevel")
            {
                HYPRE_MGRSetCoarseGridPrintLevel(
                    object, opt->second.template as<int>());
            }
            else if (key == "logging")
            {
                HYPRE_MGRSetLogging(object, opt->second.template as<int>());
            }
            else if (key == "noncpointstofpoints")
            {
                HYPRE_MGRSetNonCpointsToFpoints(object,
                                                opt->second.template as<int>());
            }
            else if (key == "truncatecoarsegridthreshold")
            {
                HYPRE_MGRSetTruncateCoarseGridThreshold(
                    object,
                    opt->second.template as<typename Matrix::DataType>());
            }
            else if (key == "maxcoarselevels")
            {
                HYPRE_MGRSetMaxCoarseLevels(object,
                                            opt->second.template as<int>());
            }
            else if (key == "coarsegridmethod")
            {
                HYPRE_MGRSetCoarseGridMethod(
                    object, opt->second.template as<std::vector<int>>().data());
            }
            else if (key == "relaxtype")
            {
                HYPRE_MGRSetRelaxType(object, opt->second.template as<int>());
            }
            else if (key == "blocksize")
            {
                HYPRE_MGRSetBlockSize(object, opt->second.template as<int>());
            }
            else if (key == "levelfrelaxmethod")
            {
                HYPRE_MGRSetLevelFRelaxMethod(
                    object, opt->second.template as<std::vector<int>>().data());
            }
            else if (key == "levelfrelaxtype")
            {
                HYPRE_MGRSetLevelFRelaxType(
                    object, opt->second.template as<std::vector<int>>().data());
            }
            else if (key == "levelfrelaxnumfunctions")
            {
                HYPRE_MGRSetLevelFRelaxNumFunctions(
                    object, opt->second.template as<std::vector<int>>().data());
            }
            else if (key == "numrestrictsweeps")
            {
                HYPRE_MGRSetNumRestrictSweeps(object,
                                              opt->second.template as<int>());
            }
            else if (key == "restricttype")
            {
                HYPRE_MGRSetRestrictType(object,
                                         opt->second.template as<int>());
            }
            else if (key == "levelrestricttype")
            {
                HYPRE_MGRSetLevelRestrictType(
                    object, opt->second.template as<std::vector<int>>().data());
            }
            else if (key == "interptype")
            {
                HYPRE_MGRSetInterpType(object, opt->second.template as<int>());
            }
            else if (key == "levelinterptype")
            {
                HYPRE_MGRSetLevelInterpType(
                    object, opt->second.template as<std::vector<int>>().data());
            }
            else if (key == "numrelaxsweeps")
            {
                HYPRE_MGRSetNumRelaxSweeps(object,
                                           opt->second.template as<int>());
            }
            else if (key == "levelnumrelaxsweeps")
            {
                HYPRE_MGRSetLevelNumRelaxSweeps(
                    object, opt->second.template as<std::vector<int>>().data());
            }
            else if (key == "numinterpsweeps")
            {
                HYPRE_MGRSetNumInterpSweeps(object,
                                            opt->second.template as<int>());
            }
            else if (key == "blockjacobiblocksize")
            {
                HYPRE_MGRSetBlockJacobiBlockSize(
                    object, opt->second.template as<int>());
            }
            else if (key == "globalsmoothcycle")
            {
                HYPRE_MGRSetGlobalSmoothCycle(object,
                                              opt->second.template as<int>());
            }
            else if (key == "globalsmoothtype")
            {
                HYPRE_MGRSetGlobalSmoothType(object,
                                             opt->second.template as<int>());
            }
            else if (key == "levelsmoothtype")
            {
                HYPRE_MGRSetLevelSmoothType(
                    object, opt->second.template as<std::vector<int>>().data());
            }
            else if (key == "maxglobalsmoothiters")
            {
                HYPRE_MGRSetMaxGlobalSmoothIters(
                    object, opt->second.template as<int>());
            }
            else if (key == "levelsmoothiters")
            {
                HYPRE_MGRSetLevelSmoothIters(
                    object, opt->second.template as<std::vector<int>>().data());
            }
            else if (key == "pmaxelmts")
            {
                HYPRE_MGRSetPMaxElmts(object, opt->second.template as<int>());
            }
            else if (key == "levelpmaxelmts")
            {
                HYPRE_MGRSetLevelPMaxElmts(
                    object, opt->second.template as<std::vector<int>>().data());
            }
        }
    }

    // set coarse grid solver (hardcoded BoomerAMG)
    setBoomerAMG_(mgr_coarse_hypre_, options);
    HYPRE_MGRSetCoarseSolver(
        object, HYPRE_BoomerAMGSolve, HYPRE_BoomerAMGSetup, mgr_coarse_hypre_);
    cb_mgr_coarse_destroy_ = HYPRE_BoomerAMGDestroy;
}

#else

template <size_t N>
class ContextHYPRE
{
};

#endif /* HAS_HYPRE */

} /* namespace linearSolver */

#undef HYPRE_CALLBACK

#endif /* CONTEXTHYPRE_H */
