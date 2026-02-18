// File       : ContextTrilinos.h
// Created    : Wed Mar 05 2025 11:04:12 (+0100)
// Author     : Mhamad Mahdi Alloush
// Description: Trilinos (Tpetra/Belos) linear solver context
// Copyright (c) 2025 CCFNUM, Lucerne University of Applied Sciences and Arts.
// SPDX-License-Identifier: BSD-3-Clause

#ifndef CONTEXTTRILINOS_H
#define CONTEXTTRILINOS_H

#include "linearSolverContext.h"
#include "memoryLayout.h"

#ifdef HAS_TRILINOS
#include <BelosLinearProblem.hpp>
#include <BelosSolverFactory.hpp>
#include <BelosTpetraAdapter.hpp>
#include <Ifpack2_Factory.hpp>
#include <Teuchos_Array.hpp>
#include <Teuchos_DefaultMpiComm.hpp>
#include <Teuchos_DefaultSerialComm.hpp>
#include <Teuchos_OpaqueWrapper.hpp>
#include <Teuchos_OrdinalTraits.hpp>
#include <Teuchos_ParameterList.hpp>
#include <Teuchos_RCP.hpp>
#include <Tpetra_CrsMatrix.hpp>
#include <Tpetra_Map.hpp>
#include <Tpetra_Vector.hpp>
#endif /* HAS_TRILINOS */

#include <algorithm>
#include <cctype>
#include <memory>
#include <string>
#include <vector>

namespace linearSolver
{

#ifdef HAS_TRILINOS
namespace details
{
inline std::string toLower(std::string value)
{
    std::transform(
        value.begin(), value.end(), value.begin(), [](unsigned char c) {
        return static_cast<char>(std::tolower(c));
    });
    return value;
}

inline bool isBoolString(const std::string& value)
{
    const std::string key = toLower(value);
    return key == "true" || key == "false" || key == "yes" || key == "no" ||
           key == "on" || key == "off";
}

inline std::string mapBelosSolver(const std::string& value)
{
    const std::string key = toLower(value);
    if (key == "cg")
        return "CG";
    if (key == "bicgstab")
        return "BiCGStab";
    if (key == "tfqmr")
        return "TFQMR";
    if (key == "lsqr")
        return "LSQR";
    return "GMRES";
}

inline std::string mapIfpackPreconditioner(const std::string& value)
{
    const std::string key = toLower(value);
    if (key == "ilu" || key == "rilu")
        return "ILUT";
    if (key == "riluk")
        return "RILUK";
    if (key == "jacobi" || key == "relaxation")
        return "RELAXATION";
    if (key == "chebyshev")
        return "CHEBYSHEV";
    return "RELAXATION";
}
} // namespace details

template <size_t N>
class ContextTrilinos : public Context<N>
{
public:
    using typename Context<N>::Coefficients;
    using Matrix = typename Coefficients::Matrix;
    using Vector = typename Coefficients::Vector;
    using Index = typename Coefficients::Index;

    using scalar_type = double;
    using local_ordinal_type = int;
    using global_ordinal_type = long;
    using node_type = Tpetra::Map<>::node_type;
    using map_type =
        Tpetra::Map<local_ordinal_type, global_ordinal_type, node_type>;
    using matrix_type = Tpetra::CrsMatrix<scalar_type,
                                          local_ordinal_type,
                                          global_ordinal_type,
                                          node_type>;
    using vector_type = Tpetra::
        Vector<scalar_type, local_ordinal_type, global_ordinal_type, node_type>;
    using multivector_type = Tpetra::MultiVector<scalar_type,
                                                 local_ordinal_type,
                                                 global_ordinal_type,
                                                 node_type>;
    using operator_type = Tpetra::Operator<scalar_type,
                                           local_ordinal_type,
                                           global_ordinal_type,
                                           node_type>;
    using prec_type = Ifpack2::Preconditioner<scalar_type,
                                              local_ordinal_type,
                                              global_ordinal_type,
                                              node_type>;

    ContextTrilinos(const std::string& family,
                    const std::string& system_name,
                    const CRSNodeGraph* graph,
                    const YAML::Node* node = nullptr)
        : Context<N>(linearSolver::ID::Trilinos,
                     family,
                     system_name,
                     graph,
                     node),
          matrix_initialized_(false), use_preconditioner_(false),
          preconditioner_initialized_(false), belos_solver_name_("GMRES"),
          preconditioner_type_("RELAXATION")
    {
        setup_(node);
    }

    ContextTrilinos(const std::string& family,
                    const std::string& system_name,
                    std::shared_ptr<Coefficients> coeffs,
                    const YAML::Node* node = nullptr)
        : Context<N>(linearSolver::ID::Trilinos,
                     family,
                     system_name,
                     coeffs,
                     node),
          matrix_initialized_(false), use_preconditioner_(false),
          preconditioner_initialized_(false), belos_solver_name_("GMRES"),
          preconditioner_type_("RELAXATION")
    {
        setup_(node);
    }

    ~ContextTrilinos()
    {
        destroySystem_();
    }

    void solvePrologue(const int solver_id, const bool preconditioner) override
    {
        copyCoefficients_();
        copyStateToTpetra_();
        Context<N>::solvePrologue(solver_id, preconditioner);
    }

    std::string info(std::ostream& os,
                     const char* prefix = nullptr) const override
    {
        const std::string pfx = Context<N>::info(os, prefix);
        os << pfx << '\t' << "Belos solver:     " << belos_solver_name_ << '\n';
        os << pfx << '\t' << "Preconditioner:   "
           << (use_preconditioner_ ? preconditioner_type_ : std::string("none"))
           << '\n';
        return pfx;
    }

    int solve()
    {
        using BelosProblem =
            Belos::LinearProblem<scalar_type, multivector_type, operator_type>;
        if (belos_problem_.is_null())
        {
            belos_problem_ =
                Teuchos::rcp(new BelosProblem(matrix_, x_vec_, b_vec_));
        }
        belos_problem_->setProblem();
        if (use_preconditioner_ && !preconditioner_.is_null())
        {
            belos_problem_->setRightPrec(preconditioner_);
        }
        else
        {
            belos_problem_->setRightPrec(Teuchos::null);
        }

        if (belos_solver_.is_null())
        {
            createBelosSolver_();
        }
        belos_solver_->setProblem(belos_problem_);
        const Belos::ReturnType result = belos_solver_->solve();
        copyStateToNative_();

        if (result != Belos::Converged && this->verbose() > 0 &&
            this->commRank() == 0)
        {
            this->cout() << "Belos solver '" << belos_solver_name_
                         << "' did not converge after "
                         << belos_solver_->getNumIters() << " iterations.\n";
        }
        return static_cast<int>(belos_solver_->getNumIters());
    }

protected:
    inline ContextTrilinos* castContextTrilinos_() override
    {
        return this;
    }

private:
    Teuchos::RCP<const Teuchos::Comm<int>> comm_;
    Teuchos::RCP<const map_type> map_;
    Teuchos::RCP<matrix_type> matrix_;
    Teuchos::RCP<vector_type> x_vec_;
    Teuchos::RCP<vector_type> b_vec_;
    Teuchos::RCP<prec_type> preconditioner_;
    Teuchos::RCP<
        Belos::LinearProblem<scalar_type, multivector_type, operator_type>>
        belos_problem_;
    Teuchos::RCP<
        Belos::SolverManager<scalar_type, multivector_type, operator_type>>
        belos_solver_;
    Teuchos::RCP<Teuchos::ParameterList> belos_params_;
    Teuchos::ParameterList preconditioner_params_;
    bool matrix_initialized_;
    bool use_preconditioner_;
    bool preconditioner_initialized_;
    std::string belos_solver_name_;
    std::string preconditioner_type_;

    void setup_(const YAML::Node* node)
    {
        buildCommunicator_();
        buildMap_();
        buildMatrix_();
        buildVectors_();
        configureFromYaml_(node);
    }

    void configureFromYaml_(const YAML::Node* node)
    {
        belos_params_ = Teuchos::parameterList("Belos");
        belos_params_->set("Maximum Iterations", this->max_iterations_);
        belos_params_->set("Convergence Tolerance", this->rtol_);

        if (node)
        {
            Context<N>::setOptions_(node);
            const YAML::Node& s = *node;
            if (s["options"])
            {
                const YAML::Node& opts = s["options"];
                if (opts["belos_solver"])
                {
                    belos_solver_name_ = details::mapBelosSolver(
                        opts["belos_solver"].template as<std::string>());
                }
                if (opts["belos_parameters"])
                {
                    const auto& belos_opt = opts["belos_parameters"];
                    for (const auto& entry : belos_opt)
                    {
                        const std::string key =
                            entry.first.template as<std::string>();
                        const std::string value =
                            entry.second.template as<std::string>();
                        belos_params_->set(key, value);
                    }
                }
                if (opts["preconditioner"])
                {
                    const std::string prec = details::toLower(
                        opts["preconditioner"].template as<std::string>());
                    if (prec == "none")
                    {
                        use_preconditioner_ = false;
                    }
                    else
                    {
                        use_preconditioner_ = true;
                        preconditioner_type_ =
                            details::mapIfpackPreconditioner(prec);
                    }
                }
                if (opts["preconditioner_parameters"])
                {
                    const auto& prec = opts["preconditioner_parameters"];
                    for (const auto& entry : prec)
                    {
                        const std::string key =
                            entry.first.template as<std::string>();
                        addPreconditionerParameter_(key, entry.second);
                    }
                }
            }
        }

        if (use_preconditioner_)
        {
            if (preconditioner_type_ == "RELAXATION" &&
                !preconditioner_params_.isParameter("relaxation: type"))
            {
                preconditioner_params_.set("relaxation: type", "Jacobi");
            }
            createPreconditioner_();
        }
    }

    void buildCommunicator_()
    {
        if (comm_.is_null())
        {
#ifdef HAVE_MPI
            auto raw = Teuchos::opaqueWrapper(this->getCommunicator());
            comm_ = Teuchos::rcp(new Teuchos::MpiComm<int>(raw));
#else
            comm_ = Teuchos::rcp(new Teuchos::SerialComm<int>());
#endif
        }
    }

    void buildMap_()
    {
        const Matrix& A = this->getAMatrix();
        const Index n_local = A.nRows();
        Teuchos::Array<global_ordinal_type> gids(n_local * N);
        for (Index i = 0; i < n_local; ++i)
        {
            const Index gid = A.localToGlobal(i);
            for (Index k = 0; k < N; ++k)
            {
                gids[i * N + k] = static_cast<global_ordinal_type>(gid * N + k);
            }
        }
        map_ = Teuchos::rcp(new map_type(
            Teuchos::OrdinalTraits<Tpetra::global_size_t>::invalid(),
            gids(),
            0,
            comm_));
    }

    void buildMatrix_()
    {
        // Compute the number of non-zeros per row from the graph
        // Each block row has nnz blocks, and each block is N x N
        const Matrix& A = this->getAMatrix();
        const Index n_local = A.nRows();

        // Create array of nnz per scalar row (each block row becomes N scalar
        // rows)
        Teuchos::ArrayRCP<size_t> nnz_per_row(n_local * N);
        for (Index i = 0; i < n_local; ++i)
        {
            // Number of block columns in this row
            const size_t block_nnz = A.getGraph()->rowGlobalIndices(i).size();
            // Each scalar row within this block row has block_nnz * N entries
            const size_t scalar_nnz = block_nnz * N;
            for (Index k = 0; k < N; ++k)
            {
                nnz_per_row[i * N + k] = scalar_nnz;
            }
        }

        matrix_ = Teuchos::rcp(new matrix_type(map_, nnz_per_row()));
        matrix_initialized_ = false;
    }

    void buildVectors_()
    {
        x_vec_ = Teuchos::rcp(new vector_type(map_));
        b_vec_ = Teuchos::rcp(new vector_type(map_));
    }

    void destroySystem_()
    {
        belos_solver_ = Teuchos::null;
        belos_problem_ = Teuchos::null;
        preconditioner_ = Teuchos::null;
        matrix_ = Teuchos::null;
        x_vec_ = Teuchos::null;
        b_vec_ = Teuchos::null;
        map_ = Teuchos::null;
        comm_ = Teuchos::null;
        preconditioner_initialized_ = false;
    }

    void createBelosSolver_()
    {
        Belos::SolverFactory<scalar_type, multivector_type, operator_type>
            factory;
        belos_solver_ = factory.create(belos_solver_name_, belos_params_);
    }

    void createPreconditioner_()
    {
        Ifpack2::Factory factory;
        const std::string factory_name =
            preconditioner_type_.empty() ? "RELAXATION" : preconditioner_type_;
        preconditioner_ = factory.create<matrix_type>(factory_name, matrix_);
        preconditioner_->setParameters(preconditioner_params_);
        preconditioner_initialized_ = false;
    }

    void addPreconditionerParameter_(const std::string& key,
                                     const YAML::Node& value)
    {
        if (!value.IsScalar())
        {
            return;
        }

        const std::string scalar = value.Scalar();
        if (details::isBoolString(scalar))
        {
            preconditioner_params_.set(key, value.template as<bool>());
            return;
        }

        try
        {
            preconditioner_params_.set(key, value.template as<int>());
            return;
        }
        catch (const YAML::BadConversion&)
        {
        }

        try
        {
            preconditioner_params_.set(key, value.template as<double>());
            return;
        }
        catch (const YAML::BadConversion&)
        {
        }

        preconditioner_params_.set(key, scalar);
    }

    void copyStateToTpetra_()
    {
        Vector& x_native = this->coeffs_->getXVector();
        Vector& b_native = this->coeffs_->getBVector();
        auto x_view = x_vec_->get1dViewNonConst();
        auto b_view = b_vec_->get1dViewNonConst();
        const size_t n_local = static_cast<size_t>(map_->getLocalNumElements());
        std::copy_n(x_native.begin(), n_local, x_view.begin());
        std::copy_n(b_native.begin(), n_local, b_view.begin());
    }

    void copyStateToNative_()
    {
        auto x_view = x_vec_->get1dView();
        Vector& x_native = this->coeffs_->getXVector();
        const size_t n_local = static_cast<size_t>(map_->getLocalNumElements());
        std::copy_n(x_view.begin(), n_local, x_native.begin());
    }

    void copyCoefficients_()
    {
        const Matrix& A = this->getAMatrix();
        std::vector<Index> row_nnz;
        std::vector<Index> row_idx;
        std::vector<Index> col_idx;
        std::vector<typename Matrix::DataType> values;
        const Index n_rows = A.nRows();

        // Resume fill mode if matrix was previously finalized
        // After resumeFill(), only replaceGlobalValues() is allowed (not
        // insert)
        if (matrix_initialized_)
        {
            matrix_->resumeFill();
        }

        for (Index i = 0; i < n_rows; ++i)
        {
            const auto cols = A.getGraph()->rowGlobalIndices(i);
            matrixLayout::blockRowToRowMajor(
                i, A, cols, row_nnz, row_idx, col_idx, values);
            const Index block_nnz = row_nnz.front();
            std::vector<global_ordinal_type> column_ids(block_nnz);
            for (Index k = 0; k < N; ++k)
            {
                const global_ordinal_type row =
                    static_cast<global_ordinal_type>(row_idx[k]);
                for (Index j = 0; j < block_nnz; ++j)
                {
                    column_ids[j] = static_cast<global_ordinal_type>(
                        col_idx[k * block_nnz + j]);
                }
                Teuchos::ArrayView<const global_ordinal_type> column_view(
                    column_ids.data(), static_cast<int>(block_nnz));
                Teuchos::ArrayView<const scalar_type> value_view(
                    &values[k * block_nnz], static_cast<int>(block_nnz));
                if (!matrix_initialized_)
                {
                    matrix_->insertGlobalValues(row, column_view, value_view);
                }
                else
                {
                    // After resumeFill(), we can only replace values in
                    // existing positions
                    matrix_->replaceGlobalValues(row, column_view, value_view);
                }
            }
        }

        // Always call fillComplete after modifying values
        // Use explicit domain and range maps for consistency across multiple
        // calls
        matrix_->fillComplete(map_, map_);

        if (!matrix_initialized_)
        {
            matrix_initialized_ = true;
            if (use_preconditioner_ && !preconditioner_.is_null())
            {
                preconditioner_->initialize();
                preconditioner_initialized_ = true;
            }
        }

        if (use_preconditioner_ && preconditioner_initialized_ &&
            !preconditioner_.is_null())
        {
            preconditioner_->compute();
        }
    }
};
#endif /* HAS_TRILINOS */

} /* namespace linearSolver */

#endif /* CONTEXTTRILINOS_H */
