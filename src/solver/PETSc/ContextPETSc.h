// File       : ContextPETSc.h
// Created    : Mon Mar 04 2024 12:33:50 (+0100)
// Author     : Fabian Wermelinger
// Description: Specialized PETSc solver context
// Copyright (c) 2024 CCFNUM, Lucerne University of Applied Sciences and Arts.
// SPDX-License-Identifier: BSD-3-Clause

#ifndef CONTEXTPETSC_H
#define CONTEXTPETSC_H

#include "linearSolverContext.h"
#include "memoryLayout.h"
#include <cassert>

#ifdef HAS_PETSC
#include <petscksp.h>
#include <petscpc.h>

#if PETSC_VERSION_LT(3, 18, 1)
#define ErrorWrapPetscCall(c) CHKERRQ(c)
#else
#define ErrorWrapPetscCall(c) PetscCall(c)
#endif
#endif /* HAS_PETSC */

namespace linearSolver
{

#ifdef HAS_PETSC
template <size_t N>
class ContextPETSc : public Context<N>
{
public:
    using Context<N>::BLOCKSIZE;
    using Coefficients = typename Context<N>::Coefficients;
    using Matrix = typename Coefficients::Matrix;
    using Vector = typename Coefficients::Vector;
    using Index = typename Coefficients::Index;
    using IndexVector = typename Coefficients::IndexVector;

    ContextPETSc() = delete;

    ContextPETSc(const std::string& family,
                 const std::string& system_name,
                 const CRSNodeGraph* graph,
                 const YAML::Node* node = nullptr)
        : Context<N>(linearSolver::ID::PETSc, family, system_name, graph, node),
          ksp_petsc_(nullptr), A_petsc_(nullptr), x_petsc_(nullptr),
          b_petsc_(nullptr)
    {
        setup_(node);
    }

    // used for preconditioning (shares coefficients)
    ContextPETSc(const std::string& family,
                 const std::string& system_name,
                 std::shared_ptr<Coefficients> coeffs,
                 const YAML::Node* node = nullptr)
        : Context<N>(linearSolver::ID::PETSc,
                     family,
                     system_name,
                     coeffs,
                     node),
          ksp_petsc_(nullptr), A_petsc_(nullptr), x_petsc_(nullptr),
          b_petsc_(nullptr)
    {
        setup_(node);
    }

    ~ContextPETSc()
    {
        destroySystem_();
        destroySolver_();
    }

    void initializeSystem()
    {
        this->createAIJ_();
        this->createVectors_();
    }

    void setSystem(Mat A, Vec b, Vec x)
    {
        A_petsc_ = A;
        b_petsc_ = b;
        x_petsc_ = x;
    }

    void destroySystem()
    {
        this->destroySystem_();
    }

    void setOperators(Mat A, Mat P)
    {
        KSPSetOperators(ksp_petsc_, A, P);
    }

    void solvePrologue(const int solver_id, const bool preconditioner) override
    {
        // default assumes matrix elements change for every iteration: matrix
        // elements must be explicitly updated and KSPSetOperators must be
        // called before every KSPSolve
        this->copyCoefficients_();
        this->setOperators(A_petsc_, A_petsc_);

        Context<N>::solvePrologue(solver_id, preconditioner);
    }

    std::string info(std::ostream& os,
                     const char* prefix = nullptr) const override
    {
        const std::string pfx = Context<N>::info(os, prefix);

        PC pc;
        PCType pc_type;
        KSPType ksp_type;
        KSPGetPC(ksp_petsc_, &pc);
        KSPGetType(ksp_petsc_, &ksp_type);
        PCGetType(pc, &pc_type);
        // clang-format off
        os << pfx << '\t' << "PETSc KSP type:   " << ksp_type << '\n';
        os << pfx << '\t' << "PETSc PC type:    " << pc_type << '\n';
        // clang-format on
        return pfx;
    }

    int solve()
    {
        ErrorWrapPetscCall(KSPSolve(ksp_petsc_, b_petsc_, x_petsc_));

        PetscInt solver_iterations;
        KSPGetIterationNumber(ksp_petsc_, &solver_iterations);
        return solver_iterations;
    }

    // clang-format off
    Mat &getAMatrixPETSc() { return A_petsc_; }
    const Mat &getAMatrixPETSc() const { return A_petsc_; }
    Vec &getXVectorPETSc() { return x_petsc_; }
    const Vec &getXVectorPETSc() const { return x_petsc_; }
    Vec &getBVectorPETSc() { return b_petsc_; }
    const Vec &getBVectorPETSc() const { return b_petsc_; }

    KSP &getKSPPETSc() { return ksp_petsc_; }
    const KSP &getKSPPETSc() const { return ksp_petsc_; }

    // clang-format on

protected:
    inline ContextPETSc* castContextPETSc_() override
    {
        return this;
    }

private:
    KSP ksp_petsc_;
    Mat A_petsc_;
    Vec x_petsc_;
    Vec b_petsc_;

    void setup_(const YAML::Node* solver_yaml); // at construction
    int setOptions_(const YAML::Node* solver_yaml) override;

    int clearOptions_()
    {
        ErrorWrapPetscCall(PetscOptionsClear(NULL));
        return 0;
    }

    int createKSP_();
    int createVectors_();
    int createAIJ_();
    int copyCoefficients_();

    int destroySystem_()
    {
        if (A_petsc_)
        {
            ErrorWrapPetscCall(MatDestroy(&A_petsc_));
        }
        if (b_petsc_)
        {
            ErrorWrapPetscCall(VecDestroy(&b_petsc_));
        }
        if (x_petsc_)
        {
            ErrorWrapPetscCall(VecDestroy(&x_petsc_));
        }
        A_petsc_ = nullptr;
        b_petsc_ = nullptr;
        x_petsc_ = nullptr;
        return 0;
    }

    int destroySolver_()
    {
        if (ksp_petsc_)
        {
            ErrorWrapPetscCall(KSPDestroy(&ksp_petsc_));
        }
        ksp_petsc_ = nullptr;
        return 0;
    }
};

template <size_t N>
void ContextPETSc<N>::setup_(const YAML::Node* solver_yaml)
{
    // start with empty options to avoid interference if multiple contexts
    // exist
    clearOptions_();
    if (solver_yaml)
    {
        setOptions_(solver_yaml);
    }

    this->initializeSystem();
    createKSP_();

#if 0
    PetscViewer viewer;
    PetscViewerASCIIOpen(this->getCommunicator(), "petsc_options.txt", &viewer);
    PetscOptionsView(NULL, viewer);
    this->cout() << *this;
#endif /* DEBUG */
}

template <size_t N>
int ContextPETSc<N>::setOptions_(const YAML::Node* solver_yaml)
{
    assert(solver_yaml);
    Context<N>::setOptions_(solver_yaml);

    using Iter = YAML::const_iterator;

    const YAML::Node& s = *solver_yaml;
    if (s["options"])
    {
        // optional PETSc specific options.  These options are identical to what
        // can be passed to the PetscOptionsSetValue() function (without the
        // leading hyphen `-`).
        const YAML::Node& petsc_options = s["options"];
        for (Iter it = petsc_options.begin(); it != petsc_options.end(); ++it)
        {
            std::string key =
                std::string("-") + it->first.template as<std::string>();
            std::string value = it->second.template as<std::string>();

            // set lower-case
            Context<N>::tolower(key);
            Context<N>::tolower(value);
            ErrorWrapPetscCall(
                PetscOptionsSetValue(NULL, key.c_str(), value.c_str()));
        }
    }
    return 0;
}

template <size_t N>
int ContextPETSc<N>::createKSP_()
{
    // Create and set the linear solver context (PETSc)
    ErrorWrapPetscCall(KSPCreate(this->getCommunicator(), &ksp_petsc_));
    ErrorWrapPetscCall(KSPSetFromOptions(ksp_petsc_));
    ErrorWrapPetscCall(KSPSetOperators(ksp_petsc_, A_petsc_, A_petsc_));
    ErrorWrapPetscCall(KSPSetTolerances(ksp_petsc_,
                                        this->rtol_,
                                        this->atol_,
                                        PETSC_DEFAULT,
                                        this->max_iterations_));
    return 0;
}

template <size_t N>
int ContextPETSc<N>::createVectors_()
{
    // storage
    Vector& x = this->getXVector();
    Vector& b = this->getBVector();

    // Create vectors x and b
    ErrorWrapPetscCall(MatCreateVecs(A_petsc_, &x_petsc_, &b_petsc_));
    ErrorWrapPetscCall(VecPlaceArray(x_petsc_, x.data()));
    ErrorWrapPetscCall(VecPlaceArray(b_petsc_, b.data()));
    return 0;
}

template <size_t N>
int ContextPETSc<N>::createAIJ_()
{
    const Matrix& A = this->getAMatrix();
    assert(A.getMemoryLayout() & GraphLayout::ColumnIndexOrder__Global);

    // Create a parallel AIJ matrix
    const Index n_global = this->coeffs_->nGlobalCoefficients();
    const Index n = A.nRows();
    std::vector<PetscInt> d_nnz(n * BLOCKSIZE, 0);
    std::vector<PetscInt> o_nnz(n * BLOCKSIZE, 0);
    for (Index i = 0; i < n; i++)
    {
        const Index row_nnz_owned = A.getGraph()->nnzOwned(i);
        const Index row_nnz_ghost = A.getGraph()->nnzGhost(i);
        for (Index k = 0; k < BLOCKSIZE; k++)
        {
            d_nnz[i * BLOCKSIZE + k] = BLOCKSIZE * row_nnz_owned;
            o_nnz[i * BLOCKSIZE + k] = BLOCKSIZE * row_nnz_ghost;
        }
    }

    // WIP: [2024-10-02] could use CSR structure identical
    // to native for BLOCKSIZE=1 that enables use of MatUpdateMPIAIJWithArray()
    // in copyCoefficients_ (not sure how much to benefit from this).
    ErrorWrapPetscCall(MatCreateAIJ(this->getCommunicator(),
                                    n * BLOCKSIZE,
                                    n * BLOCKSIZE,
                                    n_global * BLOCKSIZE,
                                    n_global * BLOCKSIZE,
                                    0,
                                    d_nnz.data(),
                                    0,
                                    o_nnz.data(),
                                    &A_petsc_));
    ErrorWrapPetscCall(MatSetFromOptions(A_petsc_));
    return 0;
}

template <size_t N>
int ContextPETSc<N>::copyCoefficients_()
{
    const Matrix& A = this->getAMatrix();
    // WIP: [2024-10-02] see comment in createAIJ_
    // ErrorWrapPetscCall(MatUpdateMPIAIJWithArray(A_petsc_, A.valuesPtr()));

    // TODO: [2024-10-14] assertion can be removed after
    // testing
    assert(A.getMemoryLayout() & GraphLayout::ColumnIndexOrder__Global);
    std::vector<Index> row_nnz;
    std::vector<Index> row_idx;
    std::vector<Index> col_idx;
    std::vector<typename Matrix::DataType> values;
    for (Index i = 0; i < A.nRows(); i++)
    {
        matrixLayout::blockRowToRowMajor(i,
                                         A,
                                         A.getGraph()->rowGlobalIndices(i),
                                         row_nnz,
                                         row_idx,
                                         col_idx,
                                         values);
        MatSetValues(A_petsc_,
                     BLOCKSIZE,
                     row_idx.data(),
                     row_nnz[0],
                     col_idx.data(),
                     values.data(),
                     INSERT_VALUES);
    }
    ErrorWrapPetscCall(MatAssemblyBegin(A_petsc_, MAT_FINAL_ASSEMBLY));
    ErrorWrapPetscCall(MatAssemblyEnd(A_petsc_, MAT_FINAL_ASSEMBLY));
    return 0;
}

#else

template <size_t N>
class ContextPETSc
{
};

#endif /* HAS_PETSC */

} /* namespace linearSolver */

#endif /* CONTEXTPETSC_H */
