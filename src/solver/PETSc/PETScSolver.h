// File       : PETScSolver.h
// Created    : Mon Jan 22 2024 13:22:13 (+0100)
// Author     : Fabian Wermelinger
// Description: PETSc based linear solver wrapper
// Copyright (c) 2024 CCFNUM, Lucerne University of Applied Sciences and Arts.
// SPDX-License-Identifier: BSD-3-Clause

#ifndef PETSC_H
#define PETSC_H

#define DEBUG_PETSC 0

#include "ContextPETSc.h"
#include "linearSolver.h"
#include <cassert>
#include <mpi.h>
#include <yaml-cpp/yaml.h>

#ifdef HAS_PETSC
#include <petscksp.h>
#include <petscpc.h>
#endif /* HAS_PETSC */

namespace linearSolver
{

template <size_t N>
class PETSc : public ::linearSolver::Base<N>
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
    PETSc(const YAML::Node& node, const bool precond = false)
        : Base(ID::PETSc, precond), yaml_conf_(node), ctx_petsc_(nullptr)
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
    ContextPETSc<N>* ctx_petsc_;

#ifdef HAS_PETSC
    int writePETScMatrix_(const MPI_Comm comm,
                          const std::string basename,
                          const Mat& A) const;
    int writePETScVector_(const MPI_Comm comm,
                          const std::string basename,
                          const Vec& v) const;
#endif /* HAS_PETSC */
};

template <size_t N>
std::shared_ptr<typename PETSc<N>::Context>
PETSc<N>::createContext(const std::string& system_name,
                        const CRSNodeGraph* graph)
{
#ifdef HAS_PETSC
    auto ctx = std::make_shared<ContextPETSc<N>>(
        "PETSc", system_name, graph, &yaml_conf_);
    ctx_ = ctx;
    ctx_petsc_ = ctx.get();
#endif /* HAS_PETSC */
    return ctx_;
}

template <size_t N>
std::shared_ptr<typename PETSc<N>::Context>
PETSc<N>::createContext(const std::string& system_name,
                        std::shared_ptr<Coefficients> coeffs)
{
#ifdef HAS_PETSC
    auto ctx = std::make_shared<ContextPETSc<N>>(
        "PETSc", system_name, coeffs, &yaml_conf_);
    ctx_ = ctx;
    ctx_petsc_ = ctx.get();
#endif /* HAS_PETSC */
    return ctx_;
}

template <size_t N>
int PETSc<N>::solve()
{
#ifdef HAS_PETSC
    assert(ctx_);
    assert(ctx_petsc_);

    if (ctx_.get() != ctx_petsc_)
    {
        ctx_petsc_ = *ctx_;
        assert(ctx_petsc_);
    }
    assert(ctx_->coeffs() == ctx_petsc_->coeffs());

    ctx_->solvePrologue(this->getID(), this->isPreconditioner());

#if DEBUG_PETSC
    writePETScMatrix_(ctx_petsc_->getCommunicator(),
                      "A_petsc",
                      ctx_petsc_->getAMatrixPETSc());
    writePETScVector_(ctx_petsc_->getCommunicator(),
                      "b_petsc",
                      ctx_petsc_->getBVectorPETSc());
#endif /* DEBUG */

    // Solve the system Ax = b
    ctx_->getControlData().n_iterations = ctx_petsc_->solve();

#if DEBUG_PETSC
    writePETScVector_(ctx_petsc_->getCommunicator(),
                      "x_petsc",
                      ctx_petsc_->getXVectorPETSc());
#endif /* DEBUG */

    ctx_->solveEpilogue(this->getID(), this->isPreconditioner());

    ++(*ctx_); // increment solver call count
#endif         /* HAS_PETSC */

    return 0;
}

#ifdef HAS_PETSC
template <size_t N>
int PETSc<N>::writePETScMatrix_(const MPI_Comm comm,
                                const std::string basename,
                                const Mat& A) const
{
    PetscViewer viewer;
    const std::string fout =
        basename + "_" + std::to_string(reinterpret_cast<size_t>(ctx_petsc_)) +
        std::string(".bin");
    ErrorWrapPetscCall(
        PetscViewerBinaryOpen(comm, fout.c_str(), FILE_MODE_WRITE, &viewer));
    ErrorWrapPetscCall(MatView(A, viewer));
    ErrorWrapPetscCall(PetscViewerDestroy(&viewer));
    return 0;
}

template <size_t N>
int PETSc<N>::writePETScVector_(const MPI_Comm comm,
                                const std::string basename,
                                const Vec& v) const
{
    PetscViewer viewer;
    const std::string fout =
        basename + "_" + std::to_string(reinterpret_cast<size_t>(ctx_petsc_)) +
        std::string(".bin");
    ErrorWrapPetscCall(
        PetscViewerBinaryOpen(comm, fout.c_str(), FILE_MODE_WRITE, &viewer));
    ErrorWrapPetscCall(VecView(v, viewer));
    ErrorWrapPetscCall(PetscViewerDestroy(&viewer));
    return 0;
}
#endif /* HAS_PETSC */

} /* namespace linearSolver */

#undef DEBUG_PETSC

#endif /* PETSC_H */
