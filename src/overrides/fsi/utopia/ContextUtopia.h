// File       : ContextUtopia.h
// Created    : Thu Oct 10 2024 09:42:32 (+0200)
// Author     : Fabian Wermelinger, Patrick Zulian
// Description: Fluya/Utopia FSI context interface
// Copyright 2024 CCFNUM HSLU T&A. All Rights Reserved.

#ifndef CONTEXTUTOPIA_H
#define CONTEXTUTOPIA_H

#ifdef HAS_UTOPIA

#include "UtopiaFSIGraph.h"
#include "coefficients.h"
#include "linearSolverContext.h"
#include "macros.h"
#include "mesh.h"
#include "simulation.h"
#include "velocity.h"

// Utopia
#include "utopia_make_unique.hpp"
#include "utopia_moonolith_stk_FETransfer.hpp"
#include "utopia_ui.hpp"

#include "utopia_stk_FunctionSpace.hpp"
#include "utopia_stk_SpaceIO.hpp"

#include "utopia_stk.hpp"
#include "utopia_stk_intrepid2.hpp"

#include "utopia_ExtractComponent.hpp"
#include "utopia_MPI.hpp"
#include "utopia_moonolith_stk_Contact.hpp"
#include "utopia_moonolith_stk_Obstacle.hpp"

#include "utopia_FEFunctionFactory.hpp"
#include "utopia_FSIHesch2014.hpp"

#include "utopia_NewmarkIntegrator.hpp"
#include "utopia_ObstacleStabilizedVelocityNewmark.hpp"
#include "utopia_ObstacleVelocityNewmark.hpp"
#include "utopia_moonolith_ObstacleFEFunctionFactory.hpp"

#include "utopia_Tracer.hpp"

#ifdef UTOPIA_ENABLE_YAML_CPP
#define FLUYA_UTOPIA_YAML_ENABLED
#include "utopia_YAMLInput.hpp"
#endif // UTOPIA_ENABLE_YAML_CPP

#ifdef HAS_PETSC
#include <petscmat.h>
#include <petscvec.h>
#endif /* HAS_PETSC */

#include <cassert>
#include <memory>
#include <string>

namespace accel
{

namespace fsi
{

template <size_t N>
class ContextUtopia : public utopia::Configurable,
                      public linearSolver::Context<N>
{
public:
    using FunctionSpace = typename utopia::stk::FunctionSpace;
    using UtopiaField = typename utopia::Field<FunctionSpace>;
    using IO = typename utopia::IO<FunctionSpace>;
    using FEFunction = typename utopia::TimeDependentFunction<FunctionSpace>;
    using UVector = typename utopia::UVector;
    using UMatrix = typename utopia::USparseMatrix;
    using Size = utopia::Traits<UVector>::SizeType;
    using Scalar_t = utopia::Traits<utopia::stk::FunctionSpace>::Scalar;
    using UtopiaTime = utopia::SimulationTime<Scalar_t>;

    using Array = typename linearSolver::Context<N>::Array;
    using Index = typename linearSolver::Context<N>::Index;
    using linearSolver::Context<N>::BLOCKSIZE;

    ContextUtopia(const std::shared_ptr<linearSolver::Context<N>> ctx);

    ~ContextUtopia();

    void init(::stk::mesh::BulkData& data,
              const YAML::Node& fsi_node,
              simulation& sim);
    void filterVelocityBoundaryNodeID(const mesh& mesh, velocity& U);
    void applyConstraints();
    void updateSolution();
    void synchronizeVelocity();
    void resetFluid();

    void preTimeStep(const Scalar_t dt);
    scalar getAlpha() const;

    inline scalar getNormIncrement() const
    {
        return norm_increment_;
    }

    inline bool isConverged() const
    {
        return norm_increment_ < stol_;
    }

    void solvePrologue(const int solver_id, const bool preconditioner) override;
    void solveEpilogue(const int solver_id, const bool preconditioner) override;

protected:
    linearSolver::ContextAMGSolver<N>* castContextAMG_() override
    {
        assert(ctx_amg_);
        return ctx_amg_;
    }

    linearSolver::ContextGMRES<N>* castContextGMRES_() override
    {
        assert(ctx_gmres_);
        return ctx_gmres_;
    }

    linearSolver::ContextPETSc<N>* castContextPETSc_() override
    {
        assert(ctx_petsc_);
        return ctx_petsc_;
    }

private:
    FunctionSpace space_structure_;
    FunctionSpace space_fluid_;
    utopia::stk::FETransfer transfer_;

    // FEM function
    std::unique_ptr<FEFunction> function_;

    // Utopia matrices
    std::unique_ptr<UMatrix> hessian_fluid_;
    UMatrix hessian_solid_;
    UMatrix hessian_fsi_;
    // local node index for which to enforce zero Dirichlet velocity BC
    std::vector<unsigned int> zero_dirichlet_vel_idx_;

    // wrapped and supported solver context
    const std::shared_ptr<linearSolver::Context<N>> ctx_wrap_;
    linearSolver::ContextAMGSolver<BLOCKSIZE>* ctx_amg_;
    linearSolver::ContextGMRES<BLOCKSIZE>* ctx_gmres_;
    linearSolver::ContextPETSc<BLOCKSIZE>* ctx_petsc_;

    // fluid and FSI augmented graphs
    const linearSolver::CRSNodeGraph* fluid_graph_;
    std::unique_ptr<UtopiaFSIGraph<BLOCKSIZE>> fsi_graph_;

    // stream for norm_increment_ output file
    const simulation* sim_{nullptr};
    std::unique_ptr<std::ostream> norm_inc_stream_{nullptr};

    // Utopia auxiliary data structures and settings
    UtopiaField velocity_;
    UtopiaField fluid_velocity_;
    UVector diff_velocity_;
    UVector gradient_solid_;
    std::unique_ptr<IO> io_;
    std::unique_ptr<IO> io_restart_;
    std::unique_ptr<IO> restart_db_;
    std::string structure_filename_;

    // used for direct IO access only
    utopia::NewmarkIntegrator<FunctionSpace>* integrator_;

    bool valid_;
    bool restart_enabled_;
    bool convert_back_; // convert back from Utopia (PETSc) matrix to native
    bool utopia_verbose_;
    bool cache_hessian_nnz_;
    bool hessian_use_solid_diagonal_;
    bool skip_hessian_solid_contribution_;
    bool import_field_data_;
    bool inverse_diagonal_scaling_;
    bool use_strong_velocity_bc_;
    bool use_initial_guess_projection_;
    bool use_initial_guess_zero_;
    bool check_matrix_;
    bool reuse_matrix_;
    bool output_relative_;

    scalar alpha_;
    scalar scale_factor_;
    scalar stol_;
    scalar norm_increment_;
    std::shared_ptr<UtopiaTime> stime_;

    // helper PETSc vectors (their utilization is dependent on wrapped ctx)
    Vec x_petsc_;
    Vec b_petsc_;
    Vec r_petsc_;

    void read(utopia::Input& in) override final;

    void updateTransfer_();
    void updateSolidAndFSI_();
    void transferFluidVelocity_();
    bool coupleStructureWithFluid_();
    bool coupleStructureWithFluid_(UVector& velocity_in);

    void checkMatrix_();

    void convertToUtopia_(UMatrix& dst);
    void convertFromUtopia_(UMatrix& src);

    void convertPETScToNative_(const Mat& mat_seq,
                               const bool diagonal_part = true);
    Array PETScRMSError_(const PetscScalar* R);
    Array PETScRMSErrorDiagNormalized_(const PetscScalar* R);

    // post transformations on linear system in Utopia PETSc representation
    void PETScApplyStrongBC_();

#ifdef DEBUG_GRAPH
    void testGraph_();
    void testAbsDiff_(UMatrix& gold);
#endif /* DEBUG_GRAPH */
};

} /* namespace fsi */
} /* namespace accel */

#include "ContextUtopia.hpp"
#include "ContextUtopiaFSI.hpp"
#include "ContextUtopiaInit.hpp"
#include "ContextUtopiaPETSc.hpp"

#undef DEBUG_GRAPH

#endif /* HAS_UTOPIA */

#endif // CONTEXTUTOPIA_H
