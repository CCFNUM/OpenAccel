// File       : ContextUtopia.hpp
// Created    : Thu Oct 10 2024 09:42:32 (+0200)
// Author     : Fabian Wermelinger, Patrick Zulian
// Description: Fluya/Utopia FSI context interface
// Copyright 2024 CCFNUM HSLU T&A. All Rights Reserved.

// no include guard: must only be included in ContextUtopia.h

#include <cassert>

namespace accel
{
namespace fsi
{

template <size_t N>
ContextUtopia<N>::ContextUtopia(
    const std::shared_ptr<linearSolver::Context<N>> ctx)
    : linearSolver::Context<N>(ctx->getFamily() + " [ContextUtopia]",
                               ctx->getSystemName() + " [FSI]"),
      ctx_wrap_(ctx), ctx_amg_(nullptr), ctx_gmres_(nullptr),
      ctx_petsc_(nullptr), velocity_("velocity_structure_utopia"),
      fluid_velocity_("velocity_fluid_utopia"), integrator_(nullptr),
      valid_(false), restart_enabled_(false), convert_back_(false),
      utopia_verbose_(false), cache_hessian_nnz_(true),
      hessian_use_solid_diagonal_(false),
      skip_hessian_solid_contribution_(false), import_field_data_(false),
      inverse_diagonal_scaling_(false), use_strong_velocity_bc_(false),
      use_initial_guess_projection_(true), use_initial_guess_zero_(false),
      check_matrix_(false), reuse_matrix_(false), output_relative_(false),
      alpha_(1.0), scale_factor_(1.0), stol_(1.0e-3), norm_increment_(0.0)
{
    this->copyProtectedSettings(ctx_wrap_);
    this->setCoeffs(ctx_wrap_->coeffs());
    fluid_graph_ = this->coeffs_->getGraph();

    // FSI augmented node graph (only required when converting back to native
    // CRS layout when AMG solver is used for example -- currently assumes local
    // column index order)
    fsi_graph_ = std::make_unique<UtopiaFSIGraph<BLOCKSIZE>>(
        ctx_wrap_->getCommunicator(), fluid_graph_->getLayout());

    assert(fluid_graph_);
    assert(fsi_graph_);

    assert(this->solver_id_ < linearSolver::ID::UNDEFINED);
    switch (this->solver_id_)
    {
        case linearSolver::ID::AMG:
            convert_back_ = true; // convert back to native matrix format
            ctx_amg_ = *ctx_wrap_;
            break;
        case linearSolver::ID::GMRES:
        case linearSolver::ID::FlexGMRES:
            convert_back_ = true; // convert back to native matrix format
            ctx_gmres_ = *ctx_wrap_;
            assert(ctx_wrap_->getPreconditioner());
            break;
        case linearSolver::ID::PETSc:
            convert_back_ = false; // use PETSc matrix directly to solve system
            ctx_petsc_ = *ctx_wrap_;
            ctx_petsc_->destroySystem(); // context manages its own system
            break;
        default:
            errorMsg("ContextUtopia: linear solver `" +
                     std::string(linearSolver::SolverName[this->solver_id_]) +
                     "` for wrapped context is not supported");
    }

    // set preconditioner
    const auto ctx_precond = ctx_wrap_->getPreconditioner();
    if (ctx_precond)
    {
        switch (ctx_precond->solverID())
        {
            case linearSolver::ID::AMG:
                assert(ctx_amg_ == nullptr);
                ctx_amg_ = *ctx_precond;
                break;
            default:
                errorMsg(
                    "ContextUtopia: preconditioner `" +
                    std::string(
                        linearSolver::SolverName[ctx_precond->solverID()]) +
                    "` for wrapped context is not supported");
        }
    }

    // auxiliary data
    x_petsc_ = nullptr;
    b_petsc_ = nullptr;
    r_petsc_ = nullptr;
}

template <size_t N>
ContextUtopia<N>::~ContextUtopia()
{
    if (ctx_petsc_)
    {
        // restore initial system in PETSc context
        ctx_petsc_->initializeSystem();
    }
    if (x_petsc_)
    {
        VecDestroy(&x_petsc_);
    }
    if (b_petsc_)
    {
        VecDestroy(&b_petsc_);
    }
    if (r_petsc_)
    {
        VecDestroy(&r_petsc_);
    }
}

template <size_t N>
void ContextUtopia<N>::solvePrologue(const int solver_id,
                                     const bool preconditioner)
{
    if (ctx_petsc_)
    {
        VecZeroEntries(x_petsc_);
        ctx_petsc_->setOperators(hessian_fluid_->raw_type(),
                                 hessian_fluid_->raw_type());
        if (!preconditioner)
        {
            auto& ctrl = this->getControlData();
            MatResidual(
                hessian_fluid_->raw_type(), b_petsc_, x_petsc_, r_petsc_);

            const PetscScalar* R;
            VecGetArrayRead(r_petsc_, &R);
            ctrl.scaled_initial_res = this->PETScRMSErrorDiagNormalized_(R);
            ctrl.solver_initial_res = this->PETScRMSError_(R);
            VecRestoreArrayRead(r_petsc_, &R);
        }
    }
    else
    {
        ctx_wrap_->solvePrologue(solver_id, preconditioner);
    }
}

template <size_t N>
void ContextUtopia<N>::solveEpilogue(const int solver_id,
                                     const bool preconditioner)
{
    auto& ctrl = this->getControlData();
    if (ctx_petsc_ && !preconditioner)
    {
        MatResidual(hessian_fluid_->raw_type(), b_petsc_, x_petsc_, r_petsc_);

        const PetscScalar* R;
        VecGetArrayRead(r_petsc_, &R);
        ctrl.scaled_final_res = this->PETScRMSErrorDiagNormalized_(R);
        ctrl.solver_final_res = this->PETScRMSError_(R);
        VecRestoreArrayRead(r_petsc_, &R);
    }
    else
    {
        ctx_wrap_->solveEpilogue(solver_id, preconditioner);

        const auto& wrap_ctrl = ctx_wrap_->getControlData();
        ctrl.scaled_initial_res = wrap_ctrl.scaled_initial_res;
        ctrl.solver_initial_res = wrap_ctrl.solver_initial_res;
        ctrl.scaled_final_res = wrap_ctrl.scaled_final_res;
        ctrl.solver_final_res = wrap_ctrl.solver_final_res;
    }
}

template <size_t N>
void ContextUtopia<N>::updateSolution()
{
    UVector increment_solid;

    auto transfer = [this](Vec& x, UVector& increment_solid) -> scalar
    {
        UVector increment_fluid; // Utopia vector
        increment_fluid.wrap(x); // Vec x is a PETSc vector
        this->transfer_.apply(increment_fluid, increment_solid);
        increment_fluid.unwrap(x);
        return utopia::norm2(increment_solid);
    };

    if (!convert_back_)
    {
        assert(x_petsc_);
        norm_increment_ = transfer(x_petsc_, increment_solid);
    }
    else
    {
        Vec x;
        auto& coeff = this->getCoefficients();
        const Index n_owned = coeff.nCoefficients();
        if (1 == this->commSize())
        {
            VecCreateSeqWithArray(PETSC_COMM_SELF,
                                  1,
                                  BLOCKSIZE * n_owned,
                                  coeff.getXVector().data(),
                                  &x);
        }
        else
        {
            // WARNING: [fab4100@posteo.net; 2024-10-24] not further checked
            // when compiled for release (since there are plenty of numerical
            // 0.0 nnz's the provided space should be far sufficient right now).
            assert(fsi_graph_->isBuilt());
            assert(BLOCKSIZE * n_owned + fsi_graph_->getPETScNumberOfGhosts() <=
                   coeff.getXVector().size());
            VecCreateGhostWithArray(
                this->getCommunicator(),
                BLOCKSIZE * n_owned,
                BLOCKSIZE * coeff.nGlobalCoefficients(),
                fsi_graph_->getPETScNumberOfGhosts(),
                &fsi_graph_->getPETScGhostLocalToGlobalMap()[0],
                coeff.getXVector().data(),
                &x);
        }
        norm_increment_ = transfer(x, increment_solid);
        VecDestroy(&x);
    }

    if (this->norm_inc_stream_)
    {
        assert(sim_);
        auto& fout = *this->norm_inc_stream_;
        fout << sim_->getGlobalIterationCount();
        fout << '\t' << sim_->getIterationCount();
        fout << '\t' << std::scientific << sim_->getSimulationTime();
        fout << '\t' << std::scientific << norm_increment_;
        fout << std::endl; // flush
    }

    if (std::isnan(norm_increment_))
    {
        space_structure_.write("NaN.e", function_->solution());
        utopia::Utopia::Abort("Detected NaN");
    }

    auto ls = function_->line_search();
    if (ls)
    {
        UVector full_increment = increment_solid - diff_velocity_;
        ls->get_alpha(*function_,
                      gradient_solid_,
                      velocity_.data(),
                      full_increment,
                      alpha_);
    }

    velocity_.data() += alpha_ * increment_solid;
    velocity_.data() -= alpha_ * diff_velocity_;

    function_->project_onto_feasibile_region(velocity_.data());
    space_structure_.apply_constraints_time_derivative(velocity_.data());
}

template <size_t N>
void ContextUtopia<N>::synchronizeVelocity()
{
    transferFluidVelocity_();
}

template <size_t N>
void ContextUtopia<N>::resetFluid()
{
    if (!convert_back_)
    {
        return;
    }

    assert(fluid_graph_->isBuilt());
    assert(fluid_graph_->isLocalColumnOrder());
    assert(fluid_graph_->nOwnedNodes() ==
           this->coeffs_->getGraph()->nOwnedNodes());
    this->coeffs_->setGraph(fluid_graph_);
}

template <size_t N>
void ContextUtopia<N>::preTimeStep(const Scalar_t dt)
{
    UTOPIA_TRACE_SCOPE("ContextUtopia::preTimeStep");

    if (utopia_verbose_)
    {
        auto& cout = this->cout();
        const scalar t = function_->time()->get();
        cout << "ContextUtopia::preTimeStep: start: " << std::scientific << t
             << '\n';
        cout << "ContextUtopia::preTimeStep: dt:    " << std::scientific << dt
             << '\n';
        cout << "ContextUtopia::preTimeStep: end:   " << std::scientific
             << t + dt << '\n';

        // clang-format off
        cout << "norm_velocity:     " << std::scientific << utopia::norm2(velocity_.data()) << '\n';
        cout << "norm_displacement: " << std::scientific << utopia::norm2(function_->solution()) << '\n';
        cout << "nnz_solid:         " << std::scientific << hessian_solid_.global_nnz() << '\n';
        // clang-format on
    }

    hessian_fsi_.clear();
    stime_->delta() = dt;

    function_->update_IVP(velocity_.data());

    updateTransfer_();

    // UNCOMMENT FOR DBG
    // output_->update_output_field(fluid_velocity_);

    transferFluidVelocity_();

    if (!hessian_use_solid_diagonal_ && !skip_hessian_solid_contribution_)
    {
        hessian_fluid_ = utopia::make_unique<UMatrix>();
    }

    if (use_initial_guess_projection_)
    {
        function_->initial_guess_for_solver(velocity_.data());
    }
    else if (use_initial_guess_zero_)
    {
        velocity_.data().set(0.);
    }

    function_->update_BVP();
}

template <size_t N>
scalar ContextUtopia<N>::getAlpha() const
{
    return alpha_;
}

template <size_t N>
void ContextUtopia<N>::updateTransfer_()
{
    space_structure_.displace(function_->solution());
    if (!transfer_.init(utopia::make_ref(space_fluid_),
                        utopia::make_ref(space_structure_)))
    {
        utopia::Utopia::Abort("Failed to compute transfer!");
    }

    space_structure_.displace(-function_->solution());
}

template <size_t N>
void ContextUtopia<N>::updateSolidAndFSI_()
{
    function_->update(velocity_.data());
    if (!function_->hessian_and_gradient(
            velocity_.data(), hessian_solid_, gradient_solid_))
    {
        assert(false);
    }

    if (gradient_solid_.has_nan_or_inf())
    {
        space_structure_.write("NaN.e", function_->solution());
        utopia::Utopia::Abort("Detected NaN in solid gradient!");
    }

    diff_velocity_ = velocity_.data() - fluid_velocity_.data();
    gradient_solid_ -= hessian_solid_ * diff_velocity_;

    if (utopia_verbose_)
    {
        double norm_diff_velocity = utopia::norm2(diff_velocity_);
        this->cout() << "norm_diff_velocity: " << norm_diff_velocity << '\n';
    }
}

template <size_t N>
void ContextUtopia<N>::checkMatrix_()
{
    auto& comm = hessian_fluid_->comm();

    std::stringstream root_ss, ss;
    root_ss << "CHECK MATRIX\n";
    root_ss << "global_nnz:\t" << hessian_fluid_->global_nnz() << "\n";

    int size = comm.size();
    int rank = comm.rank();

    std::vector<long> local_nnz(size, 0);
    std::vector<long> numerical_nnz(size, 0);
    std::vector<long> uselss_nnz(size, 0);
    std::vector<long> fluya_nnz(size, 0);

    const auto& A = this->getAMatrix();

    local_nnz[rank] = hessian_fluid_->local_nnz();
    numerical_nnz[rank] = hessian_fluid_->nnz(0.);
    fluya_nnz[rank] = A.nnzBlocks();

    const auto* block = A.valuesPtr();
    for (Index i = 0; i < A.nnzBlocks(); i++)
    {
        bool is_zero = true;
        for (Index k = 0; k < BLOCKSIZE * BLOCKSIZE && is_zero; k++)
        {
            if (block[k] != 0.)
            {
                is_zero = false;
            }
        }
        block += BLOCKSIZE * BLOCKSIZE;
        uselss_nnz[rank] += is_zero;
    }

    comm.sum(local_nnz.size(), &local_nnz[0]);
    comm.sum(numerical_nnz.size(), &numerical_nnz[0]);
    comm.sum(fluya_nnz.size(), &fluya_nnz[0]);
    comm.sum(uselss_nnz.size(), &uselss_nnz[0]);

    if (rank == 0)
    {
        root_ss << "rank:\t\t";
        for (int r = 0; r < size; ++r)
        {
            root_ss << r << "\t";
        }

        root_ss << "\n";

        root_ss << "petsc_nnz:\t";
        for (int r = 0; r < size; ++r)
        {
            root_ss << local_nnz[r] << "\t";
        }

        root_ss << "\n";

        root_ss << "numerical_nnz:\t";
        for (int r = 0; r < size; ++r)
        {
            root_ss << numerical_nnz[r] << "\t";
        }

        root_ss << "\n";

        root_ss << "fluya_nnz:\t";
        for (int r = 0; r < size; ++r)
        {
            root_ss << fluya_nnz[r] << "\t";
        }

        root_ss << "\n";

        root_ss << "uselss_nnz:\t";
        for (int r = 0; r < size; ++r)
        {
            root_ss << uselss_nnz[r] << "\t";
        }

        root_ss << "\n";
    }

    comm.root_print(root_ss.str(), utopia::out().stream());
}

#ifdef DEBUG_GRAPH

template <size_t N>
void ContextUtopia<N>::testGraph_()
{
    constexpr Index n_blocks = 2;

    // 1.) setup test matrix
    const auto& coeff = this->getCoefficients();
    const Index n_owned = BLOCKSIZE * n_blocks; // every rank owns `n_blocks`
    assert(coeff.getGraph() == this->fluid_graph_);

    UMatrix gold;
    const auto ml = utopia::square_matrix_layout(
        utopia::layout(coeff.getCommunicator(),
                       n_owned,
                       utopia::Traits<UMatrix>::determine()));
    std::vector<Size> d_nnz(n_owned, 0);
    std::vector<Size> o_nnz(n_owned, 0);
    for (Index i = 0; i < n_owned; i++)
    {
        d_nnz[i] = 1; // one diagonal entry
        o_nnz[i] = 1; // one off-diagonal entry
    }
    gold.sparse(ml, d_nnz, o_nnz);

    {
        utopia::Write<UMatrix> scoped_write_lock(gold);

        const PetscInt n_global = n_owned * this->commSize();
        const PetscInt my_start = n_owned * this->commRank();
        std::vector<PetscInt> rows(1, 0);
        std::vector<PetscInt> cols;
        std::vector<PetscScalar> vals;
        for (Index i = 0; i < n_owned; i++)
        {
            rows[0] = i + my_start;
            cols.resize(d_nnz[i] + o_nnz[i]);
            vals.resize(d_nnz[i] + o_nnz[i]);

            // diagonal
            cols[0] = i + my_start;
            vals[0] = 2.0;

            // off-diagonal
            cols[1] = (i + my_start + n_owned) % n_global;
            vals[1] = 1.0;

            gold.set_matrix(rows, cols, vals);
        }
    }

    // Check for efficiency
    MatSetOption(gold.raw_type(), MAT_NEW_NONZERO_ALLOCATION_ERR, PETSC_TRUE);

#if 0
    // dump PETSc matrix:
    PetscViewer viewer;
    PetscViewerASCIIOpen(
        coeff.getCommunicator(), "petsc_gold.txt", &viewer);
    MatView(gold.raw_type(), viewer);
    PetscViewerDestroy(&viewer);
#endif

    // 2.) compute graph
    assert(fsi_graph_);
    fsi_graph_->buildGraph(gold);

    // 3.) convert matrix
    convertFromUtopia_(gold);

#if 0
    // dump native matrix:
    this->coeffs_->writeMatrix("native_test");
#endif

    // 4.) compute absolute difference of the matrices
    testAbsDiff_(gold);
}

template <size_t N>
void ContextUtopia<N>::testAbsDiff_(UMatrix& gold)
{
    UMatrix test;
    convertToUtopia_(test);

#if 0
    dump PETSc matrix:
    PetscViewer viewer;
    PetscViewerBinaryOpen(
        this->getCommunicator(), "petsc_gold.bin", FILE_MODE_WRITE, &viewer);
    MatView(gold.raw_type(), viewer);
    PetscViewerDestroy(&viewer);

    PetscViewerBinaryOpen(
        this->getCommunicator(), "petsc_test.bin", FILE_MODE_WRITE, &viewer);
    MatView(test.raw_type(), viewer);
    PetscViewerDestroy(&viewer);
#endif

    UMatrix diff = test - gold;
    const double diff_abs = utopia::sum(utopia::abs(diff));
    const double sum_test = utopia::sum(test);
    const double sum_gold = utopia::sum(gold);
    const double abs_sum_test = utopia::sum(utopia::abs(test));
    const double abs_sum_gold = utopia::sum(utopia::abs(gold));

    // nnz > 0
    const int nnz_local[2] = {(int)test.nnz(0.0), (int)gold.nnz(0.0)};
    int nnz_global[2] = {0};
    MPI_Reduce(
        nnz_local, nnz_global, 2, MPI_INT, MPI_SUM, 0, this->getCommunicator());

    // clang-format off
    auto& cout = this->cout();
    cout << '\n';
    cout << "abs(test - gold) = " << std::scientific << diff_abs << '\n';
    cout << "sum(test)        = " << std::scientific << sum_test << '\n';
    cout << "sum(gold)        = " << std::scientific << sum_gold << '\n';
    cout << "abs(test)        = " << std::scientific << abs_sum_test << '\n';
    cout << "abs(gold)        = " << std::scientific << abs_sum_gold << '\n';
    cout << "test.global_nnz  = " << test.global_nnz() << " (!=0: " << nnz_global[0] << "; ==0: " << test.global_nnz() - nnz_global[0] << ")\n";
    cout << "gold.global_nnz  = " << gold.global_nnz() << " (!=0: " << nnz_global[1] << "; ==0: " << gold.global_nnz() - nnz_global[1] << ")";
    cout << std::endl;
    MPI_Barrier(this->getCommunicator());
    for (int r = 0; r < this->commSize(); r++)
    {
        if (r == this->commRank())
        {
            this->coutAll() << r << ":test.local_nnz = " << test.local_nnz() << " (!=0: " << nnz_local[0] << "; ==0: " << test.local_nnz() - nnz_local[0] << ")\n";
            this->coutAll() << r << ":gold.local_nnz = " << gold.local_nnz() << " (!=0: " << nnz_local[1] << "; ==0: " << gold.local_nnz() - nnz_local[1] << ")";
            this->coutAll() << std::endl;
        }
        MPI_Barrier(this->getCommunicator());
    }
    MPI_Barrier(this->getCommunicator());
    cout << std::endl;
    // clang-format on

    MPI_Abort(this->getCommunicator(), 1);
}

#endif /* DEBUG_GRAPH */

} /* namespace fsi */
} /* namespace accel */
