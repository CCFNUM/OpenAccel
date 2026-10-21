// File       : ContextUtopiaFSI.hpp
// Created    : Mon Oct 14 2024 10:07:48 (+0200)
// Author     : Fabian Wermelinger
// Description: Utopia FSI interface implementation
// Copyright 2024 CCFNUM HSLU T&A. All Rights Reserved.

// no include guard: must only be included in ContextUtopia.h

#include "memoryLayout.h"
#include <cassert>
#include <vector>

#ifdef HAS_PETSC
#include <petscmat.h>
#endif /* HAS_PETSC */

namespace accel
{
namespace fsi
{

template <size_t N>
void ContextUtopia<N>::convertToUtopia_(UMatrix& dst)
{
    UTOPIA_TRACE_REGION_BEGIN("ContextUtopia::convertToUtopia_");

    const auto& coeff = this->getCoefficients();
    const Index n_owned = coeff.nCoefficients();

#ifndef DEBUG_GRAPH
    assert(coeff.getGraph() == this->fluid_graph_);
#endif // DEBUG_GRAPH
    if (dst.empty())
    {
        const auto ml = utopia::square_matrix_layout(
            utopia::layout(coeff.getCommunicator(),
                           BLOCKSIZE * n_owned,
                           utopia::Traits<UMatrix>::determine()));
        std::vector<Size> d_nnz(BLOCKSIZE * n_owned);
        std::vector<Size> o_nnz(BLOCKSIZE * n_owned);
        for (Index i = 0; i < n_owned; i++)
        {
            const Index nnz_owned = coeff.getGraph()->nnzOwned(i);
            const Index nnz_ghost = coeff.getGraph()->nnzGhost(i);
            for (Index k = 0; k < BLOCKSIZE; k++)
            {
                d_nnz[BLOCKSIZE * i + k] = BLOCKSIZE * nnz_owned;
                o_nnz[BLOCKSIZE * i + k] = BLOCKSIZE * nnz_ghost;
            }
        }
        dst.sparse(ml, d_nnz, o_nnz);

        // Check for efficiency
        MatSetOption(
            dst.raw_type(), MAT_NEW_NONZERO_ALLOCATION_ERR, PETSC_TRUE);
    }
    else
    {
        dst *= 0.0;
    }

    {
        // calls MatAssemblyBegin/End on destructor invocation
        utopia::Write<UMatrix> scoped_write_lock(dst);

        const auto& A = coeff.getAMatrix();
        std::vector<Index> row_nnz;
        std::vector<Index> row_idx;
        std::vector<Index> col_idx;
        std::vector<scalar> values;
        for (Index i = 0; i < n_owned; i++)
        {
            linearSolver::matrixLayout::blockRowToRowMajor(
                i,
                A,
                A.getGraph()->rowGlobalIndices(i),
                row_nnz,
                row_idx,
                col_idx,
                values);
            // TODO: [fab4100@posteo.net; 2024-10-15] could use
            // dst.set_matrix() calls -- col_idx must be one row of
            // column indices only! (Does that work with HYPRE?)
            // FIXME: [fab4100@posteo.net; 2024-10-18] use direct copy approach
            // similar to one used in copy back operation
            MatSetValues(dst.raw_type(),
                         BLOCKSIZE,
                         row_idx.data(),
                         row_nnz[0],
                         col_idx.data(),
                         values.data(),
                         INSERT_VALUES);
        }
    }

    UTOPIA_TRACE_REGION_END("ContextUtopia::convertToUtopia_");
}

template <size_t N>
void ContextUtopia<N>::convertFromUtopia_(UMatrix& src)
{
    if (!convert_back_)
    {
        return;
    }

    UTOPIA_TRACE_REGION_BEGIN("ContextUtopia::convertFromUtopia_");

    assert(fsi_graph_->isBuilt());
    assert(fsi_graph_->isLocalColumnOrder());
#ifndef DEBUG_GRAPH
    assert(fsi_graph_->nOwnedNodes() ==
           this->coeffs_->getGraph()->nOwnedNodes());
#endif // DEBUG_GRAPH
    this->coeffs_->setGraph(fsi_graph_.get());
    this->coeffs_->zeroLHS();
    this->coeffs_->zeroSOL();

    Mat d_petsc, o_petsc;
    const PetscInt* o_local_to_global = nullptr;
    if (1 == this->commSize())
    {
        d_petsc = src.raw_type();
    }
    else
    {
        MatMPIAIJGetSeqAIJ(
            src.raw_type(), &d_petsc, &o_petsc, &o_local_to_global);
    }

    convertPETScToNative_(d_petsc);
    if (o_local_to_global)
    {
        convertPETScToNative_(o_petsc, false);
    }

    UTOPIA_TRACE_REGION_END("ContextUtopia::convertFromUtopia_");
}

template <size_t N>
void ContextUtopia<N>::applyConstraints()
{
    updateSolidAndFSI_();

    UVector inv_diag; // inverse diagonal scaling

    ////////////////////////////////////////////////////////////////////////////
    // convert native matrix to Utopia representation for transfer operator
    if (!hessian_fluid_)
    {
        hessian_fluid_ = utopia::make_unique<UMatrix>();
    }
    this->convertToUtopia_(*hessian_fluid_);
    if (check_matrix_)
    {
        this->checkMatrix_();
    }

    UTOPIA_TRACE_REGION_BEGIN("ContextUtopia::applyConstraints Hessian");
    {
        if (!hessian_fsi_.empty() && reuse_matrix_)
        {
            transfer_.apply(hessian_solid_, hessian_fsi_, true);
            hessian_fluid_->axpy_subset(1.0, hessian_fsi_);
        }
        else
        {
            transfer_.apply(hessian_solid_, hessian_fsi_);
            *hessian_fluid_ += hessian_fsi_;
            if (ctx_petsc_)
            {
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
                MatCreateVecs(hessian_fluid_->raw_type(), &x_petsc_, &b_petsc_);
                MatCreateVecs(hessian_fluid_->raw_type(), nullptr, &r_petsc_);
                VecPlaceArray(x_petsc_, this->getXVector().data());
                VecPlaceArray(b_petsc_, this->getBVector().data());
                VecPlaceArray(r_petsc_, this->getRVector().data());
                ctx_petsc_->setSystem(
                    hessian_fluid_->raw_type(), b_petsc_, x_petsc_);
            }
            if (convert_back_)
            {
                if (b_petsc_)
                {
                    VecDestroy(&b_petsc_);
                }
                MatCreateVecs(hessian_fluid_->raw_type(), nullptr, &b_petsc_);
                VecPlaceArray(b_petsc_, this->getBVector().data());
                fsi_graph_->buildGraph(*hessian_fluid_);
#if defined(DEBUG_GRAPH) && 0
                testGraph_();
#endif // DEBUG_GRAPH
            }
        }

        if (!reuse_matrix_)
        {
            // If not reused clear space
            hessian_fsi_.clear();
        }

        if (inverse_diagonal_scaling_)
        {
            inv_diag = diag(*hessian_fluid_);
            inv_diag = 1. / inv_diag;
            hessian_fluid_->diag_scale_left(inv_diag);
        }
    }
    UTOPIA_TRACE_REGION_END("ContextUtopia::applyConstraints Hessian");

    UTOPIA_TRACE_REGION_BEGIN("ContextUtopia::applyConstraints gradient");
    {
        UVector gradient_fsi;
        UVector negative_gradient_fluid;
        negative_gradient_fluid.wrap(b_petsc_);

        if (negative_gradient_fluid.has_nan_or_inf())
        {
            space_structure_.write("NaN.e", function_->solution());
            utopia::Utopia::Abort("Detected NaN in fluid gradient!");
        }

        transfer_.apply_transpose(gradient_solid_, gradient_fsi);
        negative_gradient_fluid -= gradient_fsi;

        if (utopia_verbose_)
        {
            auto& cout = this->cout();
            cout << "norm_grad: " << utopia::norm2(gradient_solid_) << '\n';
            cout << "norm_grad_fsi: " << utopia::norm2(gradient_fsi) << '\n';
        }

        if (inverse_diagonal_scaling_)
        {
            negative_gradient_fluid = e_mul(negative_gradient_fluid, inv_diag);
        }

        negative_gradient_fluid.unwrap(b_petsc_);
    }
    UTOPIA_TRACE_REGION_END("ContextUtopia::applyConstraints gradient");

    // apply post linear system transformations in Utopia representation here
    PETScApplyStrongBC_();

    ////////////////////////////////////////////////////////////////////////////
    // convert back from Utopia space to local representation
    this->convertFromUtopia_(*hessian_fluid_);
#ifdef DEBUG_GRAPH
    testAbsDiff_(*hessian_fluid_);
#endif // DEBUG_GRAPH
    if (check_matrix_)
    {
        this->checkMatrix_();
    }
    if (!cache_hessian_nnz_)
    {
        hessian_fluid_ = nullptr;
    }
}

template <size_t N>
void ContextUtopia<N>::transferFluidVelocity_()
{
    UVector temp;
    space_fluid_.create_vector(temp);
    space_fluid_.nodal_field_to_global_vector(temp);
    transfer_.apply(temp, fluid_velocity_.data());

    // Since we have zeros where the fluid is not coupled with the solid
    // due the Dirichlet conditions we need to set the constrain value here
    space_structure_.apply_constraints_time_derivative(fluid_velocity_.data());
}

template <size_t N>
bool ContextUtopia<N>::coupleStructureWithFluid_()
{
    UVector indicator_structure;
    UtopiaField indicator("structure_indicator");
    space_fluid_.create_nodal_vector_field(1, indicator);

    space_structure_.create_vector(indicator_structure);
    indicator_structure.set(1.0);

    transfer_.apply(indicator_structure, indicator.data());

    std::vector<int> indices;
    indices.clear();
    indices.reserve(indicator.data().local_size());

    const int n_var = space_fluid_.n_var();

    if (n_var != 1)
    {
        UVector indicator_x;
        utopia::extract_component(indicator.data(), n_var, 0, indicator_x);
        indicator.data() = std::move(indicator_x);
    }

    auto rr = range(indicator.data());
    utopia::Read<UVector> r(indicator.data());
    for (auto i = rr.begin(); i < rr.end(); ++i)
    {
        if (utopia::approxeq(indicator.data().get(i), 1.0, 1e-3))
        {
            indices.push_back(i);
        }
    }

    UtopiaField velocity("structure_velocity");
    space_fluid_.create_nodal_vector_field(
        space_fluid_.mesh().spatial_dimension(), velocity);

    velocity.data().set(0.0);
    space_fluid_.backend_set_nodal_field(velocity);
    space_fluid_.backend_set_nodal_field(indicator);
    return true;
}

template <size_t N>
bool ContextUtopia<N>::coupleStructureWithFluid_(UVector& velocity_in)
{
    UVector indicator_structure;
    UtopiaField indicator("structure_indicator");
    space_fluid_.create_nodal_vector_field(1, indicator);

    space_structure_.create_vector(indicator_structure);
    indicator_structure.set(1.0);

    transfer_.apply(indicator_structure, indicator.data());

    std::vector<int> indices;
    indices.clear();
    indices.reserve(indicator.data().local_size());

    const int n_var = space_fluid_.n_var();
    assert(n_var == 3);

    UVector indicator_x;
    utopia::extract_component(indicator.data(), n_var, 0, indicator_x);
    indicator.data() = std::move(indicator_x);

    auto rr = range(indicator.data());
    utopia::Read<UVector> r(indicator.data());
    for (auto i = rr.begin(); i < rr.end(); ++i)
    {
        if (utopia::approxeq(indicator.data().get(i), 1.0, 1e-3))
        {
            indices.push_back(i);
        }
    }

    UtopiaField velocity("structure_velocity");
    space_fluid_.create_nodal_vector_field(
        space_fluid_.mesh().spatial_dimension(), velocity);

    transfer_.apply(velocity_in, velocity.data());
    space_fluid_.backend_set_nodal_field(velocity);

    if (utopia_verbose_)
    {
        std::vector<double> norms_structure, norms_fluid;
        UtopiaField s_velocity(
            "v", make_ref(space_structure_), make_ref(velocity_in));
        s_velocity.set_tensor_size(3);

        assert(velocity.tensor_size() == 3);

        // double s_s = sum(velocity_in);
        // double s_f = sum(velocity.data());

        utopia::l2_norm(s_velocity, norms_structure);
        utopia::l2_norm(velocity, norms_fluid);

        auto& cout = this->cout();
        // utopia::out() << s_s << " -> " << s_f << '\n';
        cout << "ObstacleBase: l2 norm velocity: (" << norms_structure[0]
             << " -> " << norms_fluid[0] << ", " << norms_structure[1] << " -> "
             << norms_fluid[1] << ", " << norms_structure[2] << " -> "
             << norms_fluid[2] << ")\n";
    }

    space_fluid_.backend_set_nodal_field(indicator);
    return true;
}

} /* namespace fsi */
} /* namespace accel */
