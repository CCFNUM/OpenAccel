// File : solidDisplacementAssemblerElemTerms.cpp
// Created : Thu Dec 04 2025 10:15:23 (+0100)
// Author : Mhamad Mahdi Alloush
// Description: Element interior terms for linear elastic stress
// Copyright (c) 2025 CCFNUM, Lucerne University of Applied Sciences and Arts.
// SPDX-License-Identifier: BSD-3-Clause

#include "solidDisplacementAssembler.h"

namespace accel
{

void solidDisplacementAssembler::assembleElemTermsInterior_(
    const domain* domain,
    Context* ctx)
{
#ifdef USE_CVFEM_SOLID_MECHANICS
    const auto& mesh = field_broker_->meshRef();
    Matrix& A = ctx->getAMatrix();
    Vector& b = ctx->getBVector();

    const stk::mesh::BulkData& bulkData = mesh.bulkDataRef();
    const stk::mesh::MetaData& metaData = mesh.metaDataRef();

    // Workspace for LHS/RHS
    std::vector<scalar> lhs;
    std::vector<scalar> rhs;
    std::vector<label> scratchIds;
    std::vector<scalar> scratchVals;
    std::vector<stk::mesh::Entity> connectedNodes;

    // Nodal fields to gather
    std::vector<scalar> ws_coordinates;
    std::vector<scalar> ws_D;
    std::vector<scalar> ws_E;
    std::vector<scalar> ws_nu;

    // Geometry related
    std::vector<scalar> ws_scs_areav;
    std::vector<scalar> ws_dndx;
    std::vector<scalar> ws_deriv;
    std::vector<scalar> ws_det_j;
    std::vector<scalar> ws_shape_function;

    // Workspace for neo-hookean model (only used if simplifiedNeoHookean)
    std::vector<scalar> ws_gradTens;

    // const-size containers
    std::vector<scalar> kd(SPATIAL_DIM * SPATIAL_DIM);

    // fill Kronecker delta
    for (label i = 0; i < SPATIAL_DIM; ++i)
    {
        for (label j = 0; j < SPATIAL_DIM; ++j)
        {
            kd[i * SPATIAL_DIM + j] = (i == j) ? 1.0 : 0.0;
        }
    }

    // pointers to const-size arrays
    scalar* p_kd = &kd[0];

    // Get fields
    const auto& DSTKFieldRef = phi_->stkFieldRef();
    const auto& ESTKFieldRef = model_->ERef().stkFieldRef();
    const auto& nuSTKFieldRef = model_->nuRef().stkFieldRef();
    const auto& coordinatesRef = *metaData.get_field<scalar>(
        stk::topology::NODE_RANK, this->getCoordinatesID_(domain));

    // Get interior parts
    const stk::mesh::PartVector& partVec = domain->zonePtr()->interiorParts();
    const stk::mesh::Selector selAllElements =
        metaData.universal_part() & stk::mesh::selectUnion(partVec);

    // Check if plane stress or plane strain
    const bool planeStress = domain->solidMechanics_.planeStress_;

    // Get the solid mechanics option (linear elastic vs simplified neo-hookean)
    const solidMechanicsOption solidMechOption =
        domain->solidMechanics_.option_;

    // Shifted integration points?
    const bool isShifted = phi_->isShifted();
    const bool isGradientShifted = phi_->isGradientShifted();

    stk::mesh::BucketVector const& elementBuckets =
        bulkData.get_buckets(stk::topology::ELEMENT_RANK, selAllElements);

    for (stk::mesh::BucketVector::const_iterator ib = elementBuckets.begin();
         ib != elementBuckets.end();
         ++ib)
    {
        stk::mesh::Bucket& elementBucket = **ib;
        const stk::mesh::Bucket::size_type nElementsPerBucket =
            elementBucket.size();

        // Extract master element
        MasterElement* meSCS =
            accel::MasterElementRepo::get_surface_master_element(
                elementBucket.topology());

        const label nodesPerElement = meSCS->nodesPerElement_;
        const label numScsIp = meSCS->numIntPoints_;
        const label* lrscv = meSCS->adjacentNodes();

        // Resize arrays
        const label lhsSize =
            nodesPerElement * SPATIAL_DIM * nodesPerElement * SPATIAL_DIM;
        const label rhsSize = nodesPerElement * SPATIAL_DIM;
        lhs.resize(lhsSize);
        rhs.resize(rhsSize);
        scratchIds.resize(rhsSize);
        scratchVals.resize(rhsSize);
        connectedNodes.resize(nodesPerElement);

        ws_coordinates.resize(nodesPerElement * SPATIAL_DIM);
        ws_D.resize(nodesPerElement * SPATIAL_DIM);
        ws_E.resize(nodesPerElement);
        ws_nu.resize(nodesPerElement);
        ws_scs_areav.resize(numScsIp * SPATIAL_DIM);
        ws_dndx.resize(SPATIAL_DIM * numScsIp * nodesPerElement);
        ws_deriv.resize(SPATIAL_DIM * numScsIp * nodesPerElement);
        ws_det_j.resize(numScsIp);
        ws_shape_function.resize(numScsIp * nodesPerElement);

        // Resize workspace for neo-hookean model if needed
        if (solidMechOption == solidMechanicsOption::simplifiedNeoHookean)
        {
            ws_gradTens.resize(SPATIAL_DIM * SPATIAL_DIM);
        }

        // Pointers
        scalar* p_lhs = &lhs[0];
        scalar* p_rhs = &rhs[0];
        scalar* p_coordinates = &ws_coordinates[0];
        scalar* p_D = &ws_D[0];
        scalar* p_E = &ws_E[0];
        scalar* p_nu = &ws_nu[0];
        scalar* p_scs_areav = &ws_scs_areav[0];
        scalar* p_dndx = &ws_dndx[0];
        scalar* p_shape_function = &ws_shape_function[0];
        scalar* p_gradTens =
            solidMechOption == solidMechanicsOption::simplifiedNeoHookean
                ? &ws_gradTens[0]
                : nullptr;

        // Extract shape functions for stress terms
        if (isShifted)
        {
            meSCS->shifted_shape_fcn(p_shape_function);
        }
        else
        {
            meSCS->shape_fcn(p_shape_function);
        }

        for (stk::mesh::Bucket::size_type iElement = 0;
             iElement < nElementsPerBucket;
             ++iElement)
        {
            stk::mesh::Entity elem = elementBucket[iElement];

            // Zero lhs/rhs
            for (label p = 0; p < lhsSize; ++p)
                p_lhs[p] = 0.0;
            for (label p = 0; p < rhsSize; ++p)
                p_rhs[p] = 0.0;

            // Gather nodal data
            stk::mesh::Entity const* nodeRels = bulkData.begin_nodes(elem);
            label numNodes = bulkData.num_nodes(elem);
            STK_ThrowAssert(numNodes == nodesPerElement);

            for (label ni = 0; ni < numNodes; ++ni)
            {
                stk::mesh::Entity node = nodeRels[ni];
                connectedNodes[ni] = node;

                const scalar* coords =
                    stk::mesh::field_data(coordinatesRef, node);
                const scalar* D = stk::mesh::field_data(DSTKFieldRef, node);

                p_E[ni] = *stk::mesh::field_data(ESTKFieldRef, node);
                p_nu[ni] = *stk::mesh::field_data(nuSTKFieldRef, node);

                for (label i = 0; i < SPATIAL_DIM; ++i)
                {
                    p_coordinates[ni * SPATIAL_DIM + i] = coords[i];
                    p_D[ni * SPATIAL_DIM + i] = D[i];
                }
            }

            // Compute geometry
            scalar scs_error = 0.0;
            meSCS->determinant(1, p_coordinates, p_scs_areav, &scs_error);

            // compute dndx
            if (isGradientShifted)
            {
                meSCS->shifted_grad_op(1,
                                       &ws_coordinates[0],
                                       &ws_dndx[0],
                                       &ws_deriv[0],
                                       &ws_det_j[0],
                                       &scs_error);
            }
            else
            {
                meSCS->grad_op(1,
                               &ws_coordinates[0],
                               &ws_dndx[0],
                               &ws_deriv[0],
                               &ws_det_j[0],
                               &scs_error);
            }

            // Loop over integration points
            for (label ip = 0; ip < numScsIp; ++ip)
            {
                const label il = lrscv[2 * ip];
                const label ir = lrscv[2 * ip + 1];

                // Interpolate Lame parameters to integration point
                scalar muIp = 0.0;
                scalar lambdaIp = 0.0;
                const label offSetSF = ip * nodesPerElement;

                for (label ic = 0; ic < nodesPerElement; ++ic)
                {
                    const scalar r = p_shape_function[offSetSF + ic];
                    const scalar E = p_E[ic];
                    const scalar nu = p_nu[ic];

                    // Compute Lame parameters from E and nu
                    const scalar mu = E / (2.0 * (1.0 + nu));
                    scalar lambda;
                    if (planeStress)
                    {
                        lambda = nu * E / ((1.0 + nu) * (1.0 - nu));
                    }
                    else
                    {
                        lambda = nu * E / ((1.0 + nu) * (1.0 - 2.0 * nu));
                    }

                    muIp += r * mu;
                    lambdaIp += r * lambda;
                }

                // ================================================================
                // Choose model based on solidMechanicsOption
                // ================================================================
                if (solidMechOption == solidMechanicsOption::linearElastic)
                {
                    for (label ic = 0; ic < nodesPerElement; ++ic)
                    {
                        const label icNdim = ic * SPATIAL_DIM;

                        for (label i = 0; i < SPATIAL_DIM; ++i)
                        {
                            const label indexL = il * SPATIAL_DIM + i;
                            const label indexR = ir * SPATIAL_DIM + i;

                            const label rowL =
                                indexL * nodesPerElement * SPATIAL_DIM;
                            const label rowR =
                                indexR * nodesPerElement * SPATIAL_DIM;

                            scalar lhs_riC_i = 0.0;
                            for (label j = 0; j < SPATIAL_DIM; ++j)
                            {
                                const scalar axj =
                                    p_scs_areav[ip * SPATIAL_DIM + j];
                                const scalar Dxj = p_D[ic * SPATIAL_DIM + j];

                                const label offSetDnDx =
                                    SPATIAL_DIM * nodesPerElement * ip +
                                    ic * SPATIAL_DIM;

                                // First mu term: -mu * dN/dx_j * A_j (diagonal
                                // contribution)
                                const scalar lhsfacDiff_i =
                                    -muIp * p_dndx[offSetDnDx + j] * axj;
                                lhs_riC_i += lhsfacDiff_i;

                                // Second mu term: -mu * dN/dx_i * A_j
                                // (off-diagonal contribution)
                                const scalar lhsfacDiff_j =
                                    -muIp * p_dndx[offSetDnDx + i] * axj;

                                p_lhs[rowL + icNdim + j] += lhsfacDiff_j;
                                p_lhs[rowR + icNdim + j] -= lhsfacDiff_j;

                                p_rhs[indexL] -= lhsfacDiff_j * Dxj;
                                p_rhs[indexR] += lhsfacDiff_j * Dxj;

                                // Lambda divergence term: -lambda * delta_ij *
                                // dN/dx_l * A_j
                                for (label l = 0; l < SPATIAL_DIM; ++l)
                                {
                                    const scalar lhsfacDiv =
                                        -lambdaIp * p_kd[i * SPATIAL_DIM + j] *
                                        p_dndx[offSetDnDx + l] * axj;

                                    p_lhs[rowL + icNdim + l] += lhsfacDiv;
                                    p_lhs[rowR + icNdim + l] -= lhsfacDiv;

                                    const scalar dxl =
                                        p_D[ic * SPATIAL_DIM + l];
                                    p_rhs[indexL] -= lhsfacDiv * dxl;
                                    p_rhs[indexR] += lhsfacDiv * dxl;
                                }
                            }

                            // Accumulated diagonal term
                            p_lhs[rowL + icNdim + i] += lhs_riC_i;
                            p_lhs[rowR + icNdim + i] -= lhs_riC_i;

                            const scalar dxi = p_D[ic * SPATIAL_DIM + i];
                            p_rhs[indexL] -= lhs_riC_i * dxi;
                            p_rhs[indexR] += lhs_riC_i * dxi;
                        }
                    }
                }
                else if (solidMechOption ==
                         solidMechanicsOption::simplifiedNeoHookean)
                {
                    // Zero out gradient tensor
                    for (label i = 0; i < SPATIAL_DIM * SPATIAL_DIM; ++i)
                    {
                        p_gradTens[i] = 0.0;
                    }

                    // Compute full lagged gradient tensor: F = I + grad(u)
                    for (label ic = 0; ic < nodesPerElement; ++ic)
                    {
                        for (label i = 0; i < SPATIAL_DIM; ++i)
                        {
                            const scalar dxi = p_D[ic * SPATIAL_DIM + i];
                            for (label j = 0; j < SPATIAL_DIM; ++j)
                            {
                                const label offSetDnDx =
                                    SPATIAL_DIM * nodesPerElement * ip +
                                    ic * SPATIAL_DIM;
                                p_gradTens[i * SPATIAL_DIM + j] +=
                                    p_dndx[offSetDnDx + j] * dxi;
                            }
                        }
                    }

                    // Add identity
                    for (label i = 0; i < SPATIAL_DIM; ++i)
                    {
                        for (label j = 0; j < SPATIAL_DIM; ++j)
                        {
                            p_gradTens[i * SPATIAL_DIM + j] +=
                                p_kd[i * SPATIAL_DIM + j];
                        }
                    }

                    // Compute I = tr(F^T * F)
                    scalar Iconst = 0.0;
                    if (SPATIAL_DIM == 3)
                    {
                        for (label i = 0; i < SPATIAL_DIM; ++i)
                        {
                            Iconst += p_gradTens[i] * p_gradTens[i] +
                                      p_gradTens[i + 3] * p_gradTens[i + 3] +
                                      p_gradTens[i + 6] * p_gradTens[i + 6];
                        }
                    }
                    else
                    {
                        for (label i = 0; i < SPATIAL_DIM; ++i)
                        {
                            Iconst += p_gradTens[i] * p_gradTens[i] +
                                      p_gradTens[i + 2] * p_gradTens[i + 2];
                        }
                    }

                    // Beta parameter
                    const scalar beta = 1.0 - 1.0 / (1.0 + Iconst);

                    // Jacobian J = det(F)
                    scalar J;
                    if (SPATIAL_DIM == 3)
                    {
                        J = p_gradTens[0] * (p_gradTens[4] * p_gradTens[8] -
                                             p_gradTens[7] * p_gradTens[5]) -
                            p_gradTens[1] * (p_gradTens[3] * p_gradTens[8] -
                                             p_gradTens[6] * p_gradTens[5]) +
                            p_gradTens[2] * (p_gradTens[3] * p_gradTens[7] -
                                             p_gradTens[6] * p_gradTens[4]);
                    }
                    else
                    {
                        J = p_gradTens[0] * p_gradTens[3] -
                            p_gradTens[2] * p_gradTens[1];
                    }

                    // Assemble stress contributions
                    for (label ic = 0; ic < nodesPerElement; ++ic)
                    {
                        const label icNdim = ic * SPATIAL_DIM;

                        for (label i = 0; i < SPATIAL_DIM; ++i)
                        {
                            const label indexL = il * SPATIAL_DIM + i;
                            const label indexR = ir * SPATIAL_DIM + i;

                            const scalar dxi = p_D[ic * SPATIAL_DIM + i];

                            scalar lhs_riC_i = 0.0;

                            for (label j = 0; j < SPATIAL_DIM; ++j)
                            {
                                const scalar axj =
                                    p_scs_areav[ip * SPATIAL_DIM + j];
                                const scalar dxj = p_D[ic * SPATIAL_DIM + j];

                                const label offSetDnDx =
                                    SPATIAL_DIM * nodesPerElement * ip +
                                    ic * SPATIAL_DIM;

                                // First set of terms: mu*beta/J * F_ik * dN/dxk
                                // * Aj
                                for (label k = 0; k < SPATIAL_DIM; ++k)
                                {
                                    const scalar factor =
                                        p_gradTens[SPATIAL_DIM * i + k] / J;
                                    const scalar lhsfacDiff_k =
                                        -(muIp * beta) * factor *
                                        p_dndx[offSetDnDx + k] * axj;

                                    p_lhs[indexL * rhsSize + icNdim + j] +=
                                        lhsfacDiff_k;
                                    p_lhs[indexR * rhsSize + icNdim + j] -=
                                        lhsfacDiff_k;

                                    p_rhs[indexL] -= lhsfacDiff_k * dxj;
                                    p_rhs[indexR] += lhsfacDiff_k * dxj;
                                }

                                // Second set of terms: mu*beta/J * dN/dxj * Aj
                                const scalar lhsfacDiff_j =
                                    -((muIp * beta) / J) *
                                    p_dndx[offSetDnDx + j] * axj;
                                lhs_riC_i += lhsfacDiff_j;

                                // Third set of terms: mu*beta/J * dN/dxi * Aj
                                const scalar lhsfacDiff_i =
                                    -((muIp * beta) / J) *
                                    p_dndx[offSetDnDx + i] * axj;
                                p_lhs[indexL * rhsSize + icNdim + j] +=
                                    lhsfacDiff_i;
                                p_lhs[indexR * rhsSize + icNdim + j] -=
                                    lhsfacDiff_i;
                                p_rhs[indexL] -= lhsfacDiff_i * dxj;
                                p_rhs[indexR] += lhsfacDiff_i * dxj;
                            }

                            // Accumulated diagonal term
                            p_lhs[indexL * rhsSize + icNdim + i] += lhs_riC_i;
                            p_lhs[indexR * rhsSize + icNdim + i] -= lhs_riC_i;
                            p_rhs[indexL] -= lhs_riC_i * dxi;
                            p_rhs[indexR] += lhs_riC_i * dxi;
                        }
                    }

                    // Additional RHS terms (volumetric contribution)
                    for (label i = 0; i < SPATIAL_DIM; ++i)
                    {
                        const label indexL = il * SPATIAL_DIM + i;
                        const label indexR = ir * SPATIAL_DIM + i;

                        for (label j = 0; j < SPATIAL_DIM; ++j)
                        {
                            const scalar axj =
                                p_scs_areav[ip * SPATIAL_DIM + j];
                            const scalar rhsTerm =
                                lambdaIp * p_kd[i * SPATIAL_DIM + j] * axj *
                                    (J - (1.0 + 0.75 * muIp / lambdaIp)) +
                                (muIp * beta / J) * p_kd[i * SPATIAL_DIM + j] *
                                    axj;

                            p_rhs[indexL] += rhsTerm;
                            p_rhs[indexR] -= rhsTerm;
                        }
                    }
                }
            }

            // Apply coefficients to system (stress terms only, mass in
            // NodeTerms)
            Base::applyCoeff_(
                A, b, connectedNodes, scratchIds, scratchVals, rhs, lhs);
        }
    }
#else
    errorMsg("FEM solid mechanics not implemented yet");
#endif
}

} // namespace accel
