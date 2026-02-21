// File : pressureCorrectionAssemblerElemTerms.cpp
// Created : Wed Jan 03 2024 13:38:51 (+0100)
// Author : Mhamad Mahdi Alloush
// Description:
// Copyright (c) 2024 CCFNUM, Lucerne University of Applied Sciences and Arts.
// SPDX-License-Identifier: BSD-3-Clause

// code
#include "flowModel.h"
#include "mesh.h"
#include "pressureCorrectionAssembler.h"
#include "zoneTransformation.h"

namespace accel
{

void pressureCorrectionAssembler::assembleElemTermsInterior_(
    const domain* domain,
    Context* ctx)
{
    auto& mesh = model_->meshRef();
    Matrix& A = ctx->getAMatrix();
    Vector& b = ctx->getBVector();

    const bool compressible = domain->isMaterialCompressible();
    const scalar comp = compressible ? 1.0 : 0.0;

    stk::mesh::BulkData& bulkData = mesh.bulkDataRef();
    stk::mesh::MetaData& metaData = mesh.metaDataRef();

    // required for frame motion (MFR)
    const auto coriolisMatrix =
        domain->zonePtr()->frameRotating()
            ? domain->zonePtr()->transformationRef().rotation().coriolisMatrix_
            : utils::matrix::Zero();
    const scalar* p_mat = coriolisMatrix.data();

    const auto origin =
        domain->zonePtr()->frameRotating()
            ? domain->zonePtr()->transformationRef().rotation().origin_
            : utils::vector::Zero();
    const scalar* p_ori = origin.data();

    // mesh motion
    const bool meshMoving = domain->zonePtr()->meshMoving();

    // space for LHS/RHS; nodesPerElem*nodesPerElem and nodesPerElem
    std::vector<scalar> lhs;
    std::vector<scalar> rhs;
    std::vector<label> scratchIds;
    std::vector<scalar> scratchVals;
    std::vector<stk::mesh::Entity> connectedNodes;

    // nodal fields to gather
    std::vector<scalar> ws_U;
    std::vector<scalar> ws_Um;
    std::vector<scalar> ws_Gpdx;
    std::vector<scalar> ws_coordinates;
    std::vector<scalar> ws_p;
    std::vector<scalar> ws_rho;
    std::vector<scalar> ws_betaRho;
    std::vector<scalar> ws_gradRho;
    std::vector<scalar> ws_psi;
    std::vector<scalar> ws_du;
    std::vector<scalar> ws_F;
    std::vector<scalar> ws_FOrig;
    std::vector<scalar> ws_scv_volume;
    std::vector<scalar> ws_scv_weight;

    // geometry related to populate
    std::vector<scalar> ws_scs_areav;
    std::vector<scalar> ws_dndx;
    std::vector<scalar> ws_deriv;
    std::vector<scalar> ws_det_j;
    std::vector<scalar> ws_velocity_shape_function;
    std::vector<scalar> ws_coordinate_shape_function;

    // integration point data that depends on size
    std::vector<scalar> coordIp(SPATIAL_DIM);
    std::vector<scalar> uIp(SPATIAL_DIM);
    std::vector<scalar> umIp(SPATIAL_DIM);
    std::vector<scalar> GpdxIp(SPATIAL_DIM);
    std::vector<scalar> dpdxIp(SPATIAL_DIM);
    std::vector<scalar> duIp(SPATIAL_DIM);
    std::vector<scalar> FIp(SPATIAL_DIM);
    std::vector<scalar> FOrigIp(SPATIAL_DIM);

    // pointers to everyone...
    scalar* p_coordIp = &coordIp[0];
    scalar* p_uIp = &uIp[0];
    scalar* p_umIp = &umIp[0];
    scalar* p_GpdxIp = &GpdxIp[0];
    scalar* p_dpdxIp = &dpdxIp[0];
    scalar* p_duIp = &duIp[0];
    scalar* p_FIp = &FIp[0];
    scalar* p_FOrigIp = &FOrigIp[0];

    // Get transport fields/side fields
    const auto& rhoSTKFieldRef = model_->rhoRef().stkFieldRef();
    const auto& betaRhoSTKFieldRef =
        model_->rhoRef().blendingFactorRef().stkFieldRef();
    const auto& gradRhoSTKFieldRef = model_->rhoRef().gradRef().stkFieldRef();
    const auto& pSTKFieldRef = model_->pRef().stkFieldRef();
    const auto& gradPSTKFieldRef = model_->pRef().gradRef().stkFieldRef();
    const auto& USTKFieldRef = model_->URef().stkFieldRef();

    const auto& mDotSTKFieldRef = model_->mDotRef().stkFieldRef();

    const auto* psiSTKFieldPtr =
        compressible ? model_->psiRef().stkFieldPtr() : nullptr;

    const auto* UmSTKFieldPtr =
        meshMoving ? model_->UmRef().stkFieldPtr() : nullptr;

#ifndef NDEBUG
    // SCL check field for debug: accumulates mesh flux to verify SCL
    auto* sclCheckSTKFieldPtr =
        meshMoving ? metaData.get_field<scalar>(stk::topology::NODE_RANK,
                                                mesh::scl_check_ID)
                   : nullptr;
#endif /* NDEBUG */

    // Get pressure diffusivity coefficient field
    const auto& duSTKFieldRef =
        *metaData.get_field<scalar>(stk::topology::NODE_RANK, flowModel::du_ID);

    // Get body force fields for buoyancy pressure stabilization
    const auto* FSTKFieldPtr =
        metaData.get_field<scalar>(stk::topology::NODE_RANK, flowModel::F_ID);
    const auto* FOrigSTKFieldPtr = metaData.get_field<scalar>(
        stk::topology::NODE_RANK, flowModel::FOriginal_ID);

    // Get geometric fields
    const auto& coordinatesRef = *metaData.get_field<scalar>(
        stk::topology::NODE_RANK, this->getCoordinatesID_(domain));

    // get interior parts the domain is defined on
    const stk::mesh::PartVector& partVec = domain->zonePtr()->interiorParts();

    // define some common selectors
    stk::mesh::Selector selAllElements =
        metaData.universal_part() & stk::mesh::selectUnion(partVec);

    // shifted ip's for fields?
    const bool isUShifted = model_->URef().isShifted();

    // shifted ip's for gradients?
    const bool isPGradientShifted = model_->pRef().isGradientShifted();

    stk::mesh::BucketVector const& elementBuckets =
        bulkData.get_buckets(stk::topology::ELEMENT_RANK, selAllElements);
    for (stk::mesh::BucketVector::const_iterator ib = elementBuckets.begin();
         ib != elementBuckets.end();
         ++ib)
    {
        stk::mesh::Bucket& elementBucket = **ib;
        const stk::mesh::Bucket::size_type nElementsPerBucket =
            elementBucket.size();

        // extract master elements
        MasterElement* meSCS = MasterElementRepo::get_surface_master_element(
            elementBucket.topology());
        MasterElement* meSCV = MasterElementRepo::get_volume_master_element(
            elementBucket.topology());

        // extract master element specifics
        const label nodesPerElement = meSCS->nodesPerElement_;
        const label numScsIp = meSCS->numIntPoints_;
        const label* lrscv = meSCS->adjacentNodes();
        const label numScvIp = meSCV->numIntPoints_;
        const label* scvIpNodeMap = meSCV->ipNodeMap();

        // resize some things; matrix related
        const label lhsSize = nodesPerElement * nodesPerElement;
        const label rhsSize = nodesPerElement;
        lhs.resize(lhsSize);
        rhs.resize(rhsSize);
        scratchIds.resize(rhsSize);
        scratchVals.resize(rhsSize);
        connectedNodes.resize(nodesPerElement);

        // algorithm related
        ws_U.resize(nodesPerElement * SPATIAL_DIM);
        ws_Um.resize(nodesPerElement * SPATIAL_DIM);
        ws_Gpdx.resize(nodesPerElement * SPATIAL_DIM);
        ws_coordinates.resize(nodesPerElement * SPATIAL_DIM);
        ws_p.resize(nodesPerElement);
        ws_rho.resize(nodesPerElement);
        ws_betaRho.resize(nodesPerElement);
        ws_gradRho.resize(nodesPerElement * SPATIAL_DIM);
        ws_psi.resize(nodesPerElement);
        ws_du.resize(nodesPerElement * SPATIAL_DIM);
        ws_F.resize(nodesPerElement * SPATIAL_DIM);
        ws_FOrig.resize(nodesPerElement * SPATIAL_DIM);
        ws_scv_volume.resize(numScvIp);
        ws_scv_weight.resize(nodesPerElement);
        ws_scs_areav.resize(numScsIp * SPATIAL_DIM);
        ws_dndx.resize(SPATIAL_DIM * numScsIp * nodesPerElement);
        ws_deriv.resize(SPATIAL_DIM * numScsIp * nodesPerElement);
        ws_det_j.resize(numScsIp);
        ws_velocity_shape_function.resize(numScsIp * nodesPerElement);
        ws_coordinate_shape_function.resize(numScsIp * nodesPerElement);

        // pointers
        scalar* p_lhs = &lhs[0];
        scalar* p_rhs = &rhs[0];
        scalar* p_U = &ws_U[0];
        scalar* p_Um = &ws_Um[0];
        scalar* p_Gpdx = &ws_Gpdx[0];
        scalar* p_coordinates = &ws_coordinates[0];
        scalar* p_p = &ws_p[0];
        scalar* p_rho = &ws_rho[0];
        scalar* p_betaRho = &ws_betaRho[0];
        scalar* p_gradRho = &ws_gradRho[0];
        scalar* p_psi = &ws_psi[0];
        scalar* p_du = &ws_du[0];
        scalar* p_F = &ws_F[0];
        scalar* p_FOrig = &ws_FOrig[0];
        scalar* p_scs_areav = &ws_scs_areav[0];
        scalar* p_dndx = &ws_dndx[0];
        scalar* p_velocity_shape_function = &ws_velocity_shape_function[0];
        scalar* p_coordinate_shape_function = &ws_coordinate_shape_function[0];

        if (isUShifted)
        {
            meSCS->shifted_shape_fcn(&p_velocity_shape_function[0]);
        }
        else
        {
            meSCS->shape_fcn(&p_velocity_shape_function[0]);
        }

        // Always use trilinear (standard) shape functions for coordinates
        meSCS->shape_fcn(&p_coordinate_shape_function[0]);

        for (stk::mesh::Bucket::size_type iElement = 0;
             iElement < nElementsPerBucket;
             ++iElement)
        {
            // get elem
            stk::mesh::Entity elem = elementBucket[iElement];

            // zero lhs/rhs
            for (label p = 0; p < lhsSize; ++p)
            {
                p_lhs[p] = 0.0;
            }
            for (label p = 0; p < rhsSize; ++p)
            {
                p_rhs[p] = 0.0;
            }

            //===============================================
            // gather nodal data; this is how we do it now..
            //===============================================
            stk::mesh::Entity const* nodeRels =
                elementBucket.begin_nodes(iElement);
            label numNodes = elementBucket.num_nodes(iElement);

            // sanity check on num nodes
            STK_ThrowAssert(numNodes == nodesPerElement);

            for (label ni = 0; ni < numNodes; ++ni)
            {
                stk::mesh::Entity node = nodeRels[ni];

                // set connected nodes
                connectedNodes[ni] = node;

                // pointers to real data
                const scalar* Gjp =
                    stk::mesh::field_data(gradPSTKFieldRef, node);
                const scalar* coords =
                    stk::mesh::field_data(coordinatesRef, node);
                const scalar* U = stk::mesh::field_data(USTKFieldRef, node);
                const scalar* du = stk::mesh::field_data(duSTKFieldRef, node);
                const scalar* gradRho =
                    stk::mesh::field_data(gradRhoSTKFieldRef, node);

                // gather scalars
                p_p[ni] = *stk::mesh::field_data(pSTKFieldRef, node);
                p_rho[ni] = *stk::mesh::field_data(rhoSTKFieldRef, node);
                p_betaRho[ni] =
                    *stk::mesh::field_data(betaRhoSTKFieldRef, node);
                p_psi[ni] = compressible
                                ? *stk::mesh::field_data(*psiSTKFieldPtr, node)
                                : 0;

                // gather vectors
                const label niNdim = ni * SPATIAL_DIM;
                for (label j = 0; j < SPATIAL_DIM; ++j)
                {
                    p_U[niNdim + j] = U[j];
                    p_Gpdx[niNdim + j] = Gjp[j];
                    p_coordinates[niNdim + j] = coords[j];
                    p_du[niNdim + j] = du[j];
                    p_gradRho[niNdim + j] = gradRho[j];
                }

                if (meshMoving)
                {
                    const scalar* Um =
                        stk::mesh::field_data(*UmSTKFieldPtr, node);
                    for (label j = 0; j < SPATIAL_DIM; ++j)
                    {
                        p_Um[niNdim + j] = Um[j];
                    }
                }
                else
                {
                    for (label j = 0; j < SPATIAL_DIM; ++j)
                    {
                        p_Um[niNdim + j] = 0.0;
                    }
                }

                // gather body force fields for buoyancy stabilization
                const scalar* F = stk::mesh::field_data(*FSTKFieldPtr, node);
                const scalar* FOrig =
                    stk::mesh::field_data(*FOrigSTKFieldPtr, node);
                for (label j = 0; j < SPATIAL_DIM; ++j)
                {
                    p_F[niNdim + j] = F[j];
                    p_FOrig[niNdim + j] = FOrig[j];
                }
            }

            // compute geometry
            scalar scs_error = 0.0;
            meSCS->determinant(
                1, &p_coordinates[0], &p_scs_areav[0], &scs_error);

            // compute SCV volumes for volume-weighted element-centre averaging
            {
                scalar scv_error = 0.0;
                meSCV->determinant(
                    1, &p_coordinates[0], &ws_scv_volume[0], &scv_error);
                scalar V_el = 0.0;
                for (label ip = 0; ip < numScvIp; ++ip)
                {
                    V_el += ws_scv_volume[ip];
                }
                const scalar invV_el = 1.0 / V_el;
                for (label i = 0; i < nodesPerElement; ++i)
                {
                    ws_scv_weight[i] = 0.0;
                }
                for (label ip = 0; ip < numScvIp; ++ip)
                {
                    ws_scv_weight[scvIpNodeMap[ip]] +=
                        ws_scv_volume[ip] * invV_el;
                }
            }

            // compute dndx for residual
            if (isPGradientShifted)
            {
                meSCS->shifted_grad_op(1,
                                       &p_coordinates[0],
                                       &p_dndx[0],
                                       &ws_deriv[0],
                                       &ws_det_j[0],
                                       &scs_error);
            }
            else
            {
                meSCS->grad_op(1,
                               &p_coordinates[0],
                               &p_dndx[0],
                               &ws_deriv[0],
                               &ws_det_j[0],
                               &scs_error);
            }

            for (label ip = 0; ip < numScsIp; ++ip)
            {
                // left and right nodes for this ip
                const label il = lrscv[2 * ip];
                const label ir = lrscv[2 * ip + 1];

                // save off mDot
                const scalar tmDot =
                    (stk::mesh::field_data(mDotSTKFieldRef, elem))[ip];

                // corresponding matrix rows
                label rowL = il * nodesPerElement;
                label rowR = ir * nodesPerElement;

                // setup for ip values; sneak in geometry for possible reduced
                // sens
                for (label j = 0; j < SPATIAL_DIM; ++j)
                {
                    p_coordIp[j] = 0.0;
                    p_uIp[j] = 0.0;
                    p_umIp[j] = 0.0;
                    p_GpdxIp[j] = 0.0;
                    p_dpdxIp[j] = 0.0;
                    p_duIp[j] = 0.0;
                    p_FIp[j] = 0.0;
                    p_FOrigIp[j] = 0.0;
                }

                const label offSetSF = ip * nodesPerElement;
                for (label ic = 0; ic < nodesPerElement; ++ic)
                {
                    const scalar r_vel =
                        p_velocity_shape_function[offSetSF + ic];
                    const scalar r_coord =
                        p_coordinate_shape_function[offSetSF + ic];
                    const scalar w_scv = ws_scv_weight[ic];

                    const label offSetDnDx =
                        SPATIAL_DIM * nodesPerElement * ip + ic * SPATIAL_DIM;
                    for (label j = 0; j < SPATIAL_DIM; ++j)
                    {
                        // use velocity shape functions
                        p_uIp[j] += r_vel * p_U[SPATIAL_DIM * ic + j];
                        p_umIp[j] += r_vel * p_Um[SPATIAL_DIM * ic + j];
                        p_duIp[j] += r_vel * p_du[SPATIAL_DIM * ic + j];

                        // use coordinates shape functions
                        p_coordIp[j] +=
                            r_coord * p_coordinates[ic * SPATIAL_DIM + j];

                        // use pressure shape function derivative
                        p_dpdxIp[j] += p_dndx[offSetDnDx + j] * p_p[ic];

                        // volume-weighted average of original body force
                        // to element centre
                        p_FOrigIp[j] += r_vel * p_FOrig[SPATIAL_DIM * ic + j];

                        // volume-weighted average of pressure gradient
                        // to element centre
                        p_GpdxIp[j] += w_scv * p_Gpdx[SPATIAL_DIM * ic + j];

                        // interpolate redistributed body force to IP
                        // using velocity shape functions
                        p_FIp[j] += w_scv * p_F[SPATIAL_DIM * ic + j];
                    }
                }

                scalar dcorr = 0;
                if (tmDot > 0)
                {
                    // deferred correction
                    for (label j = 0; j < SPATIAL_DIM; ++j)
                    {
                        const scalar dxj =
                            p_coordIp[j] - p_coordinates[il * SPATIAL_DIM + j];
                        dcorr += p_betaRho[il] * dxj *
                                 p_gradRho[il * SPATIAL_DIM + j];
                    }
                }
                else
                {
                    // deferred correction
                    for (label j = 0; j < SPATIAL_DIM; ++j)
                    {
                        const scalar dxj =
                            p_coordIp[j] - p_coordinates[ir * SPATIAL_DIM + j];
                        dcorr += p_betaRho[ir] * dxj *
                                 p_gradRho[ir * SPATIAL_DIM + j];
                    }
                }

                scalar rhoUpwind, psiUpwind;
                if (tmDot > 0)
                {
                    rhoUpwind = p_rho[il];
                    psiUpwind = p_psi[il];
                }
                else
                {
                    rhoUpwind = p_rho[ir];
                    psiUpwind = p_psi[ir];
                }

                scalar rhoHR = rhoUpwind + dcorr;

                //================================
                // rhie-chow: -rhoip*Dip*(Gpn-Gp_p).Sip
                //================================
                for (label ic = 0; ic < nodesPerElement; ++ic)
                {
                    scalar lhsfac = 0.0;
                    const label offSetDnDx =
                        SPATIAL_DIM * nodesPerElement * ip + ic * SPATIAL_DIM;
                    for (label j = 0; j < SPATIAL_DIM; ++j)
                    {
                        lhsfac += -rhoHR * p_duIp[j] * p_dndx[offSetDnDx + j] *
                                  p_scs_areav[ip * SPATIAL_DIM + j];
                    }

                    // assemble to lhs; left
                    p_lhs[rowL + ic] += lhsfac;

                    // assemble to lhs; right
                    p_lhs[rowR + ic] -= lhsfac;
                }

                //================================
                // Advection: compressibility
                // Newton-Raphson linearization: ρ*U = ρ(old)*U(new) +
                // ρ(new)*U(old) - ρ(old)*U(old)
                //   where ρ(new) = ψ p(new) and ψ = ∂ρ/∂p
                // LHS coefficient = mDot / ρ_ip * ψ_upwind
                //================================

                const label rLiL_i = rowL + il;
                const label rLiR_i = rowL + ir;
                const label rRiR_i = rowR + ir;
                const label rRiL_i = rowR + il;

                // upwind advection left node: will be 0 for incompressible flow
                const scalar alhsfacL =
                    0.5 * (tmDot + std::abs(tmDot)) * psiUpwind / rhoHR * comp;
                p_lhs[rLiL_i] += alhsfacL;
                p_lhs[rRiL_i] -= alhsfacL;

                // upwind advection right node: will be 0 for incompressible
                // flow
                const scalar alhsfacR =
                    0.5 * (tmDot - std::abs(tmDot)) * psiUpwind / rhoHR * comp;
                p_lhs[rRiR_i] -= alhsfacR;
                p_lhs[rLiR_i] += alhsfacR;

                // assemble mDot
                // mDot = ρ*U·S - ρ*D*(∇p - Gp)·S + ρ*D*(F_orig - F)·S
                //        ╰───╯   ╰──────────────╯   ╰─────────────────╯
                //      divergence   pressure RC       body force stab
                //
                scalar mDot = 0.0;
                for (label j = 0; j < SPATIAL_DIM; ++j)
                {
                    // divergence: ρ*U·S
                    mDot +=
                        rhoHR * p_uIp[j] * p_scs_areav[ip * SPATIAL_DIM + j];

                    // pressure Rhie-Chow: -ρ*D*(∇p - Gp)·S
                    mDot -= rhoHR * p_duIp[j] * (p_dpdxIp[j] - p_GpdxIp[j]) *
                            p_scs_areav[ip * SPATIAL_DIM + j];

                    // body force stabilization: +ρ*D*(F_orig - F)·S
                    mDot += rhoHR * p_duIp[j] * (p_FOrigIp[j] - p_FIp[j]) *
                            p_scs_areav[ip * SPATIAL_DIM + j];
                }

                // transform mDot to relative frame

                // 1.) frame motion
                for (label i = 0; i < SPATIAL_DIM; i++)
                {
                    for (label j = 0; j < SPATIAL_DIM; j++)
                    {
                        mDot -= rhoHR * p_mat[i * SPATIAL_DIM + j] *
                                (p_coordIp[j] - p_ori[j]) *
                                p_scs_areav[ip * SPATIAL_DIM + i];
                    }
                }

                // 2.) mesh motion
                for (label i = 0; i < SPATIAL_DIM; i++)
                {
                    mDot -=
                        rhoHR * p_umIp[i] * p_scs_areav[ip * SPATIAL_DIM + i];
                }

#ifndef NDEBUG
                // SCL check: accumulate grid flux to nodes
                if (sclCheckSTKFieldPtr)
                {
                    scalar gridFlux = 0.0;
                    for (label i = 0; i < SPATIAL_DIM; i++)
                    {
                        gridFlux += rhoHR * p_umIp[i] *
                                    p_scs_areav[ip * SPATIAL_DIM + i];
                    }
                    stk::mesh::Entity nodeL = nodeRels[il];
                    stk::mesh::Entity nodeR = nodeRels[ir];
                    scalar* sclL =
                        stk::mesh::field_data(*sclCheckSTKFieldPtr, nodeL);
                    scalar* sclR =
                        stk::mesh::field_data(*sclCheckSTKFieldPtr, nodeR);
                    *sclL -= gridFlux;
                    *sclR += gridFlux;
                }
#endif /* NDEBUG */

                // residual; left and right
                p_rhs[il] -= mDot;
                p_rhs[ir] += mDot;
            }

            this->applyCoeff_(
                A, b, connectedNodes, scratchIds, scratchVals, rhs, lhs);
        }
    }
}

} // namespace accel
