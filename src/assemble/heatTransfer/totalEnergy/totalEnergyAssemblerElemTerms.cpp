// File       : totalEnergyAssemblerElemTerms.cpp
// Created    : Thu Apr 02 2025 15:40:38 (+0100)
// Author     : Mhamad Mahdi Alloush
// Description:
// Copyright 2025 CCFNUM HSLU T&A. All Rights Reserved.

#include "totalEnergyAssembler.h"

namespace accel
{

void totalEnergyAssembler::assembleElemTermsInterior_(const domain* domain,
                                                      Context* ctx)
{
    const auto& mesh = field_broker_->meshRef();
    Matrix& A = ctx->getAMatrix();
    Vector& b = ctx->getBVector();

    const stk::mesh::BulkData& bulkData = mesh.bulkDataRef();
    const stk::mesh::MetaData& metaData = mesh.metaDataRef();

    const bool includeAdv = domain->type() == domainType::fluid;
    const bool includeViscousWork =
        domain->heatTransfer_.includeViscousWork_ && includeAdv;
    // In a rotating zone, transport the algebraically equivalent rothalpy.
    // The public field is transformed back to absolute h0 after every solve.
    const bool steadyRotatingEnergyForm =
        usesSteadyRotatingEnergyForm_(domain, includeAdv);
    const bool includeRotationWork =
        includesRotatingPressureWork_(domain, includeAdv);

    // NSO active only when advection is present and enabled in expert params
    const bool NSO =
        includeAdv &&
        mesh.controlsRef().solverRef().solverControl_.expertParameters_.nso_;
    // 4th-order factor (0.0 = 2nd order, 1.0 = 4th order)
    const scalar fourthFac =
        mesh.controlsRef()
            .solverRef()
            .solverControl_.expertParameters_.nsoFourthOrderFac_;

    // TODO: Account for BLOCKSIZE space for LHS/RHS; nodesPerElem*nodesPerElem*
    // and nodesPerElem
    std::vector<scalar> lhs;
    std::vector<scalar> rhs;
    std::vector<label> scratchIds;
    std::vector<scalar> scratchVals;
    std::vector<stk::mesh::Entity> connectedNodes;

    // nodal fields to gather
    std::vector<scalar> ws_coordinates;
    std::vector<scalar> ws_h0;
    std::vector<scalar> ws_T;
    std::vector<scalar> ws_beta;
    std::vector<scalar> ws_gradH0;
    std::vector<scalar> ws_lambdaEff;
    std::vector<scalar> ws_cp;
    std::vector<scalar> ws_muEff;
    std::vector<scalar> ws_velocity;
    std::vector<scalar> ws_rho;
    std::vector<scalar> ws_p;

    // geometry related to populate
    std::vector<scalar> ws_scs_areav;
    std::vector<scalar> ws_dndx;
    std::vector<scalar> ws_deriv;
    std::vector<scalar> ws_det_j;
    std::vector<scalar> ws_shape_function;
    std::vector<scalar> ws_coordinate_shape_function;
    std::vector<scalar> ws_gijUpper;
    std::vector<scalar> ws_gijLower;
    std::vector<scalar> ws_gij_deriv;

    // ip values
    std::vector<scalar> coordIp(SPATIAL_DIM);
    std::vector<scalar> vwIp(SPATIAL_DIM);
    std::vector<scalar> velocityIp(SPATIAL_DIM);
    std::vector<scalar> gradVelocityIp(SPATIAL_DIM * SPATIAL_DIM);
    std::vector<scalar> omegaCrossR(SPATIAL_DIM);

    // const-size vectors
    std::vector<scalar> rhoUIp(SPATIAL_DIM);
    std::vector<scalar> dh0dxIp(SPATIAL_DIM);

    // pointers
    scalar* p_coordIp = &coordIp[0];
    scalar* p_vwIp = &vwIp[0];
    scalar* p_velocityIp = &velocityIp[0];
    scalar* p_gradVelocityIp = &gradVelocityIp[0];
    scalar* p_omegaCrossR = &omegaCrossR[0];
    scalar* p_rhoUIp = &rhoUIp[0];
    scalar* p_dh0dxIp = &dh0dxIp[0];

    // Get transport fields/side fields
    const auto& h0STKFieldRef = phi_->stkFieldRef();
    const auto& gradH0STKFieldRef = phi_->gradRef().stkFieldRef();
    const auto& blendingFactorSTKFieldRef =
        phi_->blendingFactorRef().stkFieldRef();
    const auto& TSTKFieldRef = model_->TRef().stkFieldRef();
    const auto& cpSTKFieldRef = model_->cpRef().stkFieldRef();
    const auto* muEffSTKFieldPtr =
        includeViscousWork ? model_->muEffRef().stkFieldPtr() : nullptr;
    const auto* USTKFieldPtr =
        includeViscousWork || NSO ? model_->URef().stkFieldPtr() : nullptr;
    const auto* rhoSTKFieldPtr = NSO ? this->rhoRef().stkFieldPtr() : nullptr;
    const auto* pSTKFieldPtr =
        includeRotationWork ? model_->pRef().stkFieldPtr() : nullptr;

    const auto rotationMatrix =
        includeRotationWork || steadyRotatingEnergyForm
            ? domain->zonePtr()->transformationRef().rotation().coriolisMatrix_
            : utils::matrix::Zero();
    const scalar* p_rotationMatrix = rotationMatrix.data();

    const auto rotationOrigin =
        includeRotationWork || steadyRotatingEnergyForm
            ? domain->zonePtr()->transformationRef().rotation().origin_
            : utils::vector::Zero();
    const scalar* p_rotationOrigin = rotationOrigin.data();
    const auto viscousWorkFrameMatrix =
        steadyRotatingEnergyForm ? rotationMatrix : utils::matrix::Zero();
    const scalar* p_viscousWorkFrameMatrix = viscousWorkFrameMatrix.data();

    // Get geometric fields
    const auto& coordinatesRef = *metaData.get_field<scalar>(
        stk::topology::NODE_RANK, this->getCoordinatesID_(domain));

    // get interior parts the domain is defined on
    const stk::mesh::PartVector& partVec = domain->zonePtr()->interiorParts();

    // define some common selectors
    const stk::mesh::Selector selAllElements =
        metaData.universal_part() & stk::mesh::selectUnion(partVec);

    // shifted ip's for field?
    const bool isShifted = phi_->isShifted();

    // shifted ip's for gradients?
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

        // extract master element
        MasterElement* meSCS = MasterElementRepo::get_surface_master_element(
            elementBucket.topology());

        // extract master element specifics
        const label nodesPerElement = meSCS->nodesPerElement_;
        const label numScsIp = meSCS->numIntPoints_;
        const label* lrscv = meSCS->adjacentNodes();

        // resize some things; matrix related
        const label lhsSize = nodesPerElement * nodesPerElement;
        const label rhsSize = nodesPerElement;
        lhs.resize(lhsSize);
        rhs.resize(rhsSize);
        scratchIds.resize(rhsSize);
        scratchVals.resize(rhsSize);
        connectedNodes.resize(nodesPerElement);

        // algorithm related
        ws_coordinates.resize(nodesPerElement * SPATIAL_DIM);
        ws_h0.resize(nodesPerElement);
        ws_T.resize(nodesPerElement);
        ws_beta.resize(nodesPerElement);
        ws_gradH0.resize(nodesPerElement * SPATIAL_DIM);
        ws_lambdaEff.resize(nodesPerElement);
        ws_cp.resize(nodesPerElement);
        ws_muEff.resize(nodesPerElement);
        ws_velocity.resize(nodesPerElement * SPATIAL_DIM);
        ws_rho.resize(nodesPerElement);
        ws_p.resize(includeRotationWork ? nodesPerElement : 0);
        ws_scs_areav.resize(numScsIp * SPATIAL_DIM);
        ws_dndx.resize(SPATIAL_DIM * numScsIp * nodesPerElement);
        ws_deriv.resize(SPATIAL_DIM * numScsIp * nodesPerElement);
        ws_det_j.resize(numScsIp);
        ws_shape_function.resize(numScsIp * nodesPerElement);
        ws_coordinate_shape_function.resize(numScsIp * nodesPerElement);
        ws_gijUpper.resize(numScsIp * SPATIAL_DIM * SPATIAL_DIM);
        ws_gijLower.resize(numScsIp * SPATIAL_DIM * SPATIAL_DIM);
        ws_gij_deriv.resize(numScsIp * nodesPerElement * SPATIAL_DIM);

        // pointer to lhs/rhs
        scalar* p_lhs = &lhs[0];
        scalar* p_rhs = &rhs[0];
        scalar* p_coordinates = &ws_coordinates[0];
        scalar* p_h0 = &ws_h0[0];
        scalar* p_T = &ws_T[0];
        scalar* p_beta = &ws_beta[0];
        scalar* p_gradH0 = &ws_gradH0[0];
        scalar* p_lambdaEff = &ws_lambdaEff[0];
        scalar* p_cp = &ws_cp[0];
        scalar* p_muEff = &ws_muEff[0];
        scalar* p_U = &ws_velocity[0];
        scalar* p_rho = &ws_rho[0];
        scalar* p_p = includeRotationWork ? &ws_p[0] : nullptr;
        scalar* p_scs_areav = &ws_scs_areav[0];
        scalar* p_dndx = &ws_dndx[0];
        scalar* p_shape_function = &ws_shape_function[0];
        scalar* p_coordinate_shape_function = &ws_coordinate_shape_function[0];
        scalar* p_gijUp = &ws_gijUpper[0];
        scalar* p_gijLow = &ws_gijLower[0];

        // extract shape function
        if (isShifted)
        {
            meSCS->shifted_shape_fcn(&p_shape_function[0]);
        }
        else
        {
            meSCS->shape_fcn(&p_shape_function[0]);
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
            stk::mesh::Entity const* nodeRels = bulkData.begin_nodes(elem);
            label numNodes = bulkData.num_nodes(elem);

            // sanity check on num nodes
            STK_ThrowAssert(numNodes == nodesPerElement);

            for (label ni = 0; ni < numNodes; ++ni)
            {
                stk::mesh::Entity node = nodeRels[ni];

                // set connected nodes
                connectedNodes[ni] = node;

                // pointers to real data
                const scalar* coords =
                    stk::mesh::field_data(coordinatesRef, node);
                const scalar* gradH0 =
                    stk::mesh::field_data(gradH0STKFieldRef, node);

                // gather scalars
                p_lambdaEff[ni] =
                    *stk::mesh::field_data(*GammaSTKFieldPtr_, node);
                p_cp[ni] = *stk::mesh::field_data(cpSTKFieldRef, node);
                p_h0[ni] = *stk::mesh::field_data(h0STKFieldRef, node);
                p_T[ni] = *stk::mesh::field_data(TSTKFieldRef, node);
                p_beta[ni] =
                    *stk::mesh::field_data(blendingFactorSTKFieldRef, node);
                if (includeRotationWork)
                {
                    p_p[ni] = *stk::mesh::field_data(*pSTKFieldPtr, node);
                }

                // gather vectors
                for (label i = 0; i < SPATIAL_DIM; ++i)
                {
                    p_gradH0[ni * SPATIAL_DIM + i] = gradH0[i];
                    p_coordinates[ni * SPATIAL_DIM + i] = coords[i];
                }

                if (includeViscousWork || NSO)
                {
                    const scalar* U =
                        stk::mesh::field_data(*USTKFieldPtr, node);
                    for (label i = 0; i < SPATIAL_DIM; ++i)
                    {
                        p_U[ni * SPATIAL_DIM + i] = U[i];
                    }
                }

                if (includeViscousWork)
                {
                    p_muEff[ni] =
                        *stk::mesh::field_data(*muEffSTKFieldPtr, node);
                }

                // NSO density field
                if (NSO)
                {
                    p_rho[ni] = *stk::mesh::field_data(*rhoSTKFieldPtr, node);
                }
            }

            // compute geometry
            scalar scs_error = 0.0;
            meSCS->determinant(
                1, &p_coordinates[0], &p_scs_areav[0], &scs_error);

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

            // compute metric tensors (only when NSO active)
            if (NSO)
            {
                meSCS->gij(&ws_coordinates[0],
                           &ws_gijUpper[0],
                           &ws_gijLower[0],
                           &ws_deriv[0]);
            }

            for (label ip = 0; ip < numScsIp; ++ip)
            {
                // left and right nodes for this ip
                const label il = lrscv[2 * ip];
                const label ir = lrscv[2 * ip + 1];

                // save off mDot
                const scalar tmDot =
                    (stk::mesh::field_data(*mDotSTKFieldPtr_, elem))[ip];

                for (label j = 0; j < SPATIAL_DIM; ++j)
                {
                    p_coordIp[j] = 0.0;
                    p_vwIp[j] = 0.0;
                    p_velocityIp[j] = 0.0;
                }
                for (label j = 0; j < SPATIAL_DIM * SPATIAL_DIM; ++j)
                {
                    p_gradVelocityIp[j] = 0.0;
                }

                // save off ip values; offset to Shape Function
                scalar lambdaEffIp = 0.0;
                scalar h0Ip = 0.0;
                scalar pIp = 0.0;
                scalar muEffIp = 0.0;
                const label offSetSF = ip * nodesPerElement;
                for (label ic = 0; ic < nodesPerElement; ++ic)
                {
                    const scalar r = p_shape_function[offSetSF + ic];
                    const scalar r_coord =
                        p_coordinate_shape_function[offSetSF + ic];

                    lambdaEffIp += r * p_lambdaEff[ic];
                    h0Ip += r * p_h0[ic];
                    if (includeRotationWork)
                    {
                        pIp += r * p_p[ic];
                    }
                    if (includeViscousWork)
                    {
                        muEffIp += r * p_muEff[ic];
                    }

                    for (label i = 0; i < SPATIAL_DIM; ++i)
                    {
                        p_coordIp[i] +=
                            r_coord * p_coordinates[ic * SPATIAL_DIM + i];
                        if (includeViscousWork)
                        {
                            p_velocityIp[i] += r * p_U[ic * SPATIAL_DIM + i];
                            const label offSetDnDx =
                                SPATIAL_DIM * nodesPerElement * ip +
                                ic * SPATIAL_DIM;
                            for (label j = 0; j < SPATIAL_DIM; ++j)
                            {
                                p_gradVelocityIp[i * SPATIAL_DIM + j] +=
                                    p_U[ic * SPATIAL_DIM + i] *
                                    p_dndx[offSetDnDx + j];
                            }
                        }
                    }
                }

                // Form tau.U from the element-consistent velocity gradient at
                // the SCS integration point.  The former nodal-gradient
                // interpolation produces false strain for exact solid-body
                // rotation, especially next to walls.
                if (includeViscousWork)
                {
                    for (label i = 0; i < SPATIAL_DIM; ++i)
                    {
                        p_omegaCrossR[i] = 0.0;
                        for (label j = 0; j < SPATIAL_DIM; ++j)
                        {
                            p_omegaCrossR[i] +=
                                p_viscousWorkFrameMatrix[i * SPATIAL_DIM + j] *
                                (p_coordIp[j] - p_rotationOrigin[j]);
                        }
                    }
                    for (label i = 0; i < SPATIAL_DIM; ++i)
                        p_velocityIp[i] -= p_omegaCrossR[i];

                    scalar divU = 0.0;
                    for (label i = 0; i < SPATIAL_DIM; ++i)
                    {
                        divU += p_gradVelocityIp[i * SPATIAL_DIM + i];
                    }
                    for (label i = 0; i < SPATIAL_DIM; ++i)
                    {
                        p_vwIp[i] =
                            -2.0 / 3.0 * muEffIp * divU * p_velocityIp[i];
                        for (label j = 0; j < SPATIAL_DIM; ++j)
                        {
                            p_vwIp[i] +=
                                muEffIp *
                                (p_gradVelocityIp[i * SPATIAL_DIM + j] +
                                 p_gradVelocityIp[j * SPATIAL_DIM + i]) *
                                p_velocityIp[j];
                        }
                    }
                }

                //================================
                // Advection
                //================================

                // assemble advection; rhs and upwind contributions
                scalar h0Upwind;
                scalar dcorr = 0;
                if (tmDot > 0)
                {
                    h0Upwind = p_h0[il];

                    // deferred correction
                    for (label j = 0; j < SPATIAL_DIM; ++j)
                    {
                        const scalar dxj =
                            p_coordIp[j] - p_coordinates[il * SPATIAL_DIM + j];
                        dcorr +=
                            p_beta[il] * dxj * p_gradH0[il * SPATIAL_DIM + j];
                    }
                }
                else
                {
                    h0Upwind = p_h0[ir];

                    // deferred correction
                    for (label j = 0; j < SPATIAL_DIM; ++j)
                    {
                        const scalar dxj =
                            p_coordIp[j] - p_coordinates[ir * SPATIAL_DIM + j];
                        dcorr +=
                            p_beta[ir] * dxj * p_gradH0[ir * SPATIAL_DIM + j];
                    }
                }

                // total upwind advection
                scalar aflux = tmDot * (h0Upwind + dcorr);

                // Conservative rotating-frame pressure work:
                // div[p (Omega x r)].
                if (includeRotationWork)
                {
                    for (label i = 0; i < SPATIAL_DIM; ++i)
                    {
                        p_omegaCrossR[i] = 0.0;
                        for (label j = 0; j < SPATIAL_DIM; ++j)
                        {
                            p_omegaCrossR[i] +=
                                p_rotationMatrix[i * SPATIAL_DIM + j] *
                                (p_coordIp[j] - p_rotationOrigin[j]);
                        }
                    }

                    for (label i = 0; i < SPATIAL_DIM; ++i)
                    {
                        aflux += pIp * p_omegaCrossR[i] *
                                 p_scs_areav[ip * SPATIAL_DIM + i];
                    }
                }

                const label rowL = il * nodesPerElement;
                const label rowR = ir * nodesPerElement;

                const label rLiL_i = rowL + il;
                const label rLiR_i = rowL + ir;
                const label rRiR_i = rowR + ir;
                const label rRiL_i = rowR + il;

                // right hand side; L and R
                p_rhs[il] -= aflux;
                p_rhs[ir] += aflux;

                // upwind advection left node
                const scalar alhsfacL = 0.5 * (tmDot + std::abs(tmDot));
                p_lhs[rLiL_i] += alhsfacL;
                p_lhs[rRiL_i] -= alhsfacL;

                // upwind advection right node
                const scalar alhsfacR = 0.5 * (tmDot - std::abs(tmDot));
                p_lhs[rRiR_i] -= alhsfacR;
                p_lhs[rLiR_i] += alhsfacR;

                //================================
                // NSO stabilization
                //================================
                if (NSO)
                {
                    constexpr scalar Cupw = 0.1;
                    constexpr scalar small = 1.0e-16;

                    const label gijOffset = ip * SPATIAL_DIM * SPATIAL_DIM;

                    for (label j = 0; j < SPATIAL_DIM; ++j)
                    {
                        p_rhoUIp[j] = 0.0;
                        p_dh0dxIp[j] = 0.0;
                    }

                    scalar dFdxAdv = 0.0;
                    scalar dFdxCont = 0.0;

                    for (label ic = 0; ic < nodesPerElement; ++ic)
                    {
                        const scalar r = p_shape_function[offSetSF + ic];
                        const scalar h0 = p_h0[ic];
                        const scalar rho = p_rho[ic];

                        const label offSetDnDx =
                            SPATIAL_DIM * nodesPerElement * ip +
                            ic * SPATIAL_DIM;
                        for (label j = 0; j < SPATIAL_DIM; ++j)
                        {
                            const scalar dnj = p_dndx[offSetDnDx + j];
                            const scalar U = p_U[ic * SPATIAL_DIM + j];
                            p_dh0dxIp[j] += h0 * dnj;
                            p_rhoUIp[j] += r * rho * U;
                            dFdxAdv += rho * U * h0 * dnj;
                            dFdxCont += rho * U * dnj;
                        }
                    }

                    // alternative residual (commutation error)
                    scalar residualIp = dFdxAdv - h0Ip * dFdxCont;
                    for (label j = 0; j < SPATIAL_DIM; ++j)
                    {
                        residualIp -= p_rhoUIp[j] * p_dh0dxIp[j];
                    }

                    // denominator for nu and upwind nu
                    scalar gUpperMagGradH0 = 0.0;
                    scalar rhoUiGLowerRhoUj = 0.0;
                    for (label j = 0; j < SPATIAL_DIM; ++j)
                    {
                        for (label k = 0; k < SPATIAL_DIM; ++k)
                        {
                            gUpperMagGradH0 +=
                                p_dh0dxIp[j] *
                                p_gijUp[gijOffset + j * SPATIAL_DIM + k] *
                                p_dh0dxIp[k];
                            rhoUiGLowerRhoUj +=
                                p_rhoUIp[j] *
                                p_gijLow[gijOffset + j * SPATIAL_DIM + k] *
                                p_rhoUIp[k];
                        }
                    }

                    // artificial viscosity
                    const scalar nuResidualIp = std::sqrt(
                        (residualIp * residualIp) / (gUpperMagGradH0 + small));
                    const scalar nuFirstOrder = std::sqrt(rhoUiGLowerRhoUj);
                    const scalar nuIp =
                        std::min(Cupw * nuFirstOrder, nuResidualIp);

                    // NSO diffusion-like term:
                    // -nu * gUpper_ij * (dh0/dxj - Gjq_j) * areav_i
                    scalar gijFac = 0.0;
                    for (label ic = 0; ic < nodesPerElement; ++ic)
                    {
                        const scalar r = p_shape_function[offSetSF + ic];
                        const scalar h0 = p_h0[ic];

                        const label offSetDnDx =
                            SPATIAL_DIM * nodesPerElement * ip +
                            ic * SPATIAL_DIM;
                        scalar lhsfac = 0.0;
                        for (label j = 0; j < SPATIAL_DIM; ++j)
                        {
                            const scalar axi =
                                p_scs_areav[ip * SPATIAL_DIM + j];
                            for (label k = 0; k < SPATIAL_DIM; ++k)
                            {
                                const scalar gUp =
                                    p_gijUp[gijOffset + j * SPATIAL_DIM + k];
                                const scalar fac =
                                    gUp * p_dndx[offSetDnDx + k] * axi;
                                const scalar facGj =
                                    r * gUp * p_gradH0[ic * SPATIAL_DIM + k] *
                                    axi;
                                // fourthFac: 0 = 2nd order, 1 = 4th order
                                gijFac += fac * h0 - facGj * fourthFac;
                                lhsfac += -fac;
                            }
                        }

                        p_lhs[rowL + ic] += nuIp * lhsfac;
                        p_lhs[rowR + ic] -= nuIp * lhsfac;
                    }

                    // NSO residual contribution
                    const scalar residualNSO = -nuIp * gijFac;
                    p_rhs[il] -= residualNSO;
                    p_rhs[ir] += residualNSO;
                }

                //================================
                // Diffusion
                //================================

                for (label ic = 0; ic < nodesPerElement; ++ic)
                {
                    const label rowL = il * nodesPerElement;
                    const label rowR = ir * nodesPerElement;

                    const label offSetDnDx =
                        SPATIAL_DIM * nodesPerElement * ip + ic * SPATIAL_DIM;
                    scalar lhs_riC_i = 0.0;
                    for (label j = 0; j < SPATIAL_DIM; ++j)
                    {
                        const scalar axj = p_scs_areav[ip * SPATIAL_DIM + j];

                        // -Gamma*dphi/dxj*A_j; fixed i over j loop; see
                        // below..
                        const scalar lhsfacDiff_i =
                            -lambdaEffIp * p_dndx[offSetDnDx + j] * axj;

                        // lhs; il then ir
                        lhs_riC_i += lhsfacDiff_i;
                    }

                    // deal with accumulated lhs and flux for
                    // -Gamma*dphi/dxj*Aj
                    p_lhs[rowL + ic] += lhs_riC_i / p_cp[il];
                    p_lhs[rowR + ic] -= lhs_riC_i / p_cp[ir];

                    const scalar T = p_T[ic];
                    p_rhs[il] -= lhs_riC_i * T;
                    p_rhs[ir] += lhs_riC_i * T;
                }

                //================================
                // Viscous Work
                //================================
                for (label j = 0; j < SPATIAL_DIM; ++j)
                {
                    const scalar axj = p_scs_areav[ip * SPATIAL_DIM + j];
                    p_rhs[il] += vwIp[j] * axj;
                    p_rhs[ir] -= vwIp[j] * axj;
                }
            }

            Base::applyCoeff_(
                A, b, connectedNodes, scratchIds, scratchVals, rhs, lhs);
        }
    }
}

} // namespace accel
