// File       : robinNeumannWall.h
// Created    : Mon Oct 26 2026 03:10:00 (+0300)
// Author     : Mohammad Nasreddine
// Description: Fluid-solid wall kernels of the single-phase and VOF
//              pressure-correction assemblers: no-slip mass flux, plus the
//              Robin-Neumann row (Tukovic et al. 2014) when requested
// Copyright 2026 CCFNUM HSLU T&A. All Rights Reserved.

#ifndef ROBINNEUMANNWALL_H
#define ROBINNEUMANNWALL_H

// code
#include "domain.h"
#include "flowModel.h"
#include "interface.h"
#include "mesh.h"
#include "robinNeumannState.h"

namespace accel
{

template <class Assembler>
void assembleWallInterfaceSide(Assembler& self,
                               const domain* domain,
                               const interfaceSideInfo* interfaceSideInfoPtr,
                               typename Assembler::Context* ctx,
                               const wallPhase& phase,
                               const bool robin)
{
    auto* model = self.model_;
    auto& mesh = model->meshRef();

    auto& A = ctx->getAMatrix();
    auto& b = ctx->getBVector();

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

    // ip values; both boundary and opposing surface
    std::vector<scalar> umBip(SPATIAL_DIM);
    std::vector<scalar> coordBip(SPATIAL_DIM);

    // pointers to fixed values
    scalar* p_umBip = &umBip[0];
    scalar* p_coordBip = &coordBip[0];

    // nodal fields to gather
    std::vector<scalar> ws_Um;
    std::vector<scalar> ws_coordinates;
    std::vector<scalar> ws_rho;
    std::vector<scalar> ws_alpha;
    std::vector<scalar> ws_rho_mix;
    std::vector<scalar> ws_du;     // pressure diffusivity (per component)
    std::vector<scalar> ws_FOrig;  // original body force (per component)
    std::vector<scalar> ws_elem_p; // element nodal p for the BIP gradient
    std::vector<scalar> ws_dndx;   // shape-function derivatives at the BIP
    std::vector<scalar> ws_det_j;  // Jacobian (general_face_grad_op)

    // master element
    std::vector<scalar> ws_velocity_face_shape_function;
    std::vector<scalar> ws_coordinate_face_shape_function;

    // Get fields
    const auto& rhoSTKFieldRef = *phase.rho;
    const auto& sideUSTKFieldRef = model->URef().sideFieldRef().stkFieldRef();
    const auto* UmSTKFieldPtr =
        meshMoving ? model->UmRef().stkFieldPtr() : nullptr;

    // Robin fields: wall acceleration a_n, flux memory Q_vol, frozen term E
    const STKScalarField* aWallBcSTKFieldPtr = nullptr;
    const STKScalarField* rnWallFluxSTKFieldPtr = nullptr;
    STKScalarField* rnWallFluxESTKFieldPtr = nullptr;
    const STKScalarField* duLhsSTKFieldPtr = nullptr;
    const STKScalarField* pSTKFieldPtr = nullptr;
    const STKScalarField* FOrigSTKFieldPtr = nullptr;
    if (robin)
    {
        aWallBcSTKFieldPtr = &model->aWallBcRef().stkFieldRef();
        rnWallFluxSTKFieldPtr = &model->rnWallFluxRef().stkFieldRef();
        rnWallFluxESTKFieldPtr = &model->rnWallFluxERef().stkFieldRef();

        // pressure diffusivity: same SIMPLEC choice (duLhs) as the interior
        const bool consistent =
            model->controlsRef()
                .solverRef()
                .solverControl_.expertParameters_.consistent_;
        duLhsSTKFieldPtr = metaData.get_field<scalar>(
            stk::topology::NODE_RANK,
            consistent ? flowModel::duTilde_ID : flowModel::du_ID);
        pSTKFieldPtr = &model->pRef().stkFieldRef();

        // interior Rhie-Chow balances grad(p) against the original body force
        FOrigSTKFieldPtr = metaData.get_field<scalar>(stk::topology::NODE_RANK,
                                                      flowModel::FOriginal_ID);
    }

#ifndef NDEBUG
    // SCL check field for debug: accumulates mesh flux to verify SCL
    auto* sclCheckSTKFieldPtr =
        meshMoving ? metaData.get_field<scalar>(stk::topology::NODE_RANK,
                                                mesh::scl_check_ID)
                   : nullptr;
#endif /* NDEBUG */

    // Get geometric fields
    const auto& coordsSTKFieldRef = *metaData.get_field<scalar>(
        stk::topology::NODE_RANK, self.getCoordinatesID_(domain));
    const auto& exposedAreaVecSTKFieldRef = *metaData.get_field<scalar>(
        metaData.side_rank(), self.getExposedAreaVectorID_(domain));

    // define vector of parent topos; should always be UNITY in size
    std::vector<stk::topology> parentTopo;

    // define some common selectors
    stk::mesh::Selector selAllSides =
        metaData.universal_part() &
        stk::mesh::selectUnion(interfaceSideInfoPtr->currentPartVec_);

    // shifted ip's for fields?
    const bool isUShifted = model->URef().isShifted();

    // h_s and rho_s: set once per timestep in preTimeStep()
    const scalar robinHs = self.robinNeumann_.robinHs();
    const scalar robinRhoSolid = self.robinNeumann_.robinRhoSolid();
    // compute and store E in this pass (first corrector) or reuse it
    const bool computeE = self.robinNeumann_.freezeE();
    const auto* rnParams =
        robin
            ? &domain->fluidSolidCouplingRef(interfaceSideInfoPtr->interfPtr())
                   .robin_
            : nullptr;

    stk::mesh::BucketVector const& sideBuckets =
        bulkData.get_buckets(metaData.side_rank(), selAllSides);
    for (stk::mesh::BucketVector::const_iterator ib = sideBuckets.begin();
         ib != sideBuckets.end();
         ++ib)
    {
        stk::mesh::Bucket& sideBucket = **ib;

        // extract connected element topology
        sideBucket.parent_topology(stk::topology::ELEMENT_RANK, parentTopo);
        STK_ThrowAssert(parentTopo.size() == 1);
        stk::topology theElemTopo = parentTopo[0];

        // volume master element
        MasterElement* meSCS =
            MasterElementRepo::get_surface_master_element(theElemTopo);
        const label nodesPerElement = meSCS->nodesPerElement_;

        // face master element
        MasterElement* meFC = MasterElementRepo::get_surface_master_element(
            sideBucket.topology());
        const label nodesPerSide = sideBucket.topology().num_nodes();
        const label numScsBip = meFC->numIntPoints_;

        // resize some things; matrix related
        const label lhsSize = nodesPerElement * nodesPerElement;
        const label rhsSize = nodesPerElement;
        lhs.resize(lhsSize);
        rhs.resize(rhsSize);
        scratchIds.resize(rhsSize);
        scratchVals.resize(rhsSize);
        connectedNodes.resize(nodesPerElement);

        // algorithm related; element
        ws_Um.resize(nodesPerSide * SPATIAL_DIM);
        ws_coordinates.resize(nodesPerElement * SPATIAL_DIM);
        ws_rho.resize(nodesPerSide);
        ws_alpha.resize(nodesPerSide);
        ws_rho_mix.resize(nodesPerSide);
        ws_du.resize(nodesPerSide * SPATIAL_DIM);
        ws_FOrig.assign(nodesPerSide * SPATIAL_DIM, 0.0);
        ws_elem_p.resize(nodesPerElement);
        ws_dndx.resize(SPATIAL_DIM * nodesPerElement);
        ws_det_j.resize(1);
        ws_velocity_face_shape_function.resize(numScsBip * nodesPerSide);
        ws_coordinate_face_shape_function.resize(numScsBip * nodesPerSide);

        // pointers
        scalar* p_lhs = &lhs[0];
        scalar* p_rhs = &rhs[0];
        scalar* p_Um = &ws_Um[0];
        scalar* p_coordinates = &ws_coordinates[0];
        scalar* p_rho = &ws_rho[0];
        scalar* p_alpha = &ws_alpha[0];
        scalar* p_du = &ws_du[0];
        scalar* p_velocity_face_shape_function =
            &ws_velocity_face_shape_function[0];
        scalar* p_coordinate_face_shape_function =
            &ws_coordinate_face_shape_function[0];

        // shape functions; boundary
        if (isUShifted)
        {
            meFC->shifted_shape_fcn(&p_velocity_face_shape_function[0]);
        }
        else
        {
            meFC->shape_fcn(&p_velocity_face_shape_function[0]);
        }

        // Always use trilinear (standard) shape functions for coordinates
        meFC->shape_fcn(&p_coordinate_face_shape_function[0]);

        const stk::mesh::Bucket::size_type nSidesPerBucket = sideBucket.size();
        for (stk::mesh::Bucket::size_type iSide = 0; iSide < nSidesPerBucket;
             ++iSide)
        {
            // zero lhs/rhs
            for (label p = 0; p < lhsSize; ++p)
            {
                p_lhs[p] = 0.0;
            }
            for (label p = 0; p < rhsSize; ++p)
            {
                p_rhs[p] = 0.0;
            }

            // get face
            stk::mesh::Entity side = sideBucket[iSide];

            //======================================
            // gather nodal data off of face
            //======================================
            stk::mesh::Entity const* sideNodeRels = bulkData.begin_nodes(side);
            label numSideNodes = bulkData.num_nodes(side);

            // sanity check on num nodes
            STK_ThrowAssert(numSideNodes == nodesPerSide);
            for (label ni = 0; ni < numSideNodes; ++ni)
            {
                stk::mesh::Entity node = sideNodeRels[ni];

                // gather scalars
                p_rho[ni] = *stk::mesh::field_data(rhoSTKFieldRef, node);
                if (phase.alpha)
                {
                    p_alpha[ni] = *stk::mesh::field_data(*phase.alpha, node);
                }

                if (robin)
                {
                    if (phase.rhoMix)
                    {
                        ws_rho_mix[ni] =
                            *stk::mesh::field_data(*phase.rhoMix, node);
                    }

                    // du has one diagonal momentum coefficient per component
                    const scalar* duNode =
                        stk::mesh::field_data(*duLhsSTKFieldPtr, node);
                    for (label j = 0; j < SPATIAL_DIM; ++j)
                    {
                        p_du[ni * SPATIAL_DIM + j] = duNode[j];
                    }

                    const scalar* FOrigNode =
                        FOrigSTKFieldPtr
                            ? stk::mesh::field_data(*FOrigSTKFieldPtr, node)
                            : nullptr;
                    for (label j = 0; j < SPATIAL_DIM; ++j)
                    {
                        ws_FOrig[ni * SPATIAL_DIM + j] =
                            FOrigNode ? FOrigNode[j] : 0.0;
                    }
                }

                if (meshMoving)
                {
                    const scalar* Um =
                        stk::mesh::field_data(*UmSTKFieldPtr, node);
                    const label offSet = ni * SPATIAL_DIM;
                    for (label j = 0; j < SPATIAL_DIM; ++j)
                    {
                        p_Um[offSet + j] = Um[j];
                    }
                }
                else
                {
                    const label offSet = ni * SPATIAL_DIM;
                    for (label j = 0; j < SPATIAL_DIM; ++j)
                    {
                        p_Um[offSet + j] = 0.0;
                    }
                }
            }

            // pointer to face data
            const scalar* areaVec =
                stk::mesh::field_data(exposedAreaVecSTKFieldRef, side);
            const scalar* UbcVec =
                stk::mesh::field_data(sideUSTKFieldRef, side);
            const scalar* aWallBcVals =
                robin ? stk::mesh::field_data(*aWallBcSTKFieldPtr, side)
                      : nullptr;
            const scalar* rnWallFluxVals =
                robin ? stk::mesh::field_data(*rnWallFluxSTKFieldPtr, side)
                      : nullptr;
            scalar* rnWallFluxEVals =
                robin ? stk::mesh::field_data(*rnWallFluxESTKFieldPtr, side)
                      : nullptr;

            // extract the connected element to this exposed face; should be
            // single in size!
            const stk::mesh::Entity* faceElemRels =
                bulkData.begin_elements(side);
            STK_ThrowAssert(bulkData.num_elements(side) == 1);

            // get element; its face ordinal number and populate
            // face_node_ordinals
            stk::mesh::Entity element = faceElemRels[0];
            const stk::mesh::ConnectivityOrdinal* face_elem_ords =
                bulkData.begin_element_ordinals(side);
            const label faceOrdinal = face_elem_ords[0];

            // mapping from ip to nodes for this ordinal
            const label* ipNodeMap = meSCS->ipNodeMap(faceOrdinal);

            // populate faceNodeOrdinals
            const label* faceNodeOrdinals =
                meSCS->side_node_ordinals(faceOrdinal);

            //======================================
            // gather nodal data off of element
            //======================================
            stk::mesh::Entity const* elemNodeRels =
                bulkData.begin_nodes(element);
            label numNodes = bulkData.num_nodes(element);

            // sanity check on num nodes
            STK_ThrowAssert(numNodes == nodesPerElement);
            for (label ni = 0; ni < numNodes; ++ni)
            {
                stk::mesh::Entity node = elemNodeRels[ni];

                // set connected nodes
                connectedNodes[ni] = node;

                // gather vectors
                const scalar* coords =
                    stk::mesh::field_data(coordsSTKFieldRef, node);
                const label offSet = ni * SPATIAL_DIM;
                for (label j = 0; j < SPATIAL_DIM; ++j)
                {
                    p_coordinates[offSet + j] = coords[j];
                }
                if (robin)
                {
                    ws_elem_p[ni] = *stk::mesh::field_data(*pSTKFieldPtr, node);
                }
            }

            std::vector<scalar> faceIsoParCoords(SPATIAL_DIM - 1);
            std::vector<scalar> elemIsoParCoords(SPATIAL_DIM);
            std::vector<scalar> pGradBip(SPATIAL_DIM);

            // loop over boundary ips
            for (label ip = 0; ip < numScsBip; ++ip)
            {
                const label offSetSF_face = ip * nodesPerSide;

                const label nearestNode = ipNodeMap[ip];

                // zero out vector quantities; form aMag
                for (label j = 0; j < SPATIAL_DIM; ++j)
                {
                    p_umBip[j] = 0.0;
                    p_coordBip[j] = 0.0;
                    pGradBip[j] = 0.0;
                }

                // grad(p) at the BIP from the element shape functions
                if (robin)
                {
                    const label faceDim = SPATIAL_DIM - 1;
                    for (label j = 0; j < faceDim; ++j)
                        faceIsoParCoords[j] = static_cast<scalar>(
                            meFC->intgLoc_[ip * faceDim + j]);
                    meSCS->sidePcoords_to_elemPcoords(faceOrdinal,
                                                      1,
                                                      &faceIsoParCoords[0],
                                                      &elemIsoParCoords[0]);
                    scalar scs_error = 0.0;
                    meSCS->general_face_grad_op(faceOrdinal,
                                                &elemIsoParCoords[0],
                                                &p_coordinates[0],
                                                &ws_dndx[0],
                                                &ws_det_j[0],
                                                &scs_error);
                    for (label ic = 0; ic < nodesPerElement; ++ic)
                    {
                        const label offSetDnDx = ic * SPATIAL_DIM;
                        for (label j = 0; j < SPATIAL_DIM; ++j)
                            pGradBip[j] +=
                                ws_dndx[offSetDnDx + j] * ws_elem_p[ic];
                    }
                }

                // interpolate to bip; single phase has alpha = 1 exactly
                scalar rhoBip = 0;
                scalar alphaBip = phase.alpha ? 0.0 : 1.0;
                scalar rhoMixBip = 0;
                scalar duBipVec[SPATIAL_DIM] = {};
                scalar FOrigBip[SPATIAL_DIM] = {};
                for (label ic = 0; ic < nodesPerSide; ++ic)
                {
                    const label inn = faceNodeOrdinals[ic];

                    const scalar r_vel =
                        p_velocity_face_shape_function[offSetSF_face + ic];
                    const scalar r_coord =
                        p_coordinate_face_shape_function[offSetSF_face + ic];

                    // use velocity shape functions
                    rhoBip += r_vel * p_rho[ic];
                    if (phase.alpha)
                    {
                        alphaBip += r_vel * p_alpha[ic];
                    }
                    if (robin)
                    {
                        if (phase.rhoMix)
                        {
                            rhoMixBip += r_vel * ws_rho_mix[ic];
                        }
                        for (label j = 0; j < SPATIAL_DIM; ++j)
                        {
                            duBipVec[j] += r_vel * p_du[ic * SPATIAL_DIM + j];
                            FOrigBip[j] +=
                                r_vel * ws_FOrig[ic * SPATIAL_DIM + j];
                        }
                    }

                    const label icNdim = ic * SPATIAL_DIM;
                    for (label j = 0; j < SPATIAL_DIM; ++j)
                    {
                        // use coordinates shape functions
                        p_coordBip[j] +=
                            r_coord * p_coordinates[inn * SPATIAL_DIM + j];

                        // use velocity shape functions
                        p_umBip[j] += r_vel * p_Um[icNdim + j];
                    }
                }

                // boundary mass flux rho*u.A (Tukovic Eq. 29)
                scalar mDot = 0.0;
                for (label j = 0; j < SPATIAL_DIM; ++j)
                {
                    const scalar axj = areaVec[ip * SPATIAL_DIM + j];
                    mDot += (rhoBip * UbcVec[ip * SPATIAL_DIM + j]) * axj;
                }

                // transform mDot to relative frame

                // 1.) frame motion
                for (label i = 0; i < SPATIAL_DIM; ++i)
                {
                    for (label j = 0; j < SPATIAL_DIM; ++j)
                    {
                        mDot -= rhoBip * p_mat[i * SPATIAL_DIM + j] *
                                (p_coordBip[j] - p_ori[j]) *
                                areaVec[ip * SPATIAL_DIM + i];
                    }
                }

                // 2.) mesh motion
                for (label i = 0; i < SPATIAL_DIM; i++)
                {
                    mDot -= rhoBip * p_umBip[i] * areaVec[ip * SPATIAL_DIM + i];
                }

#ifndef NDEBUG
                // SCL check: accumulate grid flux to nearest node
                if (sclCheckSTKFieldPtr)
                {
                    scalar gridFlux = 0.0;
                    for (label i = 0; i < SPATIAL_DIM; i++)
                    {
                        gridFlux += rhoBip * p_umBip[i] *
                                    areaVec[ip * SPATIAL_DIM + i] * alphaBip /
                                    phase.densityScale;
                    }
                    stk::mesh::Entity node = elemNodeRels[nearestNode];
                    scalar* scl =
                        stk::mesh::field_data(*sclCheckSTKFieldPtr, node);
                    *scl -= gridFlux;
                }
#endif /* NDEBUG */

                // residual
                p_rhs[nearestNode] -= mDot * alphaBip / phase.densityScale;

                if (!robin)
                {
                    continue;
                }

                // Robin flux of the previous correctors (Tukovic Eq. 30)
                p_rhs[nearestNode] -=
                    rhoBip * rnWallFluxVals[ip] * alphaBip / phase.densityScale;

                // Robin row p' + c1*dp'/dn = -rho_s*h_s*a_n (Tukovic Eq. 20)
                const scalar rhoFluid =
                    robinFluidDensity(rnParams->useReferenceFluidDensity_,
                                      rnParams->referenceFluidDensity_,
                                      phase.rhoMix ? rhoMixBip : rhoBip);

                const robinBipCoefficients rc = computeRobinBipCoefficients(
                    &areaVec[ip * SPATIAL_DIM],
                    duBipVec,
                    p_coordBip,
                    &p_coordinates[nearestNode * SPATIAL_DIM],
                    rhoFluid,
                    robinRhoSolid,
                    robinHs);

                const uint64_t vnKey =
                    (static_cast<uint64_t>(bulkData.identifier(side)) << 8) |
                    static_cast<uint64_t>(ip);

                // G = rho*du_n*|A|/(deltaBn + c1), du_n as in the interior
                const scalar gradCoeff =
                    rhoBip * rc.du_n * rc.aMag / (rc.deltaBn + rc.c1);

                // rho_s*h_s*a_n (Tukovic Eq. 20)
                const scalar accelTerm =
                    robinRhoSolid * robinHs * aWallBcVals[ip];

                self.robinNeumann_.setRobinAccelTerm(vnKey, accelTerm);
                self.robinNeumann_.recordC1(rc.c1, rc.aMag);

                // phase weight; rho_k is already inside G
                const scalar w_k = alphaBip / phase.densityScale;

                // LHS: p'_BIP = sum N_ic p'_ic (Tukovic Eq. 26)
                const label rowR = nearestNode * nodesPerElement;
                for (label ic = 0; ic < nodesPerSide; ++ic)
                {
                    const label inn = faceNodeOrdinals[ic];
                    const scalar r_p =
                        p_velocity_face_shape_function[offSetSF_face + ic];
                    p_lhs[rowR + inn] += w_k * gradCoeff * r_p;
                }

                // RHS: frozen E = rho_s*h_s*a_n + c1*n.(grad p - F_orig)
                if (computeE)
                {
                    scalar nDotGradPOld = 0.0;
                    scalar nDotFOrig = 0.0;
                    for (label j = 0; j < SPATIAL_DIM; ++j)
                    {
                        nDotGradPOld += pGradBip[j] *
                                        areaVec[ip * SPATIAL_DIM + j] *
                                        rc.invAMag;
                        nDotFOrig += FOrigBip[j] *
                                     areaVec[ip * SPATIAL_DIM + j] * rc.invAMag;
                    }
                    rnWallFluxEVals[ip] =
                        accelTerm + rc.c1 * (nDotGradPOld - nDotFOrig);
                }
                p_rhs[nearestNode] -= w_k * gradCoeff * rnWallFluxEVals[ip];
            }

            self.applyCoeff_(
                A, b, connectedNodes, scratchIds, scratchVals, rhs, lhs);
        }
    }
}

template <class Assembler>
void accumulateRnWallFluxSide(Assembler& self,
                              const domain* domain,
                              const interfaceSideInfo* interfaceSideInfoPtr)
{
    auto* model = self.model_;
    auto& mesh = model->meshRef();

    stk::mesh::BulkData& bulkData = mesh.bulkDataRef();
    stk::mesh::MetaData& metaData = mesh.metaDataRef();

    // same fields as the Robin row; rhoRef() is the mixture density for VOF
    const auto& rhoSTKFieldRef = model->rhoRef().stkFieldRef();
    const bool consistent = model->controlsRef()
                                .solverRef()
                                .solverControl_.expertParameters_.consistent_;
    const auto& duLhsSTKFieldRef = *metaData.get_field<scalar>(
        stk::topology::NODE_RANK,
        consistent ? flowModel::duTilde_ID : flowModel::du_ID);
    // full (un-relaxed) p' of the corrector just solved (correctField_)
    const auto& pCorrSTKFieldRef = model->pCorrRef().stkFieldRef();
    auto& rnWallFluxSideSTKFieldRef = model->rnWallFluxRef().stkFieldRef();

    // geometric fields
    const auto& coordsSTKFieldRef = *metaData.get_field<scalar>(
        stk::topology::NODE_RANK, self.getCoordinatesID_(domain));
    const auto& exposedAreaVecSTKFieldRef = *metaData.get_field<scalar>(
        metaData.side_rank(), self.getExposedAreaVectorID_(domain));

    const scalar robinHs = self.robinNeumann_.robinHs();
    const scalar robinRhoSolid = self.robinNeumann_.robinRhoSolid();
    const auto& rnParams =
        domain->fluidSolidCouplingRef(interfaceSideInfoPtr->interfPtr()).robin_;
    const bool useReferenceFluidDensity = rnParams.useReferenceFluidDensity_;
    const scalar referenceFluidDensity = rnParams.referenceFluidDensity_;

    // shifted ip's for fields?
    const bool isUShifted = model->URef().isShifted();

    // workspace
    std::vector<stk::topology> parentTopo;
    std::vector<scalar> ws_rho;
    std::vector<scalar> ws_du;
    std::vector<scalar> ws_pCorr;
    std::vector<scalar> ws_coordinates;
    std::vector<scalar> ws_velocity_face_shape_function;
    std::vector<scalar> ws_coordinate_face_shape_function;
    std::vector<scalar> coordBip(SPATIAL_DIM);

    stk::mesh::Selector selAllSides =
        metaData.universal_part() &
        stk::mesh::selectUnion(interfaceSideInfoPtr->currentPartVec_);

    stk::mesh::BucketVector const& sideBuckets =
        bulkData.get_buckets(metaData.side_rank(), selAllSides);
    for (stk::mesh::BucketVector::const_iterator ib = sideBuckets.begin();
         ib != sideBuckets.end();
         ++ib)
    {
        stk::mesh::Bucket& sideBucket = **ib;

        // extract connected element topology
        sideBucket.parent_topology(stk::topology::ELEMENT_RANK, parentTopo);
        STK_ThrowAssert(parentTopo.size() == 1);
        stk::topology theElemTopo = parentTopo[0];

        // volume master element
        MasterElement* meSCS =
            MasterElementRepo::get_surface_master_element(theElemTopo);
        const label nodesPerElement = meSCS->nodesPerElement_;

        // face master element
        MasterElement* meFC = MasterElementRepo::get_surface_master_element(
            sideBucket.topology());
        const label nodesPerSide = sideBucket.topology().num_nodes();
        const label numScsBip = meFC->numIntPoints_;

        ws_rho.resize(nodesPerSide);
        ws_du.resize(nodesPerSide * SPATIAL_DIM);
        ws_pCorr.resize(nodesPerSide);
        ws_coordinates.resize(nodesPerElement * SPATIAL_DIM);
        ws_velocity_face_shape_function.resize(numScsBip * nodesPerSide);
        ws_coordinate_face_shape_function.resize(numScsBip * nodesPerSide);

        // shape functions; boundary (same as the Robin row)
        if (isUShifted)
        {
            meFC->shifted_shape_fcn(&ws_velocity_face_shape_function[0]);
        }
        else
        {
            meFC->shape_fcn(&ws_velocity_face_shape_function[0]);
        }

        // Always use trilinear (standard) shape functions for coordinates
        meFC->shape_fcn(&ws_coordinate_face_shape_function[0]);

        const stk::mesh::Bucket::size_type nSidesPerBucket = sideBucket.size();
        for (stk::mesh::Bucket::size_type iSide = 0; iSide < nSidesPerBucket;
             ++iSide)
        {
            stk::mesh::Entity side = sideBucket[iSide];

            // gather nodal data off of face (same node order as the kernel)
            stk::mesh::Entity const* sideNodeRels = bulkData.begin_nodes(side);
            STK_ThrowAssert(bulkData.num_nodes(side) == nodesPerSide);
            for (label ni = 0; ni < nodesPerSide; ++ni)
            {
                stk::mesh::Entity node = sideNodeRels[ni];
                ws_rho[ni] = *stk::mesh::field_data(rhoSTKFieldRef, node);
                const scalar* duNode =
                    stk::mesh::field_data(duLhsSTKFieldRef, node);
                for (label j = 0; j < SPATIAL_DIM; ++j)
                {
                    ws_du[ni * SPATIAL_DIM + j] = duNode[j];
                }
            }

            // connected element, face ordinal and node maps
            const stk::mesh::Entity* faceElemRels =
                bulkData.begin_elements(side);
            STK_ThrowAssert(bulkData.num_elements(side) == 1);
            stk::mesh::Entity element = faceElemRels[0];
            const stk::mesh::ConnectivityOrdinal* face_elem_ords =
                bulkData.begin_element_ordinals(side);
            const label faceOrdinal = face_elem_ords[0];
            const label* ipNodeMap = meSCS->ipNodeMap(faceOrdinal);
            const label* faceNodeOrdinals =
                meSCS->side_node_ordinals(faceOrdinal);

            stk::mesh::Entity const* elemNodeRels =
                bulkData.begin_nodes(element);
            STK_ThrowAssert(bulkData.num_nodes(element) == nodesPerElement);
            for (label ni = 0; ni < nodesPerElement; ++ni)
            {
                const scalar* coords =
                    stk::mesh::field_data(coordsSTKFieldRef, elemNodeRels[ni]);
                for (label j = 0; j < SPATIAL_DIM; ++j)
                {
                    ws_coordinates[ni * SPATIAL_DIM + j] = coords[j];
                }
            }

            // p' at the element nodes of the Robin row's LHS columns
            for (label ic = 0; ic < nodesPerSide; ++ic)
            {
                ws_pCorr[ic] = *stk::mesh::field_data(
                    pCorrSTKFieldRef, elemNodeRels[faceNodeOrdinals[ic]]);
            }

            const scalar* areaVec =
                stk::mesh::field_data(exposedAreaVecSTKFieldRef, side);
            scalar* rnWallFluxVals =
                stk::mesh::field_data(rnWallFluxSideSTKFieldRef, side);

            for (label ip = 0; ip < numScsBip; ++ip)
            {
                const label offSetSF_face = ip * nodesPerSide;
                const label nearestNode = ipNodeMap[ip];

                scalar rhoBip = 0.0;
                scalar pCorrBip = 0.0;
                scalar duBipVec[SPATIAL_DIM] = {};
                for (label j = 0; j < SPATIAL_DIM; ++j)
                {
                    coordBip[j] = 0.0;
                }

                for (label ic = 0; ic < nodesPerSide; ++ic)
                {
                    const label inn = faceNodeOrdinals[ic];
                    const scalar r_vel =
                        ws_velocity_face_shape_function[offSetSF_face + ic];
                    const scalar r_coord =
                        ws_coordinate_face_shape_function[offSetSF_face + ic];

                    rhoBip += r_vel * ws_rho[ic];
                    pCorrBip += r_vel * ws_pCorr[ic];
                    for (label j = 0; j < SPATIAL_DIM; ++j)
                    {
                        duBipVec[j] += r_vel * ws_du[ic * SPATIAL_DIM + j];
                        coordBip[j] +=
                            r_coord * ws_coordinates[inn * SPATIAL_DIM + j];
                    }
                }

                const scalar rhoFluid = robinFluidDensity(
                    useReferenceFluidDensity, referenceFluidDensity, rhoBip);

                const robinBipCoefficients rc = computeRobinBipCoefficients(
                    &areaVec[ip * SPATIAL_DIM],
                    duBipVec,
                    &coordBip[0],
                    &ws_coordinates[nearestNode * SPATIAL_DIM],
                    rhoFluid,
                    robinRhoSolid,
                    robinHs);

                // Q_vol += G/rho*p'_BIP; the kernel multiplies it by rho
                rnWallFluxVals[ip] +=
                    rc.du_n * rc.aMag / (rc.deltaBn + rc.c1) * pCorrBip;
            }
        }
    }
}

// Tukovic Eq. 30: accumulate Q_vol on the Robin-Neumann sides after a corrector
template <class Assembler>
void accumulateRnWallFluxSides(Assembler& self, const domain* domain)
{
    // the Robin row is only assembled once all fluid-solid sides are coupled
    if (!self.robinNeumann_.robinParamsValid() ||
        !self.model_->meshRef().interfacesCoupled(
            self.model_->controlsRef().time))
    {
        return;
    }

    for (const interface* interf : domain->interfacesRef())
    {
        if (interf->isConformalTreatment() || interf->isInternal())
            continue;

        if (!domain->isRobinNeumann(interf))
            continue;

        accumulateRnWallFluxSide(
            self, domain, interf->interfaceSideInfoPtr(domain->index()));
    }
}

} // namespace accel

#endif // ROBINNEUMANNWALL_H
