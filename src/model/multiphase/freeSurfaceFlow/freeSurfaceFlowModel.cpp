// File : freeSurfaceFlowModel.cpp
// Created : Sun Jan 26 2025 22:06:38 (+0100)
// Author : Mhamad Mahdi Alloush
// Description:
// Copyright (c) 2025 CCFNUM, Lucerne University of Applied Sciences and Arts.
// SPDX-License-Identifier: BSD-3-Clause

#include "freeSurfaceFlowModel.h"
#include "idealGasModel.h"
#include "simulation.h"
#include "sutherlandsFormulaModel.h"

namespace accel
{

freeSurfaceFlowModel::freeSurfaceFlowModel(realm* realm)
    : multiphaseModel(realm)
{
    // collect phases
    for (const auto& domain : realm->simulationRef().domainVector())
    {
        if (domain->multiphase_.homogeneous_ &&
            domain->multiphase_.freeSurfaceModel_.option_ ==
                freeSurfaceModelOption::standard)
        {
            for (const auto& material : domain->materialVector())
            {
                std::string materialName = material.name_;

                // global material index
                label phaseIndex =
                    realm->simulationRef().materialIndex(materialName);

                // check if phase already registered
                bool registered = false;
                for (label iPhase = 0; iPhase < nPhases(); iPhase++)
                {
                    if (phaseRef(iPhase).index_ == phaseIndex)
                    {
                        registered = true;
                        break;
                    }
                }

                if (!registered)
                {
                    // add phase
                    phases_.push_back(phase(phaseIndex, materialName));
                }
            }
        }
    }

    // store fourier number: required for smoothing of alpha
    Fo_ = controlsRef()
              .solverRef()
              .solverControl_.advancedOptions_.equationControls_
              .volumeFractionSmoothing_.fourierNumber_;

    // detect primary phases: phases that need not to be solved and are usually
    // ranked last in the set of material of a domain. These phases are rather
    // deduced from volume cosnervation principle
    for (label iPhase = 0; iPhase < nPhases(); iPhase++)
    {
        label phaseIndex = phaseRef(iPhase).index_;
        std::string phaseName = phaseRef(iPhase).name_;

        for (const auto& domain : realm->simulationRef().domainVector())
        {
            if (domain->multiphase_.homogeneous_ &&
                domain->multiphase_.freeSurfaceModel_.option_ ==
                    freeSurfaceModelOption::standard)
            {
                if (domain->hasMaterial(phaseName))
                {
                    // get phase ranking in domain
                    label phaseRankingInDomain;
                    for (phaseRankingInDomain = 0;
                         phaseRankingInDomain < domain->nMaterials();
                         phaseRankingInDomain++)
                    {
                        if (phaseIndex == domain->localToGlobalMaterialIndex(
                                              phaseRankingInDomain))
                        {
                            break;
                        }
                    }

                    // if not the last one then it is a secondary phase: to be
                    // solved for
                    if (phaseRankingInDomain < domain->nMaterials() - 1)
                    {
                        phaseRef(iPhase).primaryPhase_ = false;
                    }
                }
            }
        }
    }

    // create field instances
    for (label iPhase = 0; iPhase < nPhases(); iPhase++)
    {
        label phaseIndex = phaseRef(iPhase).index_;

        this->alphaRef(phaseIndex);
        this->rhoRef(phaseIndex);
        this->muRef(phaseIndex);
        this->mDotRef(phaseIndex);
        this->nHatRef(phaseIndex);

        // create alpha smooth fields
        if (controlsRef()
                .solverRef()
                .solverControl_.advancedOptions_.equationControls_
                .volumeFractionSmoothing_.smoothVolumeFraction_ &&
            !this->phaseRef(iPhase).primaryPhase_)
        {
            alphaSmoothRef(phaseIndex);
            rhsSmoothRef(phaseIndex);
        }
    }

    // setup per-pair curvature fields for surface tension (CSF model)
    stk::mesh::MetaData& metaData = this->meshRef().metaDataRef();
    for (const auto& domain : realm->simulationRef().domainVector())
    {
        for (const auto& fpm : domain->fluidPairModels_)
        {
            if (fpm.surfaceTension_.option_ != surfaceTensionModelOption::none)
            {
                const std::string kappaName =
                    "curvature." + fpm.materialA_ + "_" + fpm.materialB_;

                // only declare once (shared across domains on same mesh)
                if (kappaSTKFieldPtrs_.find(kappaName) ==
                    kappaSTKFieldPtrs_.end())
                {
                    auto* kappaPtr = &metaData.declare_field<scalar>(
                        stk::topology::NODE_RANK, kappaName);

                    stk::io::set_field_output_type(*kappaPtr, fieldType[1]);

                    kappaSTKFieldPtrs_[kappaName] = kappaPtr;
                }

                // put field on interior parts for this domain
                STKScalarField* kappaPtr = kappaSTKFieldPtrs_[kappaName];
                const stk::mesh::PartVector& partVec =
                    domain->zonePtr()->interiorParts();
                for (const stk::mesh::Part* part : partVec)
                {
                    if (!kappaPtr->defined_on(*part))
                    {
                        stk::mesh::put_field_on_mesh(*kappaPtr, *part, nullptr);
                    }
                }
            }
        }
    }
}

void freeSurfaceFlowModel::computeCurvature_(
    const std::shared_ptr<domain> domain,
    label iPhase,
    STKScalarField* kappaFieldPtr)
{
    // Compute kappa = -div(nHat) at nodes using SCS-based divergence
    auto& mesh = this->meshRef();
    stk::mesh::BulkData& bulkData = mesh.bulkDataRef();
    stk::mesh::MetaData& metaData = mesh.metaDataRef();

    // get interior parts
    const stk::mesh::PartVector& partVec = domain->zonePtr()->interiorParts();

    // zero out curvature field
    ops::zero(kappaFieldPtr, partVec);

    // get geometry
    const auto* coordsSTKFieldPtr = metaData.get_field<scalar>(
        stk::topology::NODE_RANK, this->getCoordinatesID_(domain));

    const auto* volSTKFieldPtr = metaData.get_field<scalar>(
        stk::topology::NODE_RANK, this->getDualNodalVolumeID_(domain));

    // get nHat field for the given phase
    const STKScalarField* nHatSTKFieldPtr = this->nHatRef(iPhase).stkFieldPtr();

    // ========================================
    // Interior IP contribution
    // ========================================
    {
        // workspace
        std::vector<scalar> ws_coordinates;
        std::vector<scalar> ws_nHat;
        std::vector<scalar> ws_dualVolume;
        std::vector<scalar> ws_scsAreav;
        std::vector<scalar> ws_shape_function;

        // integration point data that depends on size
        std::vector<scalar> nHatIp(SPATIAL_DIM);

        // pointers to everyone...
        scalar* p_nHatIp = &nHatIp[0];

        stk::mesh::Selector selAllElements =
            metaData.universal_part() & stk::mesh::selectUnion(partVec);

        stk::mesh::BucketVector const& elementBuckets =
            bulkData.get_buckets(stk::topology::ELEMENT_RANK, selAllElements);

        for (auto ib = elementBuckets.begin(); ib != elementBuckets.end(); ++ib)
        {
            stk::mesh::Bucket& elementBucket = **ib;
            const auto nElemPerBucket = elementBucket.size();

            // extract master element SCS
            MasterElement* meSCS =
                MasterElementRepo::get_surface_master_element(
                    elementBucket.topology());
            const label nodesPerElement = meSCS->nodesPerElement_;
            const label numScsIp = meSCS->numIntPoints_;
            const label* lrscv = meSCS->adjacentNodes();

            // resize workspaces
            ws_coordinates.resize(nodesPerElement * SPATIAL_DIM);
            ws_nHat.resize(nodesPerElement * SPATIAL_DIM);
            ws_dualVolume.resize(nodesPerElement);
            ws_scsAreav.resize(numScsIp * SPATIAL_DIM);
            ws_shape_function.resize(numScsIp * nodesPerElement);

            // get shape functions at SCS integration points
            meSCS->shape_fcn(&ws_shape_function[0]);

            for (stk::mesh::Bucket::size_type k = 0; k < nElemPerBucket; ++k)
            {
                stk::mesh::Entity const* nodeRels =
                    elementBucket.begin_nodes(k);
                const label numNodes = elementBucket.num_nodes(k);
                STK_ThrowAssert(numNodes == nodesPerElement);

                // gather nodal data
                for (label ni = 0; ni < nodesPerElement; ++ni)
                {
                    stk::mesh::Entity node = nodeRels[ni];

                    const scalar* coords =
                        stk::mesh::field_data(*coordsSTKFieldPtr, node);
                    const scalar* nHat =
                        stk::mesh::field_data(*nHatSTKFieldPtr, node);
                    ws_dualVolume[ni] =
                        *stk::mesh::field_data(*volSTKFieldPtr, node);

                    for (label d = 0; d < SPATIAL_DIM; ++d)
                    {
                        ws_coordinates[ni * SPATIAL_DIM + d] = coords[d];
                        ws_nHat[ni * SPATIAL_DIM + d] = nHat[d];
                    }
                }

                // compute SCS area vectors
                scalar scsError = 0.0;
                meSCS->determinant(
                    1, &ws_coordinates[0], &ws_scsAreav[0], &scsError);

                // loop over SCS integration points
                for (label ip = 0; ip < numScsIp; ++ip)
                {
                    const label il = lrscv[2 * ip];
                    const label ir = lrscv[2 * ip + 1];

                    // zero-out
                    for (label j = 0; j < SPATIAL_DIM; ++j)
                    {
                        p_nHatIp[j] = 0.0;
                    }

                    // interpolate nHat to SCS ip
                    for (label ni = 0; ni < nodesPerElement; ++ni)
                    {
                        const scalar r =
                            ws_shape_function[ip * nodesPerElement + ni];
                        for (label d = 0; d < SPATIAL_DIM; ++d)
                        {
                            p_nHatIp[d] += r * ws_nHat[ni * SPATIAL_DIM + d];
                        }
                    }

                    // compute flux = nHat_ip . A_ip (dot product)
                    scalar flux = 0.0;
                    for (label d = 0; d < SPATIAL_DIM; ++d)
                    {
                        flux += p_nHatIp[d] * ws_scsAreav[ip * SPATIAL_DIM + d];
                    }

                    // kappa = -div(nHat): assemble as negative divergence
                    stk::mesh::Entity nodeL = nodeRels[il];
                    stk::mesh::Entity nodeR = nodeRels[ir];

                    scalar* kappaL =
                        stk::mesh::field_data(*kappaFieldPtr, nodeL);
                    scalar* kappaR =
                        stk::mesh::field_data(*kappaFieldPtr, nodeR);

                    *kappaL -= flux / ws_dualVolume[il];
                    *kappaR += flux / ws_dualVolume[ir];
                }
            }
        }
    }

    // ========================================
    // Boundary IP contribution
    // ========================================
    {
        // get fields
        const auto& exposedAreaVecSTKFieldRef = *metaData.get_field<scalar>(
            metaData.side_rank(), this->getExposedAreaVectorID_(domain));

        // scratch arrays
        std::vector<scalar> ws_nHat;
        std::vector<scalar> ws_shape_function;

        // fixed-size arrays
        std::vector<scalar> nHatIp(SPATIAL_DIM);

        // pointers
        scalar* p_nHatIp = &nHatIp[0];

        std::vector<stk::topology> parentTopo;

        for (label iBoundary = 0; iBoundary < domain->zonePtr()->nBoundaries();
             iBoundary++)
        {
            stk::mesh::PartVector partVec =
                domain->zonePtr()->boundaryRef(iBoundary).parts();

            stk::mesh::Selector selAllSides =
                metaData.universal_part() & stk::mesh::selectUnion(partVec);

            stk::mesh::BucketVector const& sideBuckets =
                bulkData.get_buckets(metaData.side_rank(), selAllSides);
            for (stk::mesh::BucketVector::const_iterator ib =
                     sideBuckets.begin();
                 ib != sideBuckets.end();
                 ++ib)
            {
                stk::mesh::Bucket& sideBucket = **ib;
                const stk::mesh::Bucket::size_type nSidesPerBucket =
                    sideBucket.size();

                // extract connected element topology
                sideBucket.parent_topology(stk::topology::ELEMENT_RANK,
                                           parentTopo);
                stk::topology theElemTopo = parentTopo[0];

                // volume master element
                MasterElement* meSCS =
                    MasterElementRepo::get_surface_master_element(theElemTopo);

                // face master element
                MasterElement* meFC =
                    MasterElementRepo::get_surface_master_element(
                        sideBucket.topology());

                const label nodesPerSide = sideBucket.topology().num_nodes();
                const label numScsBip = meFC->numIntPoints_;

                ws_nHat.resize(nodesPerSide * SPATIAL_DIM);
                ws_shape_function.resize(numScsBip * nodesPerSide);

                scalar* p_nHat = &ws_nHat[0];
                scalar* p_shape_function = &ws_shape_function[0];

                meFC->shape_fcn(&p_shape_function[0]);

                for (stk::mesh::Bucket::size_type iSide = 0;
                     iSide < nSidesPerBucket;
                     ++iSide)
                {
                    // get face
                    stk::mesh::Entity side = sideBucket[iSide];

                    // gather nHat from face nodes
                    stk::mesh::Entity const* sideNodeRels =
                        bulkData.begin_nodes(side);
                    const label numSideNodes = bulkData.num_nodes(side);

                    for (label ni = 0; ni < numSideNodes; ++ni)
                    {
                        stk::mesh::Entity node = sideNodeRels[ni];
                        const scalar* nHat =
                            stk::mesh::field_data(*nHatSTKFieldPtr, node);

                        const label offSet = ni * SPATIAL_DIM;
                        for (label j = 0; j < SPATIAL_DIM; ++j)
                        {
                            p_nHat[offSet + j] = nHat[j];
                        }
                    }

                    // face area vector
                    const scalar* areaVec =
                        stk::mesh::field_data(exposedAreaVecSTKFieldRef, side);

                    // get connected element and face
                    // ordinal
                    const stk::mesh::Entity* faceElemRels =
                        bulkData.begin_elements(side);
                    stk::mesh::Entity element = faceElemRels[0];
                    const stk::mesh::ConnectivityOrdinal* face_elem_ords =
                        bulkData.begin_element_ordinals(side);
                    const label faceOrdinal = face_elem_ords[0];

                    // mapping from ip to nodes for this
                    // ordinal
                    const label* ipNodeMap = meSCS->ipNodeMap(faceOrdinal);

                    // element node relations
                    stk::mesh::Entity const* elemNodeRels =
                        bulkData.begin_nodes(element);

                    // start assembly
                    for (label ip = 0; ip < numScsBip; ++ip)
                    {
                        // nearest node
                        const label nearestNode = ipNodeMap[ip];
                        stk::mesh::Entity node = elemNodeRels[nearestNode];

                        scalar* kappa =
                            stk::mesh::field_data(*kappaFieldPtr, node);
                        scalar vol =
                            *stk::mesh::field_data(*volSTKFieldPtr, node);

                        // interpolate nHat to IP
                        for (label j = 0; j < SPATIAL_DIM; ++j)
                        {
                            p_nHatIp[j] = 0.0;
                        }

                        const label offset = ip * nodesPerSide;
                        for (label ic = 0; ic < nodesPerSide; ++ic)
                        {
                            const scalar r = p_shape_function[offset + ic];
                            for (label j = 0; j < SPATIAL_DIM; ++j)
                            {
                                p_nHatIp[j] += r * p_nHat[ic * SPATIAL_DIM + j];
                            }
                        }

                        // compute flux and accumulate
                        scalar flux = 0.0;
                        for (label j = 0; j < SPATIAL_DIM; ++j)
                        {
                            flux += p_nHatIp[j] * areaVec[ip * SPATIAL_DIM + j];
                        }

                        *kappa -= flux / vol;
                    }
                }
            }
        }
    }

    // parallel communication
    if (messager::parallel())
    {
        stk::mesh::communicate_field_data(bulkData, {kappaFieldPtr});
    }
}

void freeSurfaceFlowModel::computeBodyForces(
    const std::shared_ptr<domain> domain)
{
    // Compute buoyancy, Lorentz, uniform body forces + redistribution
    flowModel::computeBodyForces(domain);

    // Add CSF surface tension forces for each fluid pair
    // (added after redistribution: surface tension is already interface-local)

    auto& mesh = this->meshRef();
    stk::mesh::BulkData& bulkData = mesh.bulkDataRef();
    stk::mesh::MetaData& metaData = mesh.metaDataRef();

    // get interior parts
    const stk::mesh::PartVector& partVec = domain->zonePtr()->interiorParts();

    for (const auto& fpm : domain->fluidPairModels_)
    {
        if (fpm.surfaceTension_.option_ !=
            surfaceTensionModelOption::continuumSurfaceForce)
            continue;

        // Use materialA as the phase for curvature/gradient
        const label phaseIndex = fpm.globalIndexA_;

        // Retrieve smoothing boolean from user input
        const bool smooth =
            controlsRef()
                .solverRef()
                .solverControl_.advancedOptions_.equationControls_
                .volumeFractionSmoothing_.smoothVolumeFraction_;

        // get geometry
        const auto* coordsSTKFieldPtr = metaData.get_field<scalar>(
            stk::topology::NODE_RANK, this->getCoordinatesID_(domain));

        const auto* volSTKFieldPtr = metaData.get_field<scalar>(
            stk::topology::NODE_RANK, this->getDualNodalVolumeID_(domain));

        // Get alpha field (smoothed or raw)
        const STKScalarField* alphaSTKFieldPtr =
            smooth ? this->alphaSmoothRef(phaseIndex).stkFieldPtr()
                   : this->alphaRef(phaseIndex).stkFieldPtr();

        const scalar sigma = fpm.surfaceTension_.coefficient_;

        // Look up the per-pair curvature field
        const std::string kappaName =
            "curvature." + fpm.materialA_ + "_" + fpm.materialB_;
        auto it = kappaSTKFieldPtrs_.find(kappaName);
        STK_ThrowAssert(it != kappaSTKFieldPtrs_.end());
        STKScalarField* kappaFieldPtr = it->second;

        // Compute curvature for this pair
        computeCurvature_(domain, phaseIndex, kappaFieldPtr);

        // workspace
        std::vector<scalar> ws_coordinates;
        std::vector<scalar> ws_alpha;
        std::vector<scalar> ws_kappa;
        std::vector<scalar> ws_dualVolume;
        std::vector<scalar> ws_scVolume;
        std::vector<scalar> ws_dndx;
        std::vector<scalar> ws_deriv;
        std::vector<scalar> ws_det_j;
        std::vector<scalar> ws_shape_function;

        // ip values
        std::vector<scalar> gradAlphaIp(SPATIAL_DIM);

        // pointers for fast access
        scalar* p_gradAlphaIp = &gradAlphaIp[0];

        stk::mesh::Selector selAllElements =
            metaData.universal_part() & stk::mesh::selectUnion(partVec);

        stk::mesh::BucketVector const& elementBuckets =
            bulkData.get_buckets(stk::topology::ELEMENT_RANK, selAllElements);

        for (auto ib = elementBuckets.begin(); ib != elementBuckets.end(); ++ib)
        {
            stk::mesh::Bucket& elementBucket = **ib;
            const auto nElemPerBucket = elementBucket.size();

            // extract master element SCV
            MasterElement* meSCV = MasterElementRepo::get_volume_master_element(
                elementBucket.topology());
            const label nodesPerElement = meSCV->nodesPerElement_;
            const label numScvIp = meSCV->numIntPoints_;
            const label* ipNodeMap = meSCV->ipNodeMap();

            // resize workspaces
            ws_coordinates.resize(nodesPerElement * SPATIAL_DIM);
            ws_alpha.resize(nodesPerElement);
            ws_kappa.resize(nodesPerElement);
            ws_dualVolume.resize(nodesPerElement);
            ws_scVolume.resize(numScvIp);
            ws_dndx.resize(SPATIAL_DIM * numScvIp * nodesPerElement);
            ws_deriv.resize(SPATIAL_DIM * numScvIp * nodesPerElement);
            ws_det_j.resize(numScvIp);
            ws_shape_function.resize(numScvIp * nodesPerElement);

            // get shape functions at SCV integration points
            meSCV->shape_fcn(&ws_shape_function[0]);

            for (stk::mesh::Bucket::size_type k = 0; k < nElemPerBucket; ++k)
            {
                stk::mesh::Entity const* nodeRels =
                    elementBucket.begin_nodes(k);
                const label numNodes = elementBucket.num_nodes(k);
                STK_ThrowAssert(numNodes == nodesPerElement);

                // gather nodal data
                for (label ni = 0; ni < nodesPerElement; ++ni)
                {
                    stk::mesh::Entity node = nodeRels[ni];

                    const scalar* coords =
                        stk::mesh::field_data(*coordsSTKFieldPtr, node);
                    ws_alpha[ni] =
                        *stk::mesh::field_data(*alphaSTKFieldPtr, node);
                    ws_kappa[ni] = *stk::mesh::field_data(*kappaFieldPtr, node);
                    ws_dualVolume[ni] =
                        *stk::mesh::field_data(*volSTKFieldPtr, node);

                    for (label d = 0; d < SPATIAL_DIM; ++d)
                    {
                        ws_coordinates[ni * SPATIAL_DIM + d] = coords[d];
                    }
                }

                // compute SCV volumes and grad operator
                scalar scv_error = 0.0;
                meSCV->determinant(
                    1, &ws_coordinates[0], &ws_scVolume[0], &scv_error);
                meSCV->grad_op(1,
                               &ws_coordinates[0],
                               &ws_dndx[0],
                               &ws_deriv[0],
                               &ws_det_j[0],
                               &scv_error);

                // loop over SCV integration points
                for (label ip = 0; ip < numScvIp; ++ip)
                {
                    const label nn = ipNodeMap[ip];

                    // interpolate kappa to ip using shape functions
                    scalar kappaIp = 0.0;
                    for (label ni = 0; ni < nodesPerElement; ++ni)
                    {
                        kappaIp +=
                            ws_shape_function[ip * nodesPerElement + ni] *
                            ws_kappa[ni];
                    }

                    // zero-out
                    for (label d = 0; d < SPATIAL_DIM; ++d)
                    {
                        p_gradAlphaIp[d] = 0;
                    }

                    // compute grad(alpha) at ip using dndx
                    for (label ni = 0; ni < nodesPerElement; ++ni)
                    {
                        const label offsetDnDx =
                            SPATIAL_DIM * nodesPerElement * ip +
                            ni * SPATIAL_DIM;
                        for (label d = 0; d < SPATIAL_DIM; ++d)
                        {
                            p_gradAlphaIp[d] +=
                                ws_dndx[offsetDnDx + d] * ws_alpha[ni];
                        }
                    }

                    // nearest node
                    stk::mesh::Entity nearestNode = nodeRels[nn];
                    scalar* F_node = stk::mesh::field_data(
                        *FOriginalSTKFieldPtr_, nearestNode);

                    // F_st = sigma * kappa * grad(alpha) * V_scv / V_dual
                    const scalar dualVol = ws_dualVolume[nn];
                    const scalar weight = ws_scVolume[ip] / dualVol;
                    for (label d = 0; d < SPATIAL_DIM; ++d)
                    {
                        F_node[d] +=
                            weight * sigma * kappaIp * p_gradAlphaIp[d];
                    }
                }
            }
        }
    }

    // Sync F from owners: This is needed because computeBodyForces accumulates
    // CSF surface tension via element-to-node scatter, and ghost nodes may not
    // have received contributions from all surrounding elements.
    if (messager::parallel())
    {
        stk::mesh::communicate_field_data(bulkData, {FOriginalSTKFieldPtr_});
    }

    // Copy FOrig values to F
    ops::copy<scalar>(FOriginalSTKFieldPtr_, FSTKFieldPtr_, partVec);
}

void freeSurfaceFlowModel::redistributeBodyForces(
    const std::shared_ptr<domain> domain)
{
    // Harmonic density-weighted body force redistribution for free surface:
    //
    // For free surface flows, standard volume-weighted averaging smears
    // the body force across the density interface, creating spurious pressure
    // gradients and parasitic currents.
    //
    // Instead, we average the specific body force (F/rho) to the element
    // centre using SCV volume weights, then multiply by the element-centre
    // density to recover the body force:
    //
    // Step 1: Average to element centre with harmonic density weighting:
    // rho_el = Σ_i (rho_i * V_scv_i) / V_el
    // B_el = rho_el * Σ_i (B_i/rho_i * V_scv_i) / V_el
    //
    // Step 2: Scatter the same element-centre value back to all nodes,
    // weighted by SCV volumes:
    // F_node += B_el * V_scv_node
    //
    // Step 3: Normalize by dual nodal volume:
    // F_node /= V_dual
    if (this->controlsRef()
            .solverRef()
            .solverControl_.expertParameters_.bodyForceRedistribution_)
    {
        auto& mesh = this->meshRef();
        stk::mesh::BulkData& bulkData = mesh.bulkDataRef();
        stk::mesh::MetaData& metaData = mesh.metaDataRef();

        // get interior parts
        const stk::mesh::PartVector& partVec =
            domain->zonePtr()->interiorParts();

        // Reset F field to zero
        ops::zero(FSTKFieldPtr_, partVec);

        // Get coordinates field
        const auto& coordinatesRef = *metaData.get_field<scalar>(
            stk::topology::NODE_RANK, this->getCoordinatesID_(domain));

        // Get dual nodal volume field for normalization
        const auto* dualNodalVolumeSTKFieldPtr = metaData.get_field<scalar>(
            stk::topology::NODE_RANK, this->getDualNodalVolumeID_(domain));

        // Get density field for harmonic weighting
        const auto* rhoSTKFieldPtr = this->rhoRef().stkFieldPtr();

        // Workspace arrays (variable size, resized per bucket)
        std::vector<scalar> ws_F;
        std::vector<scalar> ws_rho;
        std::vector<scalar> ws_coordinates;
        std::vector<scalar> ws_scv_volume;

        // Fixed-size workspace arrays (reused in inner loop)
        std::vector<scalar> F_el(SPATIAL_DIM);

        // Define selectors
        stk::mesh::Selector selAllElements =
            metaData.universal_part() & stk::mesh::selectUnion(partVec);
        stk::mesh::BucketVector const& elementBuckets =
            bulkData.get_buckets(stk::topology::ELEMENT_RANK, selAllElements);

        for (auto ib = elementBuckets.begin(); ib != elementBuckets.end(); ++ib)
        {
            stk::mesh::Bucket& elementBucket = **ib;
            const auto nElementsPerBucket = elementBucket.size();

            // Extract master element for SCV (volume element)
            MasterElement* meSCV = MasterElementRepo::get_volume_master_element(
                elementBucket.topology());
            const label nodesPerElement = meSCV->nodesPerElement_;
            const label numScvIp = meSCV->numIntPoints_;
            const label* ipNodeMap = meSCV->ipNodeMap();

            // Resize workspace arrays
            ws_F.resize(nodesPerElement * SPATIAL_DIM);
            ws_rho.resize(nodesPerElement);
            ws_coordinates.resize(nodesPerElement * SPATIAL_DIM);
            ws_scv_volume.resize(numScvIp);

            // Loop over elements in bucket
            for (stk::mesh::Bucket::size_type k = 0; k < nElementsPerBucket;
                 ++k)
            {
                stk::mesh::Entity const* nodeRels =
                    elementBucket.begin_nodes(k);
                const label numNodes = elementBucket.num_nodes(k);

                // Gather nodal data
                for (label ni = 0; ni < numNodes; ++ni)
                {
                    stk::mesh::Entity node = nodeRels[ni];

                    const scalar* coords =
                        stk::mesh::field_data(coordinatesRef, node);
                    const scalar* FOrig =
                        stk::mesh::field_data(*FOriginalSTKFieldPtr_, node);
                    for (label d = 0; d < SPATIAL_DIM; ++d)
                    {
                        ws_coordinates[ni * SPATIAL_DIM + d] = coords[d];
                        ws_F[ni * SPATIAL_DIM + d] = FOrig[d];
                    }

                    const scalar* rho =
                        stk::mesh::field_data(*rhoSTKFieldPtr, node);
                    ws_rho[ni] = rho[0];
                }

                // Compute SCV volumes
                scalar scv_error = 0.0;
                meSCV->determinant(
                    1, &ws_coordinates[0], &ws_scv_volume[0], &scv_error);

                // Compute element volume
                scalar V_el = 0.0;
                for (label ip = 0; ip < numScvIp; ++ip)
                {
                    V_el += ws_scv_volume[ip];
                }

                const scalar invV_el = 1.0 / V_el;

                // Compute volume-weighted element-centre density:
                // rho_el = Σ_i (rho_i * V_scv_i) / V_el
                scalar rho_el = 0.0;
                for (label ip = 0; ip < numScvIp; ++ip)
                {
                    const label nn = ipNodeMap[ip];
                    rho_el += ws_rho[nn] * ws_scv_volume[ip];
                }
                rho_el *= invV_el;

                // Compute harmonic density-weighted element-centre body force:
                // B_el = rho_el * Σ_i (B_i/rho_i * V_scv_i) / V_el
                for (label d = 0; d < SPATIAL_DIM; ++d)
                {
                    F_el[d] = 0.0;
                }

                for (label ip = 0; ip < numScvIp; ++ip)
                {
                    const label nn = ipNodeMap[ip];
                    const scalar rho_node = ws_rho[nn];
                    const scalar invRho = 1.0 / rho_node;
                    const scalar scV = ws_scv_volume[ip];
                    for (label d = 0; d < SPATIAL_DIM; ++d)
                    {
                        F_el[d] += ws_F[nn * SPATIAL_DIM + d] * invRho * scV;
                    }
                }

                for (label d = 0; d < SPATIAL_DIM; ++d)
                {
                    F_el[d] *= rho_el * invV_el;
                }

                // Scatter element-centre value back to nodes weighted by
                // SCV volumes: F_node += B_el * V_scv
                for (label ip = 0; ip < numScvIp; ++ip)
                {
                    const label nn = ipNodeMap[ip];
                    stk::mesh::Entity nearestNode = nodeRels[nn];
                    const scalar scV = ws_scv_volume[ip];

                    scalar* F_node =
                        stk::mesh::field_data(*FSTKFieldPtr_, nearestNode);

                    for (label d = 0; d < SPATIAL_DIM; ++d)
                    {
                        F_node[d] += F_el[d] * scV;
                    }
                }
            }
        }

        // Parallel communication
        if (messager::parallel())
        {
            stk::mesh::communicate_field_data(bulkData, {FSTKFieldPtr_});
        }

        // Normalize F by dual nodal volume: F /= dual_nodal_volume
        stk::mesh::Selector selAllNodes =
            metaData.universal_part() & stk::mesh::selectUnion(partVec);
        stk::mesh::BucketVector const& nodeBuckets =
            bulkData.get_buckets(stk::topology::NODE_RANK, selAllNodes);
        for (auto ib = nodeBuckets.begin(); ib != nodeBuckets.end(); ++ib)
        {
            stk::mesh::Bucket& nodeBucket = **ib;
            const auto nNodesPerBucket = nodeBucket.size();
            scalar* Fb = stk::mesh::field_data(*FSTKFieldPtr_, nodeBucket);
            const scalar* dualVolb =
                stk::mesh::field_data(*dualNodalVolumeSTKFieldPtr, nodeBucket);

            for (stk::mesh::Bucket::size_type iNode = 0;
                 iNode < nNodesPerBucket;
                 ++iNode)
            {
                const scalar dualVol = dualVolb[iNode];
                const scalar invDualVol = 1.0 / dualVol;
                for (label d = 0; d < SPATIAL_DIM; ++d)
                {
                    Fb[SPATIAL_DIM * iNode + d] *= invDualVol;
                }
            }
        }
    }
}

void freeSurfaceFlowModel::transformMassFlowRateToRelative(
    const std::shared_ptr<domain> domain,
    label iPhase)
{
    transformMassFlowRateToRelative_(
        domain, this->mDotRef(iPhase), this->rhoRef(iPhase));
}

void freeSurfaceFlowModel::transformMassFlowRateToAbsolute(
    const std::shared_ptr<domain> domain,
    label iPhase)
{
    transformMassFlowRateToAbsolute_(
        domain, this->mDotRef(iPhase), this->rhoRef(iPhase));
}

void freeSurfaceFlowModel::updateFlowReversalFlag(
    const std::shared_ptr<domain> domain,
    label iPhase)
{
    updateFlowReversalFlag_(domain, this->mDotRef(iPhase).sideFieldRef());
}

void freeSurfaceFlowModel::updateMassDivergenceField(
    const std::shared_ptr<domain> domain,
    label iPhase)
{
    updateMassDivergenceField_(domain,
                               this->mDotRef(iPhase),
                               this->rhoRef(iPhase),
                               this->mDotRef(iPhase).divRef());
}

void freeSurfaceFlowModel::updateSideMassFlowRateFraction(
    const std::shared_ptr<domain> domain,
    label iPhase)
{
    this->updateSideMassFlowRateFraction_(
        domain, this->rhoRef(iPhase), this->mDotRef(iPhase));
}

// Setups

void freeSurfaceFlowModel::setupVolumeFraction(
    const std::shared_ptr<domain> domain,
    label iPhase)
{
    // raw setup
    fieldBroker::setupVolumeFraction(domain, iPhase);

    // enable interface normal
    nHatRef(iPhase).setZone(domain->index());

    // also setup alpha smooth and rhs if required
    if (controlsRef()
            .solverRef()
            .solverControl_.advancedOptions_.equationControls_
            .volumeFractionSmoothing_.smoothVolumeFraction_ &&
        !this->phaseRef(iPhase).primaryPhase_)
    {
        if (alphaSmoothRef(iPhase).isZoneUnset(domain->index()))
        {
            alphaSmoothRef(iPhase).setZone(domain->index());
            rhsSmoothRef(iPhase).setZone(domain->index());
        }
    }
}

void freeSurfaceFlowModel::setupDensity(const std::shared_ptr<domain> domain)
{
    // this is the mixture (bulk) density
    if (rhoRef().isZoneUnset(domain->index()))
    {
        rhoRef().setZone(domain->index());
    }
}

void freeSurfaceFlowModel::setupDynamicViscosity(
    const std::shared_ptr<domain> domain)
{
    // this is the mixture (bulk) viscosity
    if (muRef().isZoneUnset(domain->index()))
    {
        muRef().setZone(domain->index());
    }
}

void freeSurfaceFlowModel::setupSpecificHeatCapacity(
    const std::shared_ptr<domain> domain)
{
    // this is the mixture (bulk) specific heat capacity
    if (cpRef().isZoneUnset(domain->index()))
    {
        cpRef().setZone(domain->index());
    }
}

void freeSurfaceFlowModel::setupThermalConductivity(
    const std::shared_ptr<domain> domain)
{
    // this is the mixture (bulk) thermal conductivity
    if (lambdaRef().isZoneUnset(domain->index()))
    {
        lambdaRef().setZone(domain->index());
    }
}

void freeSurfaceFlowModel::initializeDensity(
    const std::shared_ptr<domain> domain)
{
    updateDensity(domain);
}

void freeSurfaceFlowModel::initializeDensity(
    const std::shared_ptr<domain> domain,
    label iPhase)
{
    updateDensity(domain, iPhase);
}

void freeSurfaceFlowModel::initializeDynamicViscosity(
    const std::shared_ptr<domain> domain)
{
    updateDynamicViscosity(domain);
}

void freeSurfaceFlowModel::initializeDynamicViscosity(
    const std::shared_ptr<domain> domain,
    label iPhase)
{
    updateDynamicViscosity(domain, iPhase);
}

void freeSurfaceFlowModel::initializeMassFlowRate(
    const std::shared_ptr<domain> domain)
{
    // bulk mass flux must not be initialized from rho and velocity ... but
    // rather updated from accumulating phasic mass fluxes
    updateMassFlowRate(domain);
}

void freeSurfaceFlowModel::updateDensity(const std::shared_ptr<domain> domain)
{
    auto& mesh = this->meshRef();
    stk::mesh::BulkData& bulkData = mesh.bulkDataRef();
    stk::mesh::MetaData& metaData = mesh.metaDataRef();

    // Get pointer to global density field to be used
    const auto* bulkRhoSTKFieldPtr = rhoRef().stkFieldPtr();

    // get interior parts the domain is defined on
    const stk::mesh::PartVector& partVec = domain->zonePtr()->interiorParts();

    // Initialize global density to zero before every update
    rhoRef().setToValue({0}, partVec);

    // define some common selectors; select owned nodes
    stk::mesh::Selector selUniversalNodes =
        metaData.universal_part() & stk::mesh::selectUnion(partVec);

    stk::mesh::BucketVector const& nodeBuckets =
        bulkData.get_buckets(stk::topology::NODE_RANK, selUniversalNodes);

    for (label iPhase = 0; iPhase < domain->nMaterials(); iPhase++)
    {
        label phaseIndex = domain->localToGlobalMaterialIndex(iPhase);

        // Get fields for a given iPhase
        const STKScalarField* rhoSTKFieldPtr =
            this->rhoRef(phaseIndex).stkFieldPtr();
        const STKScalarField* alphaSTKFieldPtr =
            this->alphaRef(phaseIndex).stkFieldPtr();

        for (stk::mesh::BucketVector::const_iterator ib = nodeBuckets.begin();
             ib != nodeBuckets.end();
             ++ib)
        {
            stk::mesh::Bucket& nodeBucket = **ib;

            const stk::mesh::Bucket::size_type nNodesPerBucket =
                nodeBucket.size();

            // field chunks in bucket
            const scalar* rhob =
                stk::mesh::field_data(*rhoSTKFieldPtr, nodeBucket);
            const scalar* alphab =
                stk::mesh::field_data(*alphaSTKFieldPtr, nodeBucket);
            scalar* bulkRhob =
                stk::mesh::field_data(*bulkRhoSTKFieldPtr, nodeBucket);

            for (stk::mesh::Bucket::size_type iNode = 0;
                 iNode < nNodesPerBucket;
                 ++iNode)
            {
                bulkRhob[iNode] += rhob[iNode] * alphab[iNode];
            }
        }
    }
}

void freeSurfaceFlowModel::updateDensity(const std::shared_ptr<domain> domain,
                                         label iPhase)
{
    auto option =
        domain->materialRef(domain->globalToLocalMaterialIndex(iPhase))
            .thermodynamicProperties_.equationOfState_.option_;

    switch (option)
    {
        case equationOfStateOption::value:
            fieldBroker::updateDensity(domain, iPhase);
            break;

        case equationOfStateOption::idealGas:
            {
                std::unique_ptr<idealGasModel> model =
                    std::make_unique<idealGasModel>(this->realmPtr_);
                model->updateDensity(domain, iPhase);
            }
            break;
    }
}

void freeSurfaceFlowModel::updateDynamicViscosity(
    const std::shared_ptr<domain> domain)
{
    auto& mesh = this->meshRef();
    stk::mesh::BulkData& bulkData = mesh.bulkDataRef();
    stk::mesh::MetaData& metaData = mesh.metaDataRef();

    // Get pointer to global density field to be used
    const auto* bulkMuSTKFieldPtr = muRef().stkFieldPtr();

    // get interior parts the domain is defined on
    const stk::mesh::PartVector& partVec = domain->zonePtr()->interiorParts();

    // Initialize global density to zero before every update
    muRef().setToValue({0}, partVec);

    // define some common selectors; select owned nodes
    stk::mesh::Selector selUniversalNodes =
        metaData.universal_part() & stk::mesh::selectUnion(partVec);

    stk::mesh::BucketVector const& nodeBuckets =
        bulkData.get_buckets(stk::topology::NODE_RANK, selUniversalNodes);

    for (label iPhase = 0; iPhase < domain->nMaterials(); iPhase++)
    {
        label phaseIndex = domain->localToGlobalMaterialIndex(iPhase);

        // Get fields for a given iPhase
        const STKScalarField* muSTKFieldPtr =
            this->muRef(phaseIndex).stkFieldPtr();
        const STKScalarField* alphaSTKFieldPtr =
            this->alphaRef(phaseIndex).stkFieldPtr();

        for (stk::mesh::BucketVector::const_iterator ib = nodeBuckets.begin();
             ib != nodeBuckets.end();
             ++ib)
        {
            stk::mesh::Bucket& nodeBucket = **ib;

            const stk::mesh::Bucket::size_type nNodesPerBucket =
                nodeBucket.size();

            // field chunks in bucket
            const scalar* mub =
                stk::mesh::field_data(*muSTKFieldPtr, nodeBucket);
            const scalar* alphab =
                stk::mesh::field_data(*alphaSTKFieldPtr, nodeBucket);
            scalar* bulkMub =
                stk::mesh::field_data(*bulkMuSTKFieldPtr, nodeBucket);

            for (stk::mesh::Bucket::size_type iNode = 0;
                 iNode < nNodesPerBucket;
                 ++iNode)
            {
                bulkMub[iNode] += mub[iNode] * alphab[iNode];
            }
        }
    }
}

void freeSurfaceFlowModel::updateDynamicViscosity(
    const std::shared_ptr<domain> domain,
    label iPhase)
{
    auto option =
        domain->materialRef(domain->globalToLocalMaterialIndex(iPhase))
            .transportProperties_.dynamicViscosity_.option_;

    switch (option)
    {
        case dynamicViscosityOption::value:
            fieldBroker::updateDynamicViscosity(domain, iPhase);
            break;

        case dynamicViscosityOption::sutherlandsFormula:
            {
                std::unique_ptr<sutherlandsFormulaModel> model =
                    std::make_unique<sutherlandsFormulaModel>(this->realmPtr_);
                model->updateDynamicViscosity(domain, iPhase);
            }
            break;
    }
}

void freeSurfaceFlowModel::initializeMassFlowRateInterior_(
    const std::shared_ptr<domain> domain,
    label iPhase)
{
    flowModel::initializeMassFlowRateInterior_(
        domain, this->mDotRef(iPhase), this->rhoRef(iPhase));
}

#ifdef HAS_INTERFACE
void freeSurfaceFlowModel::initializeMassFlowRateInterfaceSideField_(
    const std::shared_ptr<domain> domain,
    const interfaceSideInfo* interfaceSideInfoPtr,
    label iPhase)
{
    flowModel::initializeMassFlowRateInterfaceSideField_(
        domain,
        interfaceSideInfoPtr,
        this->mDotRef(iPhase).sideFieldRef(),
        this->rhoRef(iPhase));
}
#endif /* HAS_INTERFACE */

void freeSurfaceFlowModel::initializeMassFlowRateBoundaryField_(
    const std::shared_ptr<domain> domain,
    const boundary* boundary,
    label iPhase)
{
    flowModel::initializeMassFlowRateBoundaryField_(
        domain,
        boundary,
        this->mDotRef(iPhase).sideFieldRef(),
        this->rhoRef(iPhase));
}

void freeSurfaceFlowModel::updateMassFlowRateInterior_(
    const std::shared_ptr<domain> domain)
{
    auto& mesh = this->meshRef();
    stk::mesh::BulkData& bulkData = mesh.bulkDataRef();
    stk::mesh::MetaData& metaData = mesh.metaDataRef();

    // Get pointer to bulk density field to be used
    const auto* bulkMdotfSTKFieldPtr = this->mDotRef().stkFieldPtr();

    // nodal fields to gather
    std::vector<scalar> ws_alpha;

    // geometry related to populate
    std::vector<scalar> ws_shape_function;

    // get interior parts the domain is defined on
    const stk::mesh::PartVector& partVec = domain->zonePtr()->interiorParts();

    // Initialize global density to zero before every update
    this->mDotRef().setToValue({0}, partVec);

    // define some common selectors
    const stk::mesh::Selector selAllElements =
        metaData.universal_part() & stk::mesh::selectUnion(partVec);

    stk::mesh::BucketVector const& elementBuckets =
        bulkData.get_buckets(stk::topology::ELEMENT_RANK, selAllElements);

    for (label iPhase = 0; iPhase < domain->nMaterials(); iPhase++)
    {
        label phaseIndex = domain->localToGlobalMaterialIndex(iPhase);

        // Get fields for a given iPhase
        const STKScalarField* mDotSTKFieldPtr =
            this->mDotRef(phaseIndex).stkFieldPtr();
        const STKScalarField* alphaSTKFieldPtr =
            this->alphaRef(phaseIndex).stkFieldPtr();

        // shifted ip's for field?
        bool isAlphaShifted = this->alphaRef(phaseIndex).isShifted();

        for (stk::mesh::BucketVector::const_iterator ib =
                 elementBuckets.begin();
             ib != elementBuckets.end();
             ++ib)
        {
            stk::mesh::Bucket& elementBucket = **ib;
            const stk::mesh::Bucket::size_type nElementsPerBucket =
                elementBucket.size();

            // extract master element
            MasterElement* meSCS =
                MasterElementRepo::get_surface_master_element(
                    elementBucket.topology());

            // extract master element specifics
            const label nodesPerElement = meSCS->nodesPerElement_;
            const label numScsIp = meSCS->numIntPoints_;

            // algorithm related
            ws_alpha.resize(nodesPerElement);
            ws_shape_function.resize(numScsIp * nodesPerElement);

            // pointers
            scalar* p_alpha = &ws_alpha[0];
            scalar* p_shape_function = &ws_shape_function[0];

            // extract shape function
            if (isAlphaShifted)
            {
                meSCS->shifted_shape_fcn(&p_shape_function[0]);
            }
            else
            {
                meSCS->shape_fcn(&p_shape_function[0]);
            }

            for (stk::mesh::Bucket::size_type iElement = 0;
                 iElement < nElementsPerBucket;
                 ++iElement)
            {
                // ip data for this element; scs and scv
                scalar* mDot = stk::mesh::field_data(
                    *bulkMdotfSTKFieldPtr, elementBucket, iElement);
                const scalar* pmDot = stk::mesh::field_data(
                    *mDotSTKFieldPtr, elementBucket, iElement);

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

                    // pointers to real data
                    const scalar alpha =
                        *stk::mesh::field_data(*alphaSTKFieldPtr, node);

                    // gather scalars
                    p_alpha[ni] = alpha;
                }

                for (label ip = 0; ip < numScsIp; ++ip)
                {
                    const label offSetSF = ip * nodesPerElement;

                    scalar alphaIp = 0.0;
                    for (label ic = 0; ic < nodesPerElement; ++ic)
                    {
                        const scalar r = p_shape_function[offSetSF + ic];
                        alphaIp += r * p_alpha[ic];
                    }

                    mDot[ip] += alphaIp * pmDot[ip];
                }
            }
        }
    }
}

#ifdef HAS_INTERFACE
void freeSurfaceFlowModel::updateMassFlowRateInterfaceSideField_(
    const std::shared_ptr<domain> domain,
    const interfaceSideInfo* interfaceSideInfoPtr)
{
    const auto& mesh = this->meshRef();
    const stk::mesh::MetaData& metaData = mesh.metaDataRef();
    const stk::mesh::BulkData& bulkData = mesh.bulkDataRef();

    STKScalarField& mixtureMdotfSideSTKFieldPtr =
        this->mDotRef().sideFieldRef().stkFieldRef();

    // check
    assert(this->mDotRef().sideFieldRef().definedOn(
        interfaceSideInfoPtr->currentPartVec_));

    // zero-out
    this->mDotRef().sideFieldRef().setToValue(
        {0.0}, interfaceSideInfoPtr->currentPartVec_);

    // ip values; both boundary and opposing surface
    std::vector<scalar> currentIsoParCoords(SPATIAL_DIM);
    std::vector<scalar> opposingIsoParCoords(SPATIAL_DIM);

    // interpolate nodal values to point-in-elem
    const label sizeOfScalarField = 1;

    // nodal fields to gather; face
    std::vector<scalar> ws_c_alpha;
    std::vector<scalar> ws_o_alpha;

    // extract vector of dgInfo
    const std::vector<std::vector<dgInfo*>>& dgInfoVec =
        interfaceSideInfoPtr->dgInfoVec_;

    for (label iPhase = 0; iPhase < domain->nMaterials(); iPhase++)
    {
        label phaseIndex = domain->localToGlobalMaterialIndex(iPhase);

        // Get transport fields/side fields
        const auto& alphaSTKFieldRef = this->alphaRef(phaseIndex).stkFieldRef();
        const auto& mDotSideSTKFieldPtr =
            this->mDotRef(phaseIndex).sideFieldRef().stkFieldRef();

        for (label iSide = 0; iSide < static_cast<label>(dgInfoVec.size());
             iSide++)
        {
            const std::vector<dgInfo*>& faceDgInfoVec = dgInfoVec[iSide];

            // now loop over all the DgInfo objects on this
            // particular exposed face
            for (size_t k = 0; k < faceDgInfoVec.size(); ++k)
            {
                dgInfo* dgInfo = faceDgInfoVec[k];

                // extract current/opposing face/element
                stk::mesh::Entity currentFace = dgInfo->currentFace_;
                stk::mesh::Entity opposingFace = dgInfo->opposingFace_;
                stk::mesh::Entity currentElement = dgInfo->currentElement_;
                stk::mesh::Entity opposingElement = dgInfo->opposingElement_;

                // master element; face and volume
                MasterElement* meFCCurrent = dgInfo->meFCCurrent_;
                MasterElement* meFCOpposing = dgInfo->meFCOpposing_;

                // local ip, ordinals, etc
                const label currentGaussPointId = dgInfo->currentGaussPointId_;
                currentIsoParCoords = dgInfo->currentIsoParCoords_;
                opposingIsoParCoords = dgInfo->opposingIsoParCoords_;

                // pointer to mDot
                scalar* ncmDot = stk::mesh::field_data(
                    mixtureMdotfSideSTKFieldPtr, currentFace);
                const scalar* pncmDot =
                    stk::mesh::field_data(mDotSideSTKFieldPtr, currentFace);

                // extract some master element info
                const label currentNodesPerFace = meFCCurrent->nodesPerElement_;
                const label opposingNodesPerFace =
                    meFCOpposing->nodesPerElement_;

                // algorithm related; face
                ws_c_alpha.resize(currentNodesPerFace);
                ws_o_alpha.resize(opposingNodesPerFace);

                // face
                scalar* p_c_alpha = &ws_c_alpha[0];
                scalar* p_o_alpha = &ws_o_alpha[0];

                // gather current face data
                stk::mesh::Entity const* current_face_node_rels =
                    bulkData.begin_nodes(currentFace);
                const label current_num_face_nodes =
                    bulkData.num_nodes(currentFace);
                for (label ni = 0; ni < current_num_face_nodes; ++ni)
                {
                    stk::mesh::Entity node = current_face_node_rels[ni];

                    // gather; scalar
                    p_c_alpha[ni] =
                        *stk::mesh::field_data(alphaSTKFieldRef, node);
                }

                // gather opposing face data
                stk::mesh::Entity const* opposing_face_node_rels =
                    bulkData.begin_nodes(opposingFace);
                const label opposing_num_face_nodes =
                    bulkData.num_nodes(opposingFace);
                for (label ni = 0; ni < opposing_num_face_nodes; ++ni)
                {
                    stk::mesh::Entity node = opposing_face_node_rels[ni];

                    // gather; scalar
                    p_o_alpha[ni] =
                        *stk::mesh::field_data(alphaSTKFieldRef, node);
                }

                scalar currentAlphaBip = 0.0;
                meFCCurrent->interpolatePoint(sizeOfScalarField,
                                              &currentIsoParCoords[0],
                                              &ws_c_alpha[0],
                                              &currentAlphaBip);

                scalar opposingAlphaBip = 0.0;
                meFCOpposing->interpolatePoint(sizeOfScalarField,
                                               &opposingIsoParCoords[0],
                                               &ws_o_alpha[0],
                                               &opposingAlphaBip);

                ncmDot[currentGaussPointId] +=
                    (currentAlphaBip + opposingAlphaBip) / 2.0 *
                    pncmDot[currentGaussPointId];
            }
        }
    }
}
#endif /* HAS_INTERFACE */

void freeSurfaceFlowModel::updateMassFlowRateBoundaryField_(
    const std::shared_ptr<domain> domain,
    const boundary* boundary)
{
    boundaryPhysicalType type = boundary->type();

    switch (type)
    {
        case boundaryPhysicalType::symmetry:
        case boundaryPhysicalType::wall:
            break;

        default:
            {
                using stk::mesh::Bucket;
                using stk::mesh::BucketVector;

                const auto& mesh = this->meshRef();
                const stk::mesh::MetaData& metaData = mesh.metaDataRef();
                const stk::mesh::BulkData& bulkData = mesh.bulkDataRef();

                // nodal fields to gather
                std::vector<scalar> ws_alpha;

                // master element
                std::vector<scalar> ws_face_shape_function;

                // Get fields
                const auto& mDotSideSTKFieldRef =
                    this->mDotRef().sideFieldRef().stkFieldRef();

                this->mDotRef().sideFieldRef().setToValue({0.0},
                                                          boundary->parts());

                // define vector of parent topos; should
                // always be UNITY in size
                std::vector<stk::topology> parentTopo;

                // define some common selectors
                stk::mesh::Selector selAllSides =
                    metaData.universal_part() &
                    stk::mesh::selectUnion(boundary->parts());

                for (label iPhase = 0; iPhase < domain->nMaterials(); iPhase++)
                {
                    label phaseIndex =
                        domain->localToGlobalMaterialIndex(iPhase);

                    const auto& alphaSTKFieldRef =
                        this->alphaRef(phaseIndex).stkFieldRef();
                    const auto& pmDotSideSTKFieldRef =
                        this->mDotRef(phaseIndex).sideFieldRef().stkFieldRef();

                    // shifted ip's for field?
                    bool isAlphaShifted =
                        this->alphaRef(phaseIndex).isShifted();

                    BucketVector const& sideBuckets =
                        bulkData.get_buckets(metaData.side_rank(), selAllSides);
                    for (BucketVector::const_iterator ib = sideBuckets.begin();
                         ib != sideBuckets.end();
                         ++ib)
                    {
                        Bucket& sideBucket = **ib;

                        // extract connected element
                        // topology
                        sideBucket.parent_topology(stk::topology::ELEMENT_RANK,
                                                   parentTopo);
                        STK_ThrowAssert(parentTopo.size() == 1);
                        stk::topology theElemTopo = parentTopo[0];

                        // face master element
                        MasterElement* meFC =
                            MasterElementRepo::get_surface_master_element(
                                sideBucket.topology());
                        const label nodesPerSide =
                            sideBucket.topology().num_nodes();
                        const label numScsBip = meFC->numIntPoints_;

                        // algorithm related; element
                        // (exposed face and element)
                        ws_face_shape_function.resize(numScsBip * nodesPerSide);
                        ws_alpha.resize(nodesPerSide);

                        // pointers
                        scalar* p_alpha = &ws_alpha[0];
                        scalar* p_face_shape_function =
                            &ws_face_shape_function[0];

                        // shape functions; boundary
                        if (isAlphaShifted)
                        {
                            meFC->shifted_shape_fcn(&p_face_shape_function[0]);
                        }
                        else
                        {
                            meFC->shape_fcn(&p_face_shape_function[0]);
                        }

                        const Bucket::size_type nSidesPerBucket =
                            sideBucket.size();

                        for (Bucket::size_type iSide = 0;
                             iSide < nSidesPerBucket;
                             ++iSide)
                        {
                            // get face
                            stk::mesh::Entity side = sideBucket[iSide];

                            //======================================
                            // gather nodal data off of face
                            //======================================
                            stk::mesh::Entity const* sideNodeRels =
                                bulkData.begin_nodes(side);
                            label numSideNodes = bulkData.num_nodes(side);

                            // sanity check on num nodes
                            STK_ThrowAssert(numSideNodes == nodesPerSide);
                            for (label ni = 0; ni < numSideNodes; ++ni)
                            {
                                stk::mesh::Entity node = sideNodeRels[ni];

                                // gather scalars
                                p_alpha[ni] = *stk::mesh::field_data(
                                    alphaSTKFieldRef, node);
                            }

                            // pointer to face data
                            scalar* mDot = stk::mesh::field_data(
                                mDotSideSTKFieldRef, side);
                            const scalar* pmDot = stk::mesh::field_data(
                                pmDotSideSTKFieldRef, side);

                            // loop over boundary ips
                            for (label ip = 0; ip < numScsBip; ++ip)
                            {
                                // interpolate to bip
                                scalar alphaBip = 0;
                                const label offSetSF_face = ip * nodesPerSide;
                                for (label ic = 0; ic < nodesPerSide; ++ic)
                                {
                                    const scalar r =
                                        p_face_shape_function[offSetSF_face +
                                                              ic];
                                    const scalar alpha = p_alpha[ic];

                                    alphaBip += r * alpha;
                                }

                                for (label j = 0; j < SPATIAL_DIM; ++j)
                                {
                                    mDot[ip] += alphaBip * pmDot[ip];
                                }
                            }
                        }
                    }
                }
            }
            break;
    }
}

void freeSurfaceFlowModel::updateMassFlowRateInterior_(
    const std::shared_ptr<domain> domain,
    label iPhase)
{
    flowModel::updateMassFlowRateInterior_(
        domain, this->mDotRef(iPhase), this->rhoRef(iPhase));
}

#ifdef HAS_INTERFACE
void freeSurfaceFlowModel::updateMassFlowRateInterfaceSideField_(
    const std::shared_ptr<domain> domain,
    const interfaceSideInfo* interfaceSideInfoPtr,
    label iPhase)
{
    flowModel::updateMassFlowRateInterfaceSideField_(
        domain,
        interfaceSideInfoPtr,
        this->mDotRef(iPhase).sideFieldRef(),
        this->rhoRef(iPhase));
}
#endif /* HAS_INTERFACE */

void freeSurfaceFlowModel::updateMassFlowRateBoundaryField_(
    const std::shared_ptr<domain> domain,
    const boundary* boundary,
    label iPhase)
{
    boundaryPhysicalType type = boundary->type();
    boundaryConditionType pbcType =
        this->pRef()
            .boundaryConditionRef(domain->index(), boundary->index())
            .type();

    switch (type)
    {
        case boundaryPhysicalType::symmetry:
        case boundaryPhysicalType::wall:
            break;

        case boundaryPhysicalType::inlet:
            {
                switch (pbcType)
                {
                    case boundaryConditionType::zeroGradient:
                        {
                            flowModel::
                                updateMassFlowRateBoundaryFieldInletSpecifiedVelocity_(
                                    domain,
                                    boundary,
                                    this->mDotRef(iPhase).sideFieldRef(),
                                    this->rhoRef(iPhase));
                        }
                        break;

                    case boundaryConditionType::staticPressure:
                    case boundaryConditionType::totalPressure:
                        {
                            flowModel::
                                updateMassFlowRateBoundaryFieldInletSpecifiedPressure_(
                                    domain,
                                    boundary,
                                    this->mDotRef(iPhase).sideFieldRef(),
                                    this->rhoRef(iPhase));
                        }
                        break;

                    default:
                        errorMsg("boundary condition invalid");
                }
            }
            break;

        case boundaryPhysicalType::opening:
            {
                switch (pbcType)
                {
                    case boundaryConditionType::staticPressure:
                    case boundaryConditionType::totalPressure:
                        {
                            flowModel::
                                updateMassFlowRateBoundaryFieldOpeningPressure_(
                                    domain,
                                    boundary,
                                    this->mDotRef(iPhase).sideFieldRef(),
                                    this->rhoRef(iPhase));
                        }
                        break;

                    default:
                        errorMsg("boundary condition invalid");
                }
            }
            break;

        case boundaryPhysicalType::outlet:
            {
                switch (pbcType)
                {
                    case boundaryConditionType::staticPressure:
                        {
                            flowModel::
                                updateMassFlowRateBoundaryFieldOutletSpecifiedPressure_(
                                    domain,
                                    boundary,
                                    this->mDotRef(iPhase).sideFieldRef(),
                                    this->rhoRef(iPhase));
                        }
                        break;

                    default:
                        errorMsg("boundary condition invalid");
                }
            }
            break;

        default:
            break;
    }
}

void freeSurfaceFlowModel::applyVolumeConservation(
    const std::shared_ptr<domain> domain)
{
    auto& mesh = this->meshRef();
    stk::mesh::BulkData& bulkData = mesh.bulkDataRef();
    stk::mesh::MetaData& metaData = mesh.metaDataRef();

    // get interior parts the domain is defined on
    const stk::mesh::PartVector& partVec = domain->zonePtr()->interiorParts();

    // Initialize primary volume fraction to 1.0
    this->alphaRef(domain->nMaterials() - 1).setToValue({1.0}, partVec);

    // define some common selectors; select owned nodes
    stk::mesh::Selector selUniversalNodes =
        metaData.universal_part() & stk::mesh::selectUnion(partVec);

    stk::mesh::BucketVector const& nodeBuckets =
        bulkData.get_buckets(stk::topology::NODE_RANK, selUniversalNodes);

    // get primary volume fraction field
    const STKScalarField* palphaSTKFieldPtr =
        this->alphaRef(domain->nMaterials() - 1).stkFieldPtr();

    for (label iPhase = 0; iPhase < domain->nMaterials() - 1; iPhase++)
    {
        label phaseIndex = domain->localToGlobalMaterialIndex(iPhase);

        // Get fields for a given iPhase
        const STKScalarField* alphaSTKFieldPtr =
            this->alphaRef(phaseIndex).stkFieldPtr();

        for (stk::mesh::BucketVector::const_iterator ib = nodeBuckets.begin();
             ib != nodeBuckets.end();
             ++ib)
        {
            stk::mesh::Bucket& nodeBucket = **ib;

            const stk::mesh::Bucket::size_type nNodesPerBucket =
                nodeBucket.size();

            // field chunks in bucket
            const scalar* alphab =
                stk::mesh::field_data(*alphaSTKFieldPtr, nodeBucket);
            scalar* palphab =
                stk::mesh::field_data(*palphaSTKFieldPtr, nodeBucket);

            for (stk::mesh::Bucket::size_type iNode = 0;
                 iNode < nNodesPerBucket;
                 ++iNode)
            {
                palphab[iNode] -= alphab[iNode];
            }
        }
    }
}

void freeSurfaceFlowModel::updateInterfaceNormal(
    const std::shared_ptr<domain> domain,
    label iPhase)
{
    auto& mesh = this->meshRef();
    stk::mesh::BulkData& bulkData = mesh.bulkDataRef();
    stk::mesh::MetaData& metaData = mesh.metaDataRef();

    // Integration point data that is fixed
    std::vector<scalar> d_alpha_dxIp(SPATIAL_DIM);
    scalar* p_d_alpha_dxIp = &d_alpha_dxIp[0];

    // Define scratch spaces
    std::vector<scalar> ws_coordinates;
    std::vector<scalar> ws_alpha;
    std::vector<scalar> ws_dualVolume;
    std::vector<scalar> ws_scVolume;
    std::vector<scalar> ws_dndx;
    std::vector<scalar> ws_deriv;
    std::vector<scalar> ws_det_j;

    // Retrieve smoothing boolean from user input
    bool smooth = controlsRef()
                      .solverRef()
                      .solverControl_.advancedOptions_.equationControls_
                      .volumeFractionSmoothing_.smoothVolumeFraction_;

    // get interior parts the domain is defined on
    const stk::mesh::PartVector& partVec = domain->zonePtr()->interiorParts();

    // define some common selectors; select owned nodes
    stk::mesh::Selector selUniversalElements =
        metaData.universal_part() & stk::mesh::selectUnion(partVec);

    stk::mesh::BucketVector const& elementBuckets =
        bulkData.get_buckets(stk::topology::ELEMENT_RANK, selUniversalElements);

    // get geometry
    const auto* coordsSTKFieldPtr = metaData.get_field<scalar>(
        stk::topology::NODE_RANK, this->getCoordinatesID_(domain));

    const auto* volSTKFieldPtr = metaData.get_field<scalar>(
        stk::topology::NODE_RANK, this->getDualNodalVolumeID_(domain));

    // Get fields for a given iPhase
    const STKScalarField* alphaSTKFieldPtr =
        smooth ? this->alphaSmoothRef(iPhase).stkFieldPtr()
               : this->alphaRef(iPhase).stkFieldPtr();

    STKScalarField* nHatSTKFieldPtr = this->nHatRef(iPhase).stkFieldPtr();

    // zero-out interface normal in the current domain
    ops::zero<scalar>(nHatSTKFieldPtr, partVec);

    for (stk::mesh::BucketVector::const_iterator ib = elementBuckets.begin();
         ib < elementBuckets.end();
         ib++)
    {
        stk::mesh::Bucket& elementBucket = **ib;
        const stk::mesh::Bucket::size_type nElemPerBucket =
            elementBucket.size();

        // extract master element, and its needed specifics
        MasterElement* meSCV = MasterElementRepo::get_volume_master_element(
            elementBucket.topology());
        const label nodesPerElement = meSCV->nodesPerElement_;
        const label numScvIp = meSCV->numIntPoints_;
        const label* ipNodeMap = meSCV->ipNodeMap();

        // allocate space for scratch spaces
        ws_coordinates.resize(nodesPerElement * SPATIAL_DIM);
        ws_alpha.resize(nodesPerElement);
        ws_dualVolume.resize(nodesPerElement);
        ws_scVolume.resize(numScvIp);
        ws_dndx.resize(SPATIAL_DIM * numScvIp * nodesPerElement);
        ws_deriv.resize(SPATIAL_DIM * numScvIp * nodesPerElement);
        ws_det_j.resize(numScvIp);

        // define needed pointers
        scalar* p_coordinates = &ws_coordinates[0];
        scalar* p_alpha = &ws_alpha[0];
        scalar* p_dualVolume = &ws_dualVolume[0];
        scalar* p_scVolume = &ws_scVolume[0];
        scalar* p_dndx = &ws_dndx[0];

        for (stk::mesh::Bucket::size_type iElement = 0;
             iElement < nElemPerBucket;
             iElement++)
        {
            stk::mesh::Entity const* nodeRels =
                elementBucket.begin_nodes(iElement);

            label numNodes = elementBucket.num_nodes(iElement);

            // sanity check on nb. of nodes per element before
            // proceeding
            STK_ThrowAssert(numNodes == nodesPerElement);

            for (label iNode = 0; iNode < nodesPerElement; iNode++)
            {
                stk::mesh::Entity node = nodeRels[iNode];

                p_dualVolume[iNode] =
                    *stk::mesh::field_data(*volSTKFieldPtr, node);
                p_alpha[iNode] =
                    *stk::mesh::field_data(*alphaSTKFieldPtr, node);
                const scalar* coordsb =
                    stk::mesh::field_data(*coordsSTKFieldPtr, node);

                const label offset = iNode * SPATIAL_DIM;
                for (label j = 0; j < SPATIAL_DIM; j++)
                {
                    p_coordinates[offset + j] = coordsb[j];
                }
            }

            // compute geometry and dndx
            scalar scv_error = 0.0;
            meSCV->determinant(
                1, &p_coordinates[0], &p_scVolume[0], &scv_error);
            meSCV->grad_op(1,
                           &p_coordinates[0],
                           &p_dndx[0],
                           &ws_deriv[0],
                           &ws_det_j[0],
                           &scv_error);

            for (label ip = 0; ip < numScvIp; ip++)
            {
                // zero local ip gradient
                for (label j = 0; j < SPATIAL_DIM; j++)
                {
                    p_d_alpha_dxIp[j] = 0.0;
                }

                // compute gradient
                for (label iNode = 0; iNode < nodesPerElement; iNode++)
                {
                    const label offsetDnDx =
                        SPATIAL_DIM * nodesPerElement * ip +
                        iNode * SPATIAL_DIM;
                    const scalar nodalVof = p_alpha[iNode];

                    for (label j = 0; j < SPATIAL_DIM; j++)
                    {
                        p_d_alpha_dxIp[j] += p_dndx[offsetDnDx + j] * nodalVof;
                    }
                }

                // compute magnitude
                scalar dalphaMag = SMALL;
                for (label j = 0; j < SPATIAL_DIM; j++)
                {
                    dalphaMag += p_d_alpha_dxIp[j] * p_d_alpha_dxIp[j];
                }
                dalphaMag = std::sqrt(dalphaMag);

                // nearest node for this ip
                const label nn = ipNodeMap[ip];
                stk::mesh::Entity node = nodeRels[nn];

                scalar* interNormal =
                    stk::mesh::field_data(*nHatSTKFieldPtr, node);

                const scalar volumeFrac = p_scVolume[ip] / p_dualVolume[nn];
                for (label j = 0; j < SPATIAL_DIM; j++)
                {
                    interNormal[j] +=
                        p_d_alpha_dxIp[j] / dalphaMag * volumeFrac;
                }
            }
        }
    }

    if (messager::parallel())
    {
        this->nHatRef(iPhase).synchronizeGhostedEntities(domain->index());
    }
}

void freeSurfaceFlowModel::computeSmoothRHS_(
    const std::shared_ptr<domain> domain,
    label iPhase)
{
    auto& mesh = this->meshRef();
    stk::mesh::BulkData& bulkData = mesh.bulkDataRef();
    stk::mesh::MetaData& metaData = mesh.metaDataRef();

    std::vector<scalar> ws_alpha;
    std::vector<scalar> ws_dualVolume;
    std::vector<scalar> ws_coordinates;
    std::vector<scalar> ws_scsAreav;
    std::vector<scalar> ws_dndx;
    std::vector<scalar> ws_deriv;
    std::vector<scalar> ws_detj;

    scalar l_dxMin = 1.0e16;
    scalar g_dxMin = 1.0e16;

    ops::zero(this->rhsSmoothRef(iPhase).stkFieldPtr(),
              domain->zonePtr()->interiorParts());

    // get interior parts the domain is defined on
    const stk::mesh::PartVector& partVec = domain->zonePtr()->interiorParts();

    // define some common selectors; select owned nodes
    stk::mesh::Selector selAllElements =
        metaData.universal_part() & stk::mesh::selectUnion(partVec);

    stk::mesh::BucketVector const& elementBuckets =
        bulkData.get_buckets(stk::topology::ELEMENT_RANK, selAllElements);

    // get geometry
    const auto* coordsSTKFieldPtr = metaData.get_field<scalar>(
        stk::topology::NODE_RANK, this->getCoordinatesID_(domain));

    const auto* volSTKFieldPtr = metaData.get_field<scalar>(
        stk::topology::NODE_RANK, this->getDualNodalVolumeID_(domain));

    // get fields for a given iPhase
    const STKScalarField* alphaSmoothedSTKFieldPtr =
        this->alphaSmoothRef(iPhase).stkFieldPtr();

    STKScalarField* rhsSmoothSTKFieldPtr =
        this->rhsSmoothRef(iPhase).stkFieldPtr();

    for (stk::mesh::BucketVector::const_iterator ib = elementBuckets.begin();
         ib < elementBuckets.end();
         ib++)
    {
        stk::mesh::Bucket& elementBucket = **ib;
        const stk::mesh::Bucket::size_type nElemPerBucket =
            elementBucket.size();

        // extract master element, and its needed specifics
        MasterElement* meSCS = MasterElementRepo::get_surface_master_element(
            elementBucket.topology());
        const label nodesPerElement = meSCS->nodesPerElement_;
        const label numScsIp = meSCS->numIntPoints_;
        const label* lrscv = meSCS->adjacentNodes();

        // allocate space for scratch spaces
        ws_coordinates.resize(nodesPerElement * SPATIAL_DIM);
        ws_alpha.resize(nodesPerElement);
        ws_dualVolume.resize(nodesPerElement);
        ws_scsAreav.resize(numScsIp * SPATIAL_DIM);
        ws_dndx.resize(SPATIAL_DIM * numScsIp * nodesPerElement);
        ws_deriv.resize(SPATIAL_DIM * numScsIp * nodesPerElement);
        ws_detj.resize(numScsIp);

        // define needed pointers
        scalar* p_coordinates = &ws_coordinates[0];
        scalar* p_alpha = &ws_alpha[0];
        scalar* p_dualVolume = &ws_dualVolume[0];
        scalar* p_scsAreav = &ws_scsAreav[0];
        scalar* p_dndx = &ws_dndx[0];

        for (stk::mesh::Bucket::size_type iElement = 0;
             iElement < nElemPerBucket;
             iElement++)
        {
            stk::mesh::Entity const* nodeRels =
                elementBucket.begin_nodes(iElement);

            label numNodes = elementBucket.num_nodes(iElement);

            // sanity check on nb. of nodes per element
            STK_ThrowAssert(numNodes == nodesPerElement);

            for (label iNode = 0; iNode < nodesPerElement; iNode++)
            {
                stk::mesh::Entity node = nodeRels[iNode];

                p_dualVolume[iNode] =
                    *stk::mesh::field_data(*volSTKFieldPtr, node);
                p_alpha[iNode] =
                    *stk::mesh::field_data(*alphaSmoothedSTKFieldPtr, node);
                const scalar* coordsb =
                    stk::mesh::field_data(*coordsSTKFieldPtr, node);

                const label offset = iNode * SPATIAL_DIM;
                for (label j = 0; j < SPATIAL_DIM; j++)
                {
                    p_coordinates[offset + j] = coordsb[j];
                }
            }
            // compute geometry
            scalar scsError = 0.0;
            meSCS->determinant(1, &p_coordinates[0], &p_scsAreav[0], &scsError);
            meSCS->grad_op(1,
                           &p_coordinates[0],
                           &p_dndx[0],
                           &ws_deriv[0],
                           &ws_detj[0],
                           &scsError);

            for (label ip = 0; ip < numScsIp; ++ip)
            {
                // left and right nodes for this ip
                const label il = lrscv[2 * ip];
                const label ir = lrscv[2 * ip + 1];

                // determine dx; edge distance magnitude
                scalar dx = 0.0;
                for (label j = 0; j < SPATIAL_DIM; ++j)
                {
                    const scalar dxj = p_coordinates[ir * SPATIAL_DIM + j] -
                                       p_coordinates[il * SPATIAL_DIM + j];
                    dx += dxj * dxj;
                }
                l_dxMin = std::min(l_dxMin, std::sqrt(dx));

                stk::mesh::Entity nodeL = nodeRels[il];
                stk::mesh::Entity nodeR = nodeRels[ir];

                // pointer to fields to assemble
                scalar* smoothedRhsL =
                    stk::mesh::field_data(*rhsSmoothSTKFieldPtr, nodeL);
                scalar* smoothedRhsR =
                    stk::mesh::field_data(*rhsSmoothSTKFieldPtr, nodeR);

                scalar qDiff = 0.0;
                for (label ic = 0; ic < nodesPerElement; ++ic)
                {
                    scalar lhsfacDiff = 0.0;
                    const label offSetDnDx =
                        SPATIAL_DIM * nodesPerElement * ip + ic * SPATIAL_DIM;
                    for (label j = 0; j < SPATIAL_DIM; ++j)
                    {
                        lhsfacDiff += -p_dndx[offSetDnDx + j] *
                                      p_scsAreav[ip * SPATIAL_DIM + j];
                    }
                    qDiff += lhsfacDiff * p_alpha[ic];
                }

                *smoothedRhsL -= qDiff / ws_dualVolume[il];
                *smoothedRhsR += qDiff / ws_dualVolume[ir];
            }
        }
    }

    // synchronize ghosted entities
    this->rhsSmoothRef(iPhase).synchronizeGhostedEntities(domain->index());

    // find global min dxmin
    stk::all_reduce_min(MPI_COMM_WORLD, &l_dxMin, &g_dxMin, 1);

    // store the min dxmin across domains
    dxMin_ = std::min(dxMin_, g_dxMin);
}

void freeSurfaceFlowModel::assembleSmoothingTerm_(
    const std::shared_ptr<domain> domain,
    label iPhase)
{
    auto& mesh = this->meshRef();
    stk::mesh::BulkData& bulkData = mesh.bulkDataRef();
    stk::mesh::MetaData& metaData = mesh.metaDataRef();

    // get interior parts the domain is defined on
    const stk::mesh::PartVector& partVec = domain->zonePtr()->interiorParts();

    // get fields for a given iPhase
    const STKScalarField* rhsSmoothSTKFieldPtr =
        this->rhsSmoothRef(iPhase).stkFieldPtr();

    STKScalarField* alphaSmoothSTKFieldPtr =
        this->alphaSmoothRef(iPhase).stkFieldPtr();

    // copy vof value to vof_smooth in the current domain
    stk::mesh::Selector selAllNodes =
        metaData.universal_part() & stk::mesh::selectUnion(partVec);
    stk::mesh::BucketVector const& nodeBuckets =
        bulkData.get_buckets(stk::topology::NODE_RANK, selAllNodes);
    for (stk::mesh::BucketVector::const_iterator ib = nodeBuckets.begin();
         ib != nodeBuckets.end();
         ++ib)
    {
        stk::mesh::Bucket& nodeBucket = **ib;
        const stk::mesh::Bucket::size_type nNodesPerBucket = nodeBucket.size();
        scalar* smoothedRHS =
            stk::mesh::field_data(*rhsSmoothSTKFieldPtr, nodeBucket);
        scalar* vof_smooth =
            stk::mesh::field_data(*alphaSmoothSTKFieldPtr, nodeBucket);

        for (stk::mesh::Bucket::size_type iNode = 0; iNode < nNodesPerBucket;
             ++iNode)
        {
            vof_smooth[iNode] =
                Fo_ * dxMin_ * dxMin_ * smoothedRHS[iNode] + vof_smooth[iNode];
        }
    }
}

void freeSurfaceFlowModel::applyVOFSmoothing(
    const std::shared_ptr<domain> domain,
    label iPhase)
{
    bool smooth = controlsRef()
                      .solverRef()
                      .solverControl_.advancedOptions_.equationControls_
                      .volumeFractionSmoothing_.smoothVolumeFraction_;
    if (!smooth)
        return;

    ops::copy<scalar>(this->alphaRef(iPhase).stkFieldPtr(),
                      this->alphaSmoothRef(iPhase).stkFieldPtr(),
                      domain->zonePtr()->interiorParts());

    label kmax = controlsRef()
                     .solverRef()
                     .solverControl_.advancedOptions_.equationControls_
                     .volumeFractionSmoothing_.smoothingIterations_;
    for (label k = 0; k < kmax; k++)
    {
        computeSmoothRHS_(domain, iPhase);
        assembleSmoothingTerm_(domain, iPhase);
    }
}

void freeSurfaceFlowModel::setupFCTFields(const std::shared_ptr<domain> domain,
                                          label iPhase)
{
    // Get Mesh Data
    auto& mesh = this->meshRef();
    stk::mesh::BulkData& bulkData = mesh.bulkDataRef();
    stk::mesh::MetaData& metaData = mesh.metaDataRef();

    // get interior parts the domain is defined on
    const stk::mesh::PartVector& partVec = domain->zonePtr()->interiorParts();

    // FL
    if (!FLSTKFieldPtr_)
    {
        scalar initialValue = 0.0;

        FLSTKFieldPtr_ = &metaData.declare_field<scalar>(
            stk::topology::ELEM_RANK, "FL" + std::to_string(iPhase));

        // Put the stk field on registered parts
        for (const stk::mesh::Part* part : partVec)
        {
            // Determine number of integration points in the element
            const stk::topology theTopo = part->topology();
            MasterElement* meSCS =
                MasterElementRepo::get_surface_master_element(theTopo);
            const int numScsIP = meSCS->numIntPoints_;

            // Put the field on mesh
            stk::mesh::put_field_on_mesh(
                *FLSTKFieldPtr_, *part, numScsIP, &initialValue);
        }

        sideFLSTKFieldPtr_ = &metaData.declare_field<scalar>(
            metaData.side_rank(), "side_FL" + std::to_string(iPhase));

        for (label iBoundary = 0; iBoundary < domain->zonePtr()->nBoundaries();
             iBoundary++)
        {
            const auto* boundary = domain->zonePtr()->boundaryPtr(iBoundary);

            boundaryPhysicalType type = boundary->type();

            switch (type)
            {
                case boundaryPhysicalType::inlet:
                case boundaryPhysicalType::outlet:
                case boundaryPhysicalType::opening:
                    {
                        for (const stk::mesh::Part* part : boundary->parts())
                        {
                            for (const stk::mesh::Part* subPart :
                                 part->subsets())
                            {
                                // Determine number of integration points in the
                                // face
                                MasterElement* meFC = MasterElementRepo::
                                    get_surface_master_element(
                                        subPart->topology());
                                const label numScsBip = meFC->numIntPoints_;

                                // Put the field on mesh
                                stk::mesh::put_field_on_mesh(
                                    *sideFLSTKFieldPtr_,
                                    *subPart,
                                    numScsBip,
                                    &initialValue);
                            }
                        }
                    }
                    break;

                default:
                    break;
            }
        }
    }

    // FH
    if (!FHSTKFieldPtr_)
    {
        scalar initialValue = 0.0;

        FHSTKFieldPtr_ = &metaData.declare_field<scalar>(
            stk::topology::ELEM_RANK, "FH" + std::to_string(iPhase));

        // Put the stk field on registered parts
        for (const stk::mesh::Part* part : partVec)
        {
            // Determine number of integration points in the element
            const stk::topology theTopo = part->topology();
            MasterElement* meSCS =
                MasterElementRepo::get_surface_master_element(theTopo);
            const int numScsIP = meSCS->numIntPoints_;

            // Put the field on mesh
            stk::mesh::put_field_on_mesh(
                *FHSTKFieldPtr_, *part, numScsIP, &initialValue);
        }

        sideFHSTKFieldPtr_ = &metaData.declare_field<scalar>(
            metaData.side_rank(), "side_FH" + std::to_string(iPhase));

        for (label iBoundary = 0; iBoundary < domain->zonePtr()->nBoundaries();
             iBoundary++)
        {
            const auto* boundary = domain->zonePtr()->boundaryPtr(iBoundary);

            boundaryPhysicalType type = boundary->type();

            switch (type)
            {
                case boundaryPhysicalType::inlet:
                case boundaryPhysicalType::outlet:
                case boundaryPhysicalType::opening:
                    {
                        for (const stk::mesh::Part* part : boundary->parts())
                        {
                            for (const stk::mesh::Part* subPart :
                                 part->subsets())
                            {
                                // Determine number of integration points in the
                                // face
                                MasterElement* meFC = MasterElementRepo::
                                    get_surface_master_element(
                                        subPart->topology());
                                const label numScsBip = meFC->numIntPoints_;

                                // Put the field on mesh
                                stk::mesh::put_field_on_mesh(
                                    *sideFHSTKFieldPtr_,
                                    *subPart,
                                    numScsBip,
                                    &initialValue);
                            }
                        }
                    }
                    break;

                default:
                    break;
            }
        }
    }

    // A
    if (!ASTKFieldPtr_)
    {
        scalar initialValue = 0.0;

        ASTKFieldPtr_ = &metaData.declare_field<scalar>(
            stk::topology::ELEM_RANK, "A" + std::to_string(iPhase));

        // Put the stk field on registered parts
        for (const stk::mesh::Part* part : partVec)
        {
            // Determine number of integration points in the element
            const stk::topology theTopo = part->topology();
            MasterElement* meSCS =
                MasterElementRepo::get_surface_master_element(theTopo);
            const int numScsIP = meSCS->numIntPoints_;

            // Put the field on mesh
            stk::mesh::put_field_on_mesh(
                *ASTKFieldPtr_, *part, numScsIP, &initialValue);
        }

        sideASTKFieldPtr_ = &metaData.declare_field<scalar>(
            metaData.side_rank(), "side_A" + std::to_string(iPhase));

        for (label iBoundary = 0; iBoundary < domain->zonePtr()->nBoundaries();
             iBoundary++)
        {
            const auto* boundary = domain->zonePtr()->boundaryPtr(iBoundary);

            boundaryPhysicalType type = boundary->type();

            switch (type)
            {
                case boundaryPhysicalType::inlet:
                case boundaryPhysicalType::outlet:
                case boundaryPhysicalType::opening:
                    {
                        for (const stk::mesh::Part* part : boundary->parts())
                        {
                            for (const stk::mesh::Part* subPart :
                                 part->subsets())
                            {
                                // Determine number of integration points in the
                                // face
                                MasterElement* meFC = MasterElementRepo::
                                    get_surface_master_element(
                                        subPart->topology());
                                const label numScsBip = meFC->numIntPoints_;

                                // Put the field on mesh
                                stk::mesh::put_field_on_mesh(*sideASTKFieldPtr_,
                                                             *subPart,
                                                             numScsBip,
                                                             &initialValue);
                            }
                        }
                    }
                    break;

                default:
                    break;
            }
        }
    }

    // lambda
    if (!lambdaSTKFieldPtr_)
    {
        scalar initialValue = 1.0;

        lambdaSTKFieldPtr_ = &metaData.declare_field<scalar>(
            stk::topology::ELEM_RANK, "lambda" + std::to_string(iPhase));

        // Put the stk field on registered parts
        for (const stk::mesh::Part* part : partVec)
        {
            // Determine number of integration points in the element
            const stk::topology theTopo = part->topology();
            MasterElement* meSCS =
                MasterElementRepo::get_surface_master_element(theTopo);
            const int numScsIP = meSCS->numIntPoints_;

            // Put the field on mesh
            stk::mesh::put_field_on_mesh(
                *lambdaSTKFieldPtr_, *part, numScsIP, &initialValue);
        }

        sideLambdaSTKFieldPtr_ = &metaData.declare_field<scalar>(
            metaData.side_rank(), "side_lambda" + std::to_string(iPhase));

        for (label iBoundary = 0; iBoundary < domain->zonePtr()->nBoundaries();
             iBoundary++)
        {
            const auto* boundary = domain->zonePtr()->boundaryPtr(iBoundary);

            boundaryPhysicalType type = boundary->type();

            switch (type)
            {
                case boundaryPhysicalType::inlet:
                case boundaryPhysicalType::outlet:
                case boundaryPhysicalType::opening:
                    {
                        for (const stk::mesh::Part* part : boundary->parts())
                        {
                            for (const stk::mesh::Part* subPart :
                                 part->subsets())
                            {
                                // Determine number of integration points in the
                                // face
                                MasterElement* meFC = MasterElementRepo::
                                    get_surface_master_element(
                                        subPart->topology());
                                const label numScsBip = meFC->numIntPoints_;

                                // Put the field on mesh
                                stk::mesh::put_field_on_mesh(
                                    *sideLambdaSTKFieldPtr_,
                                    *subPart,
                                    numScsBip,
                                    &initialValue);
                            }
                        }
                    }
                    break;

                default:
                    break;
            }
        }
    }

    // # Qplus
    if (!QplusSTKFieldPtr_)
    {
        scalar initialValue = 0.0;

        QplusSTKFieldPtr_ = &metaData.declare_field<scalar>(
            stk::topology::NODE_RANK, "Qplus" + std::to_string(iPhase));

        // Set field output type
        stk::io::set_field_output_type(*QplusSTKFieldPtr_, fieldType[1]);

        // Put the stk field on interior mesh parts for the current domain
        for (const stk::mesh::Part* part : partVec)
        {
            // check if already defined from a previous pass
            if (!QplusSTKFieldPtr_->defined_on(*part))
            {
                stk::mesh::put_field_on_mesh(
                    *QplusSTKFieldPtr_, *part, 1, &initialValue);
            }
        }
    }

    // # Qminus
    if (!QminusSTKFieldPtr_)
    {
        scalar initialValue = 0.0;

        QminusSTKFieldPtr_ =
            &metaData.declare_field<scalar>(stk::topology::NODE_RANK, "Qminus");

        // Set field output type
        stk::io::set_field_output_type(*QminusSTKFieldPtr_, fieldType[1]);

        // Put the stk field on interior mesh parts for the current domain
        for (const stk::mesh::Part* part : partVec)
        {
            // check if already defined from a previous pass
            if (!QminusSTKFieldPtr_->defined_on(*part))
            {
                stk::mesh::put_field_on_mesh(
                    *QminusSTKFieldPtr_, *part, 1, &initialValue);
            }
        }
    }

    // # Pplus
    if (!PplusSTKFieldPtr_)
    {
        scalar initialValue = 0.0;

        PplusSTKFieldPtr_ = &metaData.declare_field<scalar>(
            stk::topology::NODE_RANK, "Pplus" + std::to_string(iPhase));

        // Set field output type
        stk::io::set_field_output_type(*PplusSTKFieldPtr_, fieldType[1]);

        // Put the stk field on interior mesh parts for the current domain
        for (const stk::mesh::Part* part : partVec)
        {
            // check if already defined from a previous pass
            if (!PplusSTKFieldPtr_->defined_on(*part))
            {
                stk::mesh::put_field_on_mesh(
                    *PplusSTKFieldPtr_, *part, 1, &initialValue);
            }
        }
    }

    // # Pminus
    if (!PminusSTKFieldPtr_)
    {
        scalar initialValue = 0.0;

        PminusSTKFieldPtr_ = &metaData.declare_field<scalar>(
            stk::topology::NODE_RANK, "Pminus" + std::to_string(iPhase));

        // Set field output type
        stk::io::set_field_output_type(*PminusSTKFieldPtr_, fieldType[1]);

        // Put the stk field on interior mesh parts for the current domain
        for (const stk::mesh::Part* part : partVec)
        {
            // check if already defined from a previous pass
            if (!PminusSTKFieldPtr_->defined_on(*part))
            {
                stk::mesh::put_field_on_mesh(
                    *PminusSTKFieldPtr_, *part, 1, &initialValue);
            }
        }
    }

    // # sumA Plus
    if (!sumAPlusSTKFieldPtr_)
    {
        scalar initialValue = 0.0;

        sumAPlusSTKFieldPtr_ = &metaData.declare_field<scalar>(
            stk::topology::NODE_RANK, "sumAPlus" + std::to_string(iPhase));

        // Set field output type
        stk::io::set_field_output_type(*sumAPlusSTKFieldPtr_, fieldType[1]);

        // Put the stk field on interior mesh parts for the current domain
        for (const stk::mesh::Part* part : partVec)
        {
            // check if already defined from a previous pass
            if (!sumAPlusSTKFieldPtr_->defined_on(*part))
            {
                stk::mesh::put_field_on_mesh(
                    *sumAPlusSTKFieldPtr_, *part, 1, &initialValue);
            }
        }
    }

    // # sumA Minus
    if (!sumAMinusSTKFieldPtr_)
    {
        scalar initialValue = 0.0;

        sumAMinusSTKFieldPtr_ = &metaData.declare_field<scalar>(
            stk::topology::NODE_RANK, "sumAMinus" + std::to_string(iPhase));

        // Set field output type
        stk::io::set_field_output_type(*sumAMinusSTKFieldPtr_, fieldType[1]);

        // Put the stk field on interior mesh parts for the current domain
        for (const stk::mesh::Part* part : partVec)
        {
            // check if already defined from a previous pass
            if (!sumAMinusSTKFieldPtr_->defined_on(*part))
            {
                stk::mesh::put_field_on_mesh(
                    *sumAMinusSTKFieldPtr_, *part, 1, &initialValue);
            }
        }
    }

    // # Lambda Limiter
    if (!cMULESLimiterPlusSTKFieldPtr_)
    {
        scalar initialValue = 1.0;

        cMULESLimiterPlusSTKFieldPtr_ = &metaData.declare_field<scalar>(
            stk::topology::NODE_RANK, "lambdaNPlus" + std::to_string(iPhase));

        // Set field output type
        stk::io::set_field_output_type(*cMULESLimiterPlusSTKFieldPtr_,
                                       fieldType[1]);

        // Put the stk field on interior mesh parts for the current domain
        for (const stk::mesh::Part* part : partVec)
        {
            // check if already defined from a previous pass
            if (!cMULESLimiterPlusSTKFieldPtr_->defined_on(*part))
            {
                stk::mesh::put_field_on_mesh(
                    *cMULESLimiterPlusSTKFieldPtr_, *part, 1, &initialValue);
            }
        }
    }

    // # Lambda Limiter
    if (!cMULESLimiterMinusSTKFieldPtr_)
    {
        scalar initialValue = 1.0;

        cMULESLimiterMinusSTKFieldPtr_ = &metaData.declare_field<scalar>(
            stk::topology::NODE_RANK, "lambdaNMinus" + std::to_string(iPhase));

        // Set field output type
        stk::io::set_field_output_type(*cMULESLimiterMinusSTKFieldPtr_,
                                       fieldType[1]);

        // Put the stk field on interior mesh parts for the current domain
        for (const stk::mesh::Part* part : partVec)
        {
            // check if already defined from a previous pass
            if (!cMULESLimiterMinusSTKFieldPtr_->defined_on(*part))
            {
                stk::mesh::put_field_on_mesh(
                    *cMULESLimiterMinusSTKFieldPtr_, *part, 1, &initialValue);
            }
        }
    }

    // # Lambda Limiter Previous (for convergence checking)
    if (!cMULESLimiterPlusPrevSTKFieldPtr_)
    {
        scalar initialValue = 1.0;

        cMULESLimiterPlusPrevSTKFieldPtr_ = &metaData.declare_field<scalar>(
            stk::topology::NODE_RANK,
            "lambdaNPlusPrev" + std::to_string(iPhase));

        // Set field output type (no output needed for prev fields)
        stk::io::set_field_output_type(*cMULESLimiterPlusPrevSTKFieldPtr_,
                                       stk::io::FieldOutputType::SCALAR);

        // Put the stk field on interior mesh parts for the current domain
        for (const stk::mesh::Part* part : partVec)
        {
            // check if already defined from a previous pass
            if (!cMULESLimiterPlusPrevSTKFieldPtr_->defined_on(*part))
            {
                stk::mesh::put_field_on_mesh(*cMULESLimiterPlusPrevSTKFieldPtr_,
                                             *part,
                                             1,
                                             &initialValue);
            }
        }
    }

    // # Lambda Limiter Previous (for convergence checking)
    if (!cMULESLimiterMinusPrevSTKFieldPtr_)
    {
        scalar initialValue = 1.0;

        cMULESLimiterMinusPrevSTKFieldPtr_ = &metaData.declare_field<scalar>(
            stk::topology::NODE_RANK,
            "lambdaNMinusPrev" + std::to_string(iPhase));

        // Set field output type (no output needed for prev fields)
        stk::io::set_field_output_type(*cMULESLimiterMinusPrevSTKFieldPtr_,
                                       stk::io::FieldOutputType::SCALAR);

        // Put the stk field on interior mesh parts for the current domain
        for (const stk::mesh::Part* part : partVec)
        {
            // check if already defined from a previous pass
            if (!cMULESLimiterMinusPrevSTKFieldPtr_->defined_on(*part))
            {
                stk::mesh::put_field_on_mesh(
                    *cMULESLimiterMinusPrevSTKFieldPtr_,
                    *part,
                    1,
                    &initialValue);
            }
        }
    }
}

void freeSurfaceFlowModel::correctFCT(const std::shared_ptr<domain> domain,
                                      label iPhase)
{
    // Get Mesh Data
    auto& mesh = this->meshRef();
    stk::mesh::BulkData& bulkData = mesh.bulkDataRef();
    stk::mesh::MetaData& metaData = mesh.metaDataRef();

    // Step 1: Compute bounded (upwind) flux - this is our baseline
    computeFL_(domain, iPhase);

    // Corrector loop (nAlphaCorr is typically 1-2)
    const label nAlphaCorr =
        domain->multiphase_.freeSurfaceModel_.nAlphaCorrections_;
    for (label iter = 0; iter < nAlphaCorr; iter++)
    {
        // Step 2: Compute high-order flux with interface compression
        computeFH_(domain, iPhase);

        // Step 3: Compute antidiffusive (correction) flux: A = FH - FL
        computeA_(domain, iPhase);

        // Step 4: Compute limiter lambda using iterative MULES algorithm
        computeLambda_(domain, iPhase);

        // Step 5: Update alpha using the LIMITED correction flux
        // Note: updateAlpha_ now applies lambda*A directly without modifying A
        updateAlpha_(domain, iPhase, iter);

        // Sync alpha to ghost nodes for next iteration's computeFH_
        // (needed for opening boundary which reads alpha from ghost nodes)
        if (messager::parallel())
        {
            this->alphaRef(iPhase).synchronizeGhostedEntities(domain->index());
        }

        // Step 6: Update FL by adding the limited correction: FL += lambda*A
        // This makes FL approach FH over iterations
        updateFL_(domain, iPhase, iter);

        // Step 7: Update gradient and limiter for next iteration
        updateVolumeFractionGradientField(domain, iPhase);
        updateVolumeFractionBlendingFactorField(domain, iPhase);
    }
}

void freeSurfaceFlowModel::computeFL_(const std::shared_ptr<domain> domain,
                                      label iPhase)
{
    // Get Mesh Data
    auto& mesh = this->meshRef();
    stk::mesh::BulkData& bulkData = mesh.bulkDataRef();
    stk::mesh::MetaData& metaData = mesh.metaDataRef();

    // Interior
    {
        // Get fields
        const STKScalarField* mDotSTKFieldPtr =
            this->mDotRef(iPhase).stkFieldPtr();
        const STKScalarField* alphaSTKFieldPtr =
            this->alphaRef(iPhase).stkFieldPtr();

        // get interior parts the domain is defined on
        const stk::mesh::PartVector& partVec =
            domain->zonePtr()->interiorParts();

        // Define Scratch Spaces
        std::vector<scalar> ws_alpha;

        // define some common selectors; select owned nodes
        stk::mesh::Selector selUniversalElements =
            metaData.universal_part() & stk::mesh::selectUnion(partVec);

        stk::mesh::BucketVector const& elementBuckets = bulkData.get_buckets(
            stk::topology::ELEMENT_RANK, selUniversalElements);

        for (stk::mesh::BucketVector::const_iterator ib =
                 elementBuckets.begin();
             ib < elementBuckets.end();
             ib++)
        {
            stk::mesh::Bucket& elementBucket = **ib;
            const stk::mesh::Bucket::size_type nElemPerBucket =
                elementBucket.size();

            // extract master element, and its needed specifics
            MasterElement* meSCS =
                MasterElementRepo::get_surface_master_element(
                    elementBucket.topology());

            const label nodesPerElement = meSCS->nodesPerElement_;
            const label numScsIp = meSCS->numIntPoints_;
            const label* lrscv = meSCS->adjacentNodes();

            // allocate space for scratch spaces
            ws_alpha.resize(nodesPerElement);

            // get pointers
            scalar* p_alpha = &ws_alpha[0];

            // Element Loop
            for (stk::mesh::Bucket::size_type iElement = 0;
                 iElement < nElemPerBucket;
                 iElement++)
            {
                stk::mesh::Entity elem = elementBucket[iElement];
                stk::mesh::Entity const* nodeRels =
                    elementBucket.begin_nodes(iElement);

                // FCT Fields
                scalar* FL = stk::mesh::field_data(
                    *FLSTKFieldPtr_, elementBucket, iElement);

                for (label iNode = 0; iNode < nodesPerElement; iNode++)
                {
                    stk::mesh::Entity node = nodeRels[iNode];

                    p_alpha[iNode] =
                        *stk::mesh::field_data(*alphaSTKFieldPtr, node);
                }

                for (label ip = 0; ip < numScsIp; ++ip)
                {
                    // left and right nodes for this ip
                    const label il = lrscv[2 * ip];
                    const label ir = lrscv[2 * ip + 1];

                    // save off mDot
                    const scalar tmDot =
                        (stk::mesh::field_data(*mDotSTKFieldPtr, elem))[ip];

                    if (tmDot > 0.0)
                    {
                        FL[ip] = tmDot * p_alpha[il];
                    }
                    else
                    {
                        FL[ip] = tmDot * p_alpha[ir];
                    }
                }
            }
        }
    }

    // Boundary
    for (label iBoundary = 0; iBoundary < domain->zonePtr()->nBoundaries();
         iBoundary++)
    {
        const auto* boundary = domain->zonePtr()->boundaryPtr(iBoundary);

        boundaryPhysicalType type = boundary->type();

        boundaryConditionType bcType =
            this->alphaRef(iPhase)
                .boundaryConditionRef(domain->index(), iBoundary)
                .type();

        switch (type)
        {
            case boundaryPhysicalType::inlet:
                {
                    switch (bcType)
                    {
                        case boundaryConditionType::specifiedValue:
                            {
                                // Get fields
                                const auto& alphaSideSTKFieldRef =
                                    this->alphaRef(iPhase)
                                        .sideFieldRef()
                                        .stkFieldRef();

                                const STKScalarField* mDotSideSTKFieldPtr =
                                    this->mDotRef(iPhase)
                                        .sideFieldRef()
                                        .stkFieldPtr();

                                // define some common selectors
                                stk::mesh::Selector selAllSides =
                                    metaData.universal_part() &
                                    stk::mesh::selectUnion(boundary->parts());

                                stk::mesh::BucketVector const& sideBuckets =
                                    bulkData.get_buckets(metaData.side_rank(),
                                                         selAllSides);

                                for (stk::mesh::BucketVector::const_iterator
                                         ib = sideBuckets.begin();
                                     ib != sideBuckets.end();
                                     ++ib)
                                {
                                    stk::mesh::Bucket& sideBucket = **ib;

                                    const stk::mesh::Bucket::size_type
                                        nSidesPerBucket = sideBucket.size();

                                    const stk::topology theTopo =
                                        sideBucket.topology();

                                    // face master element
                                    MasterElement* meFC = MasterElementRepo::
                                        get_surface_master_element(theTopo);

                                    const label numScsBip = meFC->numIntPoints_;

                                    for (stk::mesh::Bucket::size_type iSide = 0;
                                         iSide < nSidesPerBucket;
                                         ++iSide)
                                    {
                                        const auto& side = sideBucket[iSide];

                                        const scalar* alphaBc =
                                            stk::mesh::field_data(
                                                alphaSideSTKFieldRef,
                                                sideBucket,
                                                iSide);
                                        const scalar* mDot =
                                            stk::mesh::field_data(
                                                *mDotSideSTKFieldPtr,
                                                sideBucket,
                                                iSide);

                                        // FCT Side Fields
                                        scalar* sideFL = stk::mesh::field_data(
                                            *sideFLSTKFieldPtr_,
                                            sideBucket,
                                            iSide);

                                        for (label ip = 0; ip < numScsBip; ++ip)
                                        {
                                            // save off mDot
                                            const scalar tmDot = mDot[ip];

                                            sideFL[ip] = tmDot * alphaBc[ip];
                                        }
                                    }
                                }
                            }
                            break;

                        default:
                            errorMsg("Must not reach here");
                    }
                }
                break;

            case boundaryPhysicalType::outlet:
                {
                    switch (bcType)
                    {
                        case boundaryConditionType::zeroGradient:
                            {
                                // Get fields
                                const STKScalarField* mDotSideSTKFieldPtr =
                                    this->mDotRef(iPhase)
                                        .sideFieldRef()
                                        .stkFieldPtr();
                                const STKScalarField* alphaSTKFieldPtr =
                                    this->alphaRef(iPhase).stkFieldPtr();

                                std::vector<scalar> ws_alpha;

                                // define some common selectors
                                stk::mesh::Selector selAllSides =
                                    metaData.universal_part() &
                                    stk::mesh::selectUnion(boundary->parts());

                                stk::mesh::BucketVector const& sideBuckets =
                                    bulkData.get_buckets(metaData.side_rank(),
                                                         selAllSides);

                                for (stk::mesh::BucketVector::const_iterator
                                         ib = sideBuckets.begin();
                                     ib != sideBuckets.end();
                                     ++ib)
                                {
                                    stk::mesh::Bucket& sideBucket = **ib;

                                    const stk::mesh::Bucket::size_type
                                        nSidesPerBucket = sideBucket.size();

                                    const stk::topology theTopo =
                                        sideBucket.topology();

                                    // face master element
                                    MasterElement* meFC = MasterElementRepo::
                                        get_surface_master_element(theTopo);

                                    const label numScsBip = meFC->numIntPoints_;

                                    // extract master element specifics
                                    const label numNodesPerSide =
                                        meFC->nodesPerElement_;

                                    // set sizes
                                    ws_alpha.resize(numNodesPerSide);

                                    // get pointers
                                    scalar* p_alpha = &ws_alpha[0];

                                    for (stk::mesh::Bucket::size_type iSide = 0;
                                         iSide < nSidesPerBucket;
                                         ++iSide)
                                    {
                                        const auto& side = sideBucket[iSide];

                                        // mapping from ip to nodes for this
                                        // ordinal
                                        const label* ipNodeMap =
                                            meFC->ipNodeMap();

                                        stk::mesh::Entity const* sideRels =
                                            bulkData.begin_nodes(side);

                                        const scalar* mDot =
                                            stk::mesh::field_data(
                                                *mDotSideSTKFieldPtr,
                                                sideBucket,
                                                iSide);

                                        // FCT Side Fields
                                        scalar* sideFL = stk::mesh::field_data(
                                            *sideFLSTKFieldPtr_,
                                            sideBucket,
                                            iSide);

                                        // fill with nodal values
                                        for (label iNode = 0;
                                             iNode < numNodesPerSide;
                                             iNode++)
                                        {
                                            stk::mesh::Entity node =
                                                sideRels[iNode];

                                            p_alpha[iNode] =
                                                *stk::mesh::field_data(
                                                    *alphaSTKFieldPtr, node);
                                        }

                                        for (label ip = 0; ip < numScsBip; ++ip)
                                        {
                                            const label nearestNode =
                                                ipNodeMap[ip];

                                            // save off mDot
                                            const scalar tmDot = mDot[ip];

                                            sideFL[ip] =
                                                tmDot * p_alpha[nearestNode];
                                        }
                                    }
                                }
                            }
                            break;

                        default:
                            errorMsg("Must not reach here");
                    }
                }
                break;

            case boundaryPhysicalType::opening:
                {
                    // Get fields - need both side field (for inflow) and nodal
                    // field (for outflow)
                    const auto& alphaSideSTKFieldRef =
                        this->alphaRef(iPhase).sideFieldRef().stkFieldRef();

                    const STKScalarField* mDotSideSTKFieldPtr =
                        this->mDotRef(iPhase).sideFieldRef().stkFieldPtr();
                    const STKScalarField* alphaSTKFieldPtr =
                        this->alphaRef(iPhase).stkFieldPtr();

                    std::vector<scalar> ws_alpha;

                    // define some common selectors
                    stk::mesh::Selector selAllSides =
                        metaData.universal_part() &
                        stk::mesh::selectUnion(boundary->parts());

                    stk::mesh::BucketVector const& sideBuckets =
                        bulkData.get_buckets(metaData.side_rank(), selAllSides);

                    for (stk::mesh::BucketVector::const_iterator ib =
                             sideBuckets.begin();
                         ib != sideBuckets.end();
                         ++ib)
                    {
                        stk::mesh::Bucket& sideBucket = **ib;

                        const stk::mesh::Bucket::size_type nSidesPerBucket =
                            sideBucket.size();

                        const stk::topology theTopo = sideBucket.topology();

                        // face master element
                        MasterElement* meFC =
                            MasterElementRepo::get_surface_master_element(
                                theTopo);

                        const label numScsBip = meFC->numIntPoints_;

                        // extract master element specifics
                        const label numNodesPerSide = meFC->nodesPerElement_;

                        // set sizes
                        ws_alpha.resize(numNodesPerSide);

                        // get pointers
                        scalar* p_alpha = &ws_alpha[0];

                        for (stk::mesh::Bucket::size_type iSide = 0;
                             iSide < nSidesPerBucket;
                             ++iSide)
                        {
                            const auto& side = sideBucket[iSide];

                            // mapping from ip to nodes for this ordinal
                            const label* ipNodeMap = meFC->ipNodeMap();

                            stk::mesh::Entity const* sideRels =
                                bulkData.begin_nodes(side);

                            const scalar* alphaBc = stk::mesh::field_data(
                                alphaSideSTKFieldRef, sideBucket, iSide);
                            const scalar* mDot = stk::mesh::field_data(
                                *mDotSideSTKFieldPtr, sideBucket, iSide);

                            // FCT Side Fields
                            scalar* sideFL = stk::mesh::field_data(
                                *sideFLSTKFieldPtr_, sideBucket, iSide);

                            // fill with nodal values (for outflow case)
                            for (label iNode = 0; iNode < numNodesPerSide;
                                 iNode++)
                            {
                                stk::mesh::Entity node = sideRels[iNode];

                                p_alpha[iNode] = *stk::mesh::field_data(
                                    *alphaSTKFieldPtr, node);
                            }

                            for (label ip = 0; ip < numScsBip; ++ip)
                            {
                                const label nearestNode = ipNodeMap[ip];

                                // save off mDot
                                const scalar tmDot = mDot[ip];

                                // Opening: switch based on flow direction
                                // mDot > 0: flow leaving (zero-gradient)
                                // mDot < 0: flow entering (prescribed value)
                                if (tmDot > 0.0)
                                {
                                    // Outflow: use internal value
                                    sideFL[ip] = tmDot * p_alpha[nearestNode];
                                }
                                else
                                {
                                    // Inflow: use prescribed boundary value
                                    sideFL[ip] = tmDot * alphaBc[ip];
                                }
                            }
                        }
                    }
                }
                break;

            default:
                break;
        }
    }
}

void freeSurfaceFlowModel::computeFH_(const std::shared_ptr<domain> domain,
                                      label iPhase)
{
    // Get Mesh Data
    auto& mesh = this->meshRef();
    stk::mesh::BulkData& bulkData = mesh.bulkDataRef();
    stk::mesh::MetaData& metaData = mesh.metaDataRef();

    // Interior
    {
        // Get fields
        const STKScalarField* mDotSTKFieldPtr =
            this->mDotRef(iPhase).stkFieldPtr();
        const STKScalarField* alphaSTKFieldPtr =
            this->alphaRef(iPhase).stkFieldPtr();
        const STKScalarField* gradAlphaSTKFieldPtr =
            this->alphaRef(iPhase).gradRef().stkFieldPtr();
        const STKScalarField* rhoSTKFieldPtr =
            this->rhoRef(iPhase).stkFieldPtr();
        const STKScalarField* USTKFieldPtr = this->URef().stkFieldPtr();
        const STKScalarField& betaSTKFieldRef =
            this->alphaRef(iPhase).blendingFactorRef().stkFieldRef();
        const STKScalarField* nHatSTKFieldPtr =
            this->nHatRef(iPhase).stkFieldPtr();

        // Get geometric fields
        const auto& coordsSTKFieldRef = *metaData.get_field<scalar>(
            stk::topology::NODE_RANK, this->getCoordinatesID_(domain));

        // get interior parts the domain is defined on
        const stk::mesh::PartVector& partVec =
            domain->zonePtr()->interiorParts();

        // Interface compression coefficient
        const scalar gamma = this->gamma(domain.get());

        // Define Scratch Spaces
        std::vector<scalar> ws_alpha;
        std::vector<scalar> ws_rho;
        std::vector<scalar> ws_beta;
        std::vector<scalar> ws_shape_function;
        std::vector<scalar> ws_U;
        std::vector<scalar> ws_gradAlpha;
        std::vector<scalar> ws_nHat;
        std::vector<scalar> ws_coordinates;
        std::vector<scalar> ws_scs_areav;

        // ip values
        std::vector<scalar> uIp(SPATIAL_DIM);
        std::vector<scalar> coordIp(SPATIAL_DIM);
        std::vector<scalar> gradAlphaIp(SPATIAL_DIM);
        std::vector<scalar> nHatIp(SPATIAL_DIM);

        // pointers for fast access
        scalar* p_uIp = &uIp[0];
        scalar* p_gradAlphaIp = &gradAlphaIp[0];
        scalar* p_coordIp = &coordIp[0];
        scalar* p_nHatIp = &nHatIp[0];

        // define some common selectors; select owned nodes
        stk::mesh::Selector selUniversalElements =
            metaData.universal_part() & stk::mesh::selectUnion(partVec);

        stk::mesh::BucketVector const& elementBuckets = bulkData.get_buckets(
            stk::topology::ELEMENT_RANK, selUniversalElements);

        bool isAlphaShifted = this->alphaRef(iPhase).isShifted();

        for (stk::mesh::BucketVector::const_iterator ib =
                 elementBuckets.begin();
             ib < elementBuckets.end();
             ib++)
        {
            stk::mesh::Bucket& elementBucket = **ib;
            const stk::mesh::Bucket::size_type nElemPerBucket =
                elementBucket.size();

            // extract master element, and its needed specifics
            MasterElement* meSCS =
                MasterElementRepo::get_surface_master_element(
                    elementBucket.topology());

            const label nodesPerElement = meSCS->nodesPerElement_;
            const label numScsIp = meSCS->numIntPoints_;
            const label* lrscv = meSCS->adjacentNodes();

            // allocate space for scratch spaces
            ws_alpha.resize(nodesPerElement);
            ws_rho.resize(nodesPerElement);
            ws_beta.resize(nodesPerElement);
            ws_shape_function.resize(numScsIp * nodesPerElement);
            ws_U.resize(nodesPerElement * SPATIAL_DIM);
            ws_gradAlpha.resize(nodesPerElement * SPATIAL_DIM);
            ws_nHat.resize(nodesPerElement * SPATIAL_DIM);
            ws_coordinates.resize(nodesPerElement * SPATIAL_DIM);
            ws_scs_areav.resize(numScsIp * SPATIAL_DIM);

            // get pointers
            scalar* p_alpha = &ws_alpha[0];
            scalar* p_rho = &ws_rho[0];
            scalar* p_beta = &ws_beta[0];
            scalar* p_shape_function = &ws_shape_function[0];
            scalar* p_U = &ws_U[0];
            scalar* p_gradAlpha = &ws_gradAlpha[0];
            scalar* p_nHat = &ws_nHat[0];
            scalar* p_coordinates = &ws_coordinates[0];
            scalar* p_scs_areav = &ws_scs_areav[0];

            if (isAlphaShifted)
            {
                meSCS->shifted_shape_fcn(&p_shape_function[0]);
            }
            else
            {
                meSCS->shape_fcn(&p_shape_function[0]);
            }

            // Element Loop
            for (stk::mesh::Bucket::size_type iElement = 0;
                 iElement < nElemPerBucket;
                 iElement++)
            {
                stk::mesh::Entity elem = elementBucket[iElement];

                stk::mesh::Entity const* nodeRels =
                    elementBucket.begin_nodes(iElement);

                scalar* FH = stk::mesh::field_data(
                    *FHSTKFieldPtr_, elementBucket, iElement);

                for (label iNode = 0; iNode < nodesPerElement; iNode++)
                {
                    stk::mesh::Entity node = nodeRels[iNode];

                    // gather scalars
                    p_alpha[iNode] =
                        *stk::mesh::field_data(*alphaSTKFieldPtr, node);
                    p_rho[iNode] =
                        *stk::mesh::field_data(*rhoSTKFieldPtr, node);
                    p_beta[iNode] =
                        *stk::mesh::field_data(betaSTKFieldRef, node);

                    // gather vectors
                    scalar* coords =
                        stk::mesh::field_data(coordsSTKFieldRef, node);
                    scalar* U = stk::mesh::field_data(*USTKFieldPtr, node);
                    scalar* gradAlpha =
                        stk::mesh::field_data(*gradAlphaSTKFieldPtr, node);
                    scalar* nHat =
                        stk::mesh::field_data(*nHatSTKFieldPtr, node);

                    const label offSet = iNode * SPATIAL_DIM;
                    for (label j = 0; j < SPATIAL_DIM; ++j)
                    {
                        p_coordinates[offSet + j] = coords[j];
                        p_U[offSet + j] = U[j];
                        p_gradAlpha[offSet + j] = gradAlpha[j];
                        p_nHat[offSet + j] = nHat[j];
                    }
                }

                // compute geometry
                scalar scs_error = 0.0;
                meSCS->determinant(
                    1, &p_coordinates[0], &p_scs_areav[0], &scs_error);

                for (label ip = 0; ip < numScsIp; ++ip)
                {
                    const label ipNdim = ip * SPATIAL_DIM;

                    // left and right nodes for this ip
                    const label il = lrscv[2 * ip];
                    const label ir = lrscv[2 * ip + 1];

                    // save off mDot
                    const scalar tmDot =
                        (stk::mesh::field_data(*mDotSTKFieldPtr, elem))[ip];

                    // zero out values of interest for this ip
                    for (label j = 0; j < SPATIAL_DIM; ++j)
                    {
                        p_uIp[j] = 0.0;
                        p_coordIp[j] = 0.0;
                        p_gradAlphaIp[j] = 0.0;
                        p_nHatIp[j] = 0.0;
                    }

                    // save off ip values; offset to Shape Function
                    const label offSetSF = ip * nodesPerElement;
                    scalar alphaIp = 0.0;
                    scalar rhoIp = 0.0;
                    for (label ic = 0; ic < nodesPerElement; ++ic)
                    {
                        const scalar r = p_shape_function[offSetSF + ic];

                        // compute scs ip value
                        alphaIp += r * p_alpha[ic];
                        rhoIp += r * p_rho[ic];

                        for (label j = 0; j < SPATIAL_DIM; ++j)
                        {
                            const scalar uj = p_U[ic * SPATIAL_DIM + j];
                            const scalar gradAlphaj =
                                p_gradAlpha[ic * SPATIAL_DIM + j];
                            const scalar nHatj = p_nHat[ic * SPATIAL_DIM + j];
                            p_uIp[j] += r * uj;
                            p_coordIp[j] +=
                                r * p_coordinates[ic * SPATIAL_DIM + j];
                            p_gradAlphaIp[j] += r * gradAlphaj;
                            p_nHatIp[j] += r * nHatj;
                        }
                    }

                    // assemble advection; rhs and upwind contributions
                    scalar alphaUpwind;
                    scalar dcorr = 0.0;
                    if (tmDot > 0)
                    {
                        alphaUpwind = p_alpha[il];

                        // deferred correction
                        for (label j = 0; j < SPATIAL_DIM; ++j)
                        {
                            const scalar dxj =
                                p_coordIp[j] -
                                p_coordinates[il * SPATIAL_DIM + j];
                            dcorr += p_beta[il] * dxj *
                                     p_gradAlpha[il * SPATIAL_DIM + j];
                        }
                    }
                    else
                    {
                        alphaUpwind = p_alpha[ir];

                        // deferred correction
                        for (label j = 0; j < SPATIAL_DIM; ++j)
                        {
                            const scalar dxj =
                                p_coordIp[j] -
                                p_coordinates[ir * SPATIAL_DIM + j];
                            dcorr += p_beta[ir] * dxj *
                                     p_gradAlpha[ir * SPATIAL_DIM + j];
                        }
                    }

                    // Interface compression term:
                    // vfc = cAlpha * sign(phi) * alpha * (1-alpha) * (nHat ·
                    // nSf)
                    // where nHat = interface normal (unit), nSf = face unit
                    // normal When multiplied by phi to get flux:
                    // compression_flux = cAlpha * |phi| * alpha * (1-alpha) *
                    // (nHat · nSf)
                    // Note: (nHat · nSf) is a continuous value in [-1, 1], not
                    // just sign
                    // Note: nHat is lagged by one time step to improve
                    // convergence

                    // Compute sub-control surface area magnitude
                    scalar sMag = 0.0;
                    for (label j = 0; j < SPATIAL_DIM; ++j)
                    {
                        sMag +=
                            p_scs_areav[ipNdim + j] * p_scs_areav[ipNdim + j];
                    }
                    sMag = std::sqrt(sMag);

                    // Compute (nHat · nS) = interface normal · face unit
                    // normal using stored (lagged) nHat field
                    scalar nHatDotNs = 0.0;
                    for (label j = 0; j < SPATIAL_DIM; ++j)
                    {
                        scalar nSj = p_scs_areav[ipNdim + j] / (sMag + SMALL);
                        nHatDotNs += p_nHatIp[j] * nSj;
                    }

                    // Compression flux: gamma * |mDot| * alpha * (1-alpha)
                    // * (nHat · nS)
                    scalar compression = gamma * std::abs(tmDot) * alphaIp *
                                         (1.0 - alphaIp) * nHatDotNs;

                    // calculate FH: high-order convective flux + compression
                    FH[ip] = tmDot * (alphaUpwind + dcorr) + compression;
                }
            }
        }
    }

    // Boundary
    for (label iBoundary = 0; iBoundary < domain->zonePtr()->nBoundaries();
         iBoundary++)
    {
        const auto* boundary = domain->zonePtr()->boundaryPtr(iBoundary);

        boundaryPhysicalType type = boundary->type();

        boundaryConditionType bcType =
            this->alphaRef(iPhase)
                .boundaryConditionRef(domain->index(), iBoundary)
                .type();

        switch (type)
        {
            case boundaryPhysicalType::inlet:
                {
                    switch (bcType)
                    {
                        case boundaryConditionType::specifiedValue:
                            {
                                // Get fields
                                const auto& alphaSideSTKFieldRef =
                                    this->alphaRef(iPhase)
                                        .sideFieldRef()
                                        .stkFieldRef();

                                const STKScalarField* mDotSideSTKFieldPtr =
                                    this->mDotRef(iPhase)
                                        .sideFieldRef()
                                        .stkFieldPtr();

                                // define some common selectors
                                stk::mesh::Selector selAllSides =
                                    metaData.universal_part() &
                                    stk::mesh::selectUnion(boundary->parts());

                                stk::mesh::BucketVector const& sideBuckets =
                                    bulkData.get_buckets(metaData.side_rank(),
                                                         selAllSides);

                                for (stk::mesh::BucketVector::const_iterator
                                         ib = sideBuckets.begin();
                                     ib != sideBuckets.end();
                                     ++ib)
                                {
                                    stk::mesh::Bucket& sideBucket = **ib;

                                    const stk::mesh::Bucket::size_type
                                        nSidesPerBucket = sideBucket.size();

                                    const stk::topology theTopo =
                                        sideBucket.topology();

                                    // face master element
                                    MasterElement* meFC = MasterElementRepo::
                                        get_surface_master_element(theTopo);

                                    const label numScsBip = meFC->numIntPoints_;

                                    for (stk::mesh::Bucket::size_type iSide = 0;
                                         iSide < nSidesPerBucket;
                                         ++iSide)
                                    {
                                        const auto& side = sideBucket[iSide];

                                        const scalar* alphaBc =
                                            stk::mesh::field_data(
                                                alphaSideSTKFieldRef,
                                                sideBucket,
                                                iSide);
                                        const scalar* mDot =
                                            stk::mesh::field_data(
                                                *mDotSideSTKFieldPtr,
                                                sideBucket,
                                                iSide);

                                        // FCT Side Fields
                                        scalar* sideFH = stk::mesh::field_data(
                                            *sideFHSTKFieldPtr_,
                                            sideBucket,
                                            iSide);

                                        for (label ip = 0; ip < numScsBip; ++ip)
                                        {
                                            // save off mDot
                                            const scalar tmDot = mDot[ip];

                                            // compression factor not calculated
                                            // at the boundary
                                            sideFH[ip] = tmDot * alphaBc[ip];
                                        }
                                    }
                                }
                            }
                            break;

                        default:
                            errorMsg("Must not reach here");
                    }
                }
                break;

            case boundaryPhysicalType::outlet:
                {
                    switch (bcType)
                    {
                        case boundaryConditionType::zeroGradient:
                            {
                                // Get fields
                                const STKScalarField* mDotSideSTKFieldPtr =
                                    this->mDotRef(iPhase)
                                        .sideFieldRef()
                                        .stkFieldPtr();
                                const STKScalarField* alphaSTKFieldPtr =
                                    this->alphaRef(iPhase).stkFieldPtr();

                                std::vector<scalar> ws_alpha;

                                // define some common selectors
                                stk::mesh::Selector selAllSides =
                                    metaData.universal_part() &
                                    stk::mesh::selectUnion(boundary->parts());

                                stk::mesh::BucketVector const& sideBuckets =
                                    bulkData.get_buckets(metaData.side_rank(),
                                                         selAllSides);

                                for (stk::mesh::BucketVector::const_iterator
                                         ib = sideBuckets.begin();
                                     ib != sideBuckets.end();
                                     ++ib)
                                {
                                    stk::mesh::Bucket& sideBucket = **ib;

                                    const stk::mesh::Bucket::size_type
                                        nSidesPerBucket = sideBucket.size();

                                    const stk::topology theTopo =
                                        sideBucket.topology();

                                    // face master element
                                    MasterElement* meFC = MasterElementRepo::
                                        get_surface_master_element(theTopo);

                                    const label numScsBip = meFC->numIntPoints_;

                                    // extract master element specifics
                                    const label numNodesPerSide =
                                        meFC->nodesPerElement_;

                                    // set sizes
                                    ws_alpha.resize(numNodesPerSide);

                                    // get pointers
                                    scalar* p_alpha = &ws_alpha[0];

                                    for (stk::mesh::Bucket::size_type iSide = 0;
                                         iSide < nSidesPerBucket;
                                         ++iSide)
                                    {
                                        const auto& side = sideBucket[iSide];

                                        // mapping from ip to nodes for this
                                        // ordinal
                                        const label* ipNodeMap =
                                            meFC->ipNodeMap();

                                        stk::mesh::Entity const* sideRels =
                                            bulkData.begin_nodes(side);

                                        const scalar* mDot =
                                            stk::mesh::field_data(
                                                *mDotSideSTKFieldPtr,
                                                sideBucket,
                                                iSide);

                                        // FCT Side Fields
                                        scalar* sideFH = stk::mesh::field_data(
                                            *sideFHSTKFieldPtr_,
                                            sideBucket,
                                            iSide);

                                        // fill with nodal values
                                        for (label iNode = 0;
                                             iNode < numNodesPerSide;
                                             iNode++)
                                        {
                                            stk::mesh::Entity node =
                                                sideRels[iNode];

                                            p_alpha[iNode] =
                                                *stk::mesh::field_data(
                                                    *alphaSTKFieldPtr, node);
                                        }

                                        for (label ip = 0; ip < numScsBip; ++ip)
                                        {
                                            const label nearestNode =
                                                ipNodeMap[ip];

                                            // save off mDot
                                            const scalar tmDot = mDot[ip];

                                            sideFH[ip] =
                                                tmDot * p_alpha[nearestNode];
                                        }
                                    }
                                }
                            }
                            break;

                        default:
                            errorMsg("Must not reach here");
                    }
                }
                break;

            case boundaryPhysicalType::opening:
                {
                    // Get fields - need both side field (for inflow) and nodal
                    // field (for outflow)
                    const auto& alphaSideSTKFieldRef =
                        this->alphaRef(iPhase).sideFieldRef().stkFieldRef();

                    const STKScalarField* mDotSideSTKFieldPtr =
                        this->mDotRef(iPhase).sideFieldRef().stkFieldPtr();
                    const STKScalarField* alphaSTKFieldPtr =
                        this->alphaRef(iPhase).stkFieldPtr();

                    std::vector<scalar> ws_alpha;

                    // define some common selectors
                    stk::mesh::Selector selAllSides =
                        metaData.universal_part() &
                        stk::mesh::selectUnion(boundary->parts());

                    stk::mesh::BucketVector const& sideBuckets =
                        bulkData.get_buckets(metaData.side_rank(), selAllSides);

                    for (stk::mesh::BucketVector::const_iterator ib =
                             sideBuckets.begin();
                         ib != sideBuckets.end();
                         ++ib)
                    {
                        stk::mesh::Bucket& sideBucket = **ib;

                        const stk::mesh::Bucket::size_type nSidesPerBucket =
                            sideBucket.size();

                        const stk::topology theTopo = sideBucket.topology();

                        // face master element
                        MasterElement* meFC =
                            MasterElementRepo::get_surface_master_element(
                                theTopo);

                        const label numScsBip = meFC->numIntPoints_;

                        // extract master element specifics
                        const label numNodesPerSide = meFC->nodesPerElement_;

                        // set sizes
                        ws_alpha.resize(numNodesPerSide);

                        // get pointers
                        scalar* p_alpha = &ws_alpha[0];

                        for (stk::mesh::Bucket::size_type iSide = 0;
                             iSide < nSidesPerBucket;
                             ++iSide)
                        {
                            const auto& side = sideBucket[iSide];

                            // mapping from ip to nodes for this ordinal
                            const label* ipNodeMap = meFC->ipNodeMap();

                            stk::mesh::Entity const* sideRels =
                                bulkData.begin_nodes(side);

                            const scalar* alphaBc = stk::mesh::field_data(
                                alphaSideSTKFieldRef, sideBucket, iSide);
                            const scalar* mDot = stk::mesh::field_data(
                                *mDotSideSTKFieldPtr, sideBucket, iSide);

                            // FCT Side Fields
                            scalar* sideFH = stk::mesh::field_data(
                                *sideFHSTKFieldPtr_, sideBucket, iSide);

                            // fill with nodal values (for outflow case)
                            for (label iNode = 0; iNode < numNodesPerSide;
                                 iNode++)
                            {
                                stk::mesh::Entity node = sideRels[iNode];

                                p_alpha[iNode] = *stk::mesh::field_data(
                                    *alphaSTKFieldPtr, node);
                            }

                            for (label ip = 0; ip < numScsBip; ++ip)
                            {
                                const label nearestNode = ipNodeMap[ip];

                                // save off mDot
                                const scalar tmDot = mDot[ip];

                                // Opening: switch based on flow direction
                                // mDot > 0: flow leaving (zero-gradient)
                                // mDot < 0: flow entering (prescribed value)
                                // Note: compression term not calculated at
                                // boundary
                                if (tmDot > 0.0)
                                {
                                    // Outflow: use internal value
                                    sideFH[ip] = tmDot * p_alpha[nearestNode];
                                }
                                else
                                {
                                    // Inflow: use prescribed boundary value
                                    sideFH[ip] = tmDot * alphaBc[ip];
                                }
                            }
                        }
                    }
                }
                break;

            default:
                break;
        }
    }
}

void freeSurfaceFlowModel::computeA_(const std::shared_ptr<domain> domain,
                                     label iPhase)
{
    // Get Mesh Data
    auto& mesh = this->meshRef();
    stk::mesh::BulkData& bulkData = mesh.bulkDataRef();
    stk::mesh::MetaData& metaData = mesh.metaDataRef();

    // Interior
    {
        // get interior parts the domain is defined on
        const stk::mesh::PartVector& partVec =
            domain->zonePtr()->interiorParts();
        // define some common selectors; select owned nodes
        stk::mesh::Selector selUniversalElements =
            metaData.universal_part() & stk::mesh::selectUnion(partVec);

        stk::mesh::BucketVector const& elementBuckets = bulkData.get_buckets(
            stk::topology::ELEMENT_RANK, selUniversalElements);

        for (stk::mesh::BucketVector::const_iterator ib =
                 elementBuckets.begin();
             ib < elementBuckets.end();
             ib++)
        {
            stk::mesh::Bucket& elementBucket = **ib;
            const stk::mesh::Bucket::size_type nElemPerBucket =
                elementBucket.size();

            // extract master element, and its needed specifics
            MasterElement* meSCS =
                MasterElementRepo::get_surface_master_element(
                    elementBucket.topology());

            const label numScsIp = meSCS->numIntPoints_;

            // Element Loop
            for (stk::mesh::Bucket::size_type iElement = 0;
                 iElement < nElemPerBucket;
                 iElement++)
            {
                // FCT Fields
                const scalar* FH = stk::mesh::field_data(
                    *FHSTKFieldPtr_, elementBucket, iElement);
                const scalar* FL = stk::mesh::field_data(
                    *FLSTKFieldPtr_, elementBucket, iElement);
                scalar* A = stk::mesh::field_data(
                    *ASTKFieldPtr_, elementBucket, iElement);

                for (label ip = 0; ip < numScsIp; ++ip)
                {
                    // Calculate Antidiffusive Flux A
                    A[ip] = FH[ip] - FL[ip];
                }
            }
        }
    }

    // Boundary
    for (label iBoundary = 0; iBoundary < domain->zonePtr()->nBoundaries();
         iBoundary++)
    {
        const auto* boundary = domain->zonePtr()->boundaryPtr(iBoundary);

        boundaryPhysicalType type = boundary->type();

        switch (type)
        {
            case boundaryPhysicalType::inlet:
            case boundaryPhysicalType::outlet:
            case boundaryPhysicalType::opening:
                {
                    // define some common selectors
                    stk::mesh::Selector selAllSides =
                        metaData.universal_part() &
                        stk::mesh::selectUnion(boundary->parts());

                    stk::mesh::BucketVector const& sideBuckets =
                        bulkData.get_buckets(metaData.side_rank(), selAllSides);

                    for (stk::mesh::BucketVector::const_iterator ib =
                             sideBuckets.begin();
                         ib != sideBuckets.end();
                         ++ib)
                    {
                        stk::mesh::Bucket& sideBucket = **ib;

                        const stk::mesh::Bucket::size_type nSidesPerBucket =
                            sideBucket.size();

                        const stk::topology theTopo = sideBucket.topology();

                        // face master element
                        MasterElement* meFC =
                            MasterElementRepo::get_surface_master_element(
                                theTopo);

                        const label numScsBip = meFC->numIntPoints_;

                        for (stk::mesh::Bucket::size_type iSide = 0;
                             iSide < nSidesPerBucket;
                             ++iSide)
                        {
                            const auto& side = sideBucket[iSide];

                            // FCT Side Fields
                            const scalar* sideFL = stk::mesh::field_data(
                                *sideFLSTKFieldPtr_, sideBucket, iSide);
                            const scalar* sideFH = stk::mesh::field_data(
                                *sideFHSTKFieldPtr_, sideBucket, iSide);

                            scalar* sideA = stk::mesh::field_data(
                                *sideASTKFieldPtr_, sideBucket, iSide);

                            for (label ip = 0; ip < numScsBip; ++ip)
                            {
                                // Calculate Antidiffusive Flux A
                                sideA[ip] = sideFH[ip] - sideFL[ip];
                            }
                        }
                    }
                }
                break;

            default:
                break;
        }
    }
}

void freeSurfaceFlowModel::computeLambda_(const std::shared_ptr<domain> domain,
                                          label iPhase)
{
    // Get Mesh Data
    auto& mesh = this->meshRef();
    stk::mesh::BulkData& bulkData = mesh.bulkDataRef();
    stk::mesh::MetaData& metaData = mesh.metaDataRef();

    // Step 1: Compute Q+ and Q- (max allowable change - remains constant)
    computeQ_(domain, iPhase);

    // Step 2: Compute P+ and P- from UNLIMITED fluxes (computed ONCE)
    // P+ = sum of positive incoming antidiffusive fluxes per node
    // P- = sum of negative outgoing antidiffusive fluxes per node
    computeP_(domain, iPhase);

    // Step 3: Initialize lambda fields to 1.0 (no limiting initially)
    scalar one = 1.0;
    ops::setValue(lambdaSTKFieldPtr_, &one, domain->zonePtr()->interiorParts());
    ops::setValue(
        sideLambdaSTKFieldPtr_, &one, domain->zonePtr()->interiorParts());
    ops::setValue(cMULESLimiterMinusSTKFieldPtr_,
                  &one,
                  domain->zonePtr()->interiorParts());
    ops::setValue(cMULESLimiterPlusSTKFieldPtr_,
                  &one,
                  domain->zonePtr()->interiorParts());

    // Step 4: Initialize sumA fields to zero for first iteration
    ops::zero(sumAMinusSTKFieldPtr_, domain->zonePtr()->interiorParts());
    ops::zero(sumAPlusSTKFieldPtr_, domain->zonePtr()->interiorParts());

    // Step 5: Iterative limiter computation
    scalar maxLambdaChange = 1.0;
    label iter = 0;

    while (iter < maxLambdaIterations_ && maxLambdaChange > lambdaTolerance_)
    {
        // Compute new lambda values at nodes with bounds [0,1]
        // Iteration 0: lambdaP = Q+ / P+, lambdaM = Q- / P-
        // Iteration > 0: lambdaP = (sumA+ + Q+) / P+, lambdaM = (sumA- + Q-) /
        // P-
        maxLambdaChange = computeLambdaNodeWithBounds_(domain, iPhase, iter);

        // For parallel runs: global reduction ensures all processors agree on
        // convergence, and sync limiter fields to ghost nodes
        if (messager::parallel())
        {
            scalar g_maxLambdaChange = maxLambdaChange;
            stk::all_reduce_max(
                bulkData.parallel(), &maxLambdaChange, &g_maxLambdaChange, 1);
            maxLambdaChange = g_maxLambdaChange;

            stk::mesh::communicate_field_data(bulkData,
                                              {cMULESLimiterPlusSTKFieldPtr_,
                                               cMULESLimiterMinusSTKFieldPtr_});
        }

        // Interpolate lambda to integration points using multi-dimensional
        // limiting lambda_ip = min(lambdaP[owner], lambdaM[neighbor]) or vice
        // versa based on flux direction
        computeLambdaIP_(domain, iPhase);

        // Compute sum of LIMITED fluxes for next iteration
        // sumA+ and sumA- accumulate lambda * A contributions
        computeSumA_(domain, iPhase);

        iter++;
    }
}

void freeSurfaceFlowModel::computeP_(const std::shared_ptr<domain> domain,
                                     label iPhase)
{
    // Get Mesh Data
    auto& mesh = this->meshRef();
    stk::mesh::BulkData& bulkData = mesh.bulkDataRef();
    stk::mesh::MetaData& metaData = mesh.metaDataRef();

    // Initialize P+ and P- Fields to Zero
    ops::zero(PplusSTKFieldPtr_, domain->zonePtr()->interiorParts());
    ops::zero(PminusSTKFieldPtr_, domain->zonePtr()->interiorParts());

    // Interior
    {
        // get interior parts the domain is defined on
        const stk::mesh::PartVector& partVec =
            domain->zonePtr()->interiorParts();

        // define some common selectors; select owned nodes
        stk::mesh::Selector selUniversalElements =
            metaData.universal_part() & stk::mesh::selectUnion(partVec);

        stk::mesh::BucketVector const& elementBuckets = bulkData.get_buckets(
            stk::topology::ELEMENT_RANK, selUniversalElements);

        for (stk::mesh::BucketVector::const_iterator ib =
                 elementBuckets.begin();
             ib < elementBuckets.end();
             ib++)
        {
            stk::mesh::Bucket& elementBucket = **ib;
            const stk::mesh::Bucket::size_type nElemPerBucket =
                elementBucket.size();

            // extract master element, and its needed specifics
            MasterElement* meSCS =
                MasterElementRepo::get_surface_master_element(
                    elementBucket.topology());

            const label numScsIp = meSCS->numIntPoints_;
            const label* lrscv = meSCS->adjacentNodes();

            // Element Bucket Loop
            for (stk::mesh::Bucket::size_type iElement = 0;
                 iElement < nElemPerBucket;
                 iElement++)
            {
                stk::mesh::Entity const* nodeRels =
                    elementBucket.begin_nodes(iElement);

                // FCT Fields
                const scalar* A = stk::mesh::field_data(
                    *ASTKFieldPtr_, elementBucket, iElement);

                for (label ip = 0; ip < numScsIp; ++ip)
                {
                    // left and right nodes for this ip
                    const label il = lrscv[2 * ip];
                    const label ir = lrscv[2 * ip + 1];

                    const stk::mesh::Entity nodeL = nodeRels[il];
                    const stk::mesh::Entity nodeR = nodeRels[ir];

                    // pointer to fields to assemble
                    scalar* PplusL =
                        stk::mesh::field_data(*PplusSTKFieldPtr_, nodeL);
                    scalar* PplusR =
                        stk::mesh::field_data(*PplusSTKFieldPtr_, nodeR);
                    scalar* PminusL =
                        stk::mesh::field_data(*PminusSTKFieldPtr_, nodeL);
                    scalar* PminusR =
                        stk::mesh::field_data(*PminusSTKFieldPtr_, nodeR);

                    // ##3 - Inflow and Outflow
                    if (A[ip] > 0.0)
                    {
                        PminusL[0] += A[ip]; // outflow from L
                        PplusR[0] += A[ip];  // inflow into R
                    }
                    else
                    {
                        PminusR[0] -= A[ip]; // outflow from R
                        PplusL[0] -= A[ip];  // inflow into L
                    }
                }
            }
        }
    }

    // Boundary
    for (label iBoundary = 0; iBoundary < domain->zonePtr()->nBoundaries();
         iBoundary++)
    {
        const auto* boundary = domain->zonePtr()->boundaryPtr(iBoundary);

        boundaryPhysicalType type = boundary->type();

        switch (type)
        {
            case boundaryPhysicalType::inlet:
            case boundaryPhysicalType::outlet:
            case boundaryPhysicalType::opening:
                {
                    // define some common selectors
                    stk::mesh::Selector selAllSides =
                        metaData.universal_part() &
                        stk::mesh::selectUnion(boundary->parts());

                    stk::mesh::BucketVector const& sideBuckets =
                        bulkData.get_buckets(metaData.side_rank(), selAllSides);

                    for (stk::mesh::BucketVector::const_iterator ib =
                             sideBuckets.begin();
                         ib != sideBuckets.end();
                         ++ib)
                    {
                        stk::mesh::Bucket& sideBucket = **ib;

                        const stk::mesh::Bucket::size_type nSidesPerBucket =
                            sideBucket.size();

                        const stk::topology theTopo = sideBucket.topology();

                        // face master element
                        MasterElement* meFC =
                            MasterElementRepo::get_surface_master_element(
                                theTopo);

                        const label numScsBip = meFC->numIntPoints_;

                        for (stk::mesh::Bucket::size_type iSide = 0;
                             iSide < nSidesPerBucket;
                             ++iSide)
                        {
                            const auto& side = sideBucket[iSide];

                            // mapping from ip to nodes for this ordinal
                            const label* ipNodeMap = meFC->ipNodeMap();

                            // Get node entities for this side
                            stk::mesh::Entity const* sideRels =
                                bulkData.begin_nodes(side);

                            // FCT Side Fields
                            const scalar* sideA = stk::mesh::field_data(
                                *sideASTKFieldPtr_, sideBucket, iSide);

                            for (label ip = 0; ip < numScsBip; ++ip)
                            {
                                const label nearestNode = ipNodeMap[ip];
                                stk::mesh::Entity node = sideRels[nearestNode];

                                scalar* sidePplusL = stk::mesh::field_data(
                                    *PplusSTKFieldPtr_, node);
                                scalar* sidePminusL = stk::mesh::field_data(
                                    *PminusSTKFieldPtr_, node);

                                // ##3 - Inflow and Outflow
                                if (sideA[ip] > 0.0)
                                {
                                    sidePminusL[0] += sideA[ip];
                                }
                                else
                                {
                                    sidePplusL[0] -= sideA[ip];
                                }
                            }
                        }
                    }
                }
                break;

            default:
                break;
        }
    }
}

void freeSurfaceFlowModel::computeQ_(const std::shared_ptr<domain> domain,
                                     label iPhase)
{
    // Corrector MULES approach
    //
    // Q+ and Q- represent the maximum allowable correction flux (lambda*A)
    // that keeps alpha within bounds, starting from current alpha (after
    // implicit upwind solve).
    //
    // Since the bounded flux (FL) has already been applied by the implicit
    // solve, we only need to limit the correction flux A = FH - FL.
    //
    // The correction update is:
    // alpha_new = alpha_current - (dt/(rho*V)) * sum(lambda*A)
    //
    // For alpha_new <= alphaMax:
    // Q+ = (rho*V/dt) * (alphaMax - alpha_current)
    //
    // For alpha_new >= alphaMin:
    // Q- = (rho*V/dt) * (alpha_current - alphaMin)
    //
    // Note: No sumFL term is needed because the bounded flux effect is already
    // incorporated in alpha_current from the implicit solve.

    // Get Mesh Data
    auto& mesh = this->meshRef();
    stk::mesh::BulkData& bulkData = mesh.bulkDataRef();
    stk::mesh::MetaData& metaData = mesh.metaDataRef();

    // Get fields - use CURRENT alpha (after implicit solve)
    const STKScalarField* alphaSTKFieldPtr =
        this->alphaRef(iPhase).stkFieldPtr();
    STKScalarField* alphaMaxSTKFieldPtr =
        this->alphaRef(iPhase).maxValueRef().stkFieldPtr();
    STKScalarField* alphaMinSTKFieldPtr =
        this->alphaRef(iPhase).minValueRef().stkFieldPtr();
    const STKScalarField* rhoSTKFieldPtr = this->rhoRef(iPhase).stkFieldPtr();

    // Get geometric fields
    const auto* volSTKFieldPtr = metaData.get_field<scalar>(
        stk::topology::NODE_RANK, mesh::dual_nodal_volume_ID);

    // Get user defined fields
    const bool is_transient = this->controlsRef().isTransient();
    const scalar dt = is_transient ? this->controlsRef().getTimestep()
                                   : this->controlsRef().getPhysicalTimescale();

    // get interior parts the domain is defined on
    const stk::mesh::PartVector& partVec = domain->zonePtr()->interiorParts();

    // Update local min/max fields based on current alpha
    {
        const scalar extremaCoeff = 0.0;
        const scalar smoothLimiter = 0.0;

        this->alphaRef(iPhase).updateMinMaxFields(
            domain->index(), true, extremaCoeff, smoothLimiter);
    }

    // Compute Q+ and Q- based on current alpha (corrector MULES formula)
    {
        stk::mesh::Selector selOwnedNodes =
            metaData.locally_owned_part() & stk::mesh::selectUnion(partVec);

        stk::mesh::BucketVector const& nodeBuckets =
            bulkData.get_buckets(stk::topology::NODE_RANK, selOwnedNodes);

        for (stk::mesh::BucketVector::const_iterator ib = nodeBuckets.begin();
             ib != nodeBuckets.end();
             ++ib)
        {
            stk::mesh::Bucket& nodeBucket = **ib;
            const stk::mesh::Bucket::size_type nNodesPerBucket =
                nodeBucket.size();

            for (stk::mesh::Bucket::size_type iNode = 0;
                 iNode < nNodesPerBucket;
                 ++iNode)
            {
                stk::mesh::Entity node = nodeBucket[iNode];

                // Use current alpha (after implicit upwind solve)
                const scalar* alpha =
                    stk::mesh::field_data(*alphaSTKFieldPtr, node);
                const scalar* rho =
                    stk::mesh::field_data(*rhoSTKFieldPtr, node);
                const scalar* vol =
                    stk::mesh::field_data(*volSTKFieldPtr, node);

                scalar* Qplus = stk::mesh::field_data(*QplusSTKFieldPtr_, node);
                scalar* Qminus =
                    stk::mesh::field_data(*QminusSTKFieldPtr_, node);

                // Get local extrema (computed in updateMinMaxFields)
                const scalar* alphaMax =
                    stk::mesh::field_data(*alphaMaxSTKFieldPtr, node);
                const scalar* alphaMin =
                    stk::mesh::field_data(*alphaMinSTKFieldPtr, node);

                // Clip to global bounds [0, 1]
                const scalar globalMax = 1.0;
                const scalar globalMin = 0.0;
                const scalar clippedMax =
                    std::max(globalMin, std::min(alphaMax[0], globalMax));
                const scalar clippedMin =
                    std::max(globalMin, std::min(alphaMin[0], globalMax));

                // Corrector MULES formula
                // Q+ = (rho*V/dt) * (alphaMax - alpha_current) Q- =
                // (rho*V/dt) * (alpha_current - alphaMin)
                const scalar rhoVdt = rho[0] * vol[0] / dt;
                Qplus[0] = rhoVdt * (clippedMax - alpha[0]);
                Qminus[0] = rhoVdt * (alpha[0] - clippedMin);
            }
        }
    }
}

void freeSurfaceFlowModel::computeSumA_(const std::shared_ptr<domain> domain,
                                        label iPhase)
{
    // Get Mesh Data
    auto& mesh = this->meshRef();
    stk::mesh::BulkData& bulkData = mesh.bulkDataRef();
    stk::mesh::MetaData& metaData = mesh.metaDataRef();

    // Initialize sumA+ and SumA- to zero
    ops::zero(sumAMinusSTKFieldPtr_, domain->zonePtr()->interiorParts());
    ops::zero(sumAPlusSTKFieldPtr_, domain->zonePtr()->interiorParts());

    // Interior
    {
        // get interior parts the domain is defined on
        const stk::mesh::PartVector& partVec =
            domain->zonePtr()->interiorParts();
        // define some common selectors; select owned nodes
        stk::mesh::Selector selUniversalElements =
            metaData.universal_part() & stk::mesh::selectUnion(partVec);

        stk::mesh::BucketVector const& elementBuckets = bulkData.get_buckets(
            stk::topology::ELEMENT_RANK, selUniversalElements);

        for (stk::mesh::BucketVector::const_iterator ib =
                 elementBuckets.begin();
             ib < elementBuckets.end();
             ib++)
        {
            stk::mesh::Bucket& elementBucket = **ib;
            const stk::mesh::Bucket::size_type nElemPerBucket =
                elementBucket.size();

            // extract master element, and its needed specifics
            MasterElement* meSCS =
                MasterElementRepo::get_surface_master_element(
                    elementBucket.topology());

            const label numScsIp = meSCS->numIntPoints_;
            const label* lrscv = meSCS->adjacentNodes();

            // Element Loop
            for (stk::mesh::Bucket::size_type iElement = 0;
                 iElement < nElemPerBucket;
                 iElement++)
            {
                stk::mesh::Entity const* nodeRels =
                    elementBucket.begin_nodes(iElement);

                // FCT Fields
                const scalar* A = stk::mesh::field_data(
                    *ASTKFieldPtr_, elementBucket, iElement);
                const scalar* lambda = stk::mesh::field_data(
                    *lambdaSTKFieldPtr_, elementBucket, iElement);

                for (label ip = 0; ip < numScsIp; ++ip)
                {
                    // left and right nodes for this ip
                    const label il = lrscv[2 * ip];
                    const label ir = lrscv[2 * ip + 1];

                    stk::mesh::Entity nodeL = nodeRels[il];
                    stk::mesh::Entity nodeR = nodeRels[ir];

                    // pointer to fields to assemble
                    scalar* sumAPlusL =
                        stk::mesh::field_data(*sumAPlusSTKFieldPtr_, nodeL);
                    scalar* sumAPlusR =
                        stk::mesh::field_data(*sumAPlusSTKFieldPtr_, nodeR);
                    scalar* sumAMinusL =
                        stk::mesh::field_data(*sumAMinusSTKFieldPtr_, nodeL);
                    scalar* sumAMinusR =
                        stk::mesh::field_data(*sumAMinusSTKFieldPtr_, nodeR);

                    // Sum Limited Antidiffusion flux
                    if (A[ip] > 0.0)
                    {
                        sumAMinusL[0] += lambda[ip] * A[ip]; // outflow from L
                        sumAPlusR[0] += lambda[ip] * A[ip];  // inflow into R
                    }
                    else
                    {
                        sumAPlusL[0] -= lambda[ip] * A[ip];  // outflow from R
                        sumAMinusR[0] -= lambda[ip] * A[ip]; // inflow into L
                    }
                }
            }
        }
    }

    // Boundary
    for (label iBoundary = 0; iBoundary < domain->zonePtr()->nBoundaries();
         iBoundary++)
    {
        const auto* boundary = domain->zonePtr()->boundaryPtr(iBoundary);

        boundaryPhysicalType type = boundary->type();

        switch (type)
        {
            case boundaryPhysicalType::inlet:
            case boundaryPhysicalType::outlet:
            case boundaryPhysicalType::opening:
                {
                    // define some common selectors
                    stk::mesh::Selector selAllSides =
                        metaData.universal_part() &
                        stk::mesh::selectUnion(boundary->parts());

                    stk::mesh::BucketVector const& sideBuckets =
                        bulkData.get_buckets(metaData.side_rank(), selAllSides);

                    for (stk::mesh::BucketVector::const_iterator ib =
                             sideBuckets.begin();
                         ib != sideBuckets.end();
                         ++ib)
                    {
                        stk::mesh::Bucket& sideBucket = **ib;

                        const stk::mesh::Bucket::size_type nSidesPerBucket =
                            sideBucket.size();

                        const stk::topology theTopo = sideBucket.topology();

                        // face master element
                        MasterElement* meFC =
                            MasterElementRepo::get_surface_master_element(
                                theTopo);

                        const label numScsBip = meFC->numIntPoints_;

                        for (stk::mesh::Bucket::size_type iSide = 0;
                             iSide < nSidesPerBucket;
                             ++iSide)
                        {
                            const auto& side = sideBucket[iSide];

                            const label* ipNodeMap = meFC->ipNodeMap();

                            // Get node entities for this side
                            stk::mesh::Entity const* sideRels =
                                bulkData.begin_nodes(side);

                            // FCT side Fields
                            const scalar* sideLambda = stk::mesh::field_data(
                                *sideLambdaSTKFieldPtr_, sideBucket, iSide);
                            const scalar* sideA = stk::mesh::field_data(
                                *sideASTKFieldPtr_, sideBucket, iSide);

                            for (label ip = 0; ip < numScsBip; ++ip)
                            {
                                const label nearestNode = ipNodeMap[ip];
                                stk::mesh::Entity node = sideRels[nearestNode];

                                scalar* sideSumAMinusL = stk::mesh::field_data(
                                    *sumAMinusSTKFieldPtr_, node);

                                scalar* sideSumAPlusL = stk::mesh::field_data(
                                    *sumAPlusSTKFieldPtr_, node);

                                // Sum Limited Antidiffusive Flux A
                                if (sideA[ip] > 0.0)
                                {
                                    sideSumAMinusL[0] +=
                                        sideLambda[ip] *
                                        sideA[ip]; // outflow from L
                                }
                                else
                                {
                                    sideSumAPlusL[0] -=
                                        sideLambda[ip] *
                                        sideA[ip]; // inflow into L
                                }
                            }
                        }
                    }
                }
                break;

            default:
                break;
        }
    }
}

void freeSurfaceFlowModel::computeLambdaNode_(
    const std::shared_ptr<domain> domain,
    label iPhase)
{
    // Get Mesh Data
    auto& mesh = this->meshRef();
    stk::mesh::BulkData& bulkData = mesh.bulkDataRef();
    stk::mesh::MetaData& metaData = mesh.metaDataRef();

    // Interior
    {
        // get interior parts the domain is defined on
        const stk::mesh::PartVector& partVec =
            domain->zonePtr()->interiorParts();

        // define some common selectors; select owned nodes
        stk::mesh::Selector selOwnedNodes =
            metaData.locally_owned_part() & stk::mesh::selectUnion(partVec);

        stk::mesh::BucketVector const& nodeBuckets =
            bulkData.get_buckets(stk::topology::NODE_RANK, selOwnedNodes);
        for (stk::mesh::BucketVector::const_iterator ib = nodeBuckets.begin();
             ib != nodeBuckets.end();
             ++ib)
        {
            stk::mesh::Bucket& nodeBucket = **ib;

            const stk::mesh::Bucket::size_type nNodesPerBucket =
                nodeBucket.size();

            // Node Loop
            for (stk::mesh::Bucket::size_type iNode = 0;
                 iNode < nNodesPerBucket;
                 ++iNode)
            {
                // get node
                stk::mesh::Entity node = nodeBucket[iNode];

                const scalar* Pplus =
                    stk::mesh::field_data(*PplusSTKFieldPtr_, node);
                const scalar* Pminus =
                    stk::mesh::field_data(*PminusSTKFieldPtr_, node);
                const scalar* Qplus =
                    stk::mesh::field_data(*QplusSTKFieldPtr_, node);
                const scalar* Qminus =
                    stk::mesh::field_data(*QminusSTKFieldPtr_, node);

                scalar* lambdaP =
                    stk::mesh::field_data(*cMULESLimiterPlusSTKFieldPtr_, node);
                scalar* lambdaM = stk::mesh::field_data(
                    *cMULESLimiterMinusSTKFieldPtr_, node);
                if (Pplus[0] == 0.0)
                {
                    lambdaP[0] = 1.0;
                }
                else
                {
                    lambdaP[0] = (Qplus[0] / (Pplus[0] + SMALL));
                }

                if (Pminus[0] == 0.0)
                {
                    lambdaM[0] = 1.0;
                }
                else
                {
                    lambdaM[0] = (Qminus[0] / (Pminus[0] + SMALL));
                }
            }
        }
    }
}

scalar freeSurfaceFlowModel::computeLambdaNodeWithBounds_(
    const std::shared_ptr<domain> domain,
    label iPhase,
    label iter)
{
    // Get Mesh Data
    auto& mesh = this->meshRef();
    stk::mesh::BulkData& bulkData = mesh.bulkDataRef();
    stk::mesh::MetaData& metaData = mesh.metaDataRef();

    scalar maxLambdaChange = 0.0;

    // Interior
    {
        // get interior parts the domain is defined on
        const stk::mesh::PartVector& partVec =
            domain->zonePtr()->interiorParts();

        // define some common selectors; select owned nodes
        stk::mesh::Selector selOwnedNodes =
            metaData.locally_owned_part() & stk::mesh::selectUnion(partVec);

        stk::mesh::BucketVector const& nodeBuckets =
            bulkData.get_buckets(stk::topology::NODE_RANK, selOwnedNodes);
        for (stk::mesh::BucketVector::const_iterator ib = nodeBuckets.begin();
             ib != nodeBuckets.end();
             ++ib)
        {
            stk::mesh::Bucket& nodeBucket = **ib;

            const stk::mesh::Bucket::size_type nNodesPerBucket =
                nodeBucket.size();

            // Node Loop
            for (stk::mesh::Bucket::size_type iNode = 0;
                 iNode < nNodesPerBucket;
                 ++iNode)
            {
                // get node
                stk::mesh::Entity node = nodeBucket[iNode];

                // P+ and P- are the UNLIMITED flux sums (computed once)
                const scalar* Pplus =
                    stk::mesh::field_data(*PplusSTKFieldPtr_, node);
                const scalar* Pminus =
                    stk::mesh::field_data(*PminusSTKFieldPtr_, node);

                // Q+ and Q- are the maximum allowable changes
                const scalar* Qplus =
                    stk::mesh::field_data(*QplusSTKFieldPtr_, node);
                const scalar* Qminus =
                    stk::mesh::field_data(*QminusSTKFieldPtr_, node);

                // sumA+ and sumA- are the LIMITED flux sums (from previous
                // iteration)
                const scalar* sumAPlus =
                    stk::mesh::field_data(*sumAPlusSTKFieldPtr_, node);
                const scalar* sumAMinus =
                    stk::mesh::field_data(*sumAMinusSTKFieldPtr_, node);

                scalar* lambdaP =
                    stk::mesh::field_data(*cMULESLimiterPlusSTKFieldPtr_, node);
                scalar* lambdaM = stk::mesh::field_data(
                    *cMULESLimiterMinusSTKFieldPtr_, node);

                // Store old lambda values for convergence check
                scalar lambdaPOld = lambdaP[0];
                scalar lambdaMOld = lambdaM[0];

                // Compute new lambda values with bounds [0,1]
                // Following MULES algorithm:
                // - Iteration 0: lambda = Q / P
                // - Iteration > 0: lambda = (sumA + Q) / P
                // where sumA is the sum of already-limited fluxes

                if (iter == 0)
                {
                    // First iteration: use simple Q/P ratio
                    if (Pplus[0] <= SMALL)
                    {
                        lambdaP[0] = 1.0;
                    }
                    else
                    {
                        lambdaP[0] = std::min(
                            1.0, std::max(0.0, Qplus[0] / (Pplus[0] + SMALL)));
                    }

                    if (Pminus[0] <= SMALL)
                    {
                        lambdaM[0] = 1.0;
                    }
                    else
                    {
                        lambdaM[0] = std::min(
                            1.0,
                            std::max(0.0, Qminus[0] / (Pminus[0] + SMALL)));
                    }
                }
                else
                {
                    // Subsequent iterations: account for already-limited fluxes
                    // Key insight: limited OUTFLOWS free up capacity for
                    // INFLOWS and limited INFLOWS free up capacity for OUTFLOWS
                    if (Pplus[0] <= SMALL)
                    {
                        lambdaP[0] = 1.0;
                    }
                    else
                    {
                        // Limited OUTFLOWS (sumAMinus) free up capacity for
                        // INFLOWS
                        lambdaP[0] =
                            std::min(1.0,
                                     std::max(0.0,
                                              (sumAMinus[0] + Qplus[0]) /
                                                  (Pplus[0] + SMALL)));
                    }

                    if (Pminus[0] <= SMALL)
                    {
                        lambdaM[0] = 1.0;
                    }
                    else
                    {
                        // Limited INFLOWS (sumAPlus) free up capacity for
                        // OUTFLOWS
                        lambdaM[0] =
                            std::min(1.0,
                                     std::max(0.0,
                                              (sumAPlus[0] + Qminus[0]) /
                                                  (Pminus[0] + SMALL)));
                    }
                }

                // Track maximum change for convergence
                scalar changePlus = std::abs(lambdaP[0] - lambdaPOld);
                scalar changeMinus = std::abs(lambdaM[0] - lambdaMOld);
                maxLambdaChange = std::max(maxLambdaChange,
                                           std::max(changePlus, changeMinus));
            }
        }
    }

    return maxLambdaChange;
}

void freeSurfaceFlowModel::computeLambdaIP_(
    const std::shared_ptr<domain> domain,
    label iPhase)
{
    // Get Mesh Data
    auto& mesh = this->meshRef();
    stk::mesh::BulkData& bulkData = mesh.bulkDataRef();
    stk::mesh::MetaData& metaData = mesh.metaDataRef();

    // Interior
    {
        // get interior parts the domain is defined on
        const stk::mesh::PartVector& partVec =
            domain->zonePtr()->interiorParts();

        // Define Scratch Spaces
        std::vector<scalar> ws_cMulesLimiterPlus;
        std::vector<scalar> ws_cMulesLimiterMinus;

        // define some common selectors; select universal elements
        stk::mesh::Selector selUniversalElements =
            metaData.universal_part() & stk::mesh::selectUnion(partVec);

        stk::mesh::BucketVector const& elementBuckets = bulkData.get_buckets(
            stk::topology::ELEMENT_RANK, selUniversalElements);

        for (stk::mesh::BucketVector::const_iterator ib =
                 elementBuckets.begin();
             ib < elementBuckets.end();
             ib++)
        {
            stk::mesh::Bucket& elementBucket = **ib;
            const stk::mesh::Bucket::size_type nElemPerBucket =
                elementBucket.size();

            // extract master element, and its needed specifics
            MasterElement* meSCS =
                MasterElementRepo::get_surface_master_element(
                    elementBucket.topology());

            const label nodesPerElement = meSCS->nodesPerElement_;
            const label numScsIp = meSCS->numIntPoints_;
            const label* lrscv = meSCS->adjacentNodes();

            // allocate spaces
            ws_cMulesLimiterPlus.resize(nodesPerElement);
            ws_cMulesLimiterMinus.resize(nodesPerElement);

            // get pointers
            scalar* p_lambdaP = &ws_cMulesLimiterPlus[0];
            scalar* p_lambdaM = &ws_cMulesLimiterMinus[0];

            // Element Loop
            for (stk::mesh::Bucket::size_type iElement = 0;
                 iElement < nElemPerBucket;
                 iElement++)
            {
                stk::mesh::Entity const* nodeRels =
                    elementBucket.begin_nodes(iElement);

                for (label iNode = 0; iNode < nodesPerElement; iNode++)
                {
                    stk::mesh::Entity node = nodeRels[iNode];

                    p_lambdaP[iNode] = *stk::mesh::field_data(
                        *cMULESLimiterPlusSTKFieldPtr_, node);

                    p_lambdaM[iNode] = *stk::mesh::field_data(
                        *cMULESLimiterMinusSTKFieldPtr_, node);
                }

                // FCT Fields
                const scalar* A = stk::mesh::field_data(
                    *ASTKFieldPtr_, elementBucket, iElement);

                scalar* lambdaIP = stk::mesh::field_data(
                    *lambdaSTKFieldPtr_, elementBucket, iElement);

                for (label ip = 0; ip < numScsIp; ++ip)
                {
                    // left and right nodes for this ip
                    const label il = lrscv[2 * ip];
                    const label ir = lrscv[2 * ip + 1];
                    if (A[ip] > 0.0)
                    {
                        // Flux from L to R: L loses (outflow), R gains (inflow)
                        lambdaIP[ip] = std::min(p_lambdaM[il], p_lambdaP[ir]);
                    }
                    else
                    {
                        // Flux from R to L: R loses (outflow), L gains (inflow)
                        lambdaIP[ip] = std::min(p_lambdaP[il], p_lambdaM[ir]);
                    }
                }
            }
        }
    }

    // Boundary - Only limit OUTLET faces
    // For inlet faces, keep lambda = 1.0 (no limiting on incoming flux)
    for (label iBoundary = 0; iBoundary < domain->zonePtr()->nBoundaries();
         iBoundary++)
    {
        const auto* boundary = domain->zonePtr()->boundaryPtr(iBoundary);

        boundaryPhysicalType type = boundary->type();

        switch (type)
        {
            case boundaryPhysicalType::inlet:
                {
                    // For inlet boundaries: keep lambda = 1.0 (no limiting)
                    // The bounded flux sideFL already handles the inlet
                    // correctly define some common selectors
                    stk::mesh::Selector selAllSides =
                        metaData.universal_part() &
                        stk::mesh::selectUnion(boundary->parts());

                    stk::mesh::BucketVector const& sideBuckets =
                        bulkData.get_buckets(metaData.side_rank(), selAllSides);

                    for (stk::mesh::BucketVector::const_iterator ib =
                             sideBuckets.begin();
                         ib != sideBuckets.end();
                         ++ib)
                    {
                        stk::mesh::Bucket& sideBucket = **ib;
                        const stk::mesh::Bucket::size_type nSidesPerBucket =
                            sideBucket.size();

                        const stk::topology theTopo = sideBucket.topology();
                        MasterElement* meFC =
                            MasterElementRepo::get_surface_master_element(
                                theTopo);
                        const label numScsBip = meFC->numIntPoints_;

                        for (stk::mesh::Bucket::size_type iSide = 0;
                             iSide < nSidesPerBucket;
                             ++iSide)
                        {
                            scalar* sideLambdaIP = stk::mesh::field_data(
                                *sideLambdaSTKFieldPtr_, sideBucket, iSide);

                            // Keep lambda = 1.0 for all inlet IPs
                            for (label ip = 0; ip < numScsBip; ++ip)
                            {
                                sideLambdaIP[ip] = 1.0;
                            }
                        }
                    }
                }
                break;

            case boundaryPhysicalType::outlet:
            case boundaryPhysicalType::opening:
                {
                    std::vector<scalar> ws_lambdaP;
                    std::vector<scalar> ws_lambdaM;

                    // define some common selectors
                    stk::mesh::Selector selAllSides =
                        metaData.universal_part() &
                        stk::mesh::selectUnion(boundary->parts());

                    stk::mesh::BucketVector const& sideBuckets =
                        bulkData.get_buckets(metaData.side_rank(), selAllSides);

                    for (stk::mesh::BucketVector::const_iterator ib =
                             sideBuckets.begin();
                         ib != sideBuckets.end();
                         ++ib)
                    {
                        stk::mesh::Bucket& sideBucket = **ib;

                        const stk::mesh::Bucket::size_type nSidesPerBucket =
                            sideBucket.size();

                        const stk::topology theTopo = sideBucket.topology();

                        // face master element
                        MasterElement* meFC =
                            MasterElementRepo::get_surface_master_element(
                                theTopo);

                        const label numScsBip = meFC->numIntPoints_;

                        // extract master element specifics
                        const label numNodesPerSide = meFC->nodesPerElement_;

                        // set sizes
                        ws_lambdaP.resize(numNodesPerSide);
                        ws_lambdaM.resize(numNodesPerSide);

                        // get pointers
                        scalar* p_lambdaP = &ws_lambdaP[0];
                        scalar* p_lambdaM = &ws_lambdaM[0];

                        for (stk::mesh::Bucket::size_type iSide = 0;
                             iSide < nSidesPerBucket;
                             ++iSide)
                        {
                            const auto& side = sideBucket[iSide];

                            // mapping from ip to nodes for this ordinal
                            const label* ipNodeMap = meFC->ipNodeMap();

                            stk::mesh::Entity const* sideRels =
                                bulkData.begin_nodes(side);

                            // FCT Side Fields
                            const scalar* sideA = stk::mesh::field_data(
                                *sideASTKFieldPtr_, sideBucket, iSide);
                            const scalar* sideFL = stk::mesh::field_data(
                                *sideFLSTKFieldPtr_, sideBucket, iSide);

                            scalar* sideLambdaIP = stk::mesh::field_data(
                                *sideLambdaSTKFieldPtr_, sideBucket, iSide);

                            // fill with nodal values
                            for (label iNode = 0; iNode < numNodesPerSide;
                                 iNode++)
                            {
                                stk::mesh::Entity node = sideRels[iNode];

                                p_lambdaP[iNode] = *stk::mesh::field_data(
                                    *cMULESLimiterPlusSTKFieldPtr_, node);

                                p_lambdaM[iNode] = *stk::mesh::field_data(
                                    *cMULESLimiterMinusSTKFieldPtr_, node);
                            }

                            for (label ip = 0; ip < numScsBip; ++ip)
                            {
                                const label nearestNode = ipNodeMap[ip];

                                // Only limit OUTLET faces (total flux > 0)
                                const scalar totalFlux = sideFL[ip] + sideA[ip];

                                if (totalFlux > SMALL)
                                {
                                    // This is an outlet face - apply limiting
                                    // cMULES algorithm:
                                    // - If phiCorr > 0: use lambdap (outflow
                                    // limiter = lambdaM)
                                    // - If phiCorr < 0: use lambdam (inflow
                                    // limiter = lambdaP)
                                    if (sideA[ip] > 0.0)
                                    {
                                        // Correction adds to outflow - use
                                        // outflow limiter
                                        sideLambdaIP[ip] =
                                            std::min(sideLambdaIP[ip],
                                                     p_lambdaM[nearestNode]);
                                    }
                                    else
                                    {
                                        // Correction opposes outflow - use
                                        // inflow limiter
                                        sideLambdaIP[ip] =
                                            std::min(sideLambdaIP[ip],
                                                     p_lambdaP[nearestNode]);
                                    }
                                }
                                // else: inlet flow at outlet boundary - keep
                                // lambda = 1.0
                            }
                        }
                    }
                }
                break;

            default:
                break;
        }
    }
}

void freeSurfaceFlowModel::updateA_(const std::shared_ptr<domain> domain,
                                    label iPhase)
{
    // Get Mesh Data
    auto& mesh = this->meshRef();
    stk::mesh::BulkData& bulkData = mesh.bulkDataRef();
    stk::mesh::MetaData& metaData = mesh.metaDataRef();

    // Interior
    {
        // get interior parts the domain is defined on
        const stk::mesh::PartVector& partVec =
            domain->zonePtr()->interiorParts();
        // define some common selectors; select owned nodes
        stk::mesh::Selector selUniversalElements =
            metaData.universal_part() & stk::mesh::selectUnion(partVec);

        stk::mesh::BucketVector const& elementBuckets = bulkData.get_buckets(
            stk::topology::ELEMENT_RANK, selUniversalElements);

        for (stk::mesh::BucketVector::const_iterator ib =
                 elementBuckets.begin();
             ib < elementBuckets.end();
             ib++)
        {
            stk::mesh::Bucket& elementBucket = **ib;
            const stk::mesh::Bucket::size_type nElemPerBucket =
                elementBucket.size();

            // extract master element, and its needed specifics
            MasterElement* meSCS =
                MasterElementRepo::get_surface_master_element(
                    elementBucket.topology());

            const label numScsIp = meSCS->numIntPoints_;

            // Element Loop
            for (stk::mesh::Bucket::size_type iElement = 0;
                 iElement < nElemPerBucket;
                 iElement++)
            {
                // FCT Fields
                scalar* A = stk::mesh::field_data(
                    *ASTKFieldPtr_, elementBucket, iElement);
                const scalar* lambda = stk::mesh::field_data(
                    *lambdaSTKFieldPtr_, elementBucket, iElement);

                for (label ip = 0; ip < numScsIp; ++ip)
                {
                    // Limit Antidiffusion flux A
                    A[ip] *= lambda[ip];
                }
            }
        }
    }

    // Boundary
    for (label iBoundary = 0; iBoundary < domain->zonePtr()->nBoundaries();
         iBoundary++)
    {
        const auto* boundary = domain->zonePtr()->boundaryPtr(iBoundary);

        boundaryPhysicalType type = boundary->type();

        switch (type)
        {
            case boundaryPhysicalType::inlet:
            case boundaryPhysicalType::outlet:
            case boundaryPhysicalType::opening:
                {
                    // define some common selectors
                    stk::mesh::Selector selAllSides =
                        metaData.universal_part() &
                        stk::mesh::selectUnion(boundary->parts());

                    stk::mesh::BucketVector const& sideBuckets =
                        bulkData.get_buckets(metaData.side_rank(), selAllSides);

                    for (stk::mesh::BucketVector::const_iterator ib =
                             sideBuckets.begin();
                         ib != sideBuckets.end();
                         ++ib)
                    {
                        stk::mesh::Bucket& sideBucket = **ib;

                        const stk::mesh::Bucket::size_type nSidesPerBucket =
                            sideBucket.size();

                        const stk::topology theTopo = sideBucket.topology();

                        // face master element
                        MasterElement* meFC =
                            MasterElementRepo::get_surface_master_element(
                                theTopo);

                        const label numScsBip = meFC->numIntPoints_;

                        for (stk::mesh::Bucket::size_type iSide = 0;
                             iSide < nSidesPerBucket;
                             ++iSide)
                        {
                            const auto& side = sideBucket[iSide];

                            // FCT side Fields
                            const scalar* sideLambda = stk::mesh::field_data(
                                *sideLambdaSTKFieldPtr_, sideBucket, iSide);
                            scalar* sideA = stk::mesh::field_data(
                                *sideASTKFieldPtr_, sideBucket, iSide);

                            for (label ip = 0; ip < numScsBip; ++ip)
                            {
                                // Calculate Antidiffusive Flux A
                                sideA[ip] *= sideLambda[ip];
                            }
                        }
                    }
                }
                break;

            default:
                break;
        }
    }
}

void freeSurfaceFlowModel::updateAlpha_(const std::shared_ptr<domain> domain,
                                        label iPhase,
                                        label iter)
{
    // Get Mesh Data
    auto& mesh = this->meshRef();
    stk::mesh::BulkData& bulkData = mesh.bulkDataRef();
    stk::mesh::MetaData& metaData = mesh.metaDataRef();

    // Interior
    {
        // Get fields
        STKScalarField* alphaSTKFieldPtr = this->alphaRef(iPhase).stkFieldPtr();

        const STKScalarField* rhoSTKFieldPtr =
            this->rhoRef(iPhase).stkFieldPtr();

        // Get geometric fields
        const auto* volSTKFieldPtr = metaData.get_field<scalar>(
            stk::topology::NODE_RANK, mesh::dual_nodal_volume_ID);

        // Get user defined fields
        const bool is_transient = this->controlsRef().isTransient();
        const scalar dt = is_transient
                              ? this->controlsRef().getTimestep()
                              : this->controlsRef().getPhysicalTimescale();

        // get interior parts the domain is defined on
        const stk::mesh::PartVector& partVec =
            domain->zonePtr()->interiorParts();

        // Define Scratch Spaces
        std::vector<scalar> ws_vol;
        std::vector<scalar> ws_rho;

        // define some common selectors; select owned nodes
        stk::mesh::Selector selUniversalElements =
            metaData.universal_part() & stk::mesh::selectUnion(partVec);

        stk::mesh::BucketVector const& elementBuckets = bulkData.get_buckets(
            stk::topology::ELEMENT_RANK, selUniversalElements);

        for (stk::mesh::BucketVector::const_iterator ib =
                 elementBuckets.begin();
             ib < elementBuckets.end();
             ib++)
        {
            stk::mesh::Bucket& elementBucket = **ib;
            const stk::mesh::Bucket::size_type nElemPerBucket =
                elementBucket.size();

            // extract master element, and its needed specifics
            MasterElement* meSCS =
                MasterElementRepo::get_surface_master_element(
                    elementBucket.topology());

            const label nodesPerElement = meSCS->nodesPerElement_;
            const label numScsIp = meSCS->numIntPoints_;
            const label* lrscv = meSCS->adjacentNodes();

            // allocate space for scratch spaces
            ws_vol.resize(nodesPerElement);
            ws_rho.resize(nodesPerElement);

            // get pointers
            scalar* p_vol = &ws_vol[0];
            scalar* p_rho = &ws_rho[0];

            // Element Loop
            for (stk::mesh::Bucket::size_type iElement = 0;
                 iElement < nElemPerBucket;
                 iElement++)
            {
                stk::mesh::Entity const* nodeRels =
                    elementBucket.begin_nodes(iElement);

                // FCT Fields
                const scalar* A = stk::mesh::field_data(
                    *ASTKFieldPtr_, elementBucket, iElement);

                for (label iNode = 0; iNode < nodesPerElement; iNode++)
                {
                    stk::mesh::Entity node = nodeRels[iNode];

                    p_vol[iNode] =
                        *stk::mesh::field_data(*volSTKFieldPtr, node);

                    p_rho[iNode] =
                        *stk::mesh::field_data(*rhoSTKFieldPtr, node);
                }

                // Get lambda field for this element
                const scalar* lambda = stk::mesh::field_data(
                    *lambdaSTKFieldPtr_, elementBucket, iElement);

                for (label ip = 0; ip < numScsIp; ++ip)
                {
                    // left and right nodes for this ip
                    const label il = lrscv[2 * ip];
                    const label ir = lrscv[2 * ip + 1];

                    stk::mesh::Entity nodeL = nodeRels[il];
                    stk::mesh::Entity nodeR = nodeRels[ir];

                    scalar* alphaL =
                        stk::mesh::field_data(*alphaSTKFieldPtr, nodeL);
                    scalar* alphaR =
                        stk::mesh::field_data(*alphaSTKFieldPtr, nodeR);

                    // Apply LIMITED antidiffusive flux correction
                    // No under-relaxation (kappa=1.0) - the limiter ensures
                    // boundedness Formula: alpha_new = alpha - dt/(rho*V) *
                    // sum(lambda*A)
                    const scalar limitedA = lambda[ip] * A[ip];

                    alphaL[0] -=
                        (dt / (p_rho[il] * p_vol[il] + SMALL)) * limitedA;
                    alphaR[0] +=
                        (dt / (p_rho[ir] * p_vol[ir] + SMALL)) * limitedA;
                }
            }
        }
    }

    // Boundary
    for (label iBoundary = 0; iBoundary < domain->zonePtr()->nBoundaries();
         iBoundary++)
    {
        const auto* boundary = domain->zonePtr()->boundaryPtr(iBoundary);

        boundaryPhysicalType type = boundary->type();

        boundaryConditionType bcType =
            this->alphaRef(iPhase)
                .boundaryConditionRef(domain->index(), iBoundary)
                .type();

        switch (type)
        {
            case boundaryPhysicalType::inlet:
                {
                    switch (bcType)
                    {
                        case boundaryConditionType::specifiedValue:
                            {
                                // <-- Specified Value Do not Update Alpha -->
                                errorMsg("Must not reach here");
                            }
                            break;

                        default:
                            errorMsg("Must not reach here");
                    }
                }
                break;

            case boundaryPhysicalType::outlet:
                {
                    switch (bcType)
                    {
                        case boundaryConditionType::zeroGradient:
                            {
                                // Get fields
                                STKScalarField* alphaSTKFieldPtr =
                                    this->alphaRef(iPhase).stkFieldPtr();

                                const STKScalarField* rhoSTKFieldPtr =
                                    this->rhoRef(iPhase).stkFieldPtr();

                                const auto* volSTKFieldPtr =
                                    metaData.get_field<scalar>(
                                        stk::topology::NODE_RANK,
                                        mesh::dual_nodal_volume_ID);

                                // Get user defined fields
                                const bool is_transient =
                                    this->controlsRef().isTransient();
                                const scalar dt =
                                    is_transient
                                        ? this->controlsRef().getTimestep()
                                        : this->controlsRef()
                                              .getPhysicalTimescale();

                                std::vector<scalar> ws_vol;
                                std::vector<scalar> ws_rho;

                                // define some common selectors
                                stk::mesh::Selector selAllSides =
                                    metaData.universal_part() &
                                    stk::mesh::selectUnion(boundary->parts());

                                stk::mesh::BucketVector const& sideBuckets =
                                    bulkData.get_buckets(metaData.side_rank(),
                                                         selAllSides);

                                for (stk::mesh::BucketVector::const_iterator
                                         ib = sideBuckets.begin();
                                     ib != sideBuckets.end();
                                     ++ib)
                                {
                                    stk::mesh::Bucket& sideBucket = **ib;

                                    const stk::mesh::Bucket::size_type
                                        nSidesPerBucket = sideBucket.size();

                                    const stk::topology theTopo =
                                        sideBucket.topology();

                                    // face master element
                                    MasterElement* meFC = MasterElementRepo::
                                        get_surface_master_element(theTopo);

                                    const label numScsBip = meFC->numIntPoints_;

                                    // extract master element specifics
                                    const label numNodesPerSide =
                                        meFC->nodesPerElement_;

                                    // set sizes
                                    ws_vol.resize(numNodesPerSide);
                                    ws_rho.resize(numNodesPerSide);

                                    // get pointers
                                    scalar* p_vol = &ws_vol[0];
                                    scalar* p_rho = &ws_rho[0];

                                    for (stk::mesh::Bucket::size_type iSide = 0;
                                         iSide < nSidesPerBucket;
                                         ++iSide)
                                    {
                                        const auto& side = sideBucket[iSide];

                                        // mapping from ip to nodes for this
                                        // ordinal
                                        const label* ipNodeMap =
                                            meFC->ipNodeMap();

                                        stk::mesh::Entity const* sideRels =
                                            bulkData.begin_nodes(side);

                                        // FCT Side Fields
                                        const scalar* sideA =
                                            stk::mesh::field_data(
                                                *sideASTKFieldPtr_,
                                                sideBucket,
                                                iSide);

                                        // fill with nodal values
                                        for (label iNode = 0;
                                             iNode < numNodesPerSide;
                                             iNode++)
                                        {
                                            stk::mesh::Entity node =
                                                sideRels[iNode];

                                            p_vol[iNode] =
                                                *stk::mesh::field_data(
                                                    *volSTKFieldPtr, node);

                                            p_rho[iNode] =
                                                *stk::mesh::field_data(
                                                    *rhoSTKFieldPtr, node);
                                        }

                                        // Get side lambda field
                                        const scalar* sideLambda =
                                            stk::mesh::field_data(
                                                *sideLambdaSTKFieldPtr_,
                                                sideBucket,
                                                iSide);

                                        for (label ip = 0; ip < numScsBip; ++ip)
                                        {
                                            const label nearestNode =
                                                ipNodeMap[ip];
                                            stk::mesh::Entity node =
                                                sideRels[nearestNode];

                                            scalar* sideAlphaL =
                                                stk::mesh::field_data(
                                                    *alphaSTKFieldPtr, node);

                                            // Apply LIMITED antidiffusive flux
                                            // No under-relaxation (kappa=1.0)
                                            const scalar limitedSideA =
                                                sideLambda[ip] * sideA[ip];

                                            sideAlphaL[0] -=
                                                (dt / (p_rho[nearestNode] *
                                                           p_vol[nearestNode] +
                                                       SMALL)) *
                                                limitedSideA;
                                        }
                                    }
                                }
                            }
                            break;

                        default:
                            errorMsg("Must not reach here");
                    }
                }
                break;

            case boundaryPhysicalType::opening:
                {
                    // Get fields
                    STKScalarField* alphaSTKFieldPtr =
                        this->alphaRef(iPhase).stkFieldPtr();

                    const STKScalarField* rhoSTKFieldPtr =
                        this->rhoRef(iPhase).stkFieldPtr();

                    const STKScalarField* mDotSideSTKFieldPtr =
                        this->mDotRef(iPhase).sideFieldRef().stkFieldPtr();

                    const auto* volSTKFieldPtr = metaData.get_field<scalar>(
                        stk::topology::NODE_RANK, mesh::dual_nodal_volume_ID);

                    // Get user defined fields
                    const bool is_transient = this->controlsRef().isTransient();
                    const scalar dt =
                        is_transient
                            ? this->controlsRef().getTimestep()
                            : this->controlsRef().getPhysicalTimescale();

                    std::vector<scalar> ws_vol;
                    std::vector<scalar> ws_rho;

                    // define some common selectors
                    stk::mesh::Selector selAllSides =
                        metaData.universal_part() &
                        stk::mesh::selectUnion(boundary->parts());

                    stk::mesh::BucketVector const& sideBuckets =
                        bulkData.get_buckets(metaData.side_rank(), selAllSides);

                    for (stk::mesh::BucketVector::const_iterator ib =
                             sideBuckets.begin();
                         ib != sideBuckets.end();
                         ++ib)
                    {
                        stk::mesh::Bucket& sideBucket = **ib;

                        const stk::mesh::Bucket::size_type nSidesPerBucket =
                            sideBucket.size();

                        const stk::topology theTopo = sideBucket.topology();

                        // face master element
                        MasterElement* meFC =
                            MasterElementRepo::get_surface_master_element(
                                theTopo);

                        const label numScsBip = meFC->numIntPoints_;

                        // extract master element specifics
                        const label numNodesPerSide = meFC->nodesPerElement_;

                        // set sizes
                        ws_vol.resize(numNodesPerSide);
                        ws_rho.resize(numNodesPerSide);

                        // get pointers
                        scalar* p_vol = &ws_vol[0];
                        scalar* p_rho = &ws_rho[0];

                        for (stk::mesh::Bucket::size_type iSide = 0;
                             iSide < nSidesPerBucket;
                             ++iSide)
                        {
                            const auto& side = sideBucket[iSide];

                            // mapping from ip to nodes for this ordinal
                            const label* ipNodeMap = meFC->ipNodeMap();

                            stk::mesh::Entity const* sideRels =
                                bulkData.begin_nodes(side);

                            // FCT Side Fields
                            const scalar* sideA = stk::mesh::field_data(
                                *sideASTKFieldPtr_, sideBucket, iSide);
                            const scalar* mDot = stk::mesh::field_data(
                                *mDotSideSTKFieldPtr, sideBucket, iSide);

                            // fill with nodal values
                            for (label iNode = 0; iNode < numNodesPerSide;
                                 iNode++)
                            {
                                stk::mesh::Entity node = sideRels[iNode];

                                p_vol[iNode] = *stk::mesh::field_data(
                                    *volSTKFieldPtr, node);

                                p_rho[iNode] = *stk::mesh::field_data(
                                    *rhoSTKFieldPtr, node);
                            }

                            // Get side lambda field
                            const scalar* sideLambda = stk::mesh::field_data(
                                *sideLambdaSTKFieldPtr_, sideBucket, iSide);

                            for (label ip = 0; ip < numScsBip; ++ip)
                            {
                                const label nearestNode = ipNodeMap[ip];
                                stk::mesh::Entity node = sideRels[nearestNode];

                                // save off mDot
                                const scalar tmDot = mDot[ip];

                                // Opening: only update alpha for outflow
                                // mDot > 0: flow leaving - update alpha
                                // mDot < 0: flow entering - alpha is prescribed
                                if (tmDot > 0.0)
                                {
                                    scalar* sideAlphaL = stk::mesh::field_data(
                                        *alphaSTKFieldPtr, node);

                                    // Apply LIMITED antidiffusive flux
                                    // No under-relaxation (kappa=1.0)
                                    const scalar limitedSideA =
                                        sideLambda[ip] * sideA[ip];

                                    sideAlphaL[0] -=
                                        (dt / (p_rho[nearestNode] *
                                                   p_vol[nearestNode] +
                                               SMALL)) *
                                        limitedSideA;
                                }
                                // For inflow (tmDot < 0): alpha is prescribed,
                                // do not update
                            }
                        }
                    }
                }
                break;

            default:
                break;
        }
    }
}

void freeSurfaceFlowModel::updateFL_(const std::shared_ptr<domain> domain,
                                     label iPhase,
                                     label iter)
{
    // Get Mesh Data
    auto& mesh = this->meshRef();
    stk::mesh::BulkData& bulkData = mesh.bulkDataRef();
    stk::mesh::MetaData& metaData = mesh.metaDataRef();

    // Update bounded flux by adding LIMITED correction: FL += lambda * A
    // This makes FL approach FH over iterations while maintaining boundedness

    // Interior
    {
        // get interior parts the domain is defined on
        const stk::mesh::PartVector& partVec =
            domain->zonePtr()->interiorParts();

        // define some common selectors; select owned nodes
        stk::mesh::Selector selUniversalElements =
            metaData.universal_part() & stk::mesh::selectUnion(partVec);

        stk::mesh::BucketVector const& elementBuckets = bulkData.get_buckets(
            stk::topology::ELEMENT_RANK, selUniversalElements);

        for (stk::mesh::BucketVector::const_iterator ib =
                 elementBuckets.begin();
             ib < elementBuckets.end();
             ib++)
        {
            stk::mesh::Bucket& elementBucket = **ib;
            const stk::mesh::Bucket::size_type nElemPerBucket =
                elementBucket.size();

            // extract master element, and its needed specifics
            MasterElement* meSCS =
                MasterElementRepo::get_surface_master_element(
                    elementBucket.topology());

            const label numScsIp = meSCS->numIntPoints_;

            // Element Loop
            for (stk::mesh::Bucket::size_type iElement = 0;
                 iElement < nElemPerBucket;
                 iElement++)
            {
                // FCT Fields
                scalar* FL = stk::mesh::field_data(
                    *FLSTKFieldPtr_, elementBucket, iElement);

                const scalar* A = stk::mesh::field_data(
                    *ASTKFieldPtr_, elementBucket, iElement);

                const scalar* lambda = stk::mesh::field_data(
                    *lambdaSTKFieldPtr_, elementBucket, iElement);

                // FL += lambda * A (no under-relaxation)
                for (label ip = 0; ip < numScsIp; ++ip)
                {
                    FL[ip] += lambda[ip] * A[ip];
                }
            }
        }
    }

    // Boundary
    for (label iBoundary = 0; iBoundary < domain->zonePtr()->nBoundaries();
         iBoundary++)
    {
        const auto* boundary = domain->zonePtr()->boundaryPtr(iBoundary);

        boundaryPhysicalType type = boundary->type();

        switch (type)
        {
            case boundaryPhysicalType::inlet:
            case boundaryPhysicalType::outlet:
            case boundaryPhysicalType::opening:
                {
                    // define some common selectors
                    stk::mesh::Selector selAllSides =
                        metaData.universal_part() &
                        stk::mesh::selectUnion(boundary->parts());

                    stk::mesh::BucketVector const& sideBuckets =
                        bulkData.get_buckets(metaData.side_rank(), selAllSides);

                    for (stk::mesh::BucketVector::const_iterator ib =
                             sideBuckets.begin();
                         ib != sideBuckets.end();
                         ++ib)
                    {
                        stk::mesh::Bucket& sideBucket = **ib;

                        const stk::mesh::Bucket::size_type nSidesPerBucket =
                            sideBucket.size();

                        const stk::topology theTopo = sideBucket.topology();

                        // face master element
                        MasterElement* meFC =
                            MasterElementRepo::get_surface_master_element(
                                theTopo);

                        const label numScsBip = meFC->numIntPoints_;

                        for (stk::mesh::Bucket::size_type iSide = 0;
                             iSide < nSidesPerBucket;
                             ++iSide)
                        {
                            const auto& side = sideBucket[iSide];

                            // FCT side Fields
                            scalar* sideFL = stk::mesh::field_data(
                                *sideFLSTKFieldPtr_, sideBucket, iSide);

                            const scalar* sideA = stk::mesh::field_data(
                                *sideASTKFieldPtr_, sideBucket, iSide);

                            const scalar* sideLambda = stk::mesh::field_data(
                                *sideLambdaSTKFieldPtr_, sideBucket, iSide);

                            // sideFL += sideLambda * sideA (no
                            // under-relaxation)
                            for (label ip = 0; ip < numScsBip; ++ip)
                            {
                                sideFL[ip] += sideLambda[ip] * sideA[ip];
                            }
                        }
                    }
                }
                break;

            default:
                break;
        }
    }
}

} /* namespace accel */
