// File       : freeSurfaceFlowModel.cpp
// Created    : Sun Jan 26 2025 22:06:38 (+0100)
// Author     : Mhamad Mahdi Alloush
// Description:
// Copyright 2025 CCFNUM HSLU T&A. All Rights Reserved.

#include "freeSurfaceFlowModel.h"
#include "idealGasModel.h"
#include "ipInfo.h"
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

                // balanced-force capillary potential psi (gradient-capable,
                // mirrors the pressure gradient operator/interpolation scheme)
                if (!capillaryPotentialFieldPtr_)
                {
                    capillaryPotentialFieldPtr_ =
                        std::make_unique<nodeScalarField>(&this->realmRef(),
                                                          "capillary_potential",
                                                          1,
                                                          false,
                                                          false,
                                                          true,
                                                          false);
                    const auto& basic = controlsRef()
                                            .solverRef()
                                            .solverControl_.basicSettings_
                                            .interpolationSchemeType_;
                    capillaryPotentialFieldPtr_->setInterpolationScheme(
                        basic.pressureInterpolationType_);
                    capillaryPotentialFieldPtr_->setGradientInterpolationScheme(
                        basic.pressureGradientInterpolationType_);
                    // psi is rebuilt each step: full gradient, no relaxation
                    // lag
                    capillaryPotentialFieldPtr_->setGradURF(1.0);
                }

                // sigma*kappa at nodes: the CSF force is
                // sigma*kappa*grad(alpha)
                if (!sigmaKappaSTKFieldPtr_)
                {
                    sigmaKappaSTKFieldPtr_ = &metaData.declare_field<scalar>(
                        stk::topology::NODE_RANK, "capillary_sigma_kappa");
                    stk::io::set_field_output_type(*sigmaKappaSTKFieldPtr_,
                                                   fieldType[1]);
                }
                for (const stk::mesh::Part* part : partVec)
                {
                    if (!sigmaKappaSTKFieldPtr_->defined_on(*part))
                        stk::mesh::put_field_on_mesh(
                            *sigmaKappaSTKFieldPtr_, *part, nullptr);
                }

                // scratch fields for the interface-weighted curvature extension
                if (!curvNumSTKFieldPtr_)
                {
                    curvNumSTKFieldPtr_ = &metaData.declare_field<scalar>(
                        stk::topology::NODE_RANK, "curv_num");
                    curvDenSTKFieldPtr_ = &metaData.declare_field<scalar>(
                        stk::topology::NODE_RANK, "curv_den");
                    curvAccumSTKFieldPtr_ = &metaData.declare_field<scalar>(
                        stk::topology::NODE_RANK, "curv_accum");
                    curvCntSTKFieldPtr_ = &metaData.declare_field<scalar>(
                        stk::topology::NODE_RANK, "curv_cnt");
                }
                for (const stk::mesh::Part* part : partVec)
                {
                    for (STKScalarField* fp : {curvNumSTKFieldPtr_,
                                               curvDenSTKFieldPtr_,
                                               curvAccumSTKFieldPtr_,
                                               curvCntSTKFieldPtr_})
                    {
                        if (!fp->defined_on(*part))
                            stk::mesh::put_field_on_mesh(*fp, *part, nullptr);
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

    // FOrig holds only non-ST forces (buoyancy); the CSF force enters as
    // grad(psi), added to F during redistribution (consistent with gradP).
    computeCapillaryPotential_(domain);

    // Seed F with the non-ST body forces; grad(psi) is added in redistribution
    ops::copy<scalar>(FOriginalSTKFieldPtr_,
                      FSTKFieldPtr_,
                      domain->zonePtr()->interiorParts());
}

void freeSurfaceFlowModel::computeCapillaryPotential_(
    const std::shared_ptr<domain> domain)
{
    if (!capillaryPotentialFieldPtr_)
        return;

    auto& mesh = this->meshRef();
    stk::mesh::BulkData& bulkData = mesh.bulkDataRef();
    stk::mesh::MetaData& metaData = mesh.metaDataRef();
    const stk::mesh::PartVector& partVec = domain->zonePtr()->interiorParts();

    // capillaryPotential holds alpha (its gradient = grad(alpha)); sigmaKappa
    // holds sigma*kappa. The CSF force is sigma*kappa*grad(alpha).
    STKScalarField* psiSTKFieldPtr = capillaryPotentialFieldPtr_->stkFieldPtr();
    STKScalarField* skSTKFieldPtr = sigmaKappaSTKFieldPtr_;
    ops::zero(psiSTKFieldPtr, partVec);
    ops::zero(skSTKFieldPtr, partVec);

    const bool smooth = controlsRef()
                            .solverRef()
                            .solverControl_.advancedOptions_.equationControls_
                            .volumeFractionSmoothing_.smoothVolumeFraction_;

    for (const auto& fpm : domain->fluidPairModels_)
    {
        if (fpm.surfaceTension_.option_ !=
            surfaceTensionModelOption::continuumSurfaceForce)
            continue;

        const label phaseIndex = fpm.globalIndexA_;
        const scalar sigma = fpm.surfaceTension_.coefficient_;

        const STKScalarField* alphaSTKFieldPtr =
            smooth ? this->alphaSmoothRef(phaseIndex).stkFieldPtr()
                   : this->alphaRef(phaseIndex).stkFieldPtr();

        const std::string kappaName =
            "curvature." + fpm.materialA_ + "_" + fpm.materialB_;
        auto it = kappaSTKFieldPtrs_.find(kappaName);
        STK_ThrowAssert(it != kappaSTKFieldPtrs_.end());
        STKScalarField* kappaFieldPtr = it->second;

        // curvature, then a LIGHT |grad(alpha)|-weighted band denoise (the
        // force sigma*kappa*grad(alpha) is band-localized, so no bulk extension
        // is needed; heavy smoothing would wash kappa toward uniform -> inert
        // force)
        computeCurvature_(domain, phaseIndex, kappaFieldPtr);
        extendCurvatureToBulk_(
            domain, kappaFieldPtr, this->alphaRef(phaseIndex).stkFieldPtr());

        stk::mesh::Selector selNodes =
            metaData.universal_part() & stk::mesh::selectUnion(partVec);
        stk::mesh::BucketVector const& nodeBuckets =
            bulkData.get_buckets(stk::topology::NODE_RANK, selNodes);
        for (auto ib = nodeBuckets.begin(); ib != nodeBuckets.end(); ++ib)
        {
            stk::mesh::Bucket& nodeBucket = **ib;
            const auto nNodes = nodeBucket.size();
            const scalar* alphab =
                stk::mesh::field_data(*alphaSTKFieldPtr, nodeBucket);
            const scalar* kappab =
                stk::mesh::field_data(*kappaFieldPtr, nodeBucket);
            scalar* psib = stk::mesh::field_data(*psiSTKFieldPtr, nodeBucket);
            scalar* skb = stk::mesh::field_data(*skSTKFieldPtr, nodeBucket);
            for (stk::mesh::Bucket::size_type in = 0; in < nNodes; ++in)
            {
                psib[in] = alphab[in];
                skb[in] += sigma * kappab[in];
            }
        }
    }

    // alpha is pointwise; make ghosts consistent then take grad(alpha) with the
    // pressure's gradient operator. sigma*kappa ghosts also synced.
    capillaryPotentialFieldPtr_->synchronizeGhostedEntities(domain->index());
    capillaryPotentialFieldPtr_->updateGradientField(domain->index());
    if (messager::parallel())
        stk::mesh::communicate_field_data(bulkData, {skSTKFieldPtr});
}

void freeSurfaceFlowModel::extendCurvatureToBulk_(
    const std::shared_ptr<domain> domain,
    STKScalarField* kappaFieldPtr,
    const STKScalarField* alphaFieldPtr)
{
    auto& mesh = this->meshRef();
    stk::mesh::BulkData& bulkData = mesh.bulkDataRef();
    stk::mesh::MetaData& metaData = mesh.metaDataRef();
    const stk::mesh::PartVector& partVec = domain->zonePtr()->interiorParts();

    const auto& vofSmooth = controlsRef()
                                .solverRef()
                                .solverControl_.advancedOptions_
                                .equationControls_.volumeFractionSmoothing_;
    const label nIter = vofSmooth.curvatureSmoothingIterations_;
    const bool laplacian = vofSmooth.curvatureSmootherLaplacian_;

    stk::mesh::Selector selNodes =
        metaData.universal_part() & stk::mesh::selectUnion(partVec);
    stk::mesh::BucketVector const& nodeBuckets =
        bulkData.get_buckets(stk::topology::NODE_RANK, selNodes);

    // interface weight w = |grad(alpha)| -> curvDen, then N = w*kappa ->
    // curvNum
    ops::zero(curvDenSTKFieldPtr_, partVec);
    {
        const auto* coordsPtr = metaData.get_field<scalar>(
            stk::topology::NODE_RANK, this->getCoordinatesID_(domain));
        std::vector<scalar> ws_coords, ws_a, ws_dndx, ws_deriv, ws_detj, ws_scv;
        std::vector<scalar> gIp(SPATIAL_DIM);
        stk::mesh::Selector selE =
            metaData.universal_part() & stk::mesh::selectUnion(partVec);
        stk::mesh::BucketVector const& eb =
            bulkData.get_buckets(stk::topology::ELEMENT_RANK, selE);
        for (auto ib = eb.begin(); ib != eb.end(); ++ib)
        {
            stk::mesh::Bucket& bk = **ib;
            MasterElement* meSCV =
                MasterElementRepo::get_volume_master_element(bk.topology());
            const label npe = meSCV->nodesPerElement_;
            const label nip = meSCV->numIntPoints_;
            const label* ipNodeMap = meSCV->ipNodeMap();
            ws_coords.resize(npe * SPATIAL_DIM);
            ws_a.resize(npe);
            ws_dndx.resize(SPATIAL_DIM * nip * npe);
            ws_deriv.resize(SPATIAL_DIM * nip * npe);
            ws_detj.resize(nip);
            ws_scv.resize(nip);
            for (stk::mesh::Bucket::size_type k = 0; k < bk.size(); ++k)
            {
                stk::mesh::Entity const* nodes = bk.begin_nodes(k);
                for (label ni = 0; ni < npe; ++ni)
                {
                    const scalar* c =
                        stk::mesh::field_data(*coordsPtr, nodes[ni]);
                    ws_a[ni] =
                        *stk::mesh::field_data(*alphaFieldPtr, nodes[ni]);
                    for (label d = 0; d < SPATIAL_DIM; ++d)
                        ws_coords[ni * SPATIAL_DIM + d] = c[d];
                }
                scalar err = 0.0;
                meSCV->determinant(1, &ws_coords[0], &ws_scv[0], &err);
                meSCV->grad_op(1,
                               &ws_coords[0],
                               &ws_dndx[0],
                               &ws_deriv[0],
                               &ws_detj[0],
                               &err);
                for (label ip = 0; ip < nip; ++ip)
                {
                    for (label d = 0; d < SPATIAL_DIM; ++d)
                        gIp[d] = 0.0;
                    for (label ni = 0; ni < npe; ++ni)
                    {
                        const label off =
                            SPATIAL_DIM * npe * ip + ni * SPATIAL_DIM;
                        for (label d = 0; d < SPATIAL_DIM; ++d)
                            gIp[d] += ws_dndx[off + d] * ws_a[ni];
                    }
                    scalar gmag = 0.0;
                    for (label d = 0; d < SPATIAL_DIM; ++d)
                        gmag += gIp[d] * gIp[d];
                    gmag = std::sqrt(gmag);
                    *stk::mesh::field_data(*curvDenSTKFieldPtr_,
                                           nodes[ipNodeMap[ip]]) +=
                        gmag * ws_scv[ip];
                }
            }
        }
    }
    if (messager::parallel())
        stk::mesh::communicate_field_data(bulkData, {curvDenSTKFieldPtr_});

    for (auto ib = nodeBuckets.begin(); ib != nodeBuckets.end(); ++ib)
    {
        stk::mesh::Bucket& bk = **ib;
        const auto n = bk.size();
        const scalar* kap = stk::mesh::field_data(*kappaFieldPtr, bk);
        const scalar* den = stk::mesh::field_data(*curvDenSTKFieldPtr_, bk);
        scalar* num = stk::mesh::field_data(*curvNumSTKFieldPtr_, bk);
        for (stk::mesh::Bucket::size_type in = 0; in < n; ++in)
            num[in] = den[in] * kap[in];
    }

    // normalized convolution: smooth N and D, then kappa_ext = N / D. Two
    // smoothers (curvAccum is the Laplacian rhs scratch / box-average sum).
    if (laplacian)
    {
        for (label it = 0; it < nIter; ++it)
        {
            computeSmoothRHS_(
                domain, curvNumSTKFieldPtr_, curvAccumSTKFieldPtr_);
            assembleSmoothingTerm_(
                domain, curvNumSTKFieldPtr_, curvAccumSTKFieldPtr_);
            computeSmoothRHS_(
                domain, curvDenSTKFieldPtr_, curvAccumSTKFieldPtr_);
            assembleSmoothingTerm_(
                domain, curvDenSTKFieldPtr_, curvAccumSTKFieldPtr_);
        }
    }
    else
    {
        smoothField_(domain, curvNumSTKFieldPtr_, nIter);
        smoothField_(domain, curvDenSTKFieldPtr_, nIter);
    }

    for (auto ib = nodeBuckets.begin(); ib != nodeBuckets.end(); ++ib)
    {
        stk::mesh::Bucket& bk = **ib;
        const auto n = bk.size();
        const scalar* num = stk::mesh::field_data(*curvNumSTKFieldPtr_, bk);
        const scalar* den = stk::mesh::field_data(*curvDenSTKFieldPtr_, bk);
        scalar* kap = stk::mesh::field_data(*kappaFieldPtr, bk);
        for (stk::mesh::Bucket::size_type in = 0; in < n; ++in)
            kap[in] = num[in] / (den[in] + SMALL);
    }
    if (messager::parallel())
        stk::mesh::communicate_field_data(bulkData, {kappaFieldPtr});
}

void freeSurfaceFlowModel::smoothField_(const std::shared_ptr<domain> domain,
                                        STKScalarField* fieldPtr,
                                        label nIterations)
{
    // Jacobi box-average over element neighbours (owned+aura loop is complete;
    // ghost sync per iteration refreshes the aura)
    auto& mesh = this->meshRef();
    stk::mesh::BulkData& bulkData = mesh.bulkDataRef();
    stk::mesh::MetaData& metaData = mesh.metaDataRef();
    const stk::mesh::PartVector& partVec = domain->zonePtr()->interiorParts();

    stk::mesh::Selector sel =
        metaData.universal_part() & stk::mesh::selectUnion(partVec);

    for (label it = 0; it < nIterations; ++it)
    {
        ops::zero(curvAccumSTKFieldPtr_, partVec);
        ops::zero(curvCntSTKFieldPtr_, partVec);

        stk::mesh::BucketVector const& elementBuckets =
            bulkData.get_buckets(stk::topology::ELEMENT_RANK, sel);
        for (auto ib = elementBuckets.begin(); ib != elementBuckets.end(); ++ib)
        {
            stk::mesh::Bucket& bk = **ib;
            const auto ne = bk.size();
            for (stk::mesh::Bucket::size_type k = 0; k < ne; ++k)
            {
                stk::mesh::Entity const* nodes = bk.begin_nodes(k);
                const label nn = bk.num_nodes(k);
                scalar elemSum = 0.0;
                for (label i = 0; i < nn; ++i)
                    elemSum += *stk::mesh::field_data(*fieldPtr, nodes[i]);
                for (label i = 0; i < nn; ++i)
                {
                    *stk::mesh::field_data(*curvAccumSTKFieldPtr_, nodes[i]) +=
                        elemSum;
                    *stk::mesh::field_data(*curvCntSTKFieldPtr_, nodes[i]) +=
                        static_cast<scalar>(nn);
                }
            }
        }

        stk::mesh::BucketVector const& nodeBuckets =
            bulkData.get_buckets(stk::topology::NODE_RANK, sel);
        for (auto ib = nodeBuckets.begin(); ib != nodeBuckets.end(); ++ib)
        {
            stk::mesh::Bucket& bk = **ib;
            const auto n = bk.size();
            scalar* f = stk::mesh::field_data(*fieldPtr, bk);
            const scalar* acc =
                stk::mesh::field_data(*curvAccumSTKFieldPtr_, bk);
            const scalar* cnt = stk::mesh::field_data(*curvCntSTKFieldPtr_, bk);
            for (stk::mesh::Bucket::size_type in = 0; in < n; ++in)
                if (cnt[in] > 0.0)
                    f[in] = acc[in] / cnt[in];
        }
        if (messager::parallel())
            stk::mesh::communicate_field_data(bulkData, {fieldPtr});
    }
}

void freeSurfaceFlowModel::redistributeBodyForces(
    const std::shared_ptr<domain> domain)
{
    // Vector harmonic body force redistribution for free surface:
    //
    // For free surface flows, standard volume-weighted averaging smears
    // the body force across the density interface, creating spurious pressure
    // gradients and parasitic currents.
    //
    // Instead, we use a projection-based vector harmonic average:
    //
    // Step 1: Arithmetic mean direction, unit vector
    // Step 2: Recursive harmonic mean of scalar projections onto unit vector
    // Step 3: Result = harmonic_magnitude * unit_vector
    //
    // Step 4: Scatter the element-centre value back to all nodes,
    // weighted by SCV volumes:
    // F_node += B_el * V_scv_node
    //
    // Step 5: Normalize by dual nodal volume:
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

        // Workspace arrays (variable size, resized per bucket)
        std::vector<scalar> ws_F;
        std::vector<scalar> ws_coordinates;
        std::vector<scalar> ws_scv_volume;

        // Fixed-size workspace arrays (reused in inner loop)
        std::vector<scalar> F_el(SPATIAL_DIM);
        std::vector<scalar> uhat(SPATIAL_DIM);

        // pointers to fixed-size arrays
        scalar* p_F_el = &F_el[0];
        scalar* p_uhat = &uhat[0];

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
                }

                // Compute SCV volumes
                scalar scv_error = 0.0;
                meSCV->determinant(
                    1, &ws_coordinates[0], &ws_scv_volume[0], &scv_error);

                // projection-based harmonic mean
                // Step 1: arithmetic mean direction
                for (label d = 0; d < SPATIAL_DIM; ++d)
                {
                    p_F_el[d] = 0.0;
                }
                for (label ni = 0; ni < numNodes; ++ni)
                {
                    for (label d = 0; d < SPATIAL_DIM; ++d)
                    {
                        p_F_el[d] += ws_F[ni * SPATIAL_DIM + d];
                    }
                }
                const scalar invN = 1.0 / static_cast<scalar>(numNodes);
                for (label d = 0; d < SPATIAL_DIM; ++d)
                {
                    p_F_el[d] *= invN;
                }

                // Step 2: unit direction
                scalar mag = 0.0;
                for (label d = 0; d < SPATIAL_DIM; ++d)
                {
                    mag += p_F_el[d] * p_F_el[d];
                }
                mag = std::sqrt(mag);

                if (mag < SMALL)
                {
                    for (label d = 0; d < SPATIAL_DIM; ++d)
                        p_F_el[d] = 0.0;
                }
                else
                {
                    for (label d = 0; d < SPATIAL_DIM; ++d)
                    {
                        p_uhat[d] = p_F_el[d] / mag;
                    }

                    // Step 3: recursive harmonic of scalar projections
                    scalar h = 0.0;
                    for (label d = 0; d < SPATIAL_DIM; ++d)
                    {
                        h += p_uhat[d] * ws_F[d];
                    }
                    h = std::max(h, 0.0);

                    for (label ni = 1; ni < numNodes; ++ni)
                    {
                        scalar dk = 0.0;
                        for (label d = 0; d < SPATIAL_DIM; ++d)
                        {
                            dk += p_uhat[d] * ws_F[ni * SPATIAL_DIM + d];
                        }
                        dk = std::max(dk, 0.0);
                        if (h > 0.0)
                        {
                            h = static_cast<scalar>(ni + 1) * dk * h /
                                (static_cast<scalar>(ni) * dk + h);
                        }
                        else
                        {
                            h = 0.0;
                        }
                    }

                    // Step 4: result vector
                    for (label d = 0; d < SPATIAL_DIM; ++d)
                    {
                        p_F_el[d] = h * p_uhat[d];
                    }
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
                        F_node[d] += p_F_el[d] * scV;
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

    // F += sigma*kappa*grad(alpha): the true CSF force (retains the rotational
    // part that resists deformation). For const kappa this equals grad(psi_old)
    // so the static balance/no-currents is preserved.
    if (capillaryPotentialFieldPtr_ && sigmaKappaSTKFieldPtr_)
    {
        auto& mesh = this->meshRef();
        stk::mesh::BulkData& bulkData = mesh.bulkDataRef();
        stk::mesh::MetaData& metaData = mesh.metaDataRef();
        const stk::mesh::PartVector& partVec =
            domain->zonePtr()->interiorParts();
        const STKScalarField* gradAlphaSTKFieldPtr =
            capillaryPotentialFieldPtr_->gradRef().stkFieldPtr();
        stk::mesh::Selector selNodes =
            metaData.universal_part() & stk::mesh::selectUnion(partVec);
        stk::mesh::BucketVector const& nodeBuckets =
            bulkData.get_buckets(stk::topology::NODE_RANK, selNodes);
        for (auto ib = nodeBuckets.begin(); ib != nodeBuckets.end(); ++ib)
        {
            stk::mesh::Bucket& nodeBucket = **ib;
            const auto nNodes = nodeBucket.size();
            const scalar* gradAlphab =
                stk::mesh::field_data(*gradAlphaSTKFieldPtr, nodeBucket);
            const scalar* skb =
                stk::mesh::field_data(*sigmaKappaSTKFieldPtr_, nodeBucket);
            scalar* Fb = stk::mesh::field_data(*FSTKFieldPtr_, nodeBucket);
            for (stk::mesh::Bucket::size_type in = 0; in < nNodes; ++in)
                for (label d = 0; d < SPATIAL_DIM; ++d)
                    Fb[SPATIAL_DIM * in + d] +=
                        skb[in] * gradAlphab[SPATIAL_DIM * in + d];
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

    // enable the balanced-force capillary potential on this zone (once)
    if (capillaryPotentialFieldPtr_ &&
        capillaryPotentialFieldPtr_->isZoneUnset(domain->index()))
    {
        capillaryPotentialFieldPtr_->setZone(domain->index());
    }

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

void freeSurfaceFlowModel::updateMassFlowRateInterfaceSideField_(
    const std::shared_ptr<domain> domain,
    const interfaceSideInfo* interfaceSideInfoPtr)
{
    const auto& mesh = this->meshRef();
    const stk::mesh::MetaData& metaData = mesh.metaDataRef();
    const stk::mesh::BulkData& bulkData = mesh.bulkDataRef();

    STKScalarField& mixtureMdotfSideSTKFieldRef =
        this->mDotRef().sideFieldRef().stkFieldRef();

    // check
    assert(this->mDotRef().sideFieldRef().definedOn(
        interfaceSideInfoPtr->currentPartVec_));

    // zero the mixture mass flow rate field: in thoery it is already zero. This
    // is only for safety
    ops::zero(&mixtureMdotfSideSTKFieldRef,
              interfaceSideInfoPtr->currentPartVec_);

    // ip values; both boundary and opposing surface
    std::vector<scalar> currentIsoParCoords(SPATIAL_DIM);
    std::vector<scalar> opposingIsoParCoords(SPATIAL_DIM);

    // interpolate nodal values to point-in-elem
    const label sizeOfScalarField = 1;

    // nodal fields to gather; face
    std::vector<scalar> ws_c_alpha;
    std::vector<scalar> ws_o_alpha;

    for (label iPhase = 0; iPhase < domain->nMaterials(); iPhase++)
    {
        label phaseIndex = domain->localToGlobalMaterialIndex(iPhase);

        // Get transport fields/side fields
        const auto& alphaSTKFieldRef = this->alphaRef(phaseIndex).stkFieldRef();
        const auto& mDotSideSTKFieldPtr =
            this->mDotRef(phaseIndex).sideFieldRef().stkFieldRef();

        for (const auto& faceIpInfoVec : interfaceSideInfoPtr->ipInfoVec())
        {
            for (const ipInfo* ip : faceIpInfoVec)
            {
                if (ip->isExposed_)
                    continue;

                stk::mesh::Entity currentFace = ip->currentFace_;
                stk::mesh::Entity opposingFace = ip->opposingFace_;
                MasterElement* meFCCurrent = ip->meFCCurrent_;
                MasterElement* meFCOpposing = ip->meFCOpposing_;

                const label currentGaussPointId = ip->currentGaussPointId_;
                currentIsoParCoords = ip->currentIsoParCoords_;
                opposingIsoParCoords = ip->opposingIsoParCoords_;

                // pointer to mDot
                scalar* ncmDot = stk::mesh::field_data(
                    mixtureMdotfSideSTKFieldRef, currentFace);
                const scalar* pncmDot =
                    stk::mesh::field_data(mDotSideSTKFieldPtr, currentFace);

                const label currentNodesPerFace = meFCCurrent->nodesPerElement_;
                const label opposingNodesPerFace =
                    meFCOpposing->nodesPerElement_;

                ws_c_alpha.resize(currentNodesPerFace);
                ws_o_alpha.resize(opposingNodesPerFace);

                scalar* p_c_alpha = &ws_c_alpha[0];
                scalar* p_o_alpha = &ws_o_alpha[0];

                stk::mesh::Entity const* current_face_node_rels =
                    bulkData.begin_nodes(currentFace);
                const label current_num_face_nodes =
                    bulkData.num_nodes(currentFace);
                for (label ni = 0; ni < current_num_face_nodes; ++ni)
                {
                    stk::mesh::Entity node = current_face_node_rels[ni];
                    p_c_alpha[ni] =
                        *stk::mesh::field_data(alphaSTKFieldRef, node);
                }

                stk::mesh::Entity const* opposing_face_node_rels =
                    bulkData.begin_nodes(opposingFace);
                const label opposing_num_face_nodes =
                    bulkData.num_nodes(opposingFace);
                for (label ni = 0; ni < opposing_num_face_nodes; ++ni)
                {
                    stk::mesh::Entity node = opposing_face_node_rels[ni];
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
                    ip->areaFraction_ * (currentAlphaBip + opposingAlphaBip) /
                    2.0 * pncmDot[currentGaussPointId];
            }
        }
    }
}

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
    this->updateMassFlowRateInterior_(
        domain, this->mDotRef(iPhase), this->rhoRef(iPhase));
}

void freeSurfaceFlowModel::updateMassFlowRateInterfaceSideField_(
    const std::shared_ptr<domain> domain,
    const interfaceSideInfo* interfaceSideInfoPtr,
    label iPhase)
{
    this->updateMassFlowRateInterfaceSideField_(
        domain,
        interfaceSideInfoPtr,
        this->mDotRef(iPhase).sideFieldRef(),
        this->rhoRef(iPhase));
}

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
                            this->updateMassFlowRateBoundaryFieldInletSpecifiedPressure_(
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
                            this->updateMassFlowRateBoundaryFieldOpeningPressure_(
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
                            this->updateMassFlowRateBoundaryFieldOutletSpecifiedPressure_(
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
    const STKScalarField* fieldPtr,
    STKScalarField* rhsFieldPtr)
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

    ops::zero(rhsFieldPtr, domain->zonePtr()->interiorParts());

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

    // the field being smoothed and its diffusion-RHS accumulator
    const STKScalarField* alphaSmoothedSTKFieldPtr = fieldPtr;
    STKScalarField* rhsSmoothSTKFieldPtr = rhsFieldPtr;

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
    if (messager::parallel())
        stk::mesh::communicate_field_data(bulkData, {rhsFieldPtr});

    // find global min dxmin
    stk::all_reduce_min(MPI_COMM_WORLD, &l_dxMin, &g_dxMin, 1);

    // store the min dxmin across domains
    dxMin_ = std::min(dxMin_, g_dxMin);
}

void freeSurfaceFlowModel::computeSmoothRHS_(
    const std::shared_ptr<domain> domain,
    label iPhase)
{
    computeSmoothRHS_(domain,
                      this->alphaSmoothRef(iPhase).stkFieldPtr(),
                      this->rhsSmoothRef(iPhase).stkFieldPtr());
}

void freeSurfaceFlowModel::assembleSmoothingTerm_(
    const std::shared_ptr<domain> domain,
    STKScalarField* fieldPtr,
    const STKScalarField* rhsFieldPtr)
{
    auto& mesh = this->meshRef();
    stk::mesh::BulkData& bulkData = mesh.bulkDataRef();
    stk::mesh::MetaData& metaData = mesh.metaDataRef();

    // get interior parts the domain is defined on
    const stk::mesh::PartVector& partVec = domain->zonePtr()->interiorParts();

    // explicit update: field += Fo*dx^2*lap(field)
    const STKScalarField* rhsSmoothSTKFieldPtr = rhsFieldPtr;
    STKScalarField* alphaSmoothSTKFieldPtr = fieldPtr;

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

void freeSurfaceFlowModel::assembleSmoothingTerm_(
    const std::shared_ptr<domain> domain,
    label iPhase)
{
    assembleSmoothingTerm_(domain,
                           this->alphaSmoothRef(iPhase).stkFieldPtr(),
                           this->rhsSmoothRef(iPhase).stkFieldPtr());
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

                        case boundaryConditionType::zeroGradient:
                            {
                                // zero-gradient inlet: low-order boundary flux
                                // uses the nodal (interior) alpha, mirroring
                                // the zero-gradient outlet treatment below
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

                        case boundaryConditionType::zeroGradient:
                            {
                                // zero-gradient inlet: high-order boundary
                                // flux equals the low-order one (nodal alpha,
                                // no compression at the boundary), mirroring
                                // the zero-gradient outlet treatment below
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
                                //errorMsg("Must not reach here");
                            }
                            break;

                        case boundaryConditionType::zeroGradient:
                            {
                                // zero-gradient inlet: boundary nodes are free
                                // DOFs -- apply the limited antidiffusive side
                                // correction exactly like the zero-gradient
                                // outlet treatment below
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

void freeSurfaceFlowModel::updateMassFlowRateInterior_(
    const std::shared_ptr<domain> domain,
    elementField<scalar, 1>& mDotField,
    const nodeField<1, SPATIAL_DIM>& rhoField)
{
    using stk::mesh::Bucket;
    using stk::mesh::BucketVector;

    const auto& mesh = domain->meshRef();
    const stk::mesh::MetaData& metaData = mesh.metaDataRef();
    const stk::mesh::BulkData& bulkData = mesh.bulkDataRef();

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

    // nodal fields to gather
    std::vector<scalar> ws_U;
    std::vector<scalar> ws_Gpdx;
    std::vector<scalar> ws_du;
    std::vector<scalar> ws_coordinates;
    std::vector<scalar> ws_p;
    std::vector<scalar> ws_rho;
    std::vector<scalar> ws_betaRho;
    std::vector<scalar> ws_gradRho;
    std::vector<scalar> ws_F;
    std::vector<scalar> ws_FOrig;
    std::vector<scalar> ws_capPot;
    std::vector<scalar> ws_sk;

    // geometry related to populate
    std::vector<scalar> ws_scs_areav;
    std::vector<scalar> ws_dndx;
    std::vector<scalar> ws_deriv;
    std::vector<scalar> ws_det_j;
    std::vector<scalar> ws_velocity_shape_function;
    std::vector<scalar> ws_coordinate_shape_function;

    // integration point data that depends on size
    std::vector<scalar> uIp(SPATIAL_DIM);
    std::vector<scalar> GpdxIp(SPATIAL_DIM);
    std::vector<scalar> dpdxIp(SPATIAL_DIM);
    std::vector<scalar> duIp(SPATIAL_DIM);
    std::vector<scalar> coordIp(SPATIAL_DIM);
    std::vector<scalar> FIp(SPATIAL_DIM);
    std::vector<scalar> FOrigIp(SPATIAL_DIM);
    std::vector<scalar> B_el(SPATIAL_DIM);
    std::vector<scalar> uhat(SPATIAL_DIM);

    // pointers to everyone...
    scalar* p_uIp = &uIp[0];
    scalar* p_GpdxIp = &GpdxIp[0];
    scalar* p_dpdxIp = &dpdxIp[0];
    scalar* p_duIp = &duIp[0];
    scalar* p_coordIp = &coordIp[0];
    scalar* p_FIp = &FIp[0];
    scalar* p_FOrigIp = &FOrigIp[0];
    scalar* p_uhat = &uhat[0];

    // Get pressure diffusivity coefficient field and others
    const auto& USTKFieldRef = this->URef().stkFieldRef();
    const auto& pSTKFieldRef = this->pRef().stkFieldRef();
    const auto& rhoSTKFieldRef = rhoField.stkFieldRef();
    const auto& betaRhoSTKFieldRef = rhoField.blendingFactorRef().stkFieldRef();
    const auto& gradRhoSTKFieldRef = rhoField.gradRef().stkFieldRef();
    const auto& gradPSTKFieldRef = this->pRef().gradRef().stkFieldRef();

    const auto& mDotSTKFieldRef = mDotField.stkFieldRef();
    const scalar mDotURF = mDotField.urf();

    const auto& duSTKFieldRef =
        *metaData.get_field<scalar>(stk::topology::NODE_RANK, flowModel::du_ID);

    // Get body force fields for buoyancy pressure stabilization
    const auto* FSTKFieldPtr =
        metaData.get_field<scalar>(stk::topology::NODE_RANK, flowModel::F_ID);
    const auto* FOrigSTKFieldPtr = metaData.get_field<scalar>(
        stk::topology::NODE_RANK, flowModel::FOriginal_ID);

    // balanced-force ST: FOrigIp = B_el(FOrig) + compact grad(psi), as in the
    // pressure assembler, so the advecting flux carries the same CSF balance
    const auto* capPotSTKFieldPtr = metaData.get_field<scalar>(
        stk::topology::NODE_RANK, "capillary_potential");
    const auto* sigKappaSTKFieldPtr = metaData.get_field<scalar>(
        stk::topology::NODE_RANK, "capillary_sigma_kappa");
    const bool balancedST =
        (capPotSTKFieldPtr != nullptr && sigKappaSTKFieldPtr != nullptr);

    // Geometric fields
    const auto& coordinatesRef = *metaData.get_field<scalar>(
        stk::topology::NODE_RANK, this->getCoordinatesID_(domain));

    // get interior parts the domain is defined on
    const stk::mesh::PartVector& partVec = domain->zonePtr()->interiorParts();

    // define some common selectors
    stk::mesh::Selector selAllElements =
        metaData.universal_part() & stk::mesh::selectUnion(partVec);

    // shifted ip's for fields?
    const bool isUShifted = this->URef().isShifted();

    // shifted ip's for gradients?
    const bool isPGradientShifted = this->pRef().isGradientShifted();

    // free-surface always uses pure harmonic (no comp switch)

    BucketVector const& elementBuckets =
        bulkData.get_buckets(stk::topology::ELEMENT_RANK, selAllElements);
    for (BucketVector::const_iterator ib = elementBuckets.begin();
         ib != elementBuckets.end();
         ++ib)
    {
        Bucket& elementBucket = **ib;
        const Bucket::size_type nElementsPerBucket = elementBucket.size();

        // extract master elements
        MasterElement* meSCS = MasterElementRepo::get_surface_master_element(
            elementBucket.topology());

        // extract master element specifics
        const label nodesPerElement = meSCS->nodesPerElement_;
        const label numScsIp = meSCS->numIntPoints_;
        const label* lrscv = meSCS->adjacentNodes();

        // algorithm related
        ws_U.resize(nodesPerElement * SPATIAL_DIM);
        ws_Gpdx.resize(nodesPerElement * SPATIAL_DIM);
        ws_du.resize(nodesPerElement * SPATIAL_DIM);
        ws_coordinates.resize(nodesPerElement * SPATIAL_DIM);
        ws_p.resize(nodesPerElement);
        ws_rho.resize(nodesPerElement);
        ws_betaRho.resize(nodesPerElement);
        ws_gradRho.resize(nodesPerElement * SPATIAL_DIM);
        ws_F.resize(nodesPerElement * SPATIAL_DIM);
        ws_FOrig.resize(nodesPerElement * SPATIAL_DIM);
        ws_capPot.resize(nodesPerElement);
        ws_sk.resize(nodesPerElement);
        ws_scs_areav.resize(numScsIp * SPATIAL_DIM);
        ws_dndx.resize(SPATIAL_DIM * numScsIp * nodesPerElement);
        ws_deriv.resize(SPATIAL_DIM * numScsIp * nodesPerElement);
        ws_det_j.resize(numScsIp);
        ws_velocity_shape_function.resize(numScsIp * nodesPerElement);
        ws_coordinate_shape_function.resize(numScsIp * nodesPerElement);

        // pointers
        scalar* p_U = &ws_U[0];
        scalar* p_Gpdx = &ws_Gpdx[0];
        scalar* p_du = &ws_du[0];
        scalar* p_coordinates = &ws_coordinates[0];
        scalar* p_p = &ws_p[0];
        scalar* p_rho = &ws_rho[0];
        scalar* p_betaRho = &ws_betaRho[0];
        scalar* p_gradRho = &ws_gradRho[0];
        scalar* p_F = &ws_F[0];
        scalar* p_FOrig = &ws_FOrig[0];
        scalar* p_capPot = &ws_capPot[0];
        scalar* p_sk = &ws_sk[0];
        scalar* p_scs_areav = &ws_scs_areav[0];
        scalar* p_dndx = &ws_dndx[0];
        scalar* p_velocity_shape_function = &ws_velocity_shape_function[0];
        scalar* p_coordinate_shape_function = &ws_coordinate_shape_function[0];

        // Always use trilinear (standard) shape functions for coordinates
        meSCS->shape_fcn(&p_coordinate_shape_function[0]);

        if (isUShifted)
        {
            meSCS->shifted_shape_fcn(&p_velocity_shape_function[0]);
        }
        else
        {
            meSCS->shape_fcn(&p_velocity_shape_function[0]);
        }

        for (Bucket::size_type iElement = 0; iElement < nElementsPerBucket;
             ++iElement)
        {
            // pointers to elem data
            scalar* mDot =
                stk::mesh::field_data(mDotSTKFieldRef, elementBucket, iElement);

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
                const scalar* U = stk::mesh::field_data(USTKFieldRef, node);
                const scalar* Gjp =
                    stk::mesh::field_data(gradPSTKFieldRef, node);
                const scalar* du = stk::mesh::field_data(duSTKFieldRef, node);
                const scalar* gradRho =
                    stk::mesh::field_data(gradRhoSTKFieldRef, node);
                const scalar* coords =
                    stk::mesh::field_data(coordinatesRef, node);

                // gather scalars
                p_p[ni] = *stk::mesh::field_data(pSTKFieldRef, node);
                p_rho[ni] = *stk::mesh::field_data(rhoSTKFieldRef, node);
                p_betaRho[ni] =
                    *stk::mesh::field_data(betaRhoSTKFieldRef, node);

                // gather vectors
                const label offSet = ni * SPATIAL_DIM;
                for (label j = 0; j < SPATIAL_DIM; ++j)
                {
                    p_U[offSet + j] = U[j];
                    p_Gpdx[offSet + j] = Gjp[j];
                    p_du[offSet + j] = du[j];
                    p_gradRho[offSet + j] = gradRho[j];
                    p_coordinates[offSet + j] = coords[j];
                }

                // gather body force fields for buoyancy stabilization
                const scalar* F = stk::mesh::field_data(*FSTKFieldPtr, node);
                const scalar* FOrig =
                    stk::mesh::field_data(*FOrigSTKFieldPtr, node);
                for (label j = 0; j < SPATIAL_DIM; ++j)
                {
                    p_F[offSet + j] = F[j];
                    p_FOrig[offSet + j] = FOrig[j];
                }
                if (balancedST)
                {
                    p_capPot[ni] =
                        *stk::mesh::field_data(*capPotSTKFieldPtr, node);
                    p_sk[ni] =
                        *stk::mesh::field_data(*sigKappaSTKFieldPtr, node);
                }
            }

            // compute geometry
            scalar scs_error = 0.0;
            meSCS->determinant(
                1, &p_coordinates[0], &p_scs_areav[0], &scs_error);

            // compute dndx
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

            // projection-based harmonic mean of FOrig
            for (label d = 0; d < SPATIAL_DIM; ++d)
            {
                B_el[d] = 0.0;
            }
            for (label ni = 0; ni < nodesPerElement; ++ni)
            {
                for (label d = 0; d < SPATIAL_DIM; ++d)
                {
                    B_el[d] += p_FOrig[ni * SPATIAL_DIM + d];
                }
            }

            {
                const scalar invN = 1.0 / static_cast<scalar>(nodesPerElement);
                for (label d = 0; d < SPATIAL_DIM; ++d)
                {
                    B_el[d] *= invN;
                }
            }
            scalar mag = 0.0;
            for (label d = 0; d < SPATIAL_DIM; ++d)
            {
                mag += B_el[d] * B_el[d];
            }
            mag = std::sqrt(mag);
            if (mag < SMALL)
            {
                for (label d = 0; d < SPATIAL_DIM; ++d)
                {
                    B_el[d] = 0.0;
                }
            }
            else
            {
                for (label d = 0; d < SPATIAL_DIM; ++d)
                {
                    p_uhat[d] = B_el[d] / mag;
                }
                scalar h = 0.0;
                for (label d = 0; d < SPATIAL_DIM; ++d)
                {
                    h += p_uhat[d] * p_FOrig[d];
                }
                h = std::max(h, 0.0);
                for (label ni = 1; ni < nodesPerElement; ++ni)
                {
                    scalar dk = 0.0;
                    for (label d = 0; d < SPATIAL_DIM; ++d)
                    {
                        dk += p_uhat[d] * p_FOrig[ni * SPATIAL_DIM + d];
                    }
                    dk = std::max(dk, 0.0);
                    if (h > 0.0)
                    {
                        h = static_cast<scalar>(ni + 1) * dk * h /
                            (static_cast<scalar>(ni) * dk + h);
                    }
                    else
                    {
                        h = 0.0;
                    }
                }
                for (label d = 0; d < SPATIAL_DIM; ++d)
                {
                    B_el[d] = h * p_uhat[d];
                }
            }

            for (label ip = 0; ip < numScsIp; ++ip)
            {
                // left and right nodes for this ip
                const label il = lrscv[2 * ip];
                const label ir = lrscv[2 * ip + 1];

                const label ipNdim = ip * SPATIAL_DIM;
                const label offSetSF = ip * nodesPerElement;

                // setup for ip values; p_GpdxIp and p_FIp use
                // pure harmonic from adjacent nodes (il, ir)
                for (label j = 0; j < SPATIAL_DIM; ++j)
                {
                    p_coordIp[j] = 0.0;
                    p_uIp[j] = 0.0;
                    p_dpdxIp[j] = 0.0;
                    p_duIp[j] = 0.0;
                    p_FOrigIp[j] = B_el[j];

                    // pressure gradient interpolation: harmonic for Gpdx
                    const scalar gL = p_Gpdx[il * SPATIAL_DIM + j];
                    const scalar gR = p_Gpdx[ir * SPATIAL_DIM + j];
                    p_GpdxIp[j] = (std::abs(gL) * gR + gL * std::abs(gR)) /
                                  (std::abs(gL) + std::abs(gR) + SMALL);

                    // body force interpolation: harmonic for F
                    const scalar fL = p_F[il * SPATIAL_DIM + j];
                    const scalar fR = p_F[ir * SPATIAL_DIM + j];
                    p_FIp[j] = (std::abs(fL) * fR + fL * std::abs(fR)) /
                               (std::abs(fL) + std::abs(fR) + SMALL);
                }

                scalar rhoIp = 0.0;
                scalar skIp = 0.0;
                scalar gradAlphaIp[SPATIAL_DIM] = {};
                for (label ic = 0; ic < nodesPerElement; ++ic)
                {
                    const scalar r_vel =
                        p_velocity_shape_function[offSetSF + ic];
                    const scalar r_coord =
                        p_coordinate_shape_function[offSetSF + ic];

                    // use velocity shape functions
                    rhoIp += r_vel * p_rho[ic];
                    if (balancedST)
                        skIp += r_vel * p_sk[ic];

                    const label offSetDnDx =
                        SPATIAL_DIM * nodesPerElement * ip + ic * SPATIAL_DIM;
                    for (label j = 0; j < SPATIAL_DIM; ++j)
                    {
                        // use velocity shape functions
                        p_duIp[j] += r_vel * p_du[SPATIAL_DIM * ic + j];
                        p_uIp[j] += r_vel * p_U[SPATIAL_DIM * ic + j];

                        // use pressure shape function derivative
                        p_dpdxIp[j] += p_dndx[offSetDnDx + j] * p_p[ic];

                        // compact grad(alpha) (same dndx as grad p)
                        if (balancedST)
                            gradAlphaIp[j] +=
                                p_dndx[offSetDnDx + j] * p_capPot[ic];

                        // use coordinates shape functions
                        p_coordIp[j] +=
                            r_coord * p_coordinates[SPATIAL_DIM * ic + j];
                    }
                }

                // balanced CSF: FOrigIp += sigma*kappa_ip * compact grad(alpha)
                if (balancedST)
                    for (label j = 0; j < SPATIAL_DIM; ++j)
                        p_FOrigIp[j] += skIp * gradAlphaIp[j];

                // calculate a local relative flow rate
                scalar tUDotRel = mDot[ip] / rhoIp;
                for (label i = 0; i < SPATIAL_DIM; i++)
                {
                    for (label j = 0; j < SPATIAL_DIM; j++)
                    {
                        tUDotRel -= p_mat[i * SPATIAL_DIM + j] *
                                    (p_coordIp[j] - p_ori[j]) *
                                    p_scs_areav[ipNdim + i];
                    }
                }

                scalar dcorr = 0;
                if (tUDotRel > 0)
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

                scalar rhoUpwind;
                if (tUDotRel > 0)
                {
                    rhoUpwind = p_rho[il];
                }
                else
                {
                    rhoUpwind = p_rho[ir];
                }

                scalar rhoHR = rhoUpwind + dcorr;

                // rhie-chow
                // mDot = ρ*U·S - ρ*D*(∇p - Gp)·S + ρ*D*(F_orig - F)·S
                //        ╰───╯   ╰──────────────╯   ╰─────────────────╯
                //      divergence   pressure RC       body force stab
                //
                scalar tmDot = 0.0;
                for (label j = 0; j < SPATIAL_DIM; ++j)
                {
                    // divergence: ρ*U·S
                    tmDot +=
                        rhoHR * p_uIp[j] * p_scs_areav[ip * SPATIAL_DIM + j];

                    // pressure Rhie-Chow: -ρ*D*(∇p - Gp)·S
                    tmDot -= rhoHR * p_duIp[j] * (p_dpdxIp[j] - p_GpdxIp[j]) *
                             p_scs_areav[ip * SPATIAL_DIM + j];

                    // body force stabilization: +ρ*D*(F_orig - F)·S
                    tmDot += rhoHR * p_duIp[j] * (p_FOrigIp[j] - p_FIp[j]) *
                             p_scs_areav[ip * SPATIAL_DIM + j];
                }

                // store with relaxation: mDot[ip] at this point must be in
                // absolute frame
                mDot[ip] = mDotURF * tmDot + (1.0 - mDotURF) * mDot[ip];
            }
        }
    }
}

void freeSurfaceFlowModel::updateMassFlowRateInterfaceSideField_(
    const std::shared_ptr<domain> domain,
    const interfaceSideInfo* interfaceSideInfoPtr,
    sideField<scalar, 1>& mDotSideField,
    const nodeField<1, SPATIAL_DIM>& rhoField)
{
    const auto& mesh = this->meshRef();
    const stk::mesh::MetaData& metaData = mesh.metaDataRef();
    const stk::mesh::BulkData& bulkData = mesh.bulkDataRef();

    scalar penaltyFactor = interfaceSideInfoPtr->interfPtr()->penaltyFactor();

    STKScalarField& mDotSideSTKFieldRef = mDotSideField.stkFieldRef();
    const scalar mDotURF = mDotSideField.urf();

    // check
    assert(mDotSideField.definedOn(interfaceSideInfoPtr->currentPartVec_));

    // ip values; both boundary and opposing surface
    std::vector<scalar> currentIsoParCoords(SPATIAL_DIM);
    std::vector<scalar> opposingIsoParCoords(SPATIAL_DIM);
    std::vector<scalar> cNx(SPATIAL_DIM);
    std::vector<scalar> oNx(SPATIAL_DIM);
    std::vector<scalar> cRhoUBip(SPATIAL_DIM);
    std::vector<scalar> oRhoUBip(SPATIAL_DIM);
    std::vector<scalar> currentCoordsBip(SPATIAL_DIM);
    std::vector<scalar> cDuBip(SPATIAL_DIM);
    std::vector<scalar> oDuBip(SPATIAL_DIM);

    // pressure stabilization
    std::vector<scalar> cGjpBip(SPATIAL_DIM);
    std::vector<scalar> oGjpBip(SPATIAL_DIM);
    std::vector<scalar> cDpdxBip(SPATIAL_DIM);
    std::vector<scalar> oDpdxBip(SPATIAL_DIM);

    // mapping for -1:1 -> -0.5:0.5 volume element
    std::vector<scalar> currentElementIsoParCoords(SPATIAL_DIM);
    std::vector<scalar> opposingElementIsoParCoords(SPATIAL_DIM);

    // interpolate nodal values to point-in-elem
    const label sizeOfScalarField = 1;
    const label sizeOfVectorField = SPATIAL_DIM;

    // pointers to fixed values
    scalar* p_cNx = &cNx[0];
    scalar* p_oNx = &oNx[0];

    // nodal fields to gather; face
    std::vector<scalar> ws_c_p;
    std::vector<scalar> ws_o_p;
    std::vector<scalar> ws_c_Gjp;
    std::vector<scalar> ws_o_Gjp;
    std::vector<scalar> ws_c_rhoU;
    std::vector<scalar> ws_o_rhoU;
    std::vector<scalar> ws_c_rho;
    std::vector<scalar> ws_o_rho;
    std::vector<scalar> ws_c_face_coordinates;

    // element
    std::vector<scalar> ws_c_elem_p;
    std::vector<scalar> ws_o_elem_p;
    std::vector<scalar> ws_c_elem_coordinates;
    std::vector<scalar> ws_o_elem_coordinates;
    std::vector<scalar> ws_c_du;
    std::vector<scalar> ws_o_du;

    // master element data
    std::vector<scalar> ws_c_dndx;
    std::vector<scalar> ws_o_dndx;
    std::vector<scalar> ws_c_det_j;
    std::vector<scalar> ws_o_det_j;

    // Get transport fields/side fields
    const auto& rhoSTKFieldRef = rhoField.stkFieldRef();
    const auto& pSTKFieldRef = this->pRef().stkFieldRef();
    const auto& gradPSTKFieldRef = this->pRef().gradRef().stkFieldRef();
    const auto& USTKFieldRef = this->URef().stkFieldRef();

    // Get pressure diffusivity coefficient field
    const auto& duSTKFieldRef =
        *metaData.get_field<scalar>(stk::topology::NODE_RANK, flowModel::du_ID);

    // Get geometric fields
    const auto& coordsSTKFieldRef = *metaData.get_field<scalar>(
        stk::topology::NODE_RANK, this->getCoordinatesID_(domain));
    const auto& exposedAreaVecSTKFieldRef = *metaData.get_field<scalar>(
        metaData.side_rank(), this->getExposedAreaVectorID_(domain));

    std::vector<scalar> oldMdot;
    std::vector<scalar> accumulatedMdot;
    std::vector<bool> exposedMdot;
    std::vector<bool> touchedMdot;

    // free-surface always uses pure harmonic (no comp switch)
    for (const auto& faceIpInfoVec : interfaceSideInfoPtr->ipInfoVec())
    {
        if (faceIpInfoVec.empty())
            continue;

        const ipInfo* firstIp = faceIpInfoVec.front();
        scalar* faceMdot =
            stk::mesh::field_data(mDotSideSTKFieldRef, firstIp->currentFace_);
        const label numScsBip = firstIp->meFCCurrent_->numIntPoints_;
        oldMdot.assign(faceMdot, faceMdot + numScsBip);
        accumulatedMdot.assign(numScsBip, 0.0);
        exposedMdot.assign(numScsBip, false);
        touchedMdot.assign(numScsBip, false);

        for (const ipInfo* ip : faceIpInfoVec)
        {
            // extract current/opposing face/element
            stk::mesh::Entity currentFace = ip->currentFace_;
            stk::mesh::Entity opposingFace = ip->opposingFace_;
            stk::mesh::Entity currentElement = ip->currentElement_;
            stk::mesh::Entity opposingElement = ip->opposingElement_;
            const label currentFaceOrdinal = ip->currentFaceOrdinal_;
            const label opposingFaceOrdinal = ip->opposingFaceOrdinal_;

            // master element; face and volume
            MasterElement* meFCCurrent = ip->meFCCurrent_;
            MasterElement* meFCOpposing = ip->meFCOpposing_;
            MasterElement* meSCSCurrent = ip->meSCSCurrent_;
            MasterElement* meSCSOpposing = ip->meSCSOpposing_;

            // local ip, ordinals, etc
            const label currentGaussPointId = ip->currentGaussPointId_;
            touchedMdot[currentGaussPointId] = true;
            currentIsoParCoords = ip->currentIsoParCoords_;
            opposingIsoParCoords = ip->opposingIsoParCoords_;

            // if gauss point is exposed (non-overlapping), then
            // treat as a wall
            if (ip->isExposed_)
            {
                exposedMdot[currentGaussPointId] = true;
                continue;
            }

            // extract some master element info
            const label currentNodesPerSide = meFCCurrent->nodesPerElement_;
            const label opposingNodesPerSide = meFCOpposing->nodesPerElement_;
            const label currentNodesPerElement = meSCSCurrent->nodesPerElement_;
            const label opposingNodesPerElement =
                meSCSOpposing->nodesPerElement_;

            // arithmetic weights
            const scalar f_c = 1.0 / static_cast<scalar>(currentNodesPerSide);
            const scalar f_o = 1.0 / static_cast<scalar>(opposingNodesPerSide);

            // algorithm related; face
            ws_c_p.resize(currentNodesPerSide);
            ws_o_p.resize(opposingNodesPerSide);
            ws_c_du.resize(currentNodesPerSide * SPATIAL_DIM);
            ws_o_du.resize(opposingNodesPerSide * SPATIAL_DIM);
            ws_c_Gjp.resize(currentNodesPerSide * SPATIAL_DIM);
            ws_o_Gjp.resize(opposingNodesPerSide * SPATIAL_DIM);
            ws_c_rhoU.resize(currentNodesPerSide * SPATIAL_DIM);
            ws_o_rhoU.resize(opposingNodesPerSide * SPATIAL_DIM);
            ws_c_rho.resize(currentNodesPerSide);
            ws_o_rho.resize(opposingNodesPerSide);
            ws_c_face_coordinates.resize(currentNodesPerSide * SPATIAL_DIM);

            // algorithm related; element; dndx will be at a
            // single gauss point
            ws_c_elem_p.resize(currentNodesPerElement);
            ws_o_elem_p.resize(opposingNodesPerElement);
            ws_c_elem_coordinates.resize(currentNodesPerElement * SPATIAL_DIM);
            ws_o_elem_coordinates.resize(opposingNodesPerElement * SPATIAL_DIM);
            ws_c_dndx.resize(SPATIAL_DIM * currentNodesPerElement);
            ws_o_dndx.resize(SPATIAL_DIM * opposingNodesPerElement);
            ws_c_det_j.resize(1);
            ws_o_det_j.resize(1);

            // face
            scalar* p_c_p = &ws_c_p[0];
            scalar* p_o_p = &ws_o_p[0];
            scalar* p_c_du = &ws_c_du[0];
            scalar* p_o_du = &ws_o_du[0];
            scalar* p_c_Gjp = &ws_c_Gjp[0];
            scalar* p_o_Gjp = &ws_o_Gjp[0];
            scalar* p_c_rhoU = &ws_c_rhoU[0];
            scalar* p_o_rhoU = &ws_o_rhoU[0];
            scalar* p_c_rho = &ws_c_rho[0];
            scalar* p_o_rho = &ws_o_rho[0];
            scalar* p_c_face_coordinates = &ws_c_face_coordinates[0];

            // element
            scalar* p_c_elem_p = &ws_c_elem_p[0];
            scalar* p_o_elem_p = &ws_o_elem_p[0];
            scalar* p_c_elem_coordinates = &ws_c_elem_coordinates[0];
            scalar* p_o_elem_coordinates = &ws_o_elem_coordinates[0];

            // me pointers
            scalar* p_c_dndx = &ws_c_dndx[0];
            scalar* p_o_dndx = &ws_o_dndx[0];

            // populate current face_node_ordinals
            const label* c_face_node_ordinals =
                meSCSCurrent->side_node_ordinals(currentFaceOrdinal);

            // gather current face data
            stk::mesh::Entity const* current_face_node_rels =
                bulkData.begin_nodes(currentFace);
            const label current_num_face_nodes =
                bulkData.num_nodes(currentFace);
            for (label ni = 0; ni < current_num_face_nodes; ++ni)
            {
                stk::mesh::Entity node = current_face_node_rels[ni];

                // gather; scalar
                p_c_p[ni] = *stk::mesh::field_data(pSTKFieldRef, node);
                p_c_rho[ni] = *stk::mesh::field_data(rhoSTKFieldRef, node);

                // gather; vector
                const scalar* U = stk::mesh::field_data(USTKFieldRef, node);
                const scalar* Gjp =
                    stk::mesh::field_data(gradPSTKFieldRef, node);
                const scalar* coords =
                    stk::mesh::field_data(coordsSTKFieldRef, node);
                const scalar* du = stk::mesh::field_data(duSTKFieldRef, node);
                for (label i = 0; i < SPATIAL_DIM; ++i)
                {
                    const label offSet = i * current_num_face_nodes + ni;
                    p_c_rhoU[offSet] = U[i];
                    p_c_Gjp[offSet] = Gjp[i];
                    p_c_face_coordinates[offSet] = coords[i];
                    p_c_du[offSet] = du[i];
                }
            }

            // populate opposing face_node_ordinals
            const label* o_face_node_ordinals =
                meSCSOpposing->side_node_ordinals(opposingFaceOrdinal);

            // gather opposing face data
            stk::mesh::Entity const* opposing_face_node_rels =
                bulkData.begin_nodes(opposingFace);
            const label opposing_num_face_nodes =
                bulkData.num_nodes(opposingFace);
            for (label ni = 0; ni < opposing_num_face_nodes; ++ni)
            {
                stk::mesh::Entity node = opposing_face_node_rels[ni];

                // gather; scalar
                p_o_p[ni] = *stk::mesh::field_data(pSTKFieldRef, node);
                p_o_rho[ni] = *stk::mesh::field_data(rhoSTKFieldRef, node);

                // gather; vector
                const scalar* U = stk::mesh::field_data(USTKFieldRef, node);
                const scalar* Gjp =
                    stk::mesh::field_data(gradPSTKFieldRef, node);
                const scalar* du = stk::mesh::field_data(duSTKFieldRef, node);
                for (label i = 0; i < SPATIAL_DIM; ++i)
                {
                    const label offSet = i * opposing_num_face_nodes + ni;
                    p_o_rhoU[offSet] = U[i];
                    p_o_Gjp[offSet] = Gjp[i];
                    p_o_du[offSet] = du[i];
                }
            }

            // gather current element data
            stk::mesh::Entity const* current_elem_node_rels =
                bulkData.begin_nodes(currentElement);
            const label current_num_elem_nodes =
                bulkData.num_nodes(currentElement);
            for (label ni = 0; ni < current_num_elem_nodes; ++ni)
            {
                stk::mesh::Entity node = current_elem_node_rels[ni];

                // gather; scalar
                p_c_elem_p[ni] = *stk::mesh::field_data(pSTKFieldRef, node);

                // gather; vector
                const scalar* coords =
                    stk::mesh::field_data(coordsSTKFieldRef, node);
                const label niNdim = ni * SPATIAL_DIM;
                for (label i = 0; i < SPATIAL_DIM; ++i)
                {
                    p_c_elem_coordinates[niNdim + i] = coords[i];
                }
            }

            // gather opposing element data; sneak in second
            // connected nodes
            stk::mesh::Entity const* opposing_elem_node_rels =
                bulkData.begin_nodes(opposingElement);
            const label opposing_num_elem_nodes =
                bulkData.num_nodes(opposingElement);
            for (label ni = 0; ni < opposing_num_elem_nodes; ++ni)
            {
                stk::mesh::Entity node = opposing_elem_node_rels[ni];

                // gather; scalar
                p_o_elem_p[ni] = *stk::mesh::field_data(pSTKFieldRef, node);

                // gather; vector
                const scalar* coords =
                    stk::mesh::field_data(coordsSTKFieldRef, node);
                const label niNdim = ni * SPATIAL_DIM;
                for (label i = 0; i < SPATIAL_DIM; ++i)
                {
                    p_o_elem_coordinates[niNdim + i] = coords[i];
                }
            }

            // pointer to face data
            const scalar* c_areaVec =
                stk::mesh::field_data(exposedAreaVecSTKFieldRef, currentFace);

            scalar c_amag = 0.0;
            for (label j = 0; j < SPATIAL_DIM; ++j)
            {
                const scalar c_axj =
                    c_areaVec[currentGaussPointId * SPATIAL_DIM + j];
                c_amag += c_axj * c_axj;
            }
            c_amag = std::sqrt(c_amag) * ip->areaFraction_;

            // now compute normal
            for (label i = 0; i < SPATIAL_DIM; ++i)
            {
                p_cNx[i] =
                    c_areaVec[currentGaussPointId * SPATIAL_DIM + i] / c_amag;
            }

            // compute opposing normal: in theory it is assumed
            // that the current and opposing sub-control surfaces
            // are sufficiently planar
            for (label i = 0; i < SPATIAL_DIM; ++i)
            {
                p_oNx[i] = -p_cNx[i];
            }

            // transform opposing normal back to opposing side
            interfaceSideInfoPtr->reverseRotateVector<SPATIAL_DIM>(oNx);

            // project from side to element; method deals with
            // the -1:1 isInElement range to the proper
            // underlying CVFEM range
            meSCSCurrent->sidePcoords_to_elemPcoords(
                currentFaceOrdinal,
                1,
                &currentIsoParCoords[0],
                &currentElementIsoParCoords[0]);
            meSCSOpposing->sidePcoords_to_elemPcoords(
                opposingFaceOrdinal,
                1,
                &opposingIsoParCoords[0],
                &opposingElementIsoParCoords[0]);

            // compute dndx
            scalar scs_error = 0.0;
            meSCSCurrent->general_face_grad_op(currentFaceOrdinal,
                                               &currentElementIsoParCoords[0],
                                               &p_c_elem_coordinates[0],
                                               &p_c_dndx[0],
                                               &ws_c_det_j[0],
                                               &scs_error);
            meSCSOpposing->general_face_grad_op(opposingFaceOrdinal,
                                                &opposingElementIsoParCoords[0],
                                                &p_o_elem_coordinates[0],
                                                &p_o_dndx[0],
                                                &ws_o_det_j[0],
                                                &scs_error);

            // current inverse length scale; can loop over face
            // nodes to avoid "nodesOnFace" array
            scalar currentInverseLength = 0.0;
            for (label ic = 0; ic < current_num_face_nodes; ++ic)
            {
                const label faceNodeNumber = c_face_node_ordinals[ic];
                const label offSetDnDx =
                    faceNodeNumber * SPATIAL_DIM; // single intg. point
                for (label j = 0; j < SPATIAL_DIM; ++j)
                {
                    const scalar nxj = p_cNx[j];
                    const scalar dndxj = p_c_dndx[offSetDnDx + j];
                    currentInverseLength += dndxj * nxj;
                }
            }

            // opposing inverse length scale; can loop over face
            // nodes to avoid "nodesOnFace" array
            scalar opposingInverseLength = 0.0;
            for (label ic = 0; ic < opposing_num_face_nodes; ++ic)
            {
                const label faceNodeNumber = o_face_node_ordinals[ic];
                const label offSetDnDx =
                    faceNodeNumber * SPATIAL_DIM; // single intg. point
                for (label j = 0; j < SPATIAL_DIM; ++j)
                {
                    const scalar nxj = p_oNx[j];
                    const scalar dndxj = p_o_dndx[offSetDnDx + j];
                    opposingInverseLength += dndxj * nxj;
                }
            }

            // bip gradients; zero out
            for (label j = 0; j < SPATIAL_DIM; ++j)
            {
                cDpdxBip[j] = 0.0;
                oDpdxBip[j] = 0.0;
            }

            // current pressure gradient
            for (label ic = 0; ic < currentNodesPerElement; ++ic)
            {
                const label offSetDnDx = ic * SPATIAL_DIM; // single intg. point
                const scalar pNp1 = p_c_elem_p[ic];
                for (label j = 0; j < SPATIAL_DIM; ++j)
                {
                    const scalar dndxj = p_c_dndx[offSetDnDx + j];
                    cDpdxBip[j] += dndxj * pNp1;
                }
            }

            // opposing pressure gradient
            for (label ic = 0; ic < opposingNodesPerElement; ++ic)
            {
                const label offSetDnDx = ic * SPATIAL_DIM; // single intg. point
                const scalar pNp1 = p_o_elem_p[ic];
                for (label j = 0; j < SPATIAL_DIM; ++j)
                {
                    const scalar dndxj = p_o_dndx[offSetDnDx + j];
                    oDpdxBip[j] += dndxj * pNp1;
                }
            }

            // interpolate to boundary ips
            scalar cPBip = 0.0;
            meFCCurrent->interpolatePoint(
                sizeOfScalarField, &currentIsoParCoords[0], &ws_c_p[0], &cPBip);

            scalar oPBip = 0.0;
            meFCOpposing->interpolatePoint(sizeOfScalarField,
                                           &opposingIsoParCoords[0],
                                           &ws_o_p[0],
                                           &oPBip);

            scalar cRhoBip = 0.0;
            meFCCurrent->interpolatePoint(sizeOfScalarField,
                                          &currentIsoParCoords[0],
                                          &ws_c_rho[0],
                                          &cRhoBip);

            scalar oRhoBip = 0.0;
            meFCOpposing->interpolatePoint(sizeOfScalarField,
                                           &opposingIsoParCoords[0],
                                           &ws_o_rho[0],
                                           &oRhoBip);

            // projected nodal gradient: use arithmetic
            // interpolations: zero-out first
            for (label i = 0; i < SPATIAL_DIM; i++)
            {
                cGjpBip[i] = 0;
                oGjpBip[i] = 0;
            }

            for (label ic = 0; ic < currentNodesPerSide; ++ic)
            {
                for (label i = 0; i < SPATIAL_DIM; i++)
                {
                    const label offSet = i * currentNodesPerSide + ic;
                    cGjpBip[i] += f_c * p_c_Gjp[offSet];
                }
            }

            for (label ic = 0; ic < opposingNodesPerSide; ++ic)
            {
                for (label i = 0; i < SPATIAL_DIM; i++)
                {
                    const label offSet = i * opposingNodesPerSide + ic;
                    oGjpBip[i] += f_o * p_o_Gjp[offSet];
                }
            }

            // pure harmonic interpolation between the two
            // interface sides for free-surface
            for (label i = 0; i < SPATIAL_DIM; ++i)
            {
                const scalar a = cGjpBip[i];
                const scalar b = oGjpBip[i];
                const scalar harm = (std::abs(a) * b + a * std::abs(b)) /
                                    (std::abs(a) + std::abs(b) + SMALL);
                cGjpBip[i] = harm;
                oGjpBip[i] = harm;
            }

            // product of density and velocity; current (take
            // over previous nodal value for velocity)
            for (label ni = 0; ni < current_num_face_nodes; ++ni)
            {
                const scalar rho = p_c_rho[ni];
                for (label i = 0; i < SPATIAL_DIM; ++i)
                {
                    const label offSet = i * current_num_face_nodes + ni;
                    p_c_rhoU[offSet] *= rho;
                }
            }

            // opposite
            for (label ni = 0; ni < opposing_num_face_nodes; ++ni)
            {
                const scalar rho = p_o_rho[ni];
                for (label i = 0; i < SPATIAL_DIM; ++i)
                {
                    const label offSet = i * opposing_num_face_nodes + ni;
                    p_o_rhoU[offSet] *= rho;
                }
            }

            // interpolate velocity with density scaling
            meFCCurrent->interpolatePoint(sizeOfVectorField,
                                          &currentIsoParCoords[0],
                                          &ws_c_rhoU[0],
                                          &cRhoUBip[0]);

            meFCOpposing->interpolatePoint(sizeOfVectorField,
                                           &opposingIsoParCoords[0],
                                           &ws_o_rhoU[0],
                                           &oRhoUBip[0]);

            meFCCurrent->interpolatePoint(sizeOfVectorField,
                                          &currentIsoParCoords[0],
                                          &ws_c_du[0],
                                          &cDuBip[0]);

            meFCOpposing->interpolatePoint(sizeOfVectorField,
                                           &opposingIsoParCoords[0],
                                           &ws_o_du[0],
                                           &oDuBip[0]);

            meFCCurrent->interpolatePoint(sizeOfVectorField,
                                          &currentIsoParCoords[0],
                                          &ws_c_face_coordinates[0],
                                          &currentCoordsBip[0]);

            scalar currentDiffBipmag = 0;
            scalar opposingDiffBipmag = 0;
            for (label i = 0; i < SPATIAL_DIM; i++)
            {
                currentDiffBipmag +=
                    cRhoBip * cDuBip[i] / static_cast<scalar>(SPATIAL_DIM);
                opposingDiffBipmag +=
                    oRhoBip * oDuBip[i] / static_cast<scalar>(SPATIAL_DIM);
            }

            const scalar penaltyIp =
                penaltyFactor * 0.5 *
                (currentDiffBipmag * currentInverseLength +
                 opposingDiffBipmag * opposingInverseLength);

            scalar ncFlux = 0.0;
            scalar ncPstabFlux = 0.0;
            for (label j = 0; j < SPATIAL_DIM; ++j)
            {
                const scalar cRhoU = cRhoUBip[j];
                const scalar oRhoU = oRhoUBip[j];
                ncFlux += 0.5 * (cRhoU * p_cNx[j] - oRhoU * p_oNx[j]);

                const scalar cPstab =
                    cRhoBip * cDuBip[j] * (cDpdxBip[j] - cGjpBip[j]);
                const scalar oPstab =
                    oRhoBip * oDuBip[j] * (oDpdxBip[j] - oGjpBip[j]);
                ncPstabFlux += 0.5 * (cPstab * p_cNx[j] - oPstab * p_oNx[j]);
            }

            // scatter it
            scalar tmDot =
                (ncFlux - ncPstabFlux + penaltyIp * (cPBip - oPBip)) * c_amag;

            accumulatedMdot[currentGaussPointId] += tmDot;
        }

        for (label ip = 0; ip < numScsBip; ++ip)
        {
            if (!touchedMdot[ip])
                continue;

            faceMdot[ip] = exposedMdot[ip] ? 0.0
                                           : mDotURF * accumulatedMdot[ip] +
                                                 (1.0 - mDotURF) * oldMdot[ip];
        }
    }
}

void freeSurfaceFlowModel::
    updateMassFlowRateBoundaryFieldInletSpecifiedPressure_(
        const std::shared_ptr<domain> domain,
        const boundary* boundary,
        sideField<scalar, 1>& mDotSideField,
        const nodeField<1, SPATIAL_DIM>& rhoField)
{
    using stk::mesh::Bucket;
    using stk::mesh::BucketVector;

    const auto& mesh = this->meshRef();
    const stk::mesh::MetaData& metaData = mesh.metaDataRef();
    const stk::mesh::BulkData& bulkData = mesh.bulkDataRef();

    // ip values; both boundary and opposing
    // surface
    std::vector<scalar> coordBip(SPATIAL_DIM);
    std::vector<scalar> uBip(SPATIAL_DIM);
    std::vector<scalar> GpdxBip(SPATIAL_DIM);
    std::vector<scalar> dpdxBip(SPATIAL_DIM);
    std::vector<scalar> duBip(SPATIAL_DIM);
    std::vector<scalar> FBip(SPATIAL_DIM);
    std::vector<scalar> FOrigBip(SPATIAL_DIM);
    std::vector<scalar> uhat(SPATIAL_DIM);

    // pointers to fixed values
    scalar* p_coordBip = &coordBip[0];
    scalar* p_uBip = &uBip[0];
    scalar* p_GpdxBip = &GpdxBip[0];
    scalar* p_dpdxBip = &dpdxBip[0];
    scalar* p_duBip = &duBip[0];
    scalar* p_FBip = &FBip[0];
    scalar* p_FOrigBip = &FOrigBip[0];
    scalar* p_uhat = &uhat[0];

    // nodal fields to gather
    std::vector<scalar> ws_coordinates;
    std::vector<scalar> ws_p;
    std::vector<scalar> ws_U;
    std::vector<scalar> ws_Gpdx;
    std::vector<scalar> ws_du;
    std::vector<scalar> ws_rho;
    std::vector<scalar> ws_F;
    std::vector<scalar> ws_FOrig;
    std::vector<scalar> ws_Gpdx_elem;
    std::vector<scalar> ws_F_elem;
    std::vector<scalar> ws_FOrig_elem;
    std::vector<scalar> ws_capPot;
    std::vector<scalar> ws_sk;
    std::vector<scalar> B_el(SPATIAL_DIM);

    // master element
    std::vector<scalar> ws_face_shape_function;
    std::vector<scalar> ws_coordinate_face_shape_function;
    std::vector<scalar> ws_dndx;
    std::vector<scalar> ws_det_j;

    // Get fields
    const auto& rhoSTKFieldRef = rhoField.stkFieldRef();
    const auto& pSTKFieldRef = this->pRef().stkFieldRef();
    const auto& nodalSidePSTKFieldRef =
        this->pRef().nodeSideFieldRef().stkFieldRef();
    const auto& USTKFieldRef = this->URef().stkFieldRef();
    const auto& gradPSTKFieldRef = this->pRef().gradRef().stkFieldRef();

    const auto& mDotSideSTKFieldRef = mDotSideField.stkFieldRef();
    const scalar mDotURF = mDotSideField.urf();

    const auto& duSTKFieldRef =
        *metaData.get_field<scalar>(stk::topology::NODE_RANK, flowModel::du_ID);

    // Get body force fields for buoyancy pressure stabilization
    const auto* FSTKFieldPtr =
        metaData.get_field<scalar>(stk::topology::NODE_RANK, flowModel::F_ID);
    const auto* FOrigSTKFieldPtr = metaData.get_field<scalar>(
        stk::topology::NODE_RANK, flowModel::FOriginal_ID);

    const auto* capPotSTKFieldPtr = metaData.get_field<scalar>(
        stk::topology::NODE_RANK, "capillary_potential");
    const auto* sigKappaSTKFieldPtr = metaData.get_field<scalar>(
        stk::topology::NODE_RANK, "capillary_sigma_kappa");
    const bool balancedST =
        (capPotSTKFieldPtr != nullptr && sigKappaSTKFieldPtr != nullptr);

    // Get geometric fields
    const auto& coordsSTKFieldRef = *metaData.get_field<scalar>(
        stk::topology::NODE_RANK, this->getCoordinatesID_(domain));
    const auto& exposedAreaVecSTKFieldRef = *metaData.get_field<scalar>(
        metaData.side_rank(), this->getExposedAreaVectorID_(domain));

    // define vector of parent topos; should
    // always be UNITY in size
    std::vector<stk::topology> parentTopo;

    // define some common selectors
    stk::mesh::Selector selAllSides =
        metaData.universal_part() & stk::mesh::selectUnion(boundary->parts());

    // shifted ip's for field?
    const bool isUShifted = this->URef().isShifted();

    // shifted ip's for gradients?
    const bool isPGradientShifted = this->pRef().isGradientShifted();

    // free-surface always uses pure harmonic (no comp switch)

    BucketVector const& sideBuckets =
        bulkData.get_buckets(metaData.side_rank(), selAllSides);
    for (BucketVector::const_iterator ib = sideBuckets.begin();
         ib != sideBuckets.end();
         ++ib)
    {
        Bucket& sideBucket = **ib;

        // extract connected element
        // topology
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
        const label* faceIpNodeMap = meFC->ipNodeMap();

        // algorithm related; element
        // (exposed face and element)
        ws_coordinates.resize(nodesPerElement * SPATIAL_DIM);
        ws_p.resize(nodesPerElement);
        ws_capPot.resize(nodesPerElement);
        ws_sk.resize(nodesPerElement);
        ws_U.resize(nodesPerSide * SPATIAL_DIM);
        ws_Gpdx.resize(nodesPerSide * SPATIAL_DIM);
        ws_du.resize(nodesPerSide * SPATIAL_DIM);
        ws_rho.resize(nodesPerSide);
        ws_F.resize(nodesPerSide * SPATIAL_DIM);
        ws_FOrig.resize(nodesPerSide * SPATIAL_DIM);
        ws_face_shape_function.resize(numScsBip * nodesPerSide);
        ws_coordinate_face_shape_function.resize(numScsBip * nodesPerSide);
        ws_dndx.resize(SPATIAL_DIM * numScsBip * nodesPerElement);
        ws_det_j.resize(numScsBip);
        ws_Gpdx_elem.resize(nodesPerElement * SPATIAL_DIM);
        ws_F_elem.resize(nodesPerElement * SPATIAL_DIM);
        ws_FOrig_elem.resize(nodesPerElement * SPATIAL_DIM);

        // pointers
        scalar* p_coordinates = &ws_coordinates[0];
        scalar* p_p = &ws_p[0];
        scalar* p_capPot = &ws_capPot[0];
        scalar* p_sk = &ws_sk[0];
        scalar* p_U = &ws_U[0];
        scalar* p_Gpdx = &ws_Gpdx[0];
        scalar* p_du = &ws_du[0];
        scalar* p_rho = &ws_rho[0];
        scalar* p_F = &ws_F[0];
        scalar* p_FOrig = &ws_FOrig[0];
        scalar* p_face_shape_function = &ws_face_shape_function[0];
        scalar* p_coordinate_face_shape_function =
            &ws_coordinate_face_shape_function[0];
        scalar* p_dndx = &ws_dndx[0];
        scalar* p_Gpdx_elem = &ws_Gpdx_elem[0];
        scalar* p_F_elem = &ws_F_elem[0];
        scalar* p_FOrig_elem = &ws_FOrig_elem[0];

        // Always use trilinear (standard) shape functions for coordinates
        meFC->shape_fcn(&p_coordinate_face_shape_function[0]);

        // shape functions; boundary
        if (isUShifted)
        {
            meFC->shifted_shape_fcn(&p_face_shape_function[0]);
        }
        else
        {
            meFC->shape_fcn(&p_face_shape_function[0]);
        }

        const Bucket::size_type nSidesPerBucket = sideBucket.size();

        for (Bucket::size_type iSide = 0; iSide < nSidesPerBucket; ++iSide)
        {
            // get face
            stk::mesh::Entity side = sideBucket[iSide];

            // pointer to face data
            const scalar* areaVec =
                stk::mesh::field_data(exposedAreaVecSTKFieldRef, side);
            scalar* mDot = stk::mesh::field_data(mDotSideSTKFieldRef, side);

            // extract the connected element to this exposed face; should be
            // single in size!
            const stk::mesh::Entity* faceElemRels =
                bulkData.begin_elements(side);
            STK_ThrowAssert(bulkData.num_elements(side) == 1);

            // get element; its face ordinal number and populate
            // face_node_ordinals
            stk::mesh::Entity element = faceElemRels[0];
            const label faceOrdinal = bulkData.begin_element_ordinals(side)[0];
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

                // gather scalars
                p_p[ni] = *stk::mesh::field_data(pSTKFieldRef, node);
                if (balancedST)
                {
                    p_capPot[ni] =
                        *stk::mesh::field_data(*capPotSTKFieldPtr, node);
                    p_sk[ni] =
                        *stk::mesh::field_data(*sigKappaSTKFieldPtr, node);
                }

                // gather vectors
                scalar* coords = stk::mesh::field_data(coordsSTKFieldRef, node);
                scalar* Gjp_e = stk::mesh::field_data(gradPSTKFieldRef, node);
                const scalar* F_e = stk::mesh::field_data(*FSTKFieldPtr, node);
                const scalar* FOrig_e =
                    stk::mesh::field_data(*FOrigSTKFieldPtr, node);
                const label offSet = ni * SPATIAL_DIM;
                for (label j = 0; j < SPATIAL_DIM; ++j)
                {
                    p_coordinates[offSet + j] = coords[j];
                    p_Gpdx_elem[offSet + j] = Gjp_e[j];
                    p_F_elem[offSet + j] = F_e[j];
                    p_FOrig_elem[offSet + j] = FOrig_e[j];
                }
            }

            // projection-based harmonic mean of FOrig
            for (label d = 0; d < SPATIAL_DIM; ++d)
            {
                B_el[d] = 0.0;
            }
            for (label ni = 0; ni < nodesPerElement; ++ni)
            {
                for (label d = 0; d < SPATIAL_DIM; ++d)
                {
                    B_el[d] += p_FOrig_elem[ni * SPATIAL_DIM + d];
                }
            }
            {
                const scalar invN = 1.0 / static_cast<scalar>(nodesPerElement);
                for (label d = 0; d < SPATIAL_DIM; ++d)
                {
                    B_el[d] *= invN;
                }
            }
            scalar mag = 0.0;
            for (label d = 0; d < SPATIAL_DIM; ++d)
            {
                mag += B_el[d] * B_el[d];
            }
            mag = std::sqrt(mag);
            if (mag < SMALL)
            {
                for (label d = 0; d < SPATIAL_DIM; ++d)
                {
                    B_el[d] = 0.0;
                }
            }
            else
            {
                for (label d = 0; d < SPATIAL_DIM; ++d)
                {
                    p_uhat[d] = B_el[d] / mag;
                }
                scalar h = 0.0;
                for (label d = 0; d < SPATIAL_DIM; ++d)
                {
                    h += p_uhat[d] * p_FOrig_elem[d];
                }
                h = std::max(h, 0.0);
                for (label ni = 1; ni < nodesPerElement; ++ni)
                {
                    scalar dk = 0.0;
                    for (label d = 0; d < SPATIAL_DIM; ++d)
                    {
                        dk += p_uhat[d] * p_FOrig_elem[ni * SPATIAL_DIM + d];
                    }
                    dk = std::max(dk, 0.0);
                    if (h > 0.0)
                    {
                        h = static_cast<scalar>(ni + 1) * dk * h /
                            (static_cast<scalar>(ni) * dk + h);
                    }
                    else
                    {
                        h = 0.0;
                    }
                }
                for (label d = 0; d < SPATIAL_DIM; ++d)
                {
                    B_el[d] = h * p_uhat[d];
                }
            }

            //======================================
            // gather nodal data off of face
            //======================================
            stk::mesh::Entity const* sideNodeRels = bulkData.begin_nodes(side);
            label numSideNodes = bulkData.num_nodes(side);

            // sanity check on num nodes
            STK_ThrowAssert(numSideNodes == nodesPerSide);
            for (label ni = 0; ni < numSideNodes; ++ni)
            {
                const label ic = faceNodeOrdinals[ni];

                stk::mesh::Entity node = sideNodeRels[ni];

                // gather scalars
                p_rho[ni] = *stk::mesh::field_data(rhoSTKFieldRef, node);
                p_p[ic] = *stk::mesh::field_data(nodalSidePSTKFieldRef, node);

                // gather vectors
                scalar* U = stk::mesh::field_data(USTKFieldRef, node);
                scalar* Gjp = stk::mesh::field_data(gradPSTKFieldRef, node);
                const scalar* du = stk::mesh::field_data(duSTKFieldRef, node);

                const label offSet = ni * SPATIAL_DIM;
                for (label j = 0; j < SPATIAL_DIM; ++j)
                {
                    p_U[offSet + j] = U[j];
                    p_Gpdx[offSet + j] = Gjp[j];
                    p_du[offSet + j] = du[j];
                }

                // gather body force vectors for buoyancy stabilization
                const scalar* F = stk::mesh::field_data(*FSTKFieldPtr, node);
                const scalar* FOrig =
                    stk::mesh::field_data(*FOrigSTKFieldPtr, node);
                for (label j = 0; j < SPATIAL_DIM; ++j)
                {
                    p_F[offSet + j] = F[j];
                    p_FOrig[offSet + j] = FOrig[j];
                }

                // NOTE: [2024-11-25] Correction uses
                // computed pressure values for side nodes, not actual boundary
                // condition values (after discussion with Mahdi).
            }

            // compute dndx
            scalar scs_error = 0.0;
            if (isPGradientShifted)
            {
                meSCS->shifted_face_grad_op(1,
                                            faceOrdinal,
                                            &p_coordinates[0],
                                            &p_dndx[0],
                                            &ws_det_j[0],
                                            &scs_error);
            }
            else
            {
                meSCS->face_grad_op(1,
                                    faceOrdinal,
                                    &p_coordinates[0],
                                    &p_dndx[0],
                                    &ws_det_j[0],
                                    &scs_error);
            }

            // loop over boundary ips
            for (label ip = 0; ip < numScsBip; ++ip)
            {
                // per-ip blend: local face node vs opposing interior node
                const label localFaceNode = faceIpNodeMap[ip];
                const label opposingNode =
                    meSCS->opposingNodes(faceOrdinal, ip);

                // zero out vector quantities; form aMag
                scalar aMag = 0.0;
                for (label j = 0; j < SPATIAL_DIM; ++j)
                {
                    p_coordBip[j] = 0.0;
                    p_uBip[j] = 0.0;
                    // pressure gradient: harmonic for Gpdx
                    const scalar gFace =
                        p_Gpdx[localFaceNode * SPATIAL_DIM + j];
                    const scalar gOpp =
                        p_Gpdx_elem[opposingNode * SPATIAL_DIM + j];
                    p_GpdxBip[j] =
                        (std::abs(gFace) * gOpp + gFace * std::abs(gOpp)) /
                        (std::abs(gFace) + std::abs(gOpp) + SMALL);
                    p_dpdxBip[j] = 0.0;
                    p_duBip[j] = 0.0;
                    // body force: harmonic for F
                    const scalar fFace = p_F[localFaceNode * SPATIAL_DIM + j];
                    const scalar fOpp =
                        p_F_elem[opposingNode * SPATIAL_DIM + j];
                    p_FBip[j] =
                        (std::abs(fFace) * fOpp + fFace * std::abs(fOpp)) /
                        (std::abs(fFace) + std::abs(fOpp) + SMALL);
                    p_FOrigBip[j] = B_el[j];
                }

                // interpolate to bip
                scalar rhoBip = 0;
                scalar skBip = 0.0;
                const label offSetSF_face = ip * nodesPerSide;
                for (label ic = 0; ic < nodesPerSide; ++ic)
                {
                    const label inn = faceNodeOrdinals[ic];
                    const scalar r = p_face_shape_function[offSetSF_face + ic];
                    const scalar r_coord =
                        p_coordinate_face_shape_function[offSetSF_face + ic];

                    rhoBip += r * p_rho[ic];
                    if (balancedST)
                        skBip += r * p_sk[inn];

                    const label icNdim = ic * SPATIAL_DIM;
                    for (label j = 0; j < SPATIAL_DIM; ++j)
                    {
                        p_coordBip[j] +=
                            r_coord * p_coordinates[inn * SPATIAL_DIM + j];
                        p_uBip[j] += r * p_U[icNdim + j];
                        p_duBip[j] += r * p_du[icNdim + j];
                    }
                }

                // form dpdxBip
                for (label ic = 0; ic < nodesPerElement; ++ic)
                {
                    const label offSetDnDx =
                        SPATIAL_DIM * nodesPerElement * ip + ic * SPATIAL_DIM;
                    const scalar pIc = p_p[ic];
                    for (label j = 0; j < SPATIAL_DIM; ++j)
                    {
                        p_dpdxBip[j] += p_dndx[offSetDnDx + j] * pIc;
                        // balanced CSF: sigma*kappa_ip * compact grad(alpha)
                        if (balancedST)
                            p_FOrigBip[j] +=
                                skBip * p_dndx[offSetDnDx + j] * p_capPot[ic];
                    }
                }

                // form mDot:
                // rho*uj*Aj - rho*du*(dpdxj - Gjp)*Aj + rho*du*(FOrigj - Fj)*Aj
                scalar tmDot = 0.0;
                for (label j = 0; j < SPATIAL_DIM; ++j)
                {
                    const scalar axj = areaVec[ip * SPATIAL_DIM + j];

                    // divergence + pressure Rhie-Chow
                    tmDot +=
                        (rhoBip * p_uBip[j] -
                         rhoBip * p_duBip[j] * (p_dpdxBip[j] - p_GpdxBip[j])) *
                        axj;

                    // buoyancy stabilization: +rho*D*(F_orig - F)·S
                    tmDot +=
                        rhoBip * p_duBip[j] * (p_FOrigBip[j] - p_FBip[j]) * axj;
                }

                // store with relaxation
                mDot[ip] = mDotURF * tmDot + (1.0 - mDotURF) * mDot[ip];
            }
        }
    }
}

void freeSurfaceFlowModel::
    updateMassFlowRateBoundaryFieldOutletSpecifiedPressure_(
        const std::shared_ptr<domain> domain,
        const boundary* boundary,
        sideField<scalar, 1>& mDotSideField,
        const nodeField<1, SPATIAL_DIM>& rhoField)
{
    using stk::mesh::Bucket;
    using stk::mesh::BucketVector;

    const auto& mesh = this->meshRef();
    const stk::mesh::MetaData& metaData = mesh.metaDataRef();
    const stk::mesh::BulkData& bulkData = mesh.bulkDataRef();

    // ip values; both boundary and opposing
    // surface
    std::vector<scalar> coordBip(SPATIAL_DIM);
    std::vector<scalar> uBip(SPATIAL_DIM);
    std::vector<scalar> GpdxBip(SPATIAL_DIM);
    std::vector<scalar> dpdxBip(SPATIAL_DIM);
    std::vector<scalar> duBip(SPATIAL_DIM);
    std::vector<scalar> FBip(SPATIAL_DIM);
    std::vector<scalar> FOrigBip(SPATIAL_DIM);
    std::vector<scalar> uhat(SPATIAL_DIM);

    // pointers to fixed values
    scalar* p_coordBip = &coordBip[0];
    scalar* p_uBip = &uBip[0];
    scalar* p_GpdxBip = &GpdxBip[0];
    scalar* p_dpdxBip = &dpdxBip[0];
    scalar* p_duBip = &duBip[0];
    scalar* p_FBip = &FBip[0];
    scalar* p_FOrigBip = &FOrigBip[0];
    scalar* p_uhat = &uhat[0];

    // nodal fields to gather
    std::vector<scalar> ws_coordinates;
    std::vector<scalar> ws_p;
    std::vector<scalar> ws_U;
    std::vector<scalar> ws_Gpdx;
    std::vector<scalar> ws_du;
    std::vector<scalar> ws_rho;
    std::vector<scalar> ws_F;
    std::vector<scalar> ws_FOrig;
    std::vector<scalar> ws_Gpdx_elem;
    std::vector<scalar> ws_F_elem;
    std::vector<scalar> ws_FOrig_elem;
    std::vector<scalar> ws_capPot;
    std::vector<scalar> ws_sk;
    std::vector<scalar> B_el(SPATIAL_DIM);

    // master element
    std::vector<scalar> ws_face_shape_function;
    std::vector<scalar> ws_coordinate_face_shape_function;
    std::vector<scalar> ws_dndx;
    std::vector<scalar> ws_det_j;

    // Get fields
    const auto& rhoSTKFieldRef = rhoField.stkFieldRef();
    const auto& pSTKFieldRef = this->pRef().stkFieldRef();
    const auto& USTKFieldRef = this->URef().stkFieldRef();
    const auto& gradPSTKFieldRef = this->pRef().gradRef().stkFieldRef();
    const auto& nodalSidePSTKFieldRef =
        this->pRef().nodeSideFieldRef().stkFieldRef();

    const auto& mDotSideSTKFieldRef = mDotSideField.stkFieldRef();
    const auto& reversalFlowFlagSTKFieldRef =
        this->URef().reversalFlagRef().stkFieldRef();
    const scalar mDotURF = mDotSideField.urf();

    const auto& duSTKFieldRef =
        *metaData.get_field<scalar>(stk::topology::NODE_RANK, flowModel::du_ID);

    // Get body force fields for buoyancy pressure stabilization
    const auto* FSTKFieldPtr =
        metaData.get_field<scalar>(stk::topology::NODE_RANK, flowModel::F_ID);
    const auto* FOrigSTKFieldPtr = metaData.get_field<scalar>(
        stk::topology::NODE_RANK, flowModel::FOriginal_ID);

    const auto* capPotSTKFieldPtr = metaData.get_field<scalar>(
        stk::topology::NODE_RANK, "capillary_potential");
    const auto* sigKappaSTKFieldPtr = metaData.get_field<scalar>(
        stk::topology::NODE_RANK, "capillary_sigma_kappa");
    const bool balancedST =
        (capPotSTKFieldPtr != nullptr && sigKappaSTKFieldPtr != nullptr);

    // Get geometric fields
    const auto& coordsSTKFieldRef = *metaData.get_field<scalar>(
        stk::topology::NODE_RANK, this->getCoordinatesID_(domain));
    const auto& exposedAreaVecSTKFieldRef = *metaData.get_field<scalar>(
        metaData.side_rank(), this->getExposedAreaVectorID_(domain));

    // define vector of parent topos; should
    // always be UNITY in size
    std::vector<stk::topology> parentTopo;

    // define some common selectors
    stk::mesh::Selector selAllSides =
        metaData.universal_part() & stk::mesh::selectUnion(boundary->parts());

    // shifted ip's for field?
    const bool isUShifted = this->URef().isShifted();

    // shifted ip's for gradients?
    const bool isPGradientShifted = this->pRef().isGradientShifted();

    // free-surface always uses pure harmonic (no comp switch)

    BucketVector const& sideBuckets =
        bulkData.get_buckets(metaData.side_rank(), selAllSides);
    for (BucketVector::const_iterator ib = sideBuckets.begin();
         ib != sideBuckets.end();
         ++ib)
    {
        Bucket& sideBucket = **ib;

        // extract connected element
        // topology
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

        // algorithm related; element
        // (exposed face and element)
        ws_coordinates.resize(nodesPerElement * SPATIAL_DIM);
        ws_p.resize(nodesPerElement);
        ws_capPot.resize(nodesPerElement);
        ws_sk.resize(nodesPerElement);
        ws_U.resize(nodesPerSide * SPATIAL_DIM);
        ws_Gpdx.resize(nodesPerSide * SPATIAL_DIM);
        ws_du.resize(nodesPerSide * SPATIAL_DIM);
        ws_rho.resize(nodesPerSide);
        ws_F.resize(nodesPerSide * SPATIAL_DIM);
        ws_FOrig.resize(nodesPerSide * SPATIAL_DIM);
        ws_Gpdx_elem.resize(nodesPerElement * SPATIAL_DIM);
        ws_F_elem.resize(nodesPerElement * SPATIAL_DIM);
        ws_FOrig_elem.resize(nodesPerElement * SPATIAL_DIM);
        ws_face_shape_function.resize(numScsBip * nodesPerSide);
        ws_coordinate_face_shape_function.resize(numScsBip * nodesPerSide);
        ws_dndx.resize(SPATIAL_DIM * numScsBip * nodesPerElement);
        ws_det_j.resize(numScsBip);

        // pointers
        scalar* p_coordinates = &ws_coordinates[0];
        scalar* p_p = &ws_p[0];
        scalar* p_capPot = &ws_capPot[0];
        scalar* p_sk = &ws_sk[0];
        scalar* p_U = &ws_U[0];
        scalar* p_Gpdx = &ws_Gpdx[0];
        scalar* p_du = &ws_du[0];
        scalar* p_rho = &ws_rho[0];
        scalar* p_F = &ws_F[0];
        scalar* p_FOrig = &ws_FOrig[0];
        scalar* p_Gpdx_elem = &ws_Gpdx_elem[0];
        scalar* p_F_elem = &ws_F_elem[0];
        scalar* p_FOrig_elem = &ws_FOrig_elem[0];
        scalar* p_face_shape_function = &ws_face_shape_function[0];
        scalar* p_coordinate_face_shape_function =
            &ws_coordinate_face_shape_function[0];
        scalar* p_dndx = &ws_dndx[0];

        // Always use trilinear (standard) shape functions for coordinates
        meFC->shape_fcn(&p_coordinate_face_shape_function[0]);

        const label* faceIpNodeMap = meFC->ipNodeMap();

        // shape functions; boundary
        if (isUShifted)
        {
            meFC->shifted_shape_fcn(&p_face_shape_function[0]);
        }
        else
        {
            meFC->shape_fcn(&p_face_shape_function[0]);
        }

        const Bucket::size_type nSidesPerBucket = sideBucket.size();

        for (Bucket::size_type iSide = 0; iSide < nSidesPerBucket; ++iSide)
        {
            // get face
            stk::mesh::Entity side = sideBucket[iSide];

            // pointer to face data
            const scalar* areaVec =
                stk::mesh::field_data(exposedAreaVecSTKFieldRef, side);
            scalar* mDot = stk::mesh::field_data(mDotSideSTKFieldRef, side);
            const label* rfflag =
                stk::mesh::field_data(reversalFlowFlagSTKFieldRef, side);

            // extract the connected element to this exposed face; should be
            // single in size!
            const stk::mesh::Entity* faceElemRels =
                bulkData.begin_elements(side);
            STK_ThrowAssert(bulkData.num_elements(side) == 1);

            // get element; its face ordinal number and populate
            // face_node_ordinals
            stk::mesh::Entity element = faceElemRels[0];
            const label faceOrdinal = bulkData.begin_element_ordinals(side)[0];
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

                // gather scalars
                p_p[ni] = *stk::mesh::field_data(pSTKFieldRef, node);
                if (balancedST)
                {
                    p_capPot[ni] =
                        *stk::mesh::field_data(*capPotSTKFieldPtr, node);
                    p_sk[ni] =
                        *stk::mesh::field_data(*sigKappaSTKFieldPtr, node);
                }

                // gather vectors
                scalar* coords = stk::mesh::field_data(coordsSTKFieldRef, node);
                scalar* Gjp_e = stk::mesh::field_data(gradPSTKFieldRef, node);
                const scalar* F_e = stk::mesh::field_data(*FSTKFieldPtr, node);
                const scalar* FOrig_e =
                    stk::mesh::field_data(*FOrigSTKFieldPtr, node);
                const label offSet = ni * SPATIAL_DIM;
                for (label j = 0; j < SPATIAL_DIM; ++j)
                {
                    p_coordinates[offSet + j] = coords[j];
                    p_Gpdx_elem[offSet + j] = Gjp_e[j];
                    p_F_elem[offSet + j] = F_e[j];
                    p_FOrig_elem[offSet + j] = FOrig_e[j];
                }
            }

            // projection-based harmonic mean of FOrig
            for (label d = 0; d < SPATIAL_DIM; ++d)
            {
                B_el[d] = 0.0;
            }
            for (label ni = 0; ni < nodesPerElement; ++ni)
            {
                for (label d = 0; d < SPATIAL_DIM; ++d)
                {
                    B_el[d] += p_FOrig_elem[ni * SPATIAL_DIM + d];
                }
            }
            {
                const scalar invN = 1.0 / static_cast<scalar>(nodesPerElement);
                for (label d = 0; d < SPATIAL_DIM; ++d)
                {
                    B_el[d] *= invN;
                }
            }
            scalar mag = 0.0;
            for (label d = 0; d < SPATIAL_DIM; ++d)
            {
                mag += B_el[d] * B_el[d];
            }
            mag = std::sqrt(mag);
            if (mag < SMALL)
            {
                for (label d = 0; d < SPATIAL_DIM; ++d)
                {
                    B_el[d] = 0.0;
                }
            }
            else
            {
                for (label d = 0; d < SPATIAL_DIM; ++d)
                {
                    p_uhat[d] = B_el[d] / mag;
                }
                scalar h = 0.0;
                for (label d = 0; d < SPATIAL_DIM; ++d)
                {
                    h += p_uhat[d] * p_FOrig_elem[d];
                }
                h = std::max(h, 0.0);
                for (label ni = 1; ni < nodesPerElement; ++ni)
                {
                    scalar dk = 0.0;
                    for (label d = 0; d < SPATIAL_DIM; ++d)
                    {
                        dk += p_uhat[d] * p_FOrig_elem[ni * SPATIAL_DIM + d];
                    }
                    dk = std::max(dk, 0.0);
                    if (h > 0.0)
                    {
                        h = static_cast<scalar>(ni + 1) * dk * h /
                            (static_cast<scalar>(ni) * dk + h);
                    }
                    else
                    {
                        h = 0.0;
                    }
                }
                for (label d = 0; d < SPATIAL_DIM; ++d)
                {
                    B_el[d] = h * p_uhat[d];
                }
            }

            //======================================
            // gather nodal data off of face
            //======================================
            stk::mesh::Entity const* sideNodeRels = bulkData.begin_nodes(side);
            label numSideNodes = bulkData.num_nodes(side);

            // sanity check on num nodes
            STK_ThrowAssert(numSideNodes == nodesPerSide);
            for (label ni = 0; ni < numSideNodes; ++ni)
            {
                const label ic = faceNodeOrdinals[ni];

                stk::mesh::Entity node = sideNodeRels[ni];

                // gather scalars
                p_rho[ni] = *stk::mesh::field_data(rhoSTKFieldRef, node);
                p_p[ic] = *stk::mesh::field_data(nodalSidePSTKFieldRef, node);

                // gather vectors
                scalar* U = stk::mesh::field_data(USTKFieldRef, node);
                scalar* Gjp = stk::mesh::field_data(gradPSTKFieldRef, node);
                const scalar* du = stk::mesh::field_data(duSTKFieldRef, node);

                const label offSet = ni * SPATIAL_DIM;
                for (label j = 0; j < SPATIAL_DIM; ++j)
                {
                    p_U[offSet + j] = U[j];
                    p_Gpdx[offSet + j] = Gjp[j];
                    p_du[offSet + j] = du[j];
                }

                // gather body force vectors for buoyancy stabilization
                const scalar* F = stk::mesh::field_data(*FSTKFieldPtr, node);
                const scalar* FOrig =
                    stk::mesh::field_data(*FOrigSTKFieldPtr, node);
                for (label j = 0; j < SPATIAL_DIM; ++j)
                {
                    p_F[offSet + j] = F[j];
                    p_FOrig[offSet + j] = FOrig[j];
                }

                // NOTE: Correction uses
                // computed pressure values for side nodes, not actual boundary
                // condition values (after discussion with Mahdi).
            }

            // compute dndx
            scalar scs_error = 0.0;
            if (isPGradientShifted)
            {
                meSCS->shifted_face_grad_op(1,
                                            faceOrdinal,
                                            &p_coordinates[0],
                                            &p_dndx[0],
                                            &ws_det_j[0],
                                            &scs_error);
            }
            else
            {
                meSCS->face_grad_op(1,
                                    faceOrdinal,
                                    &p_coordinates[0],
                                    &p_dndx[0],
                                    &ws_det_j[0],
                                    &scs_error);
            }

            // loop over boundary ips
            for (label ip = 0; ip < numScsBip; ++ip)
            {
                if (rfflag[ip] == 1) // slip-wall when reversed
                {
                    mDot[ip] = 0.0;
                    continue;
                }

                const label localFaceNode = faceIpNodeMap[ip];
                const label opposingNode =
                    meSCS->opposingNodes(faceOrdinal, ip);

                // zero out vector quantities; form aMag
                for (label j = 0; j < SPATIAL_DIM; ++j)
                {
                    p_coordBip[j] = 0.0;
                    p_uBip[j] = 0.0;
                    // pressure gradient: harmonic for Gpdx
                    const scalar gFace =
                        p_Gpdx[localFaceNode * SPATIAL_DIM + j];
                    const scalar gOpp =
                        p_Gpdx_elem[opposingNode * SPATIAL_DIM + j];
                    p_GpdxBip[j] =
                        (std::abs(gFace) * gOpp + gFace * std::abs(gOpp)) /
                        (std::abs(gFace) + std::abs(gOpp) + SMALL);
                    p_dpdxBip[j] = 0.0;
                    p_duBip[j] = 0.0;
                    // body force: harmonic for F
                    const scalar fFace = p_F[localFaceNode * SPATIAL_DIM + j];
                    const scalar fOpp =
                        p_F_elem[opposingNode * SPATIAL_DIM + j];
                    p_FBip[j] =
                        (std::abs(fFace) * fOpp + fFace * std::abs(fOpp)) /
                        (std::abs(fFace) + std::abs(fOpp) + SMALL);
                    p_FOrigBip[j] = B_el[j];
                }

                // interpolate to bip
                scalar rhoBip = 0;
                scalar skBip = 0.0;
                const label offSetSF_face = ip * nodesPerSide;
                for (label ic = 0; ic < nodesPerSide; ++ic)
                {
                    const label inn = faceNodeOrdinals[ic];

                    const scalar r = p_face_shape_function[offSetSF_face + ic];
                    const scalar r_coord =
                        p_coordinate_face_shape_function[offSetSF_face + ic];

                    rhoBip += r * p_rho[ic];

                    const label icNdim = ic * SPATIAL_DIM;
                    for (label j = 0; j < SPATIAL_DIM; ++j)
                    {
                        p_coordBip[j] +=
                            r_coord * p_coordinates[inn * SPATIAL_DIM + j];
                        p_uBip[j] += r * p_U[icNdim + j];
                        p_duBip[j] += r * p_du[icNdim + j];
                    }
                }

                // form dpdxBip
                for (label ic = 0; ic < nodesPerElement; ++ic)
                {
                    const label offSetDnDx =
                        SPATIAL_DIM * nodesPerElement * ip + ic * SPATIAL_DIM;
                    const scalar pIc = p_p[ic];
                    for (label j = 0; j < SPATIAL_DIM; ++j)
                    {
                        p_dpdxBip[j] += p_dndx[offSetDnDx + j] * pIc;
                        // balanced CSF: sigma*kappa_ip * compact grad(alpha)
                        if (balancedST)
                            p_FOrigBip[j] +=
                                skBip * p_dndx[offSetDnDx + j] * p_capPot[ic];
                    }
                }

                // form mDot:
                // rho*uj*Aj - rho*du*(dpdxj - Gjp)*Aj + rho*du*(FOrigj - Fj)*Aj
                scalar tmDot = 0.0;
                for (label j = 0; j < SPATIAL_DIM; ++j)
                {
                    const scalar axj = areaVec[ip * SPATIAL_DIM + j];

                    // divergence + pressure Rhie-Chow
                    tmDot +=
                        (rhoBip * p_uBip[j] -
                         rhoBip * p_duBip[j] * (p_dpdxBip[j] - p_GpdxBip[j])) *
                        axj;

                    // buoyancy stabilization: +rho*D*(F_orig - F)·S
                    tmDot +=
                        rhoBip * p_duBip[j] * (p_FOrigBip[j] - p_FBip[j]) * axj;
                }

                // store with relaxation
                mDot[ip] = mDotURF * tmDot + (1.0 - mDotURF) * mDot[ip];
            }
        }
    }
}

void freeSurfaceFlowModel::updateMassFlowRateBoundaryFieldOpeningPressure_(
    const std::shared_ptr<domain> domain,
    const boundary* boundary,
    sideField<scalar, 1>& mDotSideField,
    const nodeField<1, SPATIAL_DIM>& rhoField)
{
    using stk::mesh::Bucket;
    using stk::mesh::BucketVector;

    const auto& mesh = this->meshRef();
    const stk::mesh::MetaData& metaData = mesh.metaDataRef();
    const stk::mesh::BulkData& bulkData = mesh.bulkDataRef();

    // ip values; both boundary and opposing
    // surface
    std::vector<scalar> coordBip(SPATIAL_DIM);
    std::vector<scalar> uBip(SPATIAL_DIM);
    std::vector<scalar> GpdxBip(SPATIAL_DIM);
    std::vector<scalar> dpdxBip(SPATIAL_DIM);
    std::vector<scalar> duBip(SPATIAL_DIM);
    std::vector<scalar> FBip(SPATIAL_DIM);
    std::vector<scalar> FOrigBip(SPATIAL_DIM);
    std::vector<scalar> uhat(SPATIAL_DIM);

    // pointers to fixed values
    scalar* p_coordBip = &coordBip[0];
    scalar* p_uBip = &uBip[0];
    scalar* p_GpdxBip = &GpdxBip[0];
    scalar* p_dpdxBip = &dpdxBip[0];
    scalar* p_duBip = &duBip[0];
    scalar* p_FBip = &FBip[0];
    scalar* p_FOrigBip = &FOrigBip[0];
    scalar* p_uhat = &uhat[0];

    // nodal fields to gather
    std::vector<scalar> ws_coordinates;
    std::vector<scalar> ws_p;
    std::vector<scalar> ws_U;
    std::vector<scalar> ws_Gpdx;
    std::vector<scalar> ws_du;
    std::vector<scalar> ws_rho;
    std::vector<scalar> ws_F;
    std::vector<scalar> ws_FOrig;
    std::vector<scalar> ws_Gpdx_elem;
    std::vector<scalar> ws_F_elem;
    std::vector<scalar> ws_FOrig_elem;
    std::vector<scalar> ws_capPot;
    std::vector<scalar> ws_sk;
    std::vector<scalar> B_el(SPATIAL_DIM);

    // master element
    std::vector<scalar> ws_face_shape_function;
    std::vector<scalar> ws_coordinate_face_shape_function;
    std::vector<scalar> ws_dndx;
    std::vector<scalar> ws_det_j;

    // Get fields
    const auto& rhoSTKFieldRef = rhoField.stkFieldRef();
    const auto& pSTKFieldRef = this->pRef().stkFieldRef();
    const auto& USTKFieldRef = this->URef().stkFieldRef();
    const auto& gradPSTKFieldRef = this->pRef().gradRef().stkFieldRef();

    const auto& mDotSideSTKFieldRef = mDotSideField.stkFieldRef();
    const scalar mDotURF = mDotSideField.urf();

    const auto& duSTKFieldRef =
        *metaData.get_field<scalar>(stk::topology::NODE_RANK, flowModel::du_ID);

    // Get body force fields for buoyancy pressure stabilization
    const auto* FSTKFieldPtr =
        metaData.get_field<scalar>(stk::topology::NODE_RANK, flowModel::F_ID);
    const auto* FOrigSTKFieldPtr = metaData.get_field<scalar>(
        stk::topology::NODE_RANK, flowModel::FOriginal_ID);

    const auto* capPotSTKFieldPtr = metaData.get_field<scalar>(
        stk::topology::NODE_RANK, "capillary_potential");
    const auto* sigKappaSTKFieldPtr = metaData.get_field<scalar>(
        stk::topology::NODE_RANK, "capillary_sigma_kappa");
    const bool balancedST =
        (capPotSTKFieldPtr != nullptr && sigKappaSTKFieldPtr != nullptr);

    // Get geometric fields
    const auto& coordsSTKFieldRef = *metaData.get_field<scalar>(
        stk::topology::NODE_RANK, this->getCoordinatesID_(domain));
    const auto& exposedAreaVecSTKFieldRef = *metaData.get_field<scalar>(
        metaData.side_rank(), this->getExposedAreaVectorID_(domain));

    // define vector of parent topos; should
    // always be UNITY in size
    std::vector<stk::topology> parentTopo;

    // define some common selectors
    stk::mesh::Selector selAllSides =
        metaData.universal_part() & stk::mesh::selectUnion(boundary->parts());

    // shifted ip's for field?
    const bool isUShifted = this->URef().isShifted();

    // shifted ip's for gradients?
    const bool isPGradientShifted = this->pRef().isGradientShifted();

    // free-surface always uses pure harmonic (no comp switch)

    BucketVector const& sideBuckets =
        bulkData.get_buckets(metaData.side_rank(), selAllSides);
    for (BucketVector::const_iterator ib = sideBuckets.begin();
         ib != sideBuckets.end();
         ++ib)
    {
        Bucket& sideBucket = **ib;

        // extract connected element
        // topology
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

        // algorithm related; element
        // (exposed face and element)
        ws_coordinates.resize(nodesPerElement * SPATIAL_DIM);
        ws_p.resize(nodesPerElement);
        ws_capPot.resize(nodesPerElement);
        ws_sk.resize(nodesPerElement);
        ws_U.resize(nodesPerSide * SPATIAL_DIM);
        ws_Gpdx.resize(nodesPerSide * SPATIAL_DIM);
        ws_du.resize(nodesPerSide * SPATIAL_DIM);
        ws_rho.resize(nodesPerSide);
        ws_F.resize(nodesPerSide * SPATIAL_DIM);
        ws_FOrig.resize(nodesPerSide * SPATIAL_DIM);
        ws_Gpdx_elem.resize(nodesPerElement * SPATIAL_DIM);
        ws_F_elem.resize(nodesPerElement * SPATIAL_DIM);
        ws_FOrig_elem.resize(nodesPerElement * SPATIAL_DIM);
        ws_face_shape_function.resize(numScsBip * nodesPerSide);
        ws_coordinate_face_shape_function.resize(numScsBip * nodesPerSide);
        ws_dndx.resize(SPATIAL_DIM * numScsBip * nodesPerElement);
        ws_det_j.resize(numScsBip);

        // pointers
        scalar* p_coordinates = &ws_coordinates[0];
        scalar* p_p = &ws_p[0];
        scalar* p_capPot = &ws_capPot[0];
        scalar* p_sk = &ws_sk[0];
        scalar* p_U = &ws_U[0];
        scalar* p_Gpdx = &ws_Gpdx[0];
        scalar* p_du = &ws_du[0];
        scalar* p_rho = &ws_rho[0];
        scalar* p_F = &ws_F[0];
        scalar* p_FOrig = &ws_FOrig[0];
        scalar* p_Gpdx_elem = &ws_Gpdx_elem[0];
        scalar* p_F_elem = &ws_F_elem[0];
        scalar* p_FOrig_elem = &ws_FOrig_elem[0];
        scalar* p_face_shape_function = &ws_face_shape_function[0];
        scalar* p_coordinate_face_shape_function =
            &ws_coordinate_face_shape_function[0];
        scalar* p_dndx = &ws_dndx[0];

        // shape functions; boundary
        if (isUShifted)
        {
            meFC->shifted_shape_fcn(&p_face_shape_function[0]);
        }
        else
        {
            meFC->shape_fcn(&p_face_shape_function[0]);
        }

        // Always use trilinear (standard) shape functions for coordinates
        meFC->shape_fcn(&p_coordinate_face_shape_function[0]);

        const label* faceIpNodeMap = meFC->ipNodeMap();

        const Bucket::size_type nSidesPerBucket = sideBucket.size();

        for (Bucket::size_type iSide = 0; iSide < nSidesPerBucket; ++iSide)
        {
            // get face
            stk::mesh::Entity side = sideBucket[iSide];

            // pointer to face data
            const scalar* areaVec =
                stk::mesh::field_data(exposedAreaVecSTKFieldRef, side);
            scalar* mDot = stk::mesh::field_data(mDotSideSTKFieldRef, side);

            // extract the connected element to this exposed face; should be
            // single in size!
            const stk::mesh::Entity* faceElemRels =
                bulkData.begin_elements(side);
            STK_ThrowAssert(bulkData.num_elements(side) == 1);

            // get element; its face ordinal number and populate
            // face_node_ordinals
            stk::mesh::Entity element = faceElemRels[0];
            const label faceOrdinal = bulkData.begin_element_ordinals(side)[0];
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

                // gather scalars
                p_p[ni] = *stk::mesh::field_data(pSTKFieldRef, node);
                if (balancedST)
                {
                    p_capPot[ni] =
                        *stk::mesh::field_data(*capPotSTKFieldPtr, node);
                    p_sk[ni] =
                        *stk::mesh::field_data(*sigKappaSTKFieldPtr, node);
                }

                // gather vectors
                scalar* coords = stk::mesh::field_data(coordsSTKFieldRef, node);
                scalar* Gjp_e = stk::mesh::field_data(gradPSTKFieldRef, node);
                const scalar* F_e = stk::mesh::field_data(*FSTKFieldPtr, node);
                const scalar* FOrig_e =
                    stk::mesh::field_data(*FOrigSTKFieldPtr, node);
                const label offSet = ni * SPATIAL_DIM;
                for (label j = 0; j < SPATIAL_DIM; ++j)
                {
                    p_coordinates[offSet + j] = coords[j];
                    p_Gpdx_elem[offSet + j] = Gjp_e[j];
                    p_F_elem[offSet + j] = F_e[j];
                    p_FOrig_elem[offSet + j] = FOrig_e[j];
                }
            }

            // projection-based harmonic mean of FOrig
            for (label d = 0; d < SPATIAL_DIM; ++d)
            {
                B_el[d] = 0.0;
            }
            for (label ni = 0; ni < nodesPerElement; ++ni)
            {
                for (label d = 0; d < SPATIAL_DIM; ++d)
                {
                    B_el[d] += p_FOrig_elem[ni * SPATIAL_DIM + d];
                }
            }
            {
                const scalar invN = 1.0 / static_cast<scalar>(nodesPerElement);
                for (label d = 0; d < SPATIAL_DIM; ++d)
                {
                    B_el[d] *= invN;
                }
            }
            scalar mag = 0.0;
            for (label d = 0; d < SPATIAL_DIM; ++d)
            {
                mag += B_el[d] * B_el[d];
            }
            mag = std::sqrt(mag);
            if (mag < SMALL)
            {
                for (label d = 0; d < SPATIAL_DIM; ++d)
                {
                    B_el[d] = 0.0;
                }
            }
            else
            {
                for (label d = 0; d < SPATIAL_DIM; ++d)
                {
                    p_uhat[d] = B_el[d] / mag;
                }
                scalar h = 0.0;
                for (label d = 0; d < SPATIAL_DIM; ++d)
                {
                    h += p_uhat[d] * p_FOrig_elem[d];
                }
                h = std::max(h, 0.0);
                for (label ni = 1; ni < nodesPerElement; ++ni)
                {
                    scalar dk = 0.0;
                    for (label d = 0; d < SPATIAL_DIM; ++d)
                    {
                        dk += p_uhat[d] * p_FOrig_elem[ni * SPATIAL_DIM + d];
                    }
                    dk = std::max(dk, 0.0);
                    if (h > 0.0)
                    {
                        h = static_cast<scalar>(ni + 1) * dk * h /
                            (static_cast<scalar>(ni) * dk + h);
                    }
                    else
                    {
                        h = 0.0;
                    }
                }
                for (label d = 0; d < SPATIAL_DIM; ++d)
                {
                    B_el[d] = h * p_uhat[d];
                }
            }

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

                // gather vectors
                scalar* U = stk::mesh::field_data(USTKFieldRef, node);
                scalar* Gjp = stk::mesh::field_data(gradPSTKFieldRef, node);
                const scalar* du = stk::mesh::field_data(duSTKFieldRef, node);

                const label offSet = ni * SPATIAL_DIM;
                for (label j = 0; j < SPATIAL_DIM; ++j)
                {
                    p_U[offSet + j] = U[j];
                    p_Gpdx[offSet + j] = Gjp[j];
                    p_du[offSet + j] = du[j];
                }

                // gather body force vectors for buoyancy stabilization
                const scalar* F = stk::mesh::field_data(*FSTKFieldPtr, node);
                const scalar* FOrig =
                    stk::mesh::field_data(*FOrigSTKFieldPtr, node);
                for (label j = 0; j < SPATIAL_DIM; ++j)
                {
                    p_F[offSet + j] = F[j];
                    p_FOrig[offSet + j] = FOrig[j];
                }

                // NOTE: Correction uses
                // computed pressure values for side nodes, not actual boundary
                // condition values (after discussion with Mahdi).
            }

            // compute dndx
            scalar scs_error = 0.0;
            if (isPGradientShifted)
            {
                meSCS->shifted_face_grad_op(1,
                                            faceOrdinal,
                                            &p_coordinates[0],
                                            &p_dndx[0],
                                            &ws_det_j[0],
                                            &scs_error);
            }
            else
            {
                meSCS->face_grad_op(1,
                                    faceOrdinal,
                                    &p_coordinates[0],
                                    &p_dndx[0],
                                    &ws_det_j[0],
                                    &scs_error);
            }

            // loop over boundary ips
            for (label ip = 0; ip < numScsBip; ++ip)
            {
                const label localFaceNode = faceIpNodeMap[ip];
                const label opposingNode =
                    meSCS->opposingNodes(faceOrdinal, ip);

                // zero out vector quantities; form aMag
                for (label j = 0; j < SPATIAL_DIM; ++j)
                {
                    p_coordBip[j] = 0.0;
                    p_uBip[j] = 0.0;
                    // pressure gradient: harmonic for Gpdx
                    const scalar gFace =
                        p_Gpdx[localFaceNode * SPATIAL_DIM + j];
                    const scalar gOpp =
                        p_Gpdx_elem[opposingNode * SPATIAL_DIM + j];
                    p_GpdxBip[j] =
                        (std::abs(gFace) * gOpp + gFace * std::abs(gOpp)) /
                        (std::abs(gFace) + std::abs(gOpp) + SMALL);
                    p_dpdxBip[j] = 0.0;
                    p_duBip[j] = 0.0;
                    // body force: harmonic for F
                    const scalar fFace = p_F[localFaceNode * SPATIAL_DIM + j];
                    const scalar fOpp =
                        p_F_elem[opposingNode * SPATIAL_DIM + j];
                    p_FBip[j] =
                        (std::abs(fFace) * fOpp + fFace * std::abs(fOpp)) /
                        (std::abs(fFace) + std::abs(fOpp) + SMALL);
                    p_FOrigBip[j] = B_el[j];
                }

                // interpolate to bip
                scalar rhoBip = 0;
                scalar skBip = 0.0;
                const label offSetSF_face = ip * nodesPerSide;
                for (label ic = 0; ic < nodesPerSide; ++ic)
                {
                    const label inn = faceNodeOrdinals[ic];
                    const scalar r = p_face_shape_function[offSetSF_face + ic];
                    const scalar r_coord =
                        p_coordinate_face_shape_function[offSetSF_face + ic];

                    rhoBip += r * p_rho[ic];
                    if (balancedST)
                        skBip += r * p_sk[inn];

                    const label icNdim = ic * SPATIAL_DIM;
                    for (label j = 0; j < SPATIAL_DIM; ++j)
                    {
                        p_coordBip[j] +=
                            r_coord * p_coordinates[inn * SPATIAL_DIM + j];
                        p_uBip[j] += r * p_U[icNdim + j];
                        p_duBip[j] += r * p_du[icNdim + j];
                    }
                }

                // form dpdxBip
                for (label ic = 0; ic < nodesPerElement; ++ic)
                {
                    const label offSetDnDx =
                        SPATIAL_DIM * nodesPerElement * ip + ic * SPATIAL_DIM;
                    const scalar pIc = p_p[ic];
                    for (label j = 0; j < SPATIAL_DIM; ++j)
                    {
                        p_dpdxBip[j] += p_dndx[offSetDnDx + j] * pIc;
                        // balanced CSF: sigma*kappa_ip * compact grad(alpha)
                        if (balancedST)
                            p_FOrigBip[j] +=
                                skBip * p_dndx[offSetDnDx + j] * p_capPot[ic];
                    }
                }

                // form mDot:
                // rho*uj*Aj - rho*du*(dpdxj - Gjp)*Aj + rho*du*(FOrigj - Fj)*Aj
                scalar tmDot = 0.0;
                for (label j = 0; j < SPATIAL_DIM; ++j)
                {
                    const scalar axj = areaVec[ip * SPATIAL_DIM + j];

                    // divergence + pressure Rhie-Chow
                    tmDot +=
                        (rhoBip * p_uBip[j] -
                         rhoBip * p_duBip[j] * (p_dpdxBip[j] - p_GpdxBip[j])) *
                        axj;

                    // buoyancy stabilization: +rho*D*(F_orig - F)·S
                    tmDot +=
                        rhoBip * p_duBip[j] * (p_FOrigBip[j] - p_FBip[j]) * axj;
                }

                // store with relaxation
                mDot[ip] = mDotURF * tmDot + (1.0 - mDotURF) * mDot[ip];
            }
        }
    }
}

} /* namespace accel */
