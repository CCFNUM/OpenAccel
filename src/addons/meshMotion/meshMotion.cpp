// File       : meshMotion.cpp
// Created    : Fri Dec 13 2024 12:55:24 (+0100)
// Author     : Mhamad Mahdi Alloush
// Description:
// Copyright 2024 CCFNUM HSLU T&A. All Rights Reserved.

// code
#include "meshMotion.h"
#include "deformation.h"
#include "simulation.h"
#include "transformation.h"
#include "zoneDeformation.h"
#include "zoneTransformation.h"

// std
#include <limits>

namespace accel
{

meshMotion::meshMotion(realm* realm) : fieldBroker(realm)
{
    if (messager::master())
    {
        std::cout << "Creating mesh motion functionality" << std::endl;
    }

    // Create field instances
    DtRef();
    UmRef();

    // create deformation
    if (this->meshRef().anyZoneMeshDeforming())
    {
        deformationPtr_ = std::make_unique<deformation>(this);

        // required for GCL
        divUmRef();

        // ensure that all deforming zones are consistent in their relative
        // motion
        bool foundTotal = false;
        bool foundRelative = false;
        for (auto domain : realmPtr_->simulationRef().domainVector())
        {
            if (domain->zonePtr()->meshDeforming() &&
                !domain->zonePtr()
                     ->deformationRef()
                     .displacementRelativeToPreviousMesh())
            {
                foundTotal = true;
            }
            else
            {
                foundRelative = true;
            }
        }

        if (foundRelative && foundTotal)
        {
            errorMsg("Displacement in deforming zones must either be all "
                     "relative or all absolute");
        }
    }

    // create transformation
    if (this->meshRef().anyZoneMeshTransforming())
    {
        transformationPtr_ = std::make_unique<transformation>(this);
    }
}

void meshMotion::setup()
{
    for (auto domain : realmPtr_->simulationRef().domainVector())
    {
        if (domain->zonePtr()->meshMoving())
        {
            setupTotalDisplacement(domain);
            setupMeshVelocity(domain);

            if (domain->type() == domainType::fluid &&
                domain->zonePtr()->meshDeforming())
            {
                setupDivMeshVelocity(domain);
            }
        }
    }

    // setup deformation
    if (deformationPtr_)
    {
        deformationPtr_->setup();
    }

    // setup transformation
    if (transformationPtr_)
    {
        transformationPtr_->setup();
    }
}

void meshMotion::reset()
{
    // store prev time states of total displacement field: ineffective in the
    // first time step
    for (auto domain : realmPtr_->simulationRef().domainVector())
    {
        if (domain->zonePtr()->meshMoving())
        {
            updateTotalDisplacementPrevTimeField(domain);

            resetTotalDisplacement(domain);
            resetMeshVelocity(domain);
        }
    }
}

void meshMotion::initialize()
{
    // set total displacement and mesh velocity to zero
    reset();

    // no previous a_n after start-up or restart: the first one is unrelaxed
    hasAWallBc_ = false;

    // initialize deformation
    if (deformationPtr_)
    {
        deformationPtr_->initialize();
    }

    // initialize transformation
    if (transformationPtr_)
    {
        transformationPtr_->initialize();
    }
}

bool meshMotion::followsSolution() const
{
    for (const interface* interf : this->meshRef().interfaceVector())
    {
        if (interf->isFluidSolidType())
        {
            return true;
        }
    }
    return false;
}

void meshMotion::update()
{
    if (messager::master())
    {
        std::cout << "mesh motion: updating mesh ..." << std::endl;
    }

    // update total displacement due to deformation
    if (deformationPtr_)
    {
        deformationPtr_->update();
    }

    // update total displacement due to transformation
    if (transformationPtr_)
    {
        transformationPtr_->update();
    }

    // Update mesh coordinates from Dt (coords = original + Dt).
    updateCoordinates_();

    // calculate mesh velocity from total displacement
    updateMeshVelocityField_();

    // update div(Um) for deforming zones
    if (deformationPtr_)
    {
        updateMeshVelocityDivergenceField_();
    }

    // wall velocity from Dt on the fluid-solid sides (no-slip wall velocity)
    this->updateFsiWallVelocity_();

    // Robin-Neumann wall normal acceleration from Dt
    this->updateFsiWallAcceleration_();

    if (messager::master())
    {
        std::cout << "mesh motion: done updating mesh" << std::endl;
    }
}

void meshMotion::updateCoordinates_()
{
    // update coordinates from mesh deformation/motion
    for (label iZone = 0; iZone < this->meshRef().nZones(); iZone++)
    {
        if (this->meshRef().zoneRef(iZone).meshMoving())
        {
#ifndef NDEBUG
            if (messager::master())
            {
                std::cout << "updating coordinates of zone: "
                          << this->meshRef().zoneRef(iZone).name() << std::endl;
            }
#endif
            auto& mesh = this->meshRef();
            stk::mesh::BulkData& bulkData = mesh.bulkDataRef();
            stk::mesh::MetaData& metaData = mesh.metaDataRef();

            // get interior parts the domain is defined on
            const stk::mesh::PartVector& partVec =
                this->meshRef().zoneRef(iZone).interiorParts();

            // define some common selectors; select all nodes
            stk::mesh::Selector selUniversalNodes =
                metaData.universal_part() & stk::mesh::selectUnion(partVec);

            stk::mesh::BucketVector const& nodeBuckets = bulkData.get_buckets(
                stk::topology::NODE_RANK, selUniversalNodes);

            auto& DtSTKFieldRef = DtRef().stkFieldRef();
            const auto& originalCoordsSTKFieldRef = *metaData.get_field<scalar>(
                stk::topology::NODE_RANK, mesh::original_coordinates_ID);
            auto& coordsSTKFieldRef = *metaData.get_field<scalar>(
                stk::topology::NODE_RANK, mesh::coordinates_ID);

            for (stk::mesh::BucketVector::const_iterator ib =
                     nodeBuckets.begin();
                 ib != nodeBuckets.end();
                 ++ib)
            {
                stk::mesh::Bucket& nodeBucket = **ib;

                const stk::mesh::Bucket::size_type nNodesPerBucket =
                    nodeBucket.size();

                // field chunks in bucket
                scalar* Dtb = stk::mesh::field_data(DtSTKFieldRef, nodeBucket);
                const scalar* orgCoordsb = stk::mesh::field_data(
                    originalCoordsSTKFieldRef, nodeBucket);
                scalar* coordsb =
                    stk::mesh::field_data(coordsSTKFieldRef, nodeBucket);

                for (stk::mesh::Bucket::size_type iNode = 0;
                     iNode < nNodesPerBucket;
                     ++iNode)
                {
                    for (label i = 0; i < SPATIAL_DIM; i++)
                    {
                        // update coords
                        coordsb[SPATIAL_DIM * iNode + i] =
                            orgCoordsb[SPATIAL_DIM * iNode + i] +
                            Dtb[SPATIAL_DIM * iNode + i];
                    }
                }
            }
        }
    }

    // update mesh related structures
    this->meshRef().update();
}

void meshMotion::updateMeshVelocityField_()
{
    // calculate mesh velocity from mesh displacement: hack a mesh velocity to
    // be first order backward Euler always
    for (label iZone = 0; iZone < this->meshRef().nZones(); iZone++)
    {
        if (this->meshRef().zoneRef(iZone).meshMoving())
        {
            const auto& zone = this->meshRef().zoneRef(iZone);
            // chord velocity (D - D_old)/dt has divergence 2(cos(w dt)-1)/dt: a
            // spurious compression of the whole block that grows with w dt
            if (zone.meshTransforming() && !zone.meshDeforming() &&
                zone.transformationRef().type() == meshMotionType::rotating)
            {
                updateRigidRotationMeshVelocity_(iZone);
                continue;
            }

            auto& mesh = this->meshRef();
            stk::mesh::BulkData& bulkData = mesh.bulkDataRef();
            stk::mesh::MetaData& metaData = mesh.metaDataRef();

            // get interior parts the domain is defined on
            const stk::mesh::PartVector& partVec =
                this->meshRef().zoneRef(iZone).interiorParts();

            // define some common selectors; select owned nodes
            stk::mesh::Selector selUniversalNodes =
                metaData.universal_part() & stk::mesh::selectUnion(partVec);

            stk::mesh::BucketVector const& nodeBuckets = bulkData.get_buckets(
                stk::topology::NODE_RANK, selUniversalNodes);

            // Get fields
            const auto& DtSTKFieldRef = DtRef().stkFieldRef();
            const auto& DtSTKFieldRefOld = DtRef().prevTimeRef().stkFieldRef();
            const auto& DtSTKFieldRefOldOld =
                DtRef().prevTimeRef().prevTimeRef().stkFieldRef();
            auto& UmSTKFieldRef = UmRef().stkFieldRef();

            const scalar dt = meshRef().controlsRef().getTimestep();

            const auto c =
                BDF2::coeff(dt, meshRef().controlsRef().getTimestep(-1));

            for (stk::mesh::BucketVector::const_iterator ib =
                     nodeBuckets.begin();
                 ib != nodeBuckets.end();
                 ++ib)
            {
                stk::mesh::Bucket& nodeBucket = **ib;

                const stk::mesh::Bucket::size_type nNodesPerBucket =
                    nodeBucket.size();

                // field chunks in bucket
                const scalar* Dtb =
                    stk::mesh::field_data(DtSTKFieldRef, nodeBucket);
                const scalar* DtbOld =
                    stk::mesh::field_data(DtSTKFieldRefOld, nodeBucket);
                const scalar* DtbOldOld =
                    stk::mesh::field_data(DtSTKFieldRefOldOld, nodeBucket);
                scalar* Umb = stk::mesh::field_data(UmSTKFieldRef, nodeBucket);

                for (stk::mesh::Bucket::size_type iNode = 0;
                     iNode < nNodesPerBucket;
                     ++iNode)
                {
                    for (label i = 0; i < SPATIAL_DIM; i++)
                    {
                        Umb[SPATIAL_DIM * iNode + i] =
                            (c[0] * Dtb[SPATIAL_DIM * iNode + i] +
                             c[1] * DtbOld[SPATIAL_DIM * iNode + i] +
                             c[2] * DtbOldOld[SPATIAL_DIM * iNode + i]) /
                            dt;
                    }
                }
            }
        }
    }
}

void meshMotion::updateRigidRotationMeshVelocity_(label iZone)
{
    auto& mesh = this->meshRef();
    stk::mesh::BulkData& bulkData = mesh.bulkDataRef();
    stk::mesh::MetaData& metaData = mesh.metaDataRef();
    const auto& zone = mesh.zoneRef(iZone);

    const STKScalarField* coordsSTKFieldPtr = metaData.get_field<scalar>(
        stk::topology::NODE_RANK, metaData.coordinate_field_name());
    auto& UmSTKFieldRef = UmRef().stkFieldRef();

    const scalar* p_mat =
        zone.transformationRef().rotation().coriolisMatrix_.data();
    const scalar* p_ori = zone.transformationRef().rotation().origin_.data();

    const stk::mesh::Selector selStationary =
        stk::mesh::selectUnion(zone.stationaryParts());
    stk::mesh::Selector selNodes = metaData.universal_part() &
                                   stk::mesh::selectUnion(zone.interiorParts());

    for (const stk::mesh::Bucket* b :
         bulkData.get_buckets(stk::topology::NODE_RANK, selNodes))
    {
        const bool stationary = selStationary(*b);
        const scalar* x = stk::mesh::field_data(*coordsSTKFieldPtr, *b);
        scalar* Um = stk::mesh::field_data(UmSTKFieldRef, *b);
        for (stk::mesh::Bucket::size_type k = 0; k < b->size(); ++k)
        {
            for (label i = 0; i < SPATIAL_DIM; i++)
            {
                scalar v = 0.0;
                if (!stationary)
                {
                    for (label j = 0; j < SPATIAL_DIM; j++)
                    {
                        v += p_mat[i * SPATIAL_DIM + j] *
                             (x[SPATIAL_DIM * k + j] - p_ori[j]);
                    }
                }
                Um[SPATIAL_DIM * k + i] = v;
            }
        }
    }
}

void meshMotion::updateMeshVelocityDivergenceField_()
{
    // Update divergence of mesh velocity field using Gauss theorem:
    // div(Um)_i = sum(Um_ip . S_ip) / V_i
    for (auto domain : realmPtr_->simulationRef().domainVector())
    {
        if (domain->type() == domainType::fluid &&
            domain->zonePtr()->meshDeforming())
        {
            // reset divUm before accumulation (necessary when called
            // multiple times per time step, e.g. mesh motion per Newton
            // iteration)
            resetDivMeshVelocity(domain);

            auto& mesh = this->meshRef();
            stk::mesh::BulkData& bulkData = mesh.bulkDataRef();
            stk::mesh::MetaData& metaData = mesh.metaDataRef();

            const zone* zonePtr = domain->zonePtr();
            const label iZone = zonePtr->index();

            // shifted ip's for fields?
            const bool isUShifted = URef().isShifted();

            // Get fields
            const auto& UmSTKFieldRef = UmRef().stkFieldRef();
            auto& divUmSTKFieldRef = divUmRef().stkFieldRef();

            // Geometric fields
            STKScalarField* coordinates = metaData.get_field<scalar>(
                stk::topology::NODE_RANK, mesh::coordinates_ID);
            const STKScalarField* dualNodalVolumeSTKFieldPtr =
                metaData.get_field<scalar>(stk::topology::NODE_RANK,
                                           mesh::dual_nodal_volume_ID);

            // ========================================
            // Interior IP contribution
            // ========================================
            {
                std::vector<scalar> ws_shape_function;
                std::vector<scalar> ws_Um;
                std::vector<scalar> ws_coordinates;
                std::vector<scalar> ws_scs_areav;
                std::vector<scalar> ws_dualVolume;

                // integration point data that depends on size
                std::vector<scalar> umIp(SPATIAL_DIM);

                // pointers to everyone...
                scalar* p_umIp = &umIp[0];

                stk::mesh::Selector selAllElements =
                    metaData.universal_part() &
                    stk::mesh::selectUnion(zonePtr->interiorParts());

                stk::mesh::BucketVector const& elementBuckets =
                    bulkData.get_buckets(stk::topology::ELEMENT_RANK,
                                         selAllElements);
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
                    const label nodesPerElement = meSCS->nodesPerElement_;
                    const label numScsIp = meSCS->numIntPoints_;

                    ws_shape_function.resize(numScsIp * nodesPerElement);
                    scalar* p_shape_function = &ws_shape_function[0];
                    if (isUShifted)
                    {
                        meSCS->shifted_shape_fcn(&p_shape_function[0]);
                    }
                    else
                    {
                        meSCS->shape_fcn(&p_shape_function[0]);
                    }

                    // set sizes
                    ws_Um.resize(nodesPerElement * SPATIAL_DIM);
                    ws_coordinates.resize(nodesPerElement * SPATIAL_DIM);
                    ws_scs_areav.resize(numScsIp * SPATIAL_DIM);
                    ws_dualVolume.resize(nodesPerElement);

                    // pointers
                    scalar* p_Um = &ws_Um[0];
                    scalar* p_coordinates = &ws_coordinates[0];
                    scalar* p_scs_areav = &ws_scs_areav[0];
                    scalar* p_dualVolume = &ws_dualVolume[0];

                    const label* lrscv = meSCS->adjacentNodes();

                    for (stk::mesh::Bucket::size_type k = 0;
                         k < nElementsPerBucket;
                         ++k)
                    {
                        stk::mesh::Entity const* nodeRels =
                            elementBucket.begin_nodes(k);
                        const label numNodes = elementBucket.num_nodes(k);

                        for (label ni = 0; ni < numNodes; ++ni)
                        {
                            stk::mesh::Entity node = nodeRels[ni];

                            const scalar* Um =
                                stk::mesh::field_data(UmSTKFieldRef, node);
                            const scalar* coords =
                                stk::mesh::field_data(*coordinates, node);

                            p_dualVolume[ni] = *stk::mesh::field_data(
                                *dualNodalVolumeSTKFieldPtr, node);

                            for (label i = 0; i < SPATIAL_DIM; ++i)
                            {
                                p_Um[ni * SPATIAL_DIM + i] = Um[i];
                                p_coordinates[ni * SPATIAL_DIM + i] = coords[i];
                            }
                        }

                        // compute geometry
                        scalar scs_error = 0.0;
                        meSCS->determinant(
                            1, &p_coordinates[0], &p_scs_areav[0], &scs_error);

                        // start IP loop
                        for (label ip = 0; ip < numScsIp; ++ip)
                        {
                            // left and right nodes for this ip
                            const label il = lrscv[2 * ip];
                            const label ir = lrscv[2 * ip + 1];

                            stk::mesh::Entity nodeL = nodeRels[il];
                            stk::mesh::Entity nodeR = nodeRels[ir];

                            scalar* divUmL =
                                stk::mesh::field_data(divUmSTKFieldRef, nodeL);
                            scalar* divUmR =
                                stk::mesh::field_data(divUmSTKFieldRef, nodeR);

                            // zero-out
                            for (label j = 0; j < SPATIAL_DIM; ++j)
                            {
                                p_umIp[j] = 0.0;
                            }

                            const label offset = ip * nodesPerElement;
                            for (label ic = 0; ic < nodesPerElement; ++ic)
                            {
                                const scalar r = p_shape_function[offset + ic];
                                for (label j = 0; j < SPATIAL_DIM; ++j)
                                {
                                    p_umIp[j] += r * p_Um[ic * SPATIAL_DIM + j];
                                }
                            }

                            // compute flux: Um_ip . S_ip
                            scalar flux = 0.0;
                            for (label j = 0; j < SPATIAL_DIM; ++j)
                            {
                                flux += p_umIp[j] *
                                        p_scs_areav[ip * SPATIAL_DIM + j];
                            }

                            // accumulate to left/right nodes
                            scalar inv_volL = 1.0 / p_dualVolume[il];
                            scalar inv_volR = 1.0 / p_dualVolume[ir];

                            *divUmL += flux * inv_volL;
                            *divUmR -= flux * inv_volR;
                        }
                    }
                }
            }

            // ========================================
            // Scope 2: Interface IP contribution
            // ========================================
            {
                const auto& exposedAreaVecSTKFieldRef =
                    *metaData.get_field<scalar>(metaData.side_rank(),
                                                mesh::exposed_area_vector_ID);

                std::vector<scalar> ws_Um;
                std::vector<scalar> ws_shape_function;
                std::vector<scalar> currentIsoParCoords(SPATIAL_DIM);
                std::vector<scalar> currentUmBip(SPATIAL_DIM);
                scalar* p_currentUmBip = &currentUmBip[0];

                std::vector<const interfaceSideInfo*> sidesToProcess;

                for (const interface* interf : zonePtr->interfacesRef())
                {
                    if (interf->isConformalTreatment())
                        continue;

                    if (interf->isInternal())
                    {
                        // Internal interface: both master and slave sides
                        // are in this domain. Use ipInfo iteration.
                        sidesToProcess.clear();
                        sidesToProcess.push_back(interf->masterInfoPtr());
                        sidesToProcess.push_back(interf->slaveInfoPtr());

                        for (const auto* sideInfo : sidesToProcess)
                        {
                            const auto& ipInfoVec = sideInfo->ipInfoVec();

                            for (label iSide = 0;
                                 iSide < static_cast<label>(ipInfoVec.size());
                                 iSide++)
                            {
                                const auto& faceIpInfoVec = ipInfoVec[iSide];

                                for (label k = 0; k < static_cast<label>(
                                                          faceIpInfoVec.size());
                                     ++k)
                                {
                                    ipInfo* dg = faceIpInfoVec[k];

                                    stk::mesh::Entity currentFace =
                                        dg->currentFace_;
                                    MasterElement* meFCCurrent =
                                        dg->meFCCurrent_;
                                    const label currentGaussPointId =
                                        dg->currentGaussPointId_;
                                    currentIsoParCoords =
                                        dg->currentIsoParCoords_;

                                    const label* faceIpNodeMap =
                                        meFCCurrent->ipNodeMap();
                                    const label currentNodesPerFace =
                                        meFCCurrent->nodesPerElement_;

                                    ws_Um.resize(currentNodesPerFace *
                                                 SPATIAL_DIM);
                                    scalar* p_Um = &ws_Um[0];

                                    // gather Um on current face
                                    // nodes
                                    stk::mesh::Entity const*
                                        current_face_node_rels =
                                            bulkData.begin_nodes(currentFace);
                                    const label current_num_face_nodes =
                                        bulkData.num_nodes(currentFace);

                                    for (label ni = 0;
                                         ni < current_num_face_nodes;
                                         ++ni)
                                    {
                                        stk::mesh::Entity node =
                                            current_face_node_rels[ni];
                                        const scalar* Um =
                                            stk::mesh::field_data(UmSTKFieldRef,
                                                                  node);
                                        for (label i = 0; i < SPATIAL_DIM; ++i)
                                        {
                                            p_Um[i * current_num_face_nodes +
                                                 ni] = Um[i];
                                        }
                                    }

                                    // interpolate Um to IP
                                    meFCCurrent->interpolatePoint(
                                        SPATIAL_DIM,
                                        &currentIsoParCoords[0],
                                        &ws_Um[0],
                                        &currentUmBip[0]);

                                    // nearest node
                                    const label nn =
                                        faceIpNodeMap[currentGaussPointId];
                                    stk::mesh::Entity nNode =
                                        current_face_node_rels[nn];

                                    scalar vol = *stk::mesh::field_data(
                                        *dualNodalVolumeSTKFieldPtr, nNode);
                                    scalar* divUm = stk::mesh::field_data(
                                        divUmSTKFieldRef, nNode);

                                    // face area vector
                                    const scalar* c_areaVec =
                                        stk::mesh::field_data(
                                            exposedAreaVecSTKFieldRef,
                                            currentFace);

                                    // compute flux and accumulate
                                    scalar flux = 0.0;
                                    for (label j = 0; j < SPATIAL_DIM; ++j)
                                    {
                                        flux += currentUmBip[j] *
                                                c_areaVec[currentGaussPointId *
                                                              SPATIAL_DIM +
                                                          j];
                                    }

                                    *divUm += flux / vol;
                                }
                            }
                        }
                    }
                    else if (interf->isFluidSolidType())
                    {
                        // External fluid-solid interface: treat as no-slip
                        // wall using side bucket iteration to match the
                        // assembler's no-slip wall interface treatment
                        const auto* sideInfo =
                            interf->interfaceSideInfoPtr(iZone);

                        std::vector<stk::topology> parentTopo;

                        stk::mesh::Selector selAllSides =
                            metaData.universal_part() &
                            stk::mesh::selectUnion(sideInfo->currentPartVec_);

                        stk::mesh::BucketVector const& sideBuckets =
                            bulkData.get_buckets(metaData.side_rank(),
                                                 selAllSides);
                        for (stk::mesh::BucketVector::const_iterator ib =
                                 sideBuckets.begin();
                             ib != sideBuckets.end();
                             ++ib)
                        {
                            stk::mesh::Bucket& sideBucket = **ib;
                            const stk::mesh::Bucket::size_type nSidesPerBucket =
                                sideBucket.size();

                            // extract connected element topology
                            sideBucket.parent_topology(
                                stk::topology::ELEMENT_RANK, parentTopo);
                            stk::topology theElemTopo = parentTopo[0];

                            // volume master element
                            MasterElement* meSCS =
                                MasterElementRepo::get_surface_master_element(
                                    theElemTopo);

                            // face master element
                            MasterElement* meFC =
                                MasterElementRepo::get_surface_master_element(
                                    sideBucket.topology());

                            const label nodesPerSide =
                                sideBucket.topology().num_nodes();
                            const label numScsBip = meFC->numIntPoints_;

                            ws_Um.resize(nodesPerSide * SPATIAL_DIM);
                            ws_shape_function.resize(numScsBip * nodesPerSide);

                            scalar* p_Um = &ws_Um[0];
                            scalar* p_shape_function = &ws_shape_function[0];

                            if (isUShifted)
                            {
                                meFC->shifted_shape_fcn(&p_shape_function[0]);
                            }
                            else
                            {
                                meFC->shape_fcn(&p_shape_function[0]);
                            }

                            for (stk::mesh::Bucket::size_type iSide = 0;
                                 iSide < nSidesPerBucket;
                                 ++iSide)
                            {
                                // get face
                                stk::mesh::Entity side = sideBucket[iSide];

                                // gather Um from face nodes
                                stk::mesh::Entity const* sideNodeRels =
                                    bulkData.begin_nodes(side);
                                const label numSideNodes =
                                    bulkData.num_nodes(side);

                                for (label ni = 0; ni < numSideNodes; ++ni)
                                {
                                    stk::mesh::Entity node = sideNodeRels[ni];
                                    const scalar* Um = stk::mesh::field_data(
                                        UmSTKFieldRef, node);

                                    const label offSet = ni * SPATIAL_DIM;
                                    for (label j = 0; j < SPATIAL_DIM; ++j)
                                    {
                                        p_Um[offSet + j] = Um[j];
                                    }
                                }

                                // face area vector
                                const scalar* areaVec = stk::mesh::field_data(
                                    exposedAreaVecSTKFieldRef, side);

                                // get connected element and face ordinal
                                const stk::mesh::Entity* faceElemRels =
                                    bulkData.begin_elements(side);
                                stk::mesh::Entity element = faceElemRels[0];
                                const stk::mesh::ConnectivityOrdinal*
                                    face_elem_ords =
                                        bulkData.begin_element_ordinals(side);
                                const label faceOrdinal = face_elem_ords[0];

                                // mapping from ip to nodes for this ordinal
                                const label* ipNodeMap =
                                    meSCS->ipNodeMap(faceOrdinal);

                                // element node relations
                                stk::mesh::Entity const* elemNodeRels =
                                    bulkData.begin_nodes(element);

                                for (label ip = 0; ip < numScsBip; ++ip)
                                {
                                    const label nearestNode = ipNodeMap[ip];
                                    stk::mesh::Entity node =
                                        elemNodeRels[nearestNode];

                                    scalar* divUm = stk::mesh::field_data(
                                        divUmSTKFieldRef, node);
                                    scalar vol = *stk::mesh::field_data(
                                        *dualNodalVolumeSTKFieldPtr, node);

                                    // interpolate Um to IP
                                    for (label j = 0; j < SPATIAL_DIM; ++j)
                                    {
                                        p_currentUmBip[j] = 0.0;
                                    }

                                    const label offset = ip * nodesPerSide;
                                    for (label ic = 0; ic < nodesPerSide; ++ic)
                                    {
                                        const scalar r =
                                            p_shape_function[offset + ic];
                                        for (label j = 0; j < SPATIAL_DIM; ++j)
                                        {
                                            p_currentUmBip[j] +=
                                                r * p_Um[ic * SPATIAL_DIM + j];
                                        }
                                    }

                                    scalar flux = 0.0;
                                    for (label j = 0; j < SPATIAL_DIM; ++j)
                                    {
                                        flux += p_currentUmBip[j] *
                                                areaVec[ip * SPATIAL_DIM + j];
                                    }

                                    *divUm += flux / vol;
                                }
                            }
                        }
                    }
                    else
                    {
                        // External non-fluid-solid interface: use ipInfo
                        // iteration for our side only
                        const auto* sideInfo =
                            interf->interfaceSideInfoPtr(iZone);

                        const auto& ipInfoVec = sideInfo->ipInfoVec();

                        for (label iSide = 0;
                             iSide < static_cast<label>(ipInfoVec.size());
                             iSide++)
                        {
                            const auto& faceIpInfoVec = ipInfoVec[iSide];

                            for (label k = 0;
                                 k < static_cast<label>(faceIpInfoVec.size());
                                 ++k)
                            {
                                ipInfo* dg = faceIpInfoVec[k];

                                stk::mesh::Entity currentFace =
                                    dg->currentFace_;
                                MasterElement* meFCCurrent = dg->meFCCurrent_;
                                const label currentGaussPointId =
                                    dg->currentGaussPointId_;
                                currentIsoParCoords = dg->currentIsoParCoords_;

                                const label* faceIpNodeMap =
                                    meFCCurrent->ipNodeMap();
                                const label currentNodesPerFace =
                                    meFCCurrent->nodesPerElement_;

                                ws_Um.resize(currentNodesPerFace * SPATIAL_DIM);
                                scalar* p_Um = &ws_Um[0];

                                // gather Um on current face nodes
                                stk::mesh::Entity const*
                                    current_face_node_rels =
                                        bulkData.begin_nodes(currentFace);
                                const label current_num_face_nodes =
                                    bulkData.num_nodes(currentFace);

                                for (label ni = 0; ni < current_num_face_nodes;
                                     ++ni)
                                {
                                    stk::mesh::Entity node =
                                        current_face_node_rels[ni];
                                    const scalar* Um = stk::mesh::field_data(
                                        UmSTKFieldRef, node);
                                    for (label i = 0; i < SPATIAL_DIM; ++i)
                                    {
                                        p_Um[i * current_num_face_nodes + ni] =
                                            Um[i];
                                    }
                                }

                                // interpolate Um to IP
                                meFCCurrent->interpolatePoint(
                                    SPATIAL_DIM,
                                    &currentIsoParCoords[0],
                                    &ws_Um[0],
                                    &currentUmBip[0]);

                                // nearest node
                                const label nn =
                                    faceIpNodeMap[currentGaussPointId];
                                stk::mesh::Entity nNode =
                                    current_face_node_rels[nn];

                                scalar vol = *stk::mesh::field_data(
                                    *dualNodalVolumeSTKFieldPtr, nNode);
                                scalar* divUm = stk::mesh::field_data(
                                    divUmSTKFieldRef, nNode);

                                // face area vector
                                const scalar* c_areaVec = stk::mesh::field_data(
                                    exposedAreaVecSTKFieldRef, currentFace);

                                // compute flux and accumulate
                                scalar flux = 0.0;
                                for (label j = 0; j < SPATIAL_DIM; ++j)
                                {
                                    flux += currentUmBip[j] *
                                            c_areaVec[currentGaussPointId *
                                                          SPATIAL_DIM +
                                                      j];
                                }

                                *divUm += flux / vol;
                            }
                        }
                    }
                }
            }

            // ========================================
            // Boundary IP contribution
            // ========================================
            {
                // get fields
                const auto& exposedAreaVecSTKFieldRef =
                    *metaData.get_field<scalar>(metaData.side_rank(),
                                                mesh::exposed_area_vector_ID);

                // scratch arrays
                std::vector<scalar> ws_Um;
                std::vector<scalar> ws_shape_function;

                // fixed-size arrays
                std::vector<scalar> UmIp(SPATIAL_DIM);

                // pointers
                scalar* p_UmIp = &UmIp[0];

                std::vector<stk::topology> parentTopo;

                for (label iBoundary = 0; iBoundary < zonePtr->nBoundaries();
                     iBoundary++)
                {
                    stk::mesh::PartVector partVec =
                        zonePtr->boundaryRef(iBoundary).parts();

                    stk::mesh::Selector selAllSides =
                        metaData.universal_part() &
                        stk::mesh::selectUnion(partVec);

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
                            MasterElementRepo::get_surface_master_element(
                                theElemTopo);

                        // face master element
                        MasterElement* meFC =
                            MasterElementRepo::get_surface_master_element(
                                sideBucket.topology());

                        const label nodesPerSide =
                            sideBucket.topology().num_nodes();
                        const label numScsBip = meFC->numIntPoints_;

                        ws_Um.resize(nodesPerSide * SPATIAL_DIM);
                        ws_shape_function.resize(numScsBip * nodesPerSide);

                        scalar* p_Um = &ws_Um[0];
                        scalar* p_shape_function = &ws_shape_function[0];

                        if (isUShifted)
                        {
                            meFC->shifted_shape_fcn(&p_shape_function[0]);
                        }
                        else
                        {
                            meFC->shape_fcn(&p_shape_function[0]);
                        }

                        for (stk::mesh::Bucket::size_type iSide = 0;
                             iSide < nSidesPerBucket;
                             ++iSide)
                        {
                            // get face
                            stk::mesh::Entity side = sideBucket[iSide];

                            // gather Um from face nodes
                            stk::mesh::Entity const* sideNodeRels =
                                bulkData.begin_nodes(side);
                            const label numSideNodes = bulkData.num_nodes(side);

                            for (label ni = 0; ni < numSideNodes; ++ni)
                            {
                                stk::mesh::Entity node = sideNodeRels[ni];
                                const scalar* Um =
                                    stk::mesh::field_data(UmSTKFieldRef, node);

                                const label offSet = ni * SPATIAL_DIM;
                                for (label j = 0; j < SPATIAL_DIM; ++j)
                                {
                                    p_Um[offSet + j] = Um[j];
                                }
                            }

                            // face area vector
                            const scalar* areaVec = stk::mesh::field_data(
                                exposedAreaVecSTKFieldRef, side);

                            // get connected element and face
                            // ordinal
                            const stk::mesh::Entity* faceElemRels =
                                bulkData.begin_elements(side);
                            stk::mesh::Entity element = faceElemRels[0];
                            const stk::mesh::ConnectivityOrdinal*
                                face_elem_ords =
                                    bulkData.begin_element_ordinals(side);
                            const label faceOrdinal = face_elem_ords[0];

                            // mapping from ip to nodes for this
                            // ordinal
                            const label* ipNodeMap =
                                meSCS->ipNodeMap(faceOrdinal);

                            // element node relations
                            stk::mesh::Entity const* elemNodeRels =
                                bulkData.begin_nodes(element);

                            // start assembly
                            for (label ip = 0; ip < numScsBip; ++ip)
                            {
                                // nearest node
                                const label nearestNode = ipNodeMap[ip];
                                stk::mesh::Entity node =
                                    elemNodeRels[nearestNode];

                                scalar* divUm = stk::mesh::field_data(
                                    divUmSTKFieldRef, node);
                                scalar vol = *stk::mesh::field_data(
                                    *dualNodalVolumeSTKFieldPtr, node);

                                // interpolate Um to IP
                                for (label j = 0; j < SPATIAL_DIM; ++j)
                                {
                                    p_UmIp[j] = 0.0;
                                }

                                const label offset = ip * nodesPerSide;
                                for (label ic = 0; ic < nodesPerSide; ++ic)
                                {
                                    const scalar r =
                                        p_shape_function[offset + ic];
                                    for (label j = 0; j < SPATIAL_DIM; ++j)
                                    {
                                        p_UmIp[j] +=
                                            r * p_Um[ic * SPATIAL_DIM + j];
                                    }
                                }

                                // compute flux and accumulate
                                scalar flux = 0.0;
                                for (label j = 0; j < SPATIAL_DIM; ++j)
                                {
                                    flux += p_UmIp[j] *
                                            areaVec[ip * SPATIAL_DIM + j];
                                }

                                *divUm += flux / vol;
                            }
                        }
                    }
                }
            }

            // sync for parallel
            if (messager::parallel())
            {
                divUmRef().synchronizeGhostedEntities(domain->index());
            }
        }
    }
}

void meshMotion::updateFsiWallVelocity_()
{
    // wall velocity from Dt at the fluid-solid BIPs (BDF2 for backward)

    // only fluid-solid interfaces carry a wall velocity
    bool anyFsi = false;
    for (auto& zonePtr : realmPtr_->simulationRef().domainVector())
    {
        if (zonePtr->type() != domainType::fluid)
            continue;

        for (const interface* interf : zonePtr->interfacesRef())
        {
            if (interf->isFluidSolidType())
            {
                anyFsi = true;
                break;
            }
        }
        if (anyFsi)
            break;
    }
    if (!anyFsi)
        return;

    auto& mesh = meshRef();
    stk::mesh::BulkData& bulkData = mesh.bulkDataRef();
    stk::mesh::MetaData& metaData = mesh.metaDataRef();

    const scalar dt = mesh.controlsRef().getTimestep();
    const scalar invDt = 1.0 / (dt + SMALL);

    const auto& DtSTKFieldRef = DtRef().stkFieldRef();
    const auto& DtOldSTKFieldRef = DtRef().prevTimeRef().stkFieldRef();
    auto& sideUSTKFieldRef = URef().sideFieldRef().stkFieldRef();

    // BDF2 wall velocity only with solid_transient_scheme: backward
    const label tsCount = mesh.controlsRef().getTimeStepCount();
    const bool useVelocityBdf2 =
        (mesh.controlsRef().getSolidTransientScheme() ==
         solidTransientSchemeType::backwardD2dt2) &&
        (tsCount >= 2);

    // variable-step BDF2 coefficients (as in updateFsiWallAcceleration_())
    const scalar dt0 = mesh.controlsRef().getTimestep(-1);
    const scalar safedt0 = (dt0 > 0.0) ? dt0 : dt;
    const scalar coefft = 1.0 + dt / (dt + safedt0);
    const scalar coefft00 = dt * dt / (safedt0 * (dt + safedt0));
    const scalar coefft0 = coefft + coefft00;

    // DtOld2 field: only accessed when BDF2 is active
    const auto* DtOld2STKFieldPtr =
        useVelocityBdf2 ? &DtRef().prevTimeRef().prevTimeRef().stkFieldRef()
                        : nullptr;

    // write UbcVec on the fluid-side faces of each fluid-solid interface
    for (auto& zonePtr : realmPtr_->simulationRef().domainVector())
    {
        if (zonePtr->type() != domainType::fluid)
            continue;

        const label iZone = zonePtr->index();

        for (const interface* interf : zonePtr->interfacesRef())
        {
            if (!interf->isFluidSolidType())
                continue;

            // Only write on the fluid side (this zone's own side info).
            const auto* fluidSideInfoPtr = interf->interfaceSideInfoPtr(iZone);

            stk::mesh::Selector selAllSides =
                metaData.universal_part() &
                stk::mesh::selectUnion(fluidSideInfoPtr->currentPartVec_);

            stk::mesh::BucketVector const& sideBuckets =
                bulkData.get_buckets(metaData.side_rank(), selAllSides);

            for (auto ib = sideBuckets.begin(); ib != sideBuckets.end(); ++ib)
            {
                stk::mesh::Bucket& sideBucket = **ib;

                MasterElement* meFC =
                    MasterElementRepo::get_surface_master_element(
                        sideBucket.topology());
                const label numScsBip = meFC->numIntPoints_;
                const label nodesPerSide = sideBucket.topology().num_nodes();

                std::vector<scalar> ws_shape(numScsBip * nodesPerSide);
                meFC->shape_fcn(&ws_shape[0]);

                const stk::mesh::Bucket::size_type nSides = sideBucket.size();
                for (stk::mesh::Bucket::size_type iSide = 0; iSide < nSides;
                     ++iSide)
                {
                    stk::mesh::Entity side = sideBucket[iSide];
                    scalar* UbcVec =
                        stk::mesh::field_data(sideUSTKFieldRef, side);

                    stk::mesh::Entity const* sideNodeRels =
                        bulkData.begin_nodes(side);
                    const label numSideNodes = bulkData.num_nodes(side);

                    // nodal wall velocity
                    std::vector<scalar> vNodal(numSideNodes * SPATIAL_DIM);
                    for (label ni = 0; ni < numSideNodes; ++ni)
                    {
                        stk::mesh::Entity node = sideNodeRels[ni];
                        const scalar* Dt =
                            stk::mesh::field_data(DtSTKFieldRef, node);
                        const scalar* DtOld =
                            stk::mesh::field_data(DtOldSTKFieldRef, node);
                        const scalar* DtOld2 = nullptr;
                        if (useVelocityBdf2)
                            DtOld2 =
                                stk::mesh::field_data(*DtOld2STKFieldPtr, node);

                        for (label j = 0; j < SPATIAL_DIM; ++j)
                        {
                            if (useVelocityBdf2)
                            {
                                vNodal[ni * SPATIAL_DIM + j] =
                                    (coefft * Dt[j] - coefft0 * DtOld[j] +
                                     coefft00 * DtOld2[j]) *
                                    invDt;
                            }
                            else
                            {
                                vNodal[ni * SPATIAL_DIM + j] =
                                    (Dt[j] - DtOld[j]) * invDt;
                            }
                        }
                    }

                    // Interpolate to BIPs and write UbcVec.
                    for (label ip = 0; ip < numScsBip; ++ip)
                    {
                        const label offSF = ip * nodesPerSide;
                        for (label j = 0; j < SPATIAL_DIM; ++j)
                        {
                            scalar vip = 0.0;
                            for (label ic = 0; ic < nodesPerSide; ++ic)
                            {
                                vip += ws_shape[offSF + ic] *
                                       vNodal[ic * SPATIAL_DIM + j];
                            }
                            UbcVec[ip * SPATIAL_DIM + j] = vip;
                        }
                    }
                }
            }
        }
    }
}

void meshMotion::updateFsiWallAcceleration_()
{
    // RN wall normal acceleration from Dt (solid d2dt2 operator)

    // the aWallBc side field exists on Robin-Neumann interfaces only
    bool anyRnFsi = false;
    for (auto& zonePtr : realmPtr_->simulationRef().domainVector())
    {
        if (zonePtr->type() != domainType::fluid)
            continue;

        for (const interface* interf : zonePtr->interfacesRef())
        {
            if (interf->isFluidSolidType() && interf->isRobinNeumann())
            {
                anyRnFsi = true;
                break;
            }
        }
        if (anyRnFsi)
            break;
    }
    if (!anyRnFsi)
        return;

    auto& mesh = meshRef();
    stk::mesh::BulkData& bulkData = mesh.bulkDataRef();
    stk::mesh::MetaData& metaData = mesh.metaDataRef();

    const scalar dt = mesh.controlsRef().getTimestep();
    const scalar rDeltaT2 = 1.0 / (dt * dt);

    const bool isBackward = (mesh.controlsRef().getSolidTransientScheme() ==
                             solidTransientSchemeType::backwardD2dt2);

    // Dt3/Dt4 exist only with the backward scheme
    const auto& DtSTKFieldRef = DtRef().stkFieldRef();
    const auto& Dt1STKFieldRef = DtRef().prevTimeRef().stkFieldRef();
    const auto& Dt2STKFieldRef =
        DtRef().prevTimeRef().prevTimeRef().stkFieldRef();

    const stk::mesh::Field<scalar>* Dt3STKFieldPtr = nullptr;
    const stk::mesh::Field<scalar>* Dt4STKFieldPtr = nullptr;
    if (isBackward)
    {
        Dt3STKFieldPtr =
            &DtRef().prevTimeRef().prevTimeRef().prevTimeRef().stkFieldRef();
        Dt4STKFieldPtr = &DtRef()
                              .prevTimeRef()
                              .prevTimeRef()
                              .prevTimeRef()
                              .prevTimeRef()
                              .stkFieldRef();
    }

    // outer BDF2 coefficients as in assembleBackwardD2dt2NodeLoop_()
    scalar coefft = 0.0, coefft0 = 0.0, coefft00 = 0.0;
    const scalar inner_c0 = 1.5;
    const scalar inner_c1 = -2.0;
    const scalar inner_c2 = 0.5;
    if (isBackward)
    {
        const scalar dt0 = mesh.controlsRef().getTimestep(-1);
        const scalar safedt0 = (dt0 > 0.0) ? dt0 : dt;
        coefft = 1.0 + dt / (dt + safedt0);
        coefft00 = dt * dt / (safedt0 * (dt + safedt0));
        coefft0 = coefft + coefft00;
    }

    // exposed area vector: same field and normalisation as the RN kernel
    const auto& exposedAreaVecSTKFieldRef = *metaData.get_field<scalar>(
        metaData.side_rank(), mesh::exposed_area_vector_ID);

    auto& aWallBcSideSTKFieldRef = aWallBcRef().stkFieldRef();

    std::vector<scalar> aBip(SPATIAL_DIM);
    scalar* p_aBip = &aBip[0];

    for (auto& zonePtr : realmPtr_->simulationRef().domainVector())
    {
        if (zonePtr->type() != domainType::fluid)
            continue;

        const label iZone = zonePtr->index();

        for (const interface* interf : zonePtr->interfacesRef())
        {
            if (!interf->isFluidSolidType())
                continue;

            if (!interf->isRobinNeumann())
                continue;

            const auto* fluidSideInfoPtr = interf->interfaceSideInfoPtr(iZone);

            // RN a_n under-relaxation factor of this interface (1 = none)
            const scalar rnAccelRelax = interf->rnAccelRelaxation();

            stk::mesh::Selector selAllSides =
                metaData.universal_part() &
                stk::mesh::selectUnion(fluidSideInfoPtr->currentPartVec_);

            stk::mesh::BucketVector const& sideBuckets =
                bulkData.get_buckets(metaData.side_rank(), selAllSides);

            for (auto ib = sideBuckets.begin(); ib != sideBuckets.end(); ++ib)
            {
                stk::mesh::Bucket& sideBucket = **ib;

                MasterElement* meFC =
                    MasterElementRepo::get_surface_master_element(
                        sideBucket.topology());
                const label numScsBip = meFC->numIntPoints_;
                const label nodesPerSide = sideBucket.topology().num_nodes();

                std::vector<scalar> ws_shape(numScsBip * nodesPerSide);
                meFC->shape_fcn(&ws_shape[0]);

                const stk::mesh::Bucket::size_type nSides = sideBucket.size();
                for (stk::mesh::Bucket::size_type iSide = 0; iSide < nSides;
                     ++iSide)
                {
                    stk::mesh::Entity side = sideBucket[iSide];
                    scalar* aWallBcVals =
                        stk::mesh::field_data(aWallBcSideSTKFieldRef, side);
                    const scalar* areaVec =
                        stk::mesh::field_data(exposedAreaVecSTKFieldRef, side);

                    stk::mesh::Entity const* sideNodeRels =
                        bulkData.begin_nodes(side);
                    const label numSideNodes = bulkData.num_nodes(side);

                    // nodal wall acceleration
                    std::vector<scalar> aNodal(numSideNodes * SPATIAL_DIM);
                    for (label ni = 0; ni < numSideNodes; ++ni)
                    {
                        stk::mesh::Entity node = sideNodeRels[ni];
                        const scalar* Dt =
                            stk::mesh::field_data(DtSTKFieldRef, node);
                        const scalar* Dt1 =
                            stk::mesh::field_data(Dt1STKFieldRef, node);
                        const scalar* Dt2 =
                            stk::mesh::field_data(Dt2STKFieldRef, node);

                        if (isBackward)
                        {
                            const scalar* Dt3 =
                                stk::mesh::field_data(*Dt3STKFieldPtr, node);
                            const scalar* Dt4 =
                                stk::mesh::field_data(*Dt4STKFieldPtr, node);

                            for (label j = 0; j < SPATIAL_DIM; ++j)
                            {
                                const scalar bdf2_n = inner_c0 * Dt[j] +
                                                      inner_c1 * Dt1[j] +
                                                      inner_c2 * Dt2[j];
                                const scalar bdf2_n1 = inner_c0 * Dt1[j] +
                                                       inner_c1 * Dt2[j] +
                                                       inner_c2 * Dt3[j];
                                const scalar bdf2_n2 = inner_c0 * Dt2[j] +
                                                       inner_c1 * Dt3[j] +
                                                       inner_c2 * Dt4[j];

                                aNodal[ni * SPATIAL_DIM + j] =
                                    rDeltaT2 *
                                    (coefft * bdf2_n - coefft0 * bdf2_n1 +
                                     coefft00 * bdf2_n2);
                            }
                        }
                        else
                        {
                            // Euler d2dt2: (D - 2*D_old + D_old2)/dt^2
                            for (label j = 0; j < SPATIAL_DIM; ++j)
                            {
                                aNodal[ni * SPATIAL_DIM + j] =
                                    rDeltaT2 * (Dt[j] - 2.0 * Dt1[j] + Dt2[j]);
                            }
                        }
                    }

                    // BIP value projected on the exposed normal
                    for (label ip = 0; ip < numScsBip; ++ip)
                    {
                        const label offSF = ip * nodesPerSide;
                        for (label j = 0; j < SPATIAL_DIM; ++j)
                        {
                            p_aBip[j] = 0.0;
                        }
                        for (label ic = 0; ic < nodesPerSide; ++ic)
                        {
                            const scalar r = ws_shape[offSF + ic];
                            for (label j = 0; j < SPATIAL_DIM; ++j)
                            {
                                p_aBip[j] += r * aNodal[ic * SPATIAL_DIM + j];
                            }
                        }

                        scalar aMagSq = 0.0;
                        for (label j = 0; j < SPATIAL_DIM; ++j)
                        {
                            const scalar axj = areaVec[ip * SPATIAL_DIM + j];
                            aMagSq += axj * axj;
                        }
                        const scalar invAMag =
                            1.0 / (std::sqrt(aMagSq) + SMALL);

                        scalar an = 0.0;
                        for (label j = 0; j < SPATIAL_DIM; ++j)
                        {
                            an += p_aBip[j] * areaVec[ip * SPATIAL_DIM + j] *
                                  invAMag;
                        }

                        // a = beta*a_new + (1 - beta)*a_prev
                        aWallBcVals[ip] =
                            (hasAWallBc_ && rnAccelRelax < 1.0)
                                ? rnAccelRelax * an +
                                      (1.0 - rnAccelRelax) * aWallBcVals[ip]
                                : an;
                    }
                }
            }
        }
    }

    hasAWallBc_ = true;
}

} // namespace accel
