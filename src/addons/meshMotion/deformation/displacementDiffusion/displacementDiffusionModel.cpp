// File : displacementDiffusionModel.cpp
// Created : Fri Feb 14 2025 12:55:24 (+0100)
// Author : Mhamad Mahdi Alloush
// Description:
// Copyright (c) 2025 CCFNUM, Lucerne University of Applied Sciences and Arts.
// SPDX-License-Identifier: BSD-3-Clause

#include "displacementDiffusionModel.h"
#include "meshMotion.h"
#include "simulation.h"

namespace accel
{

void crossProduct(scalar* force, scalar* cross, scalar* rad)
{
    cross[0] = rad[1] * force[2] - rad[2] * force[1];
    cross[1] = -(rad[0] * force[2] - rad[2] * force[0]);
    cross[2] = rad[0] * force[1] - rad[1] * force[0];
}

displacementDiffusionModel::displacementDiffusionModel(realm* realm)
    : model(realm)
{
    // create field instances
    DRef();
}

void displacementDiffusionModel::initializeDisplacement(
    const std::shared_ptr<domain> domain)
{
    DRef().initializeField(domain->index());

    // model-based side updates
    updateDisplacementSideFields_(domain);
}

void displacementDiffusionModel::updateDisplacement(
    const std::shared_ptr<domain> domain)
{
    // model-based side updates
    updateDisplacementSideFields_(domain);
}

void displacementDiffusionModel::updateDisplacementSideFields_(
    const std::shared_ptr<domain> domain)
{
#ifdef HAS_INTERFACE
    // Interface
    for (const interface* interf : domain->interfacesRef())
    {
        updateDisplacementInterfaceSideFieldDeformation_(
            domain, interf->interfaceSideInfoPtr(domain->index()));
    }
#endif /* HAS_INTERFACE */

    // Boundary
    for (label iBoundary = 0;
         iBoundary < this->meshRef().zonePtr(domain->index())->nBoundaries();
         iBoundary++)
    {
        const auto* boundary =
            this->meshRef().zonePtr(domain->index())->boundaryPtr(iBoundary);

        boundaryPhysicalType physicalType = boundary->type();

        const auto& bcType =
            this->DRef()
                .boundaryConditionRef(domain->index(), iBoundary)
                .type();

        switch (physicalType)
        {
            case boundaryPhysicalType::symmetry:
                break;

            default:
                {
                    switch (bcType)
                    {
                        case boundaryConditionType::specifiedValue:
                            {
                                updateDisplacementBoundarySideFieldSpecifiedValue_(
                                    domain, boundary);
                            }
                            break;

                        case boundaryConditionType::periodicDisplacement:
                            {
                                updateDisplacementBoundarySideFieldPeriodicDisplacement_(
                                    domain, boundary);
                            }
                            break;

                        case boundaryConditionType::rigidBodySolution:
                            {
                                updateDisplacementBoundarySideFieldRigidBodySolution_(
                                    domain, boundary);
                            }
                            break;

                        default:
                            break;
                    }
                }
                break;
        }
    }
}

void displacementDiffusionModel::
    updateDisplacementBoundarySideFieldSpecifiedValue_(
        const std::shared_ptr<domain> domain,
        const boundary* boundary)
{
    label relativeDisp = domain->zonePtr()
                             ->deformationRef()
                             .displacementRelativeToPreviousMesh();

    auto& bc =
        this->DRef().boundaryConditionRef(domain->index(), boundary->index());
    auto& data = bc.template data<SPATIAL_DIM>("value");

    switch (data.type())
    {
        case inputDataType::null:
            break;

        case inputDataType::constant:
            {
                std::vector<scalar> constValue(SPATIAL_DIM);
                std::copy(data.value(),
                          data.value() + SPATIAL_DIM,
                          constValue.begin());

                // Get fields
                auto& DSTKFieldRef = this->DRef().stkFieldRef();
                auto& nodeSideDSTKFieldRef =
                    this->DRef().nodeSideFieldRef().stkFieldRef();

                const auto& DtSTKFieldRefOld =
                    this->DtRef().prevTimeRef().stkFieldRef();

                bool correctedBoundaryNodeValues =
                    this->DRef().correctedBoundaryNodeValues();

                // select all nodes relevant to
                // the node side field
                stk::mesh::Selector selAllNodes =
                    this->meshRef().metaDataRef().universal_part() &
                    stk::mesh::selectUnion(boundary->parts());
                stk::mesh::BucketVector const& sideNodeBuckets =
                    this->meshRef().bulkDataRef().get_buckets(
                        stk::topology::NODE_RANK, selAllNodes);
                for (stk::mesh::BucketVector::const_iterator ib =
                         sideNodeBuckets.begin();
                     ib != sideNodeBuckets.end();
                     ++ib)
                {
                    stk::mesh::Bucket& sideNodeBucket = **ib;
                    const stk::mesh::Bucket::size_type nSideNodesPerBucket =
                        sideNodeBucket.size();
                    scalar* snDb = stk::mesh::field_data(nodeSideDSTKFieldRef,
                                                         sideNodeBucket);
                    scalar* Db =
                        stk::mesh::field_data(DSTKFieldRef, sideNodeBucket);
                    const scalar* DtbOld =
                        stk::mesh::field_data(DtSTKFieldRefOld, sideNodeBucket);

                    for (stk::mesh::Bucket::size_type iSideNode = 0;
                         iSideNode < nSideNodesPerBucket;
                         ++iSideNode)
                    {
                        for (label i = 0; i < SPATIAL_DIM; i++)
                        {
                            snDb[SPATIAL_DIM * iSideNode + i] =
                                constValue[i] -
                                relativeDisp *
                                    DtbOld[SPATIAL_DIM * iSideNode + i];
                        }

                        // override internal node values
                        if (correctedBoundaryNodeValues)
                        {
                            for (label i = 0; i < SPATIAL_DIM; i++)
                            {
                                Db[SPATIAL_DIM * iSideNode + i] =
                                    snDb[SPATIAL_DIM * iSideNode + i];
                            }
                        }
                    }
                }

                // Update side field on the
                // current boundary
                this->DRef().sideFieldRef().interpolate(
                    this->DRef().nodeSideFieldRef(),
                    domain->index(),
                    boundary->index(),
                    this->DRef().isShifted());
            }
            break;

        case inputDataType::expression:
            {
                typedef exprtk::symbol_table<scalar> symbol_table_t;
                typedef exprtk::expression<scalar> expression_t;
                typedef exprtk::parser<scalar> parser_t;
                std::vector<expression_t> expression_list;

                symbol_table_t symbol_table;
                symbol_table.add_constants();

                // declare function vars
                scalar t, x, y, z;
                symbol_table.add_variable("t", t);
                symbol_table.add_variable("x", x);
                symbol_table.add_variable("y", y);
                symbol_table.add_variable("z", z);

                expression_t componentExpression;
                componentExpression.register_symbol_table(symbol_table);

                parser_t parser;

                for (label i = 0; i < static_cast<label>(SPATIAL_DIM); i++)
                {
                    if (parser.compile(data.expression()[i],
                                       componentExpression))
                    {
                        expression_list.push_back(componentExpression);
                    }
                    else
                    {
                        errorMsg("Error in the "
                                 "expression "
                                 "provided for "
                                 "field " +
                                 this->DRef().name() + ": " +
                                 data.expression()[i]);
                    }
                }

                const auto& coordsSTKFieldRef =
                    *this->meshRef().metaDataRef().get_field<scalar>(
                        stk::topology::NODE_RANK,
                        domain->zonePtr()
                                ->deformationRef()
                                .displacementRelativeToPreviousMesh()
                            ? mesh::coordinates_ID
                            : mesh::original_coordinates_ID);

                bool correctedBoundaryNodeValues =
                    this->DRef().correctedBoundaryNodeValues();

                // Get fields
                auto& DSTKFieldRef = this->DRef().stkFieldRef();
                auto& nodeSideDSTKFieldRef =
                    this->DRef().nodeSideFieldRef().stkFieldRef();

                const auto& DSTKFieldRefOld =
                    this->DRef().prevTimeRef().stkFieldRef();

                // select all nodes relevant to
                // the node side field
                stk::mesh::Selector selAllNodes =
                    this->meshRef().metaDataRef().universal_part() &
                    stk::mesh::selectUnion(boundary->parts());
                stk::mesh::BucketVector const& sideNodeBuckets =
                    this->meshRef().bulkDataRef().get_buckets(
                        stk::topology::NODE_RANK, selAllNodes);
                for (stk::mesh::BucketVector::const_iterator ib =
                         sideNodeBuckets.begin();
                     ib != sideNodeBuckets.end();
                     ++ib)
                {
                    stk::mesh::Bucket& sideNodeBucket = **ib;
                    const stk::mesh::Bucket::size_type nSideNodesPerBucket =
                        sideNodeBucket.size();
                    scalar* snDb = stk::mesh::field_data(nodeSideDSTKFieldRef,
                                                         sideNodeBucket);
                    scalar* Db =
                        stk::mesh::field_data(DSTKFieldRef, sideNodeBucket);
                    scalar* DbOld =
                        stk::mesh::field_data(DSTKFieldRefOld, sideNodeBucket);
                    for (stk::mesh::Bucket::size_type iSideNode = 0;
                         iSideNode < nSideNodesPerBucket;
                         ++iSideNode)
                    {
                        for (label i = 0; i < SPATIAL_DIM; i++)
                        {
                            snDb[SPATIAL_DIM * iSideNode + i] =
                                expression_list[i].value() -
                                relativeDisp *
                                    DbOld[SPATIAL_DIM * iSideNode + i];
                        }

                        // override internal
                        // node values
                        if (correctedBoundaryNodeValues)
                        {
                            for (label i = 0; i < SPATIAL_DIM; i++)
                            {
                                Db[SPATIAL_DIM * iSideNode + i] =
                                    snDb[SPATIAL_DIM * iSideNode + i];
                            }
                        }
                    }
                }

                // Update side field on the
                // current boundary
                this->DRef().sideFieldRef().interpolate(
                    this->DRef().nodeSideFieldRef(),
                    domain->index(),
                    boundary->index(),
                    this->DRef().isShifted());
            }
            break;

        case inputDataType::timeTable:
            {
                // Set current value to that of
                // the starting time of
                // simulation
                auto interpolatedValue =
                    data.interpolate(this->meshRef().controlsRef().time);

                bool correctedBoundaryNodeValues =
                    this->DRef().correctedBoundaryNodeValues();

                // Get fields
                auto& DSTKFieldRef = this->DRef().stkFieldRef();
                auto& nodeSideDSTKFieldRef =
                    this->DRef().nodeSideFieldRef().stkFieldRef();

                const auto& DSTKFieldRefOld =
                    this->DRef().prevTimeRef().stkFieldRef();

                // select all nodes relevant to
                // the node side field
                stk::mesh::Selector selAllNodes =
                    this->meshRef().metaDataRef().universal_part() &
                    stk::mesh::selectUnion(boundary->parts());
                stk::mesh::BucketVector const& sideNodeBuckets =
                    this->meshRef().bulkDataRef().get_buckets(
                        stk::topology::NODE_RANK, selAllNodes);
                for (stk::mesh::BucketVector::const_iterator ib =
                         sideNodeBuckets.begin();
                     ib != sideNodeBuckets.end();
                     ++ib)
                {
                    stk::mesh::Bucket& sideNodeBucket = **ib;
                    const stk::mesh::Bucket::size_type nSideNodesPerBucket =
                        sideNodeBucket.size();
                    scalar* snDb = stk::mesh::field_data(nodeSideDSTKFieldRef,
                                                         sideNodeBucket);
                    scalar* Db =
                        stk::mesh::field_data(DSTKFieldRef, sideNodeBucket);
                    scalar* DbOld =
                        stk::mesh::field_data(DSTKFieldRefOld, sideNodeBucket);
                    for (stk::mesh::Bucket::size_type iSideNode = 0;
                         iSideNode < nSideNodesPerBucket;
                         ++iSideNode)
                    {
                        for (label i = 0; i < SPATIAL_DIM; i++)
                        {
                            snDb[SPATIAL_DIM * iSideNode + i] =
                                interpolatedValue[i] -
                                relativeDisp *
                                    DbOld[SPATIAL_DIM * iSideNode + i];
                        }

                        // override internal
                        // node values
                        if (correctedBoundaryNodeValues)
                        {
                            for (label i = 0; i < SPATIAL_DIM; i++)
                            {
                                Db[SPATIAL_DIM * iSideNode + i] =
                                    snDb[SPATIAL_DIM * iSideNode + i];
                            }
                        }
                    }
                }

                // Update side field on the
                // current boundary
                this->DRef().sideFieldRef().interpolate(
                    this->DRef().nodeSideFieldRef(),
                    domain->index(),
                    boundary->index(),
                    this->DRef().isShifted());
            }
            break;

        case inputDataType::profileData:
            {
                errorMsg("profile data not "
                         "provided yet");
            }
            break;
    }
}

void displacementDiffusionModel::
    updateDisplacementBoundarySideFieldPeriodicDisplacement_(
        const std::shared_ptr<domain> domain,
        const boundary* boundary)
{
    label relativeDisp = domain->zonePtr()
                             ->deformationRef()
                             .displacementRelativeToPreviousMesh();

    // get boundary info
    const auto& bc =
        this->DRef().boundaryConditionRef(domain->index(), boundary->index());
    const auto& frequencyData = bc.data<1>("frequency");
    const auto& valueData = bc.data<SPATIAL_DIM>("value");

    auto& DSTKFieldRef = this->DRef().stkFieldRef();
    auto& DtSTKFieldRefOld = this->DtRef().prevTimeRef().stkFieldRef();
    auto& nodeSideDSTKFieldRef = this->DRef().nodeSideFieldRef();

    bool correctBoundaryNodeValues = this->DRef().correctedBoundaryNodeValues();

    // get time info
    scalar t = this->meshRef().controlsRef().time;

    std::vector<scalar> specValue(SPATIAL_DIM);

    scalar f = *frequencyData.value();
    for (label i = 0; i < SPATIAL_DIM; i++)
    {
        specValue[i] = valueData.value()[i] * std::sin(2.0 * M_PI * f * t);
    }

    // select all nodes relevant to the node
    // side field
    stk::mesh::Selector selAllNodes =
        this->meshRef().metaDataRef().universal_part() &
        stk::mesh::selectUnion(boundary->parts());
    stk::mesh::BucketVector const& sideNodeBuckets =
        this->meshRef().bulkDataRef().get_buckets(stk::topology::NODE_RANK,
                                                  selAllNodes);
    for (stk::mesh::BucketVector::const_iterator ib = sideNodeBuckets.begin();
         ib != sideNodeBuckets.end();
         ++ib)
    {
        stk::mesh::Bucket& sideNodeBucket = **ib;
        const stk::mesh::Bucket::size_type nSideNodesPerBucket =
            sideNodeBucket.size();

        scalar* nsDb = stk::mesh::field_data(nodeSideDSTKFieldRef.stkFieldRef(),
                                             sideNodeBucket);
        scalar* Db = stk::mesh::field_data(DSTKFieldRef, sideNodeBucket);
        scalar* DtbOld =
            stk::mesh::field_data(DtSTKFieldRefOld, sideNodeBucket);

        for (stk::mesh::Bucket::size_type iSideNode = 0;
             iSideNode < nSideNodesPerBucket;
             ++iSideNode)
        {
            for (label i = 0; i < SPATIAL_DIM; i++)
            {
                nsDb[SPATIAL_DIM * iSideNode + i] =
                    specValue[i] -
                    relativeDisp * DtbOld[SPATIAL_DIM * iSideNode + i];
            }

            // override internal node values
            if (correctBoundaryNodeValues)
            {
                for (label i = 0; i < SPATIAL_DIM; i++)
                {
                    Db[SPATIAL_DIM * iSideNode + i] =
                        nsDb[SPATIAL_DIM * iSideNode + i];
                }
            }
        }
    }

    // Update side field on the current
    // boundary
    this->DRef().sideFieldRef().interpolate(this->DRef().nodeSideFieldRef(),
                                            domain->index(),
                                            boundary->index(),
                                            this->DRef().isShifted());
}

void displacementDiffusionModel::
    updateDisplacementBoundarySideFieldRigidBodySolution_(
        const std::shared_ptr<domain> domain,
        const boundary* boundary)
{
    label relativeDisp = domain->zonePtr()
                             ->deformationRef()
                             .displacementRelativeToPreviousMesh();

    // get fields
    const auto& orgCoordsSTKFieldRef =
        *this->meshRef().metaDataRef().get_field<scalar>(
            stk::topology::NODE_RANK, mesh::original_coordinates_ID);
    auto& DSTKFieldRef = this->DRef().stkFieldRef();
    auto& DtSTKFieldRefOld = this->DtRef().prevTimeRef().stkFieldRef();
    auto& nodeSideDSTKFieldRef = this->DRef().nodeSideFieldRef().stkFieldRef();

    // get boundary info
    const auto& bc =
        this->DRef().boundaryConditionRef(domain->index(), boundary->index());
    const std::string rigidBodyName = bc.rawStringValue("rigid_body");

    rigidBody* rigidBodyPtr =
        realmRef().simulationRef().rigidBodyPtr(rigidBodyName);

    // get initial centroid position of the
    // boundary
    utils::vector centroid0 = boundary->stats().centroid0_;

    // get total centroid displacement
    utils::vector dx = rigidBodyPtr->dx();

    // to be filled
    utils::vector F;
    utils::vector M;

    // calculate force and moment
    calculateSurfaceForceAndMoment_(boundary, centroid0 + dx, F, M);

    // get time info
    scalar dt = this->meshRef().controlsRef().getTimestep();

    // calculate body motion based on the force
    // and moment
    rigidBodyPtr->update(F, M, dt);

    // Calculate full Quaternion
    std::vector<scalar> q(4, 0.0);
    utils::matrix rotMat;

    q = rigidBody::getQuatFromEulerxyz(rigidBodyPtr->theta());
    rotMat = rigidBody::getRotmatFromQuat(q);
    scalar* p_rotMat = rotMat.data();

    // local space; current coords and
    // rotated coords; generalized for 2D and 3D
    scalar mcX[3] = {0.0, 0.0, 0.0};
    scalar rcX[3] = {0.0, 0.0, 0.0};

    bool correctBoundaryNodeValues = this->DRef().correctedBoundaryNodeValues();

    // select all nodes relevant to the node
    // side field
    stk::mesh::Selector selAllNodes =
        this->meshRef().metaDataRef().universal_part() &
        stk::mesh::selectUnion(boundary->parts());
    stk::mesh::BucketVector const& sideNodeBuckets =
        this->meshRef().bulkDataRef().get_buckets(stk::topology::NODE_RANK,
                                                  selAllNodes);
    for (stk::mesh::BucketVector::const_iterator ib = sideNodeBuckets.begin();
         ib != sideNodeBuckets.end();
         ++ib)
    {
        stk::mesh::Bucket& sideNodeBucket = **ib;
        const stk::mesh::Bucket::size_type nSideNodesPerBucket =
            sideNodeBucket.size();

        scalar* nsDb =
            stk::mesh::field_data(nodeSideDSTKFieldRef, sideNodeBucket);
        scalar* Db = stk::mesh::field_data(DSTKFieldRef, sideNodeBucket);
        scalar* DtbOld =
            stk::mesh::field_data(DtSTKFieldRefOld, sideNodeBucket);

        for (stk::mesh::Bucket::size_type iSideNode = 0;
             iSideNode < nSideNodesPerBucket;
             ++iSideNode)
        {
            const scalar* orgCoords = stk::mesh::field_data(
                orgCoordsSTKFieldRef, sideNodeBucket, iSideNode);

            // load the current and model coords
            for (label i = 0; i < SPATIAL_DIM; ++i)
            {
                mcX[i] = orgCoords[i];
            }

            const scalar cX = mcX[0] - centroid0[0];
            const scalar cY = mcX[1] - centroid0[1];
            const scalar cZ = mcX[2] - centroid0[2];

            // rotated model coordinates;
            // converted to displacement; add
            // back in centroid

            rcX[0] = p_rotMat[0] * cX + p_rotMat[1] * cY + p_rotMat[2] * cZ -
                     mcX[0] + centroid0[0] + dx[0];
            rcX[1] = p_rotMat[3] * cX + p_rotMat[4] * cY + p_rotMat[5] * cZ -
                     mcX[1] + centroid0[1] + dx[1];
            rcX[2] = p_rotMat[6] * cX + p_rotMat[7] * cY + p_rotMat[8] * cZ -
                     mcX[2] + centroid0[2] + dx[2];

            for (label i = 0; i < SPATIAL_DIM; i++)
            {
                nsDb[SPATIAL_DIM * iSideNode + i] =
                    rcX[i] - relativeDisp * DtbOld[SPATIAL_DIM * iSideNode + i];
            }

            // override internal node values
            if (correctBoundaryNodeValues)
            {
                for (label i = 0; i < SPATIAL_DIM; i++)
                {
                    Db[SPATIAL_DIM * iSideNode + i] =
                        nsDb[SPATIAL_DIM * iSideNode + i];
                }
            }
        }
    }

    // Update side field on the current boundary
    this->DRef().sideFieldRef().interpolate(this->DRef().nodeSideFieldRef(),
                                            domain->index(),
                                            boundary->index(),
                                            this->DRef().isShifted());
}

#ifdef HAS_INTERFACE
void displacementDiffusionModel::
    updateDisplacementInterfaceSideFieldDeformation_(
        const std::shared_ptr<domain> domain,
        const interfaceSideInfo* interfaceSideInfoPtr)
{
    if (!interfaceSideInfoPtr->interfPtr()->isFluidSolidType())
        return;

    label relativeDisp = domain->zonePtr()
                             ->deformationRef()
                             .displacementRelativeToPreviousMesh();

    // Get fields
    auto& DSTKFieldRef = this->DRef().stkFieldRef();
    const auto& DtSTKFieldRefOld = this->DtRef().prevTimeRef().stkFieldRef();

    // select all nodes relevant to
    // the node side field
    stk::mesh::Selector selAllSideNodes =
        this->meshRef().metaDataRef().universal_part() &
        stk::mesh::selectUnion(interfaceSideInfoPtr->currentPartVec_);
    stk::mesh::BucketVector const& sideNodeBuckets =
        this->meshRef().bulkDataRef().get_buckets(stk::topology::NODE_RANK,
                                                  selAllSideNodes);
    for (stk::mesh::BucketVector::const_iterator ib = sideNodeBuckets.begin();
         ib != sideNodeBuckets.end();
         ++ib)
    {
        stk::mesh::Bucket& sideNodeBucket = **ib;
        const stk::mesh::Bucket::size_type nSideNodesPerBucket =
            sideNodeBucket.size();
        scalar* Db = stk::mesh::field_data(DSTKFieldRef, sideNodeBucket);
        const scalar* DtbOld =
            stk::mesh::field_data(DtSTKFieldRefOld, sideNodeBucket);

        for (stk::mesh::Bucket::size_type iSideNode = 0;
             iSideNode < nSideNodesPerBucket;
             ++iSideNode)
        {
            for (label i = 0; i < SPATIAL_DIM; i++)
            {
                Db[SPATIAL_DIM * iSideNode + i] -=
                    relativeDisp * DtbOld[SPATIAL_DIM * iSideNode + i];
            }
        }
    }
}
#endif /* HAS_INTERFACE */

void displacementDiffusionModel::calculateSurfaceForceAndMoment_(
    const boundary* boundary,
    const utils::vector& center,
    utils::vector& force,
    utils::vector& moment)
{
    const auto& mesh = boundary->zonePtr()->meshRef();
    const stk::mesh::MetaData& metaData = mesh.metaDataRef();
    const stk::mesh::BulkData& bulkData = mesh.bulkDataRef();

    // bip values
    std::vector<scalar> nx(SPATIAL_DIM);

    // pointers to fixed values
    scalar* p_nx = &nx[0];

    // nodal fields to gather
    std::vector<scalar> ws_p;

    // master element
    std::vector<scalar> ws_face_shape_function;

    // get fields
    const auto& pSTKFieldRef = this->pRef().stkFieldRef();
    const auto& wallShearStressSTKFieldRef =
        this->wallShearStressRef().stkFieldRef();

    // Get geometric fields
    const auto& exposedAreaVecSTKFieldRef = *metaData.get_field<scalar>(
        metaData.side_rank(),
        boundary->zonePtr()
                ->deformationRef()
                .displacementRelativeToPreviousMesh()
            ? mesh::exposed_area_vector_ID
            : mesh::original_exposed_area_vector_ID);
    const auto& coordsSTKFieldRef =
        *this->meshRef().metaDataRef().get_field<scalar>(
            stk::topology::NODE_RANK,
            boundary->zonePtr()
                    ->deformationRef()
                    .displacementRelativeToPreviousMesh()
                ? mesh::coordinates_ID
                : mesh::original_coordinates_ID);

    // local force
    scalar l_force_moment[9] = {0.0};

    // work force
    scalar ws_p_force[3] = {};
    scalar ws_v_force[3] = {};
    scalar ws_t_force[3] = {};
    scalar ws_moment[3] = {};
    scalar ws_radius[3] = {};

    // define some common selectors
    stk::mesh::Selector selOwnedSides =
        metaData.locally_owned_part() &
        stk::mesh::selectUnion(boundary->parts());

    // shifted ip's for field?
    const bool isPShifted = this->pRef().isShifted();

    stk::mesh::BucketVector const& sideBuckets =
        bulkData.get_buckets(metaData.side_rank(), selOwnedSides);
    for (stk::mesh::BucketVector::const_iterator ib = sideBuckets.begin();
         ib != sideBuckets.end();
         ++ib)
    {
        stk::mesh::Bucket& sideBucket = **ib;

        // face master element
        MasterElement* meFC = MasterElementRepo::get_surface_master_element(
            sideBucket.topology());
        const label nodesPerFace = meFC->nodesPerElement_;
        const label numScsBip = meFC->numIntPoints_;

        // mapping from ip to nodes for this ordinal
        const label* faceIpNodeMap = meFC->ipNodeMap();

        // algorithm related; element
        ws_p.resize(nodesPerFace);
        ws_face_shape_function.resize(numScsBip * nodesPerFace);

        // pointers
        scalar* p_p = &ws_p[0];
        scalar* p_face_shape_function = &ws_face_shape_function[0];

        // shape functions
        if (isPShifted)
        {
            meFC->shifted_shape_fcn(&p_face_shape_function[0]);
        }
        else
        {
            meFC->shape_fcn(&p_face_shape_function[0]);
        }

        const stk::mesh::Bucket::size_type nBoundaryFaces = sideBucket.size();

        for (stk::mesh::Bucket::size_type iSide = 0; iSide < nBoundaryFaces;
             ++iSide)
        {
            // get face
            stk::mesh::Entity side = sideBucket[iSide];

            // face node relations
            stk::mesh::Entity const* sideNodeRels = bulkData.begin_nodes(side);

            //======================================
            // gather nodal data off of face
            //======================================
            for (label ni = 0; ni < nodesPerFace; ++ni)
            {
                stk::mesh::Entity node = sideNodeRels[ni];

                // gather scalars
                p_p[ni] = *stk::mesh::field_data(pSTKFieldRef, node);
            }

            // pointer to face data
            const scalar* areaVec =
                stk::mesh::field_data(exposedAreaVecSTKFieldRef, side);
            const scalar* wallShearStressBip =
                stk::mesh::field_data(wallShearStressSTKFieldRef, side);

            for (label ip = 0; ip < numScsBip; ++ip)
            {
                // offsets
                const label offSetAveraVec = ip * SPATIAL_DIM;
                const label offSetSF_face = ip * nodesPerFace;
                const label localFaceNode = faceIpNodeMap[ip];

                // zero out vector quantities; squeeze in aMag
                scalar aMag = 0.0;
                for (label j = 0; j < SPATIAL_DIM; ++j)
                {
                    const scalar axj = areaVec[offSetAveraVec + j];
                    aMag += axj * axj;
                }
                aMag = std::sqrt(aMag);

                // interpolate to bip
                scalar pBip = 0.0;
                for (label ic = 0; ic < nodesPerFace; ++ic)
                {
                    const scalar r = p_face_shape_function[offSetSF_face + ic];
                    pBip += r * p_p[ic];
                }

                // form unit normal
                for (label j = 0; j < SPATIAL_DIM; ++j)
                {
                    p_nx[j] = areaVec[offSetAveraVec + j] / aMag;
                }

                // extract nodal fields
                stk::mesh::Entity node = sideNodeRels[localFaceNode];
                const scalar* coord =
                    stk::mesh::field_data(coordsSTKFieldRef, node);

                // load radius; assemble force -sigma_ij*njdS
                for (label i = 0; i < SPATIAL_DIM; ++i)
                {
                    const scalar ai = areaVec[offSetAveraVec + i];
                    ws_radius[i] = coord[i] - center[i];
                    ws_p_force[i] = pBip * ai;
                    ws_v_force[i] =
                        wallShearStressBip[SPATIAL_DIM * ip + i] * aMag;
                    ws_t_force[i] = ws_p_force[i] + ws_v_force[i];
                }

                crossProduct(&ws_t_force[0], &ws_moment[0], &ws_radius[0]);

                // assemble for and moment
                for (label j = 0; j < 3; ++j)
                {
                    l_force_moment[j] += ws_p_force[j];
                    l_force_moment[j + 3] += ws_v_force[j];
                    l_force_moment[j + 6] += ws_moment[j];
                }
            }
        }
    }

    // parallel assemble and output
    scalar g_force_moment[9] = {};

    // Parallel assembly of L2
    stk::all_reduce_sum(
        bulkData.parallel(), &l_force_moment[0], &g_force_moment[0], 9);

    // Copy to arg force and moment
    for (label j = 0; j < 3; ++j)
    {
        force[j] = g_force_moment[j] + g_force_moment[j + 3];
        moment[j] = g_force_moment[j + 6];
    }
}

} /* namespace accel */
