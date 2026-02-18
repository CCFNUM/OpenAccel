// File : velocity.cpp
// Created : Tue Apr 20 2024 12:55:24 (+0100)
// Author : Mhamad Mahdi Alloush
// Description:
// Copyright (c) 2024 CCFNUM, Lucerne University of Applied Sciences and Arts.
// SPDX-License-Identifier: BSD-3-Clause

// code
#include "velocity.h"
#include "controls.h"
#include "mesh.h"
#include "realm.h"
#include "simulation.h"
#include "zoneTransformation.h"

namespace accel
{

velocity::velocity(realm* realmPtr,
                   const std::string name,
                   unsigned numberOfStates,
                   bool highResolution)
    : nodeVectorField(realmPtr,
                      name,
                      numberOfStates,
                      true,
                      highResolution,
                      true,
                      false)
{
    realmPtr->registerRestartField(name);

    interpolationScheme_ =
        meshRef()
            .controlsRef()
            .solverRef()
            .solverControl_.basicSettings_.interpolationSchemeType_
            .velocityInterpolationType_;

    gradientInterpolationScheme_ =
        meshRef()
            .controlsRef()
            .solverRef()
            .solverControl_.basicSettings_.interpolationSchemeType_
            .velocityGradientInterpolationType_;

    // velocity may only apply to fluid domains
    mediumIndependent_ = false;

    // force correct gradient for velocity: remove symmetric contributions to
    // gradient
    correctGradient_ = true;

    // instantiate the reversal flag field and put on pressure boundaries
    for (label iZone = 0; iZone < this->meshPtr()->nZones(); iZone++)
    {
        for (label iBoundary = 0;
             iBoundary < this->meshPtr()->zonePtr(iZone)->nBoundaries();
             iBoundary++)
        {
            boundaryPhysicalType physicalType =
                this->meshPtr()->zonePtr(iZone)->boundaryRef(iBoundary).type();

            switch (physicalType)
            {
                case boundaryPhysicalType::inlet:
                case boundaryPhysicalType::outlet:
                case boundaryPhysicalType::opening:
                    {
                        if (!reversalFlagPtr_)
                        {
                            reversalFlagPtr_ =
                                std::make_unique<sideField<label, 1>>(
                                    realmPtr->meshPtr(), "reversal_flag", 1);
                        }
                        // Put the side field on the corresponding boundary
                        // part
                        for (auto* part : this->meshPtr()
                                              ->zonePtr(iZone)
                                              ->boundaryRef(iBoundary)
                                              .parts())
                        {
                            this->reversalFlagRef().putFieldOnPart(*part);
                        }
                    }
                    break;

                default:
                    break;
            }
        }
    }
}

void velocity::updateBoundarySideField(label iZone, label iBoundary)
{
    assert(this->isZoneSet(iZone));

    boundaryPhysicalType physicalType =
        this->meshPtr()->zonePtr(iZone)->boundaryRef(iBoundary).type();
    const auto& bcType = this->boundaryConditionRef(iZone, iBoundary).type();

    switch (physicalType)
    {
        case boundaryPhysicalType::inlet:
            {
                switch (bcType)
                {
                    case boundaryConditionType::specifiedValue:
                        updateBoundarySideFieldSpecifiedValue(iZone, iBoundary);
                        break;

                    case boundaryConditionType::normalSpeed:
                        updateBoundarySideFieldNormalSpeed(iZone, iBoundary);
                        break;

                    case boundaryConditionType::massFlowRate:
                    case boundaryConditionType::specifiedDirection:
                        updateBoundarySideDirectionFields(iZone, iBoundary);
                        break;

                    default:
                        nodeVectorField::updateBoundarySideField(iZone,
                                                                 iBoundary);
                        break;
                }
            }
            break;

        case boundaryPhysicalType::outlet:
            {
                switch (bcType)
                {
                    case boundaryConditionType::massFlowRate:
                        {
                            // model-related
                        }
                        break;

                    default:
                        nodeVectorField::updateBoundarySideField(iZone,
                                                                 iBoundary);
                        break;
                }
            }
            break;

        case boundaryPhysicalType::opening:
            {
                switch (bcType)
                {
                    case boundaryConditionType::specifiedDirection:
                        updateBoundarySideDirectionFields(iZone, iBoundary);
                        break;

                    default:
                        nodeVectorField::updateBoundarySideField(iZone,
                                                                 iBoundary);
                        break;
                }
            }
            break;

        case boundaryPhysicalType::wall:
            {
                switch (bcType)
                {
                    case boundaryConditionType::noSlip:
                        updateBoundarySideFieldNoSlipWall(iZone, iBoundary);
                        break;

                    case boundaryConditionType::slip:
                        break;

                    default:
                        nodeVectorField::updateBoundarySideField(iZone,
                                                                 iBoundary);
                        break;
                }
            }
            break;

        default:
            nodeVectorField::updateBoundarySideField(iZone, iBoundary);
            break;
    }
}

void velocity::updateBoundarySideFieldSpecifiedValue(label iZone,
                                                     label iBoundary)
{
    const auto& mesh = this->meshRef();
    const stk::mesh::BulkData& bulkData = mesh.bulkDataRef();
    const stk::mesh::MetaData& metaData = mesh.metaDataRef();

    const auto* zone = mesh.zonePtr(iZone);
    const auto* boundary = zone->boundaryPtr(iBoundary);

    if (zone->frameRotating())
    {
        switch (zone->transformationRef().type())
        {
            case meshMotionType::stationary:
            case meshMotionType::translating:
                errorMsg("Must not reach here");
                break;

            case meshMotionType::rotating:
                {
                    if (boundary->frameType() ==
                        boundaryRelativeFrameType::absolute)
                    {
                        nodeVectorField::updateBoundarySideFieldSpecifiedValue(
                            iZone, iBoundary);
                    }
                    else
                    {
                        // the velocity at the boundary is specified in the
                        // zone's rotating frame of reference
                        auto& bc = this->boundaryConditionRef(iZone, iBoundary);
                        auto& data = bc.template data<SPATIAL_DIM>("value");

                        // pointers to rotation data
                        const scalar* p_mat = zone->transformationRef()
                                                  .rotation()
                                                  .coriolisMatrix_.data();
                        const scalar* p_ori =
                            zone->transformationRef().rotation().origin_.data();

                        switch (data.type())
                        {
                            case inputDataType::null:
                                break;

                            case inputDataType::constant:
                            case inputDataType::timeTable:
                                {
                                    std::array<scalar, SPATIAL_DIM> inputValue;
                                    if (data.type() == inputDataType::constant)
                                    {
                                        std::copy(data.value(),
                                                  data.value() + SPATIAL_DIM,
                                                  inputValue.begin());
                                    }
                                    else
                                    {
                                        inputValue = data.interpolate(
                                            this->meshRef().controlsRef().time);
                                    }

                                    // select all nodes relevant to the node
                                    // side field
                                    stk::mesh::Selector selAllNodes =
                                        metaData.universal_part() &
                                        stk::mesh::selectUnion(
                                            boundary->parts());

                                    // get fields
                                    auto& stkFieldRef = this->stkFieldRef();
                                    auto& nodeSideSTKFieldRef =
                                        this->nodeSideFieldRef().stkFieldRef();

                                    // extract geometric fields
                                    const STKScalarField* coordsSTKFieldPtr =
                                        metaData.template get_field<scalar>(
                                            stk::topology::NODE_RANK,
                                            this->getCoordinatesID_(iZone));

                                    stk::mesh::BucketVector const&
                                        sideNodeBuckets = bulkData.get_buckets(
                                            stk::topology::NODE_RANK,
                                            selAllNodes);
                                    for (stk::mesh::BucketVector::const_iterator
                                             ib = sideNodeBuckets.begin();
                                         ib != sideNodeBuckets.end();
                                         ++ib)
                                    {
                                        stk::mesh::Bucket& sideNodeBucket =
                                            **ib;
                                        const stk::mesh::Bucket::size_type
                                            nSideNodesPerBucket =
                                                sideNodeBucket.size();
                                        scalar* snvalueb =
                                            stk::mesh::field_data(
                                                nodeSideSTKFieldRef,
                                                sideNodeBucket);
                                        scalar* valueb = stk::mesh::field_data(
                                            stkFieldRef, sideNodeBucket);
                                        for (stk::mesh::Bucket::size_type
                                                 iSideNode = 0;
                                             iSideNode < nSideNodesPerBucket;
                                             ++iSideNode)
                                        {
                                            const scalar* coords =
                                                stk::mesh::field_data(
                                                    *coordsSTKFieldPtr,
                                                    sideNodeBucket,
                                                    iSideNode);

                                            for (label i = 0;
                                                 i < static_cast<label>(
                                                         SPATIAL_DIM);
                                                 i++)
                                            {
                                                snvalueb[SPATIAL_DIM *
                                                             iSideNode +
                                                         i] = inputValue[i];
                                            }

                                            // add frame velocity
                                            for (label i = 0; i < SPATIAL_DIM;
                                                 i++)
                                            {
                                                for (label j = 0;
                                                     j < SPATIAL_DIM;
                                                     j++)
                                                {
                                                    snvalueb[SPATIAL_DIM *
                                                                 iSideNode +
                                                             i] +=
                                                        p_mat[i * SPATIAL_DIM +
                                                              j] *
                                                        (coords[j] - p_ori[j]);
                                                }
                                            }

                                            // override internal node values
                                            if (correctedBoundaryNodeValues_)
                                            {
                                                for (label i = 0;
                                                     i < SPATIAL_DIM;
                                                     i++)
                                                {
                                                    valueb[SPATIAL_DIM *
                                                               iSideNode +
                                                           i] =
                                                        snvalueb[SPATIAL_DIM *
                                                                     iSideNode +
                                                                 i];
                                                }
                                            }
                                        }
                                    }

                                    // Update side field on the current boundary
                                    this->sideFieldRef().interpolate(
                                        this->nodeSideFieldRef(),
                                        iZone,
                                        iBoundary,
                                        this->isShifted());
                                }
                                break;

                            case inputDataType::expression:
                                {
                                    typedef exprtk::symbol_table<scalar>
                                        symbol_table_t;
                                    typedef exprtk::expression<scalar>
                                        expression_t;
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
                                    componentExpression.register_symbol_table(
                                        symbol_table);

                                    parser_t parser;

                                    for (label i = 0; i < SPATIAL_DIM; i++)
                                    {
                                        if (parser.compile(data.expression()[i],
                                                           componentExpression))
                                        {
                                            expression_list.push_back(
                                                componentExpression);
                                        }
                                        else
                                        {
                                            errorMsg("Error in the expression "
                                                     "provided "
                                                     "for field " +
                                                     this->name() + ": " +
                                                     data.expression()[i]);
                                        }
                                    }

                                    // set time value
                                    t = this->meshRef().controlsRef().time;

                                    // select all nodes relevant to the node
                                    // side field
                                    stk::mesh::Selector selAllNodes =
                                        metaData.universal_part() &
                                        stk::mesh::selectUnion(
                                            boundary->parts());

                                    // get fields
                                    auto& stkFieldRef = this->stkFieldRef();
                                    auto& nodeSideSTKFieldRef =
                                        this->nodeSideFieldRef().stkFieldRef();

                                    // extract geometric fields
                                    const STKScalarField* coordsSTKFieldPtr =
                                        metaData.template get_field<scalar>(
                                            stk::topology::NODE_RANK,
                                            this->getCoordinatesID_(iZone));

                                    stk::mesh::BucketVector const&
                                        sideNodeBuckets = bulkData.get_buckets(
                                            stk::topology::NODE_RANK,
                                            selAllNodes);
                                    for (stk::mesh::BucketVector::const_iterator
                                             ib = sideNodeBuckets.begin();
                                         ib != sideNodeBuckets.end();
                                         ++ib)
                                    {
                                        stk::mesh::Bucket& sideNodeBucket =
                                            **ib;
                                        const stk::mesh::Bucket::size_type
                                            nSideNodesPerBucket =
                                                sideNodeBucket.size();
                                        scalar* snvalueb =
                                            stk::mesh::field_data(
                                                nodeSideSTKFieldRef,
                                                sideNodeBucket);
                                        scalar* valueb = stk::mesh::field_data(
                                            stkFieldRef, sideNodeBucket);
                                        for (stk::mesh::Bucket::size_type
                                                 iSideNode = 0;
                                             iSideNode < nSideNodesPerBucket;
                                             ++iSideNode)
                                        {
                                            const scalar* coords =
                                                stk::mesh::field_data(
                                                    *coordsSTKFieldPtr,
                                                    sideNodeBucket,
                                                    iSideNode);

#if SPATIAL_DIM == 3
                                            x = coords[0];
                                            y = coords[1];
                                            z = coords[2];
#elif SPATIAL_DIM == 2
                                            x = coords[0];
                                            y = coords[1];
#endif

                                            for (label i = 0; i < SPATIAL_DIM;
                                                 i++)
                                            {
                                                snvalueb[SPATIAL_DIM *
                                                             iSideNode +
                                                         i] =
                                                    expression_list[i].value();
                                            }

                                            // add frame velocity
                                            for (label i = 0; i < SPATIAL_DIM;
                                                 i++)
                                            {
                                                for (label j = 0;
                                                     j < SPATIAL_DIM;
                                                     j++)
                                                {
                                                    snvalueb[SPATIAL_DIM *
                                                                 iSideNode +
                                                             i] +=
                                                        p_mat[i * SPATIAL_DIM +
                                                              j] *
                                                        (coords[j] - p_ori[j]);
                                                }
                                            }

                                            // override internal node values
                                            if (correctedBoundaryNodeValues_)
                                            {
                                                for (label i = 0;
                                                     i < SPATIAL_DIM;
                                                     i++)
                                                {
                                                    valueb[SPATIAL_DIM *
                                                               iSideNode +
                                                           i] =
                                                        snvalueb[SPATIAL_DIM *
                                                                     iSideNode +
                                                                 i];
                                                }
                                            }
                                        }
                                    }

                                    // Update side field on the current boundary
                                    this->sideFieldRef().interpolate(
                                        this->nodeSideFieldRef(),
                                        iZone,
                                        iBoundary,
                                        this->isShifted());
                                }
                                break;

                            case inputDataType::profileData:
                                {
                                    errorMsg("profile data not provided yet");
                                }
                                break;
                        }
                    }
                }
                break;
        }
    }
    else if (zone->meshTransforming())
    {
        switch (zone->transformationRef().type())
        {
            case meshMotionType::translating:
            case meshMotionType::rotating:
                {
                    if (boundary->frameType() ==
                        boundaryRelativeFrameType::absolute)
                    {
                        nodeVectorField::updateBoundarySideFieldSpecifiedValue(
                            iZone, iBoundary);
                    }
                    else
                    {
                        auto& bc = this->boundaryConditionRef(iZone, iBoundary);
                        auto& data = bc.template data<SPATIAL_DIM>("value");

                        // Get mesh velocity field
                        const auto& UmSTKFieldRef = this->simulationRef()
                                                        .meshMotionRef()
                                                        .UmRef()
                                                        .stkFieldRef();

                        switch (data.type())
                        {
                            case inputDataType::null:
                                break;

                            case inputDataType::constant:
                            case inputDataType::timeTable:
                                {
                                    std::array<scalar, SPATIAL_DIM> inputValue;
                                    if (data.type() == inputDataType::constant)
                                    {
                                        std::copy(data.value(),
                                                  data.value() + SPATIAL_DIM,
                                                  inputValue.begin());
                                    }
                                    else
                                    {
                                        inputValue = data.interpolate(
                                            this->meshRef().controlsRef().time);
                                    }

                                    // select all nodes relevant to the node
                                    // side field
                                    stk::mesh::Selector selAllNodes =
                                        metaData.universal_part() &
                                        stk::mesh::selectUnion(
                                            boundary->parts());

                                    // get fields
                                    auto& stkFieldRef = this->stkFieldRef();
                                    auto& nodeSideSTKFieldRef =
                                        this->nodeSideFieldRef().stkFieldRef();

                                    stk::mesh::BucketVector const&
                                        sideNodeBuckets = bulkData.get_buckets(
                                            stk::topology::NODE_RANK,
                                            selAllNodes);
                                    for (stk::mesh::BucketVector::const_iterator
                                             ib = sideNodeBuckets.begin();
                                         ib != sideNodeBuckets.end();
                                         ++ib)
                                    {
                                        stk::mesh::Bucket& sideNodeBucket =
                                            **ib;
                                        const stk::mesh::Bucket::size_type
                                            nSideNodesPerBucket =
                                                sideNodeBucket.size();
                                        scalar* snvalueb =
                                            stk::mesh::field_data(
                                                nodeSideSTKFieldRef,
                                                sideNodeBucket);
                                        scalar* valueb = stk::mesh::field_data(
                                            stkFieldRef, sideNodeBucket);
                                        for (stk::mesh::Bucket::size_type
                                                 iSideNode = 0;
                                             iSideNode < nSideNodesPerBucket;
                                             ++iSideNode)
                                        {
                                            const scalar* Um =
                                                stk::mesh::field_data(
                                                    UmSTKFieldRef,
                                                    sideNodeBucket,
                                                    iSideNode);

                                            for (label i = 0; i < SPATIAL_DIM;
                                                 i++)
                                            {
                                                snvalueb[SPATIAL_DIM *
                                                             iSideNode +
                                                         i] = inputValue[i];
                                            }

                                            // add mesh velocity
                                            for (label i = 0; i < SPATIAL_DIM;
                                                 i++)
                                            {
                                                snvalueb[SPATIAL_DIM *
                                                             iSideNode +
                                                         i] += Um[i];
                                            }

                                            // override internal node values
                                            if (correctedBoundaryNodeValues_)
                                            {
                                                for (label i = 0;
                                                     i < SPATIAL_DIM;
                                                     i++)
                                                {
                                                    valueb[SPATIAL_DIM *
                                                               iSideNode +
                                                           i] =
                                                        snvalueb[SPATIAL_DIM *
                                                                     iSideNode +
                                                                 i];
                                                }
                                            }
                                        }
                                    }

                                    // Update side field on the current boundary
                                    this->sideFieldRef().interpolate(
                                        this->nodeSideFieldRef(),
                                        iZone,
                                        iBoundary,
                                        this->isShifted());
                                }
                                break;

                            case inputDataType::expression:
                                {
                                    typedef exprtk::symbol_table<scalar>
                                        symbol_table_t;
                                    typedef exprtk::expression<scalar>
                                        expression_t;
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
                                    componentExpression.register_symbol_table(
                                        symbol_table);

                                    parser_t parser;

                                    for (label i = 0; i < SPATIAL_DIM; i++)
                                    {
                                        if (parser.compile(data.expression()[i],
                                                           componentExpression))
                                        {
                                            expression_list.push_back(
                                                componentExpression);
                                        }
                                        else
                                        {
                                            errorMsg("Error in the expression "
                                                     "provided "
                                                     "for field " +
                                                     this->name() + ": " +
                                                     data.expression()[i]);
                                        }
                                    }

                                    // set time value
                                    t = this->meshRef().controlsRef().time;

                                    // select all nodes relevant to the node
                                    // side field
                                    stk::mesh::Selector selAllNodes =
                                        metaData.universal_part() &
                                        stk::mesh::selectUnion(
                                            boundary->parts());

                                    // get fields
                                    auto& stkFieldRef = this->stkFieldRef();
                                    auto& nodeSideSTKFieldRef =
                                        this->nodeSideFieldRef().stkFieldRef();

                                    // extract geometric fields
                                    const STKScalarField* coordsSTKFieldPtr =
                                        metaData.template get_field<scalar>(
                                            stk::topology::NODE_RANK,
                                            this->getCoordinatesID_(iZone));

                                    stk::mesh::BucketVector const&
                                        sideNodeBuckets = bulkData.get_buckets(
                                            stk::topology::NODE_RANK,
                                            selAllNodes);
                                    for (stk::mesh::BucketVector::const_iterator
                                             ib = sideNodeBuckets.begin();
                                         ib != sideNodeBuckets.end();
                                         ++ib)
                                    {
                                        stk::mesh::Bucket& sideNodeBucket =
                                            **ib;
                                        const stk::mesh::Bucket::size_type
                                            nSideNodesPerBucket =
                                                sideNodeBucket.size();
                                        scalar* snvalueb =
                                            stk::mesh::field_data(
                                                nodeSideSTKFieldRef,
                                                sideNodeBucket);
                                        scalar* valueb = stk::mesh::field_data(
                                            stkFieldRef, sideNodeBucket);
                                        for (stk::mesh::Bucket::size_type
                                                 iSideNode = 0;
                                             iSideNode < nSideNodesPerBucket;
                                             ++iSideNode)
                                        {
                                            const scalar* coords =
                                                stk::mesh::field_data(
                                                    *coordsSTKFieldPtr,
                                                    sideNodeBucket,
                                                    iSideNode);

#if SPATIAL_DIM == 3
                                            x = coords[0];
                                            y = coords[1];
                                            z = coords[2];
#elif SPATIAL_DIM == 2
                                            x = coords[0];
                                            y = coords[1];
#endif

                                            for (label i = 0; i < SPATIAL_DIM;
                                                 i++)
                                            {
                                                snvalueb[SPATIAL_DIM *
                                                             iSideNode +
                                                         i] =
                                                    expression_list[i].value();
                                            }

                                            const scalar* Um =
                                                stk::mesh::field_data(
                                                    UmSTKFieldRef,
                                                    sideNodeBucket,
                                                    iSideNode);

                                            // add mesh velocity
                                            for (label i = 0; i < SPATIAL_DIM;
                                                 i++)
                                            {
                                                snvalueb[SPATIAL_DIM *
                                                             iSideNode +
                                                         i] += Um[i];
                                            }

                                            // override internal node values
                                            if (correctedBoundaryNodeValues_)
                                            {
                                                for (label i = 0;
                                                     i < SPATIAL_DIM;
                                                     i++)
                                                {
                                                    valueb[SPATIAL_DIM *
                                                               iSideNode +
                                                           i] =
                                                        snvalueb[SPATIAL_DIM *
                                                                     iSideNode +
                                                                 i];
                                                }
                                            }
                                        }
                                    }

                                    // Update side field on the current boundary
                                    this->sideFieldRef().interpolate(
                                        this->nodeSideFieldRef(),
                                        iZone,
                                        iBoundary,
                                        this->isShifted());
                                }
                                break;

                            case inputDataType::profileData:
                                {
                                    errorMsg("profile data not provided yet");
                                }
                                break;
                        }
                    }
                }
                break;

            default:
                errorMsg("Must not reach here");
                break;
        }
    }
    else if (zone->meshDeforming())
    {
        // FIXME: Inlet boundary may be subject to mesh deformation
        assert(boundary->frameType() == boundaryRelativeFrameType::absolute);

        nodeVectorField::updateBoundarySideFieldSpecifiedValue(iZone,
                                                               iBoundary);
    }
    else
    {
        nodeVectorField::updateBoundarySideFieldSpecifiedValue(iZone,
                                                               iBoundary);
    }
}

void velocity::updateBoundarySideFieldNoSlipWall(label iZone, label iBoundary)
{
    const auto& mesh = this->meshRef();
    const stk::mesh::BulkData& bulkData = mesh.bulkDataRef();
    const stk::mesh::MetaData& metaData = mesh.metaDataRef();

    const auto* zone = mesh.zonePtr(iZone);
    const auto* wall =
        dynamic_cast<const class wall*>(zone->boundaryPtr(iBoundary));

    if (zone->frameRotating())
    {
        assert(wall->frameType() == boundaryRelativeFrameType::relative);

        if (wall->counterRotating())
        {
            // This wall is moving with the MFR, but has the opposite
            // velocity vector, that is, 0 absolute velocity
            nodeVectorField::updateBoundarySideFieldSpecifiedValue(iZone,
                                                                   iBoundary);
        }
        else
        {
            // This wall is moving with the MFR, and it has an absolute velocity
            // that is a summation of the rotational speed v = Omega x r and
            // the relative velocity provided in the input
            auto& bc = this->boundaryConditionRef(iZone, iBoundary);
            auto& data = bc.template data<SPATIAL_DIM>("value");

            // rotation data pointers
            const scalar* p_mat =
                zone->transformationRef().rotation().coriolisMatrix_.data();
            const scalar* p_ori =
                zone->transformationRef().rotation().origin_.data();

            switch (data.type())
            {
                case inputDataType::null:
                    break;

                case inputDataType::constant:
                case inputDataType::timeTable:
                    {
                        std::array<scalar, SPATIAL_DIM> inputValue;
                        if (data.type() == inputDataType::constant)
                        {
                            std::copy(data.value(),
                                      data.value() + SPATIAL_DIM,
                                      inputValue.begin());
                        }
                        else
                        {
                            inputValue = data.interpolate(
                                this->meshRef().controlsRef().time);
                        }

                        // select all nodes relevant to the node side field
                        stk::mesh::Selector selAllNodes =
                            metaData.universal_part() &
                            stk::mesh::selectUnion(wall->parts());

                        // get fields
                        auto& stkFieldRef = this->stkFieldRef();
                        auto& nodeSideSTKFieldRef =
                            this->nodeSideFieldRef().stkFieldRef();

                        // extract geometric fields
                        const STKScalarField* coordsSTKFieldPtr =
                            metaData.template get_field<scalar>(
                                stk::topology::NODE_RANK,
                                this->getCoordinatesID_(iZone));

                        stk::mesh::BucketVector const& sideNodeBuckets =
                            bulkData.get_buckets(stk::topology::NODE_RANK,
                                                 selAllNodes);
                        for (stk::mesh::BucketVector::const_iterator ib =
                                 sideNodeBuckets.begin();
                             ib != sideNodeBuckets.end();
                             ++ib)
                        {
                            stk::mesh::Bucket& sideNodeBucket = **ib;
                            const stk::mesh::Bucket::size_type
                                nSideNodesPerBucket = sideNodeBucket.size();
                            scalar* snvalueb = stk::mesh::field_data(
                                nodeSideSTKFieldRef, sideNodeBucket);
                            scalar* valueb = stk::mesh::field_data(
                                stkFieldRef, sideNodeBucket);
                            for (stk::mesh::Bucket::size_type iSideNode = 0;
                                 iSideNode < nSideNodesPerBucket;
                                 ++iSideNode)
                            {
                                const scalar* coords =
                                    stk::mesh::field_data(*coordsSTKFieldPtr,
                                                          sideNodeBucket,
                                                          iSideNode);

                                for (label i = 0; i < SPATIAL_DIM; i++)
                                {
                                    snvalueb[SPATIAL_DIM * iSideNode + i] =
                                        inputValue[i];
                                }

                                // add frame velocity
                                for (label i = 0; i < SPATIAL_DIM; i++)
                                {
                                    for (label j = 0; j < SPATIAL_DIM; j++)
                                    {
                                        snvalueb[SPATIAL_DIM * iSideNode + i] +=
                                            p_mat[i * SPATIAL_DIM + j] *
                                            (coords[j] - p_ori[j]);
                                    }
                                }

                                // override internal node values
                                if (correctedBoundaryNodeValues_)
                                {
                                    for (label i = 0; i < SPATIAL_DIM; i++)
                                    {
                                        valueb[SPATIAL_DIM * iSideNode + i] =
                                            snvalueb[SPATIAL_DIM * iSideNode +
                                                     i];
                                    }
                                }
                            }
                        }

                        // Update side field on the current boundary
                        this->sideFieldRef().interpolate(
                            this->nodeSideFieldRef(),
                            iZone,
                            iBoundary,
                            this->isShifted());
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

                        for (label i = 0; i < SPATIAL_DIM; i++)
                        {
                            if (parser.compile(data.expression()[i],
                                               componentExpression))
                            {
                                expression_list.push_back(componentExpression);
                            }
                            else
                            {
                                errorMsg("Error in the expression provided for "
                                         "field " +
                                         this->name() + ": " +
                                         data.expression()[i]);
                            }
                        }

                        // set time value
                        t = this->meshRef().controlsRef().time;

                        // select all nodes relevant to the node side field
                        stk::mesh::Selector selAllNodes =
                            metaData.universal_part() &
                            stk::mesh::selectUnion(wall->parts());

                        // get fields
                        auto& stkFieldRef = this->stkFieldRef();
                        auto& nodeSideSTKFieldRef =
                            this->nodeSideFieldRef().stkFieldRef();

                        // extract geometric fields
                        const STKScalarField* coordsSTKFieldPtr =
                            metaData.template get_field<scalar>(
                                stk::topology::NODE_RANK,
                                this->getCoordinatesID_(iZone));

                        stk::mesh::BucketVector const& sideNodeBuckets =
                            bulkData.get_buckets(stk::topology::NODE_RANK,
                                                 selAllNodes);
                        for (stk::mesh::BucketVector::const_iterator ib =
                                 sideNodeBuckets.begin();
                             ib != sideNodeBuckets.end();
                             ++ib)
                        {
                            stk::mesh::Bucket& sideNodeBucket = **ib;
                            const stk::mesh::Bucket::size_type
                                nSideNodesPerBucket = sideNodeBucket.size();
                            scalar* snvalueb = stk::mesh::field_data(
                                nodeSideSTKFieldRef, sideNodeBucket);
                            scalar* valueb = stk::mesh::field_data(
                                stkFieldRef, sideNodeBucket);
                            for (stk::mesh::Bucket::size_type iSideNode = 0;
                                 iSideNode < nSideNodesPerBucket;
                                 ++iSideNode)
                            {
                                const scalar* coords =
                                    stk::mesh::field_data(*coordsSTKFieldPtr,
                                                          sideNodeBucket,
                                                          iSideNode);

#if SPATIAL_DIM == 3
                                x = coords[0];
                                y = coords[1];
                                z = coords[2];
#elif SPATIAL_DIM == 2
                                x = coords[0];
                                y = coords[1];
#endif

                                for (label i = 0; i < SPATIAL_DIM; i++)
                                {
                                    snvalueb[SPATIAL_DIM * iSideNode + i] =
                                        expression_list[i].value();
                                }

                                // add frame velocity
                                for (label i = 0; i < SPATIAL_DIM; i++)
                                {
                                    for (label j = 0; j < SPATIAL_DIM; j++)
                                    {
                                        snvalueb[SPATIAL_DIM * iSideNode + i] +=
                                            p_mat[i * SPATIAL_DIM + j] *
                                            (coords[j] - p_ori[j]);
                                    }
                                }

                                // override internal node values
                                if (correctedBoundaryNodeValues_)
                                {
                                    for (label i = 0; i < SPATIAL_DIM; i++)
                                    {
                                        valueb[SPATIAL_DIM * iSideNode + i] =
                                            snvalueb[SPATIAL_DIM * iSideNode +
                                                     i];
                                    }
                                }
                            }
                        }

                        // Update side field on the current boundary
                        this->sideFieldRef().interpolate(
                            this->nodeSideFieldRef(),
                            iZone,
                            iBoundary,
                            this->isShifted());
                    }
                    break;

                case inputDataType::profileData:
                    {
                        errorMsg("profile data not provided yet");
                    }
                    break;
            }
        }
    }
    else if (zone->meshTransforming())
    {
        switch (zone->transformationRef().type())
        {
            case meshMotionType::stationary:
                errorMsg("Must not reach here");
                break;

            case meshMotionType::translating:
                {
                    assert(wall->frameType() ==
                           boundaryRelativeFrameType::relative);

                    auto& bc = this->boundaryConditionRef(iZone, iBoundary);
                    auto& data = bc.template data<SPATIAL_DIM>("value");

                    // Get mesh velocity field
                    const auto& UmSTKFieldRef = this->simulationRef()
                                                    .meshMotionRef()
                                                    .UmRef()
                                                    .stkFieldRef();

                    switch (data.type())
                    {
                        case inputDataType::null:
                            break;

                        case inputDataType::constant:
                        case inputDataType::timeTable:
                            {
                                std::array<scalar, SPATIAL_DIM> inputValue;
                                if (data.type() == inputDataType::constant)
                                {
                                    std::copy(data.value(),
                                              data.value() + SPATIAL_DIM,
                                              inputValue.begin());
                                }
                                else
                                {
                                    inputValue = data.interpolate(
                                        this->meshRef().controlsRef().time);
                                }

                                // select all nodes relevant to the node
                                // side field
                                stk::mesh::Selector selAllNodes =
                                    metaData.universal_part() &
                                    stk::mesh::selectUnion(wall->parts());

                                // get fields
                                auto& stkFieldRef = this->stkFieldRef();
                                auto& nodeSideSTKFieldRef =
                                    this->nodeSideFieldRef().stkFieldRef();

                                stk::mesh::BucketVector const& sideNodeBuckets =
                                    bulkData.get_buckets(
                                        stk::topology::NODE_RANK, selAllNodes);
                                for (stk::mesh::BucketVector::const_iterator
                                         ib = sideNodeBuckets.begin();
                                     ib != sideNodeBuckets.end();
                                     ++ib)
                                {
                                    stk::mesh::Bucket& sideNodeBucket = **ib;
                                    const stk::mesh::Bucket::size_type
                                        nSideNodesPerBucket =
                                            sideNodeBucket.size();
                                    scalar* snvalue = stk::mesh::field_data(
                                        nodeSideSTKFieldRef, sideNodeBucket);
                                    scalar* value = stk::mesh::field_data(
                                        stkFieldRef, sideNodeBucket);
                                    for (stk::mesh::Bucket::size_type
                                             iSideNode = 0;
                                         iSideNode < nSideNodesPerBucket;
                                         ++iSideNode)
                                    {
                                        const scalar* Um =
                                            stk::mesh::field_data(
                                                UmSTKFieldRef,
                                                sideNodeBucket,
                                                iSideNode);

                                        for (label i = 0; i < SPATIAL_DIM; i++)
                                        {
                                            snvalue[SPATIAL_DIM * iSideNode +
                                                    i] = inputValue[i];
                                        }

                                        // add mesh velocity
                                        for (label i = 0; i < SPATIAL_DIM; i++)
                                        {
                                            snvalue[SPATIAL_DIM * iSideNode +
                                                    i] += Um[i];
                                        }

                                        // override internal node values
                                        if (correctedBoundaryNodeValues_)
                                        {
                                            for (label i = 0; i < SPATIAL_DIM;
                                                 i++)
                                            {
                                                value[SPATIAL_DIM * iSideNode +
                                                      i] =
                                                    snvalue[SPATIAL_DIM *
                                                                iSideNode +
                                                            i];
                                            }
                                        }
                                    }
                                }

                                // Update side field on the current boundary
                                this->sideFieldRef().interpolate(
                                    this->nodeSideFieldRef(),
                                    iZone,
                                    iBoundary,
                                    this->isShifted());
                            }
                            break;

                        case inputDataType::expression:
                            {
                                typedef exprtk::symbol_table<scalar>
                                    symbol_table_t;
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
                                componentExpression.register_symbol_table(
                                    symbol_table);

                                parser_t parser;

                                for (label i = 0; i < SPATIAL_DIM; i++)
                                {
                                    if (parser.compile(data.expression()[i],
                                                       componentExpression))
                                    {
                                        expression_list.push_back(
                                            componentExpression);
                                    }
                                    else
                                    {
                                        errorMsg("Error in the expression "
                                                 "provided "
                                                 "for field " +
                                                 this->name() + ": " +
                                                 data.expression()[i]);
                                    }
                                }

                                // set time value
                                t = this->meshRef().controlsRef().time;

                                // select all nodes relevant to the node
                                // side field
                                stk::mesh::Selector selAllNodes =
                                    metaData.universal_part() &
                                    stk::mesh::selectUnion(wall->parts());

                                // get fields
                                auto& stkFieldRef = this->stkFieldRef();
                                auto& nodeSideSTKFieldRef =
                                    this->nodeSideFieldRef().stkFieldRef();

                                // extract geometric fields
                                const STKScalarField* coordsSTKFieldPtr =
                                    metaData.template get_field<scalar>(
                                        stk::topology::NODE_RANK,
                                        this->getCoordinatesID_(iZone));

                                stk::mesh::BucketVector const& sideNodeBuckets =
                                    bulkData.get_buckets(
                                        stk::topology::NODE_RANK, selAllNodes);
                                for (stk::mesh::BucketVector::const_iterator
                                         ib = sideNodeBuckets.begin();
                                     ib != sideNodeBuckets.end();
                                     ++ib)
                                {
                                    stk::mesh::Bucket& sideNodeBucket = **ib;
                                    const stk::mesh::Bucket::size_type
                                        nSideNodesPerBucket =
                                            sideNodeBucket.size();
                                    scalar* snvalueb = stk::mesh::field_data(
                                        nodeSideSTKFieldRef, sideNodeBucket);
                                    scalar* valueb = stk::mesh::field_data(
                                        stkFieldRef, sideNodeBucket);
                                    for (stk::mesh::Bucket::size_type
                                             iSideNode = 0;
                                         iSideNode < nSideNodesPerBucket;
                                         ++iSideNode)
                                    {
                                        const scalar* coords =
                                            stk::mesh::field_data(
                                                *coordsSTKFieldPtr,
                                                sideNodeBucket,
                                                iSideNode);

#if SPATIAL_DIM == 3
                                        x = coords[0];
                                        y = coords[1];
                                        z = coords[2];
#elif SPATIAL_DIM == 2
                                        x = coords[0];
                                        y = coords[1];
#endif

                                        for (label i = 0; i < SPATIAL_DIM; i++)
                                        {
                                            snvalueb[SPATIAL_DIM * iSideNode +
                                                     i] =
                                                expression_list[i].value();
                                        }

                                        const scalar* Um =
                                            stk::mesh::field_data(
                                                UmSTKFieldRef,
                                                sideNodeBucket,
                                                iSideNode);

                                        // add mesh velocity
                                        for (label i = 0; i < SPATIAL_DIM; i++)
                                        {
                                            snvalueb[SPATIAL_DIM * iSideNode +
                                                     i] += Um[i];
                                        }

                                        // override internal node values
                                        if (correctedBoundaryNodeValues_)
                                        {
                                            for (label i = 0; i < SPATIAL_DIM;
                                                 i++)
                                            {
                                                valueb[SPATIAL_DIM * iSideNode +
                                                       i] =
                                                    snvalueb[SPATIAL_DIM *
                                                                 iSideNode +
                                                             i];
                                            }
                                        }
                                    }
                                }

                                // Update side field on the current boundary
                                this->sideFieldRef().interpolate(
                                    this->nodeSideFieldRef(),
                                    iZone,
                                    iBoundary,
                                    this->isShifted());
                            }
                            break;

                        case inputDataType::profileData:
                            {
                                errorMsg("profile data not provided yet");
                            }
                            break;
                    }
                }
                break;

            case meshMotionType::rotating:
                {
                    assert(wall->frameType() ==
                           boundaryRelativeFrameType::relative);

                    if (wall->counterRotating())
                    {
                        // This wall is moving with the zone, but has the
                        // opposite velocity vector, that is, 0 absolute
                        // velocity
                        nodeVectorField::updateBoundarySideFieldSpecifiedValue(
                            iZone, iBoundary);
                    }
                    else
                    {
                        auto& bc = this->boundaryConditionRef(iZone, iBoundary);
                        auto& data = bc.template data<SPATIAL_DIM>("value");

                        // Get mesh velocity field
                        const auto& UmSTKFieldRef = this->simulationRef()
                                                        .meshMotionRef()
                                                        .UmRef()
                                                        .stkFieldRef();

                        switch (data.type())
                        {
                            case inputDataType::null:
                                break;

                            case inputDataType::constant:
                            case inputDataType::timeTable:
                                {
                                    std::array<scalar, SPATIAL_DIM> inputValue;
                                    if (data.type() == inputDataType::constant)
                                    {
                                        std::copy(data.value(),
                                                  data.value() + SPATIAL_DIM,
                                                  inputValue.begin());
                                    }
                                    else
                                    {
                                        inputValue = data.interpolate(
                                            this->meshRef().controlsRef().time);
                                    }

                                    // select all nodes relevant to the node
                                    // side field
                                    stk::mesh::Selector selAllNodes =
                                        metaData.universal_part() &
                                        stk::mesh::selectUnion(wall->parts());

                                    // get fields
                                    auto& stkFieldRef = this->stkFieldRef();
                                    auto& nodeSideSTKFieldRef =
                                        this->nodeSideFieldRef().stkFieldRef();

                                    stk::mesh::BucketVector const&
                                        sideNodeBuckets = bulkData.get_buckets(
                                            stk::topology::NODE_RANK,
                                            selAllNodes);
                                    for (stk::mesh::BucketVector::const_iterator
                                             ib = sideNodeBuckets.begin();
                                         ib != sideNodeBuckets.end();
                                         ++ib)
                                    {
                                        stk::mesh::Bucket& sideNodeBucket =
                                            **ib;
                                        const stk::mesh::Bucket::size_type
                                            nSideNodesPerBucket =
                                                sideNodeBucket.size();
                                        scalar* snvalueb =
                                            stk::mesh::field_data(
                                                nodeSideSTKFieldRef,
                                                sideNodeBucket);
                                        scalar* valueb = stk::mesh::field_data(
                                            stkFieldRef, sideNodeBucket);
                                        for (stk::mesh::Bucket::size_type
                                                 iSideNode = 0;
                                             iSideNode < nSideNodesPerBucket;
                                             ++iSideNode)
                                        {
                                            const scalar* Um =
                                                stk::mesh::field_data(
                                                    UmSTKFieldRef,
                                                    sideNodeBucket,
                                                    iSideNode);

                                            for (label i = 0; i < SPATIAL_DIM;
                                                 i++)
                                            {
                                                snvalueb[SPATIAL_DIM *
                                                             iSideNode +
                                                         i] = inputValue[i];
                                            }

                                            // add mesh velocity
                                            for (label i = 0; i < SPATIAL_DIM;
                                                 i++)
                                            {
                                                snvalueb[SPATIAL_DIM *
                                                             iSideNode +
                                                         i] += Um[i];
                                            }

                                            // override internal node values
                                            if (correctedBoundaryNodeValues_)
                                            {
                                                for (label i = 0;
                                                     i < SPATIAL_DIM;
                                                     i++)
                                                {
                                                    valueb[SPATIAL_DIM *
                                                               iSideNode +
                                                           i] =
                                                        snvalueb[SPATIAL_DIM *
                                                                     iSideNode +
                                                                 i];
                                                }
                                            }
                                        }
                                    }

                                    // Update side field on the current boundary
                                    this->sideFieldRef().interpolate(
                                        this->nodeSideFieldRef(),
                                        iZone,
                                        iBoundary,
                                        this->isShifted());
                                }
                                break;

                            case inputDataType::expression:
                                {
                                    typedef exprtk::symbol_table<scalar>
                                        symbol_table_t;
                                    typedef exprtk::expression<scalar>
                                        expression_t;
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
                                    componentExpression.register_symbol_table(
                                        symbol_table);

                                    parser_t parser;

                                    for (label i = 0; i < SPATIAL_DIM; i++)
                                    {
                                        if (parser.compile(data.expression()[i],
                                                           componentExpression))
                                        {
                                            expression_list.push_back(
                                                componentExpression);
                                        }
                                        else
                                        {
                                            errorMsg("Error in the expression "
                                                     "provided "
                                                     "for field " +
                                                     this->name() + ": " +
                                                     data.expression()[i]);
                                        }
                                    }

                                    // set time value
                                    t = this->meshRef().controlsRef().time;

                                    // select all nodes relevant to the node
                                    // side field
                                    stk::mesh::Selector selAllNodes =
                                        metaData.universal_part() &
                                        stk::mesh::selectUnion(wall->parts());

                                    // get fields
                                    auto& stkFieldRef = this->stkFieldRef();
                                    auto& nodeSideSTKFieldRef =
                                        this->nodeSideFieldRef().stkFieldRef();

                                    // extract geometric fields
                                    const STKScalarField* coordsSTKFieldPtr =
                                        metaData.template get_field<scalar>(
                                            stk::topology::NODE_RANK,
                                            this->getCoordinatesID_(iZone));

                                    stk::mesh::BucketVector const&
                                        sideNodeBuckets = bulkData.get_buckets(
                                            stk::topology::NODE_RANK,
                                            selAllNodes);
                                    for (stk::mesh::BucketVector::const_iterator
                                             ib = sideNodeBuckets.begin();
                                         ib != sideNodeBuckets.end();
                                         ++ib)
                                    {
                                        stk::mesh::Bucket& sideNodeBucket =
                                            **ib;
                                        const stk::mesh::Bucket::size_type
                                            nSideNodesPerBucket =
                                                sideNodeBucket.size();
                                        scalar* snvalueb =
                                            stk::mesh::field_data(
                                                nodeSideSTKFieldRef,
                                                sideNodeBucket);
                                        scalar* valueb = stk::mesh::field_data(
                                            stkFieldRef, sideNodeBucket);
                                        for (stk::mesh::Bucket::size_type
                                                 iSideNode = 0;
                                             iSideNode < nSideNodesPerBucket;
                                             ++iSideNode)
                                        {
                                            const scalar* coords =
                                                stk::mesh::field_data(
                                                    *coordsSTKFieldPtr,
                                                    sideNodeBucket,
                                                    iSideNode);

#if SPATIAL_DIM == 3
                                            x = coords[0];
                                            y = coords[1];
                                            z = coords[2];
#elif SPATIAL_DIM == 2
                                            x = coords[0];
                                            y = coords[1];
#endif

                                            for (label i = 0; i < SPATIAL_DIM;
                                                 i++)
                                            {
                                                snvalueb[SPATIAL_DIM *
                                                             iSideNode +
                                                         i] =
                                                    expression_list[i].value();
                                            }

                                            const scalar* Um =
                                                stk::mesh::field_data(
                                                    UmSTKFieldRef,
                                                    sideNodeBucket,
                                                    iSideNode);

                                            // add mesh velocity
                                            for (label i = 0; i < SPATIAL_DIM;
                                                 i++)
                                            {
                                                snvalueb[SPATIAL_DIM *
                                                             iSideNode +
                                                         i] += Um[i];
                                            }

                                            // override internal node values
                                            if (correctedBoundaryNodeValues_)
                                            {
                                                for (label i = 0;
                                                     i < SPATIAL_DIM;
                                                     i++)
                                                {
                                                    valueb[SPATIAL_DIM *
                                                               iSideNode +
                                                           i] =
                                                        snvalueb[SPATIAL_DIM *
                                                                     iSideNode +
                                                                 i];
                                                }
                                            }
                                        }
                                    }

                                    // Update side field on the current boundary
                                    this->sideFieldRef().interpolate(
                                        this->nodeSideFieldRef(),
                                        iZone,
                                        iBoundary,
                                        this->isShifted());
                                }
                                break;

                            case inputDataType::profileData:
                                {
                                    errorMsg("profile data not provided yet");
                                }
                                break;
                        }
                    }
                }
                break;
        }
    }
    else if (zone->meshDeforming())
    {
        assert(wall->frameType() == boundaryRelativeFrameType::absolute);

        auto& bc = this->boundaryConditionRef(iZone, iBoundary);
        auto& data = bc.template data<SPATIAL_DIM>("value");

        bool wallVelRelativeToMeshMotion =
            wall->wallVelocityRelativeToMeshMotion();

        // Get mesh velocity field
        const auto& UmSTKFieldRef =
            this->simulationRef().meshMotionRef().UmRef().stkFieldRef();

        switch (data.type())
        {
            case inputDataType::null:
                break;

            case inputDataType::constant:
            case inputDataType::timeTable:
                {
                    std::array<scalar, SPATIAL_DIM> inputValue;
                    if (data.type() == inputDataType::constant)
                    {
                        std::copy(data.value(),
                                  data.value() + SPATIAL_DIM,
                                  inputValue.begin());
                    }
                    else
                    {
                        inputValue = data.interpolate(
                            this->meshRef().controlsRef().time);
                    }

                    // select all nodes relevant to the node side field
                    stk::mesh::Selector selAllNodes =
                        metaData.universal_part() &
                        stk::mesh::selectUnion(wall->parts());

                    // get fields
                    auto& stkFieldRef = this->stkFieldRef();
                    auto& nodeSideSTKFieldRef =
                        this->nodeSideFieldRef().stkFieldRef();

                    stk::mesh::BucketVector const& sideNodeBuckets =
                        bulkData.get_buckets(stk::topology::NODE_RANK,
                                             selAllNodes);
                    for (stk::mesh::BucketVector::const_iterator ib =
                             sideNodeBuckets.begin();
                         ib != sideNodeBuckets.end();
                         ++ib)
                    {
                        stk::mesh::Bucket& sideNodeBucket = **ib;
                        const stk::mesh::Bucket::size_type nSideNodesPerBucket =
                            sideNodeBucket.size();
                        scalar* snvalue = stk::mesh::field_data(
                            nodeSideSTKFieldRef, sideNodeBucket);
                        scalar* value =
                            stk::mesh::field_data(stkFieldRef, sideNodeBucket);
                        for (stk::mesh::Bucket::size_type iSideNode = 0;
                             iSideNode < nSideNodesPerBucket;
                             ++iSideNode)
                        {
                            const scalar* Um = stk::mesh::field_data(
                                UmSTKFieldRef, sideNodeBucket, iSideNode);

                            for (label i = 0; i < SPATIAL_DIM; i++)
                            {
                                snvalue[SPATIAL_DIM * iSideNode + i] =
                                    inputValue[i];
                            }

                            // add mesh velocity
                            if (wallVelRelativeToMeshMotion)
                            {
                                for (label i = 0; i < SPATIAL_DIM; i++)
                                {
                                    snvalue[SPATIAL_DIM * iSideNode + i] +=
                                        Um[i];
                                }
                            }

                            // override internal node values
                            if (correctedBoundaryNodeValues_)
                            {
                                for (label i = 0; i < SPATIAL_DIM; i++)
                                {
                                    value[SPATIAL_DIM * iSideNode + i] =
                                        snvalue[SPATIAL_DIM * iSideNode + i];
                                }
                            }
                        }
                    }

                    // Update side field on the current boundary
                    this->sideFieldRef().interpolate(this->nodeSideFieldRef(),
                                                     iZone,
                                                     iBoundary,
                                                     this->isShifted());
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

                    for (label i = 0; i < SPATIAL_DIM; i++)
                    {
                        if (parser.compile(data.expression()[i],
                                           componentExpression))
                        {
                            expression_list.push_back(componentExpression);
                        }
                        else
                        {
                            errorMsg("Error in the expression "
                                     "provided "
                                     "for field " +
                                     this->name() + ": " +
                                     data.expression()[i]);
                        }
                    }

                    // set time value
                    t = this->meshRef().controlsRef().time;

                    // select all nodes relevant to the node side field
                    stk::mesh::Selector selAllNodes =
                        metaData.universal_part() &
                        stk::mesh::selectUnion(wall->parts());

                    // get fields
                    auto& stkFieldRef = this->stkFieldRef();
                    auto& nodeSideSTKFieldRef =
                        this->nodeSideFieldRef().stkFieldRef();

                    // extract geometric fields
                    const STKScalarField* coordsSTKFieldPtr =
                        metaData.template get_field<scalar>(
                            stk::topology::NODE_RANK,
                            this->getCoordinatesID_(iZone));

                    stk::mesh::BucketVector const& sideNodeBuckets =
                        bulkData.get_buckets(stk::topology::NODE_RANK,
                                             selAllNodes);
                    for (stk::mesh::BucketVector::const_iterator ib =
                             sideNodeBuckets.begin();
                         ib != sideNodeBuckets.end();
                         ++ib)
                    {
                        stk::mesh::Bucket& sideNodeBucket = **ib;
                        const stk::mesh::Bucket::size_type nSideNodesPerBucket =
                            sideNodeBucket.size();
                        scalar* snvalueb = stk::mesh::field_data(
                            nodeSideSTKFieldRef, sideNodeBucket);
                        scalar* valueb =
                            stk::mesh::field_data(stkFieldRef, sideNodeBucket);
                        for (stk::mesh::Bucket::size_type iSideNode = 0;
                             iSideNode < nSideNodesPerBucket;
                             ++iSideNode)
                        {
                            const scalar* coords = stk::mesh::field_data(
                                *coordsSTKFieldPtr, sideNodeBucket, iSideNode);

#if SPATIAL_DIM == 3
                            x = coords[0];
                            y = coords[1];
                            z = coords[2];
#elif SPATIAL_DIM == 2
                            x = coords[0];
                            y = coords[1];
#endif

                            for (label i = 0; i < SPATIAL_DIM; i++)
                            {
                                snvalueb[SPATIAL_DIM * iSideNode + i] =
                                    expression_list[i].value();
                            }

                            const scalar* Um = stk::mesh::field_data(
                                UmSTKFieldRef, sideNodeBucket, iSideNode);

                            // add mesh velocity
                            if (wallVelRelativeToMeshMotion)
                            {
                                for (label i = 0; i < SPATIAL_DIM; i++)
                                {
                                    snvalueb[SPATIAL_DIM * iSideNode + i] +=
                                        Um[i];
                                }
                            }

                            // override internal node values
                            if (correctedBoundaryNodeValues_)
                            {
                                for (label i = 0; i < SPATIAL_DIM; i++)
                                {
                                    valueb[SPATIAL_DIM * iSideNode + i] =
                                        snvalueb[SPATIAL_DIM * iSideNode + i];
                                }
                            }
                        }
                    }

                    // Update side field on the current boundary
                    this->sideFieldRef().interpolate(this->nodeSideFieldRef(),
                                                     iZone,
                                                     iBoundary,
                                                     this->isShifted());
                }
                break;

            case inputDataType::profileData:
                {
                    errorMsg("profile data not provided yet");
                }
                break;
        }
    }
    else
    {
        nodeVectorField::updateBoundarySideFieldSpecifiedValue(iZone,
                                                               iBoundary);
    }
}

void velocity::updateBoundarySideFieldNormalSpeed(label iZone, label iBoundary)
{
    const auto& mesh = this->meshRef();
    const stk::mesh::BulkData& bulkData = mesh.bulkDataRef();
    const stk::mesh::MetaData& metaData = mesh.metaDataRef();

    const auto* zone = mesh.zonePtr(iZone);
    const auto* boundary = zone->boundaryPtr(iBoundary);

    auto& bc = this->boundaryConditionRef(iZone, iBoundary);
    auto& data = bc.template data<1>("value");

    if (zone->frameRotating())
    {
        switch (zone->transformationRef().type())
        {
            case meshMotionType::stationary:
            case meshMotionType::translating:
                errorMsg("Must not reach here");
                break;

            case meshMotionType::rotating:
                {
                    if (boundary->frameType() ==
                        boundaryRelativeFrameType::absolute)
                    {
                        // Reset node side field
                        this->nodeSideFieldRef().setToValue(
                            std::vector<scalar>(SPATIAL_DIM, 0.0),
                            boundary->parts());

                        switch (data.type())
                        {
                            case inputDataType::null:
                                break;

                            case inputDataType::constant:
                            case inputDataType::timeTable:
                                {
                                    scalar normalSpeed =
                                        data.type() == inputDataType::constant
                                            ? *data.value()
                                            : data.interpolate(
                                                  this->meshRef()
                                                      .controlsRef()
                                                      .time)[0];

                                    // get fields
                                    auto& stkFieldRef = this->stkFieldRef();
                                    auto& nodeSideSTKFieldRef =
                                        this->nodeSideFieldRef().stkFieldRef();

                                    // Find exposed node area vector (sum of
                                    // surrounding ip's)
                                    {
                                        // define some common selectors
                                        stk::mesh::Selector selAllSides =
                                            metaData.universal_part() &
                                            stk::mesh::selectUnion(
                                                boundary->parts());

                                        // Get fields
                                        STKScalarField&
                                            exposedAreaVecSTKFieldRef =
                                                *metaData.get_field<scalar>(
                                                    metaData.side_rank(),
                                                    this->getExposedAreaVectorID_(
                                                        iZone));

                                        // define vector of parent topos; should
                                        // always be UNITY in size
                                        std::vector<stk::topology> parentTopo;

                                        // setup for buckets; union parts and
                                        // ask for universal part
                                        stk::mesh::BucketVector const&
                                            sideBuckets = bulkData.get_buckets(
                                                metaData.side_rank(),
                                                selAllSides);

                                        for (stk::mesh::BucketVector::
                                                 const_iterator ib =
                                                     sideBuckets.begin();
                                             ib != sideBuckets.end();
                                             ++ib)
                                        {
                                            stk::mesh::Bucket& sideBucket =
                                                **ib;

                                            // extract master element
                                            MasterElement* meFC =
                                                MasterElementRepo::
                                                    get_surface_master_element(
                                                        sideBucket.topology());

                                            // mapping from ip to nodes for this
                                            // ordinal
                                            const label* ipNodeMap =
                                                meFC->ipNodeMap();

                                            // extract master element specifics
                                            const label numScsIp =
                                                meFC->numIntPoints_;

                                            const stk::mesh::Bucket::size_type
                                                nSidesPerBucket =
                                                    sideBucket.size();
                                            for (stk::mesh::Bucket::size_type
                                                     iSide = 0;
                                                 iSide < nSidesPerBucket;
                                                 ++iSide)
                                            {
                                                // get face
                                                stk::mesh::Entity side =
                                                    sideBucket[iSide];

                                                stk::mesh::Entity const*
                                                    sideNodeRels =
                                                        bulkData.begin_nodes(
                                                            side);

                                                // face data
                                                scalar* areaVec =
                                                    stk::mesh::field_data(
                                                        exposedAreaVecSTKFieldRef,
                                                        sideBucket,
                                                        iSide);

                                                // accumulate ip areas around
                                                // node
                                                for (label ip = 0;
                                                     ip < numScsIp;
                                                     ++ip)
                                                {
                                                    const label nearestNode =
                                                        ipNodeMap[ip];
                                                    stk::mesh::Entity node =
                                                        sideNodeRels
                                                            [nearestNode];

                                                    scalar* snvalue =
                                                        stk::mesh::field_data(
                                                            nodeSideSTKFieldRef,
                                                            node);

                                                    for (label i = 0;
                                                         i < SPATIAL_DIM;
                                                         i++)
                                                    {
                                                        snvalue[i] -= areaVec
                                                            [ip * SPATIAL_DIM +
                                                             i];
                                                    }
                                                }
                                            }
                                        }
                                    }

                                    // multiply by speed and normalize
                                    {
                                        // select all nodes relevant to the node
                                        // side field
                                        stk::mesh::Selector selAllNodes =
                                            metaData.universal_part() &
                                            stk::mesh::selectUnion(
                                                boundary->parts());

                                        stk::mesh::BucketVector const&
                                            sideNodeBuckets =
                                                bulkData.get_buckets(
                                                    stk::topology::NODE_RANK,
                                                    selAllNodes);
                                        for (stk::mesh::BucketVector::
                                                 const_iterator ib =
                                                     sideNodeBuckets.begin();
                                             ib != sideNodeBuckets.end();
                                             ++ib)
                                        {
                                            stk::mesh::Bucket& sideNodeBucket =
                                                **ib;
                                            const stk::mesh::Bucket::size_type
                                                nSideNodesPerBucket =
                                                    sideNodeBucket.size();
                                            scalar* snvalue =
                                                stk::mesh::field_data(
                                                    nodeSideSTKFieldRef,
                                                    sideNodeBucket);
                                            scalar* value =
                                                stk::mesh::field_data(
                                                    stkFieldRef,
                                                    sideNodeBucket);

                                            for (stk::mesh::Bucket::size_type
                                                     iSideNode = 0;
                                                 iSideNode <
                                                 nSideNodesPerBucket;
                                                 ++iSideNode)
                                            {
                                                scalar asq = 0.0;
                                                for (label i = 0;
                                                     i < SPATIAL_DIM;
                                                     ++i)
                                                {
                                                    const scalar axi =
                                                        snvalue[SPATIAL_DIM *
                                                                    iSideNode +
                                                                i];
                                                    asq += axi * axi;
                                                }
                                                const scalar amag =
                                                    std::sqrt(asq);

                                                for (label i = 0;
                                                     i < SPATIAL_DIM;
                                                     i++)
                                                {
                                                    snvalue[SPATIAL_DIM *
                                                                iSideNode +
                                                            i] *=
                                                        normalSpeed / amag;
                                                }

                                                // override internal node values
                                                if (correctedBoundaryNodeValues_)
                                                {
                                                    for (label i = 0;
                                                         i < SPATIAL_DIM;
                                                         i++)
                                                    {
                                                        value[SPATIAL_DIM *
                                                                  iSideNode +
                                                              i] = snvalue
                                                            [SPATIAL_DIM *
                                                                 iSideNode +
                                                             i];
                                                    }
                                                }
                                            }
                                        }
                                    }

                                    // Update side field on the current boundary
                                    this->sideFieldRef().interpolate(
                                        this->nodeSideFieldRef(),
                                        iZone,
                                        iBoundary,
                                        this->isShifted());
                                }
                                break;

                            case inputDataType::expression:
                                {
                                    typedef exprtk::symbol_table<scalar>
                                        symbol_table_t;
                                    typedef exprtk::expression<scalar>
                                        expression_t;
                                    typedef exprtk::parser<scalar> parser_t;

                                    symbol_table_t symbol_table;
                                    symbol_table.add_constants();

                                    // declare function vars
                                    scalar t, x, y, z;
                                    symbol_table.add_variable("t", t);
                                    symbol_table.add_variable("x", x);
                                    symbol_table.add_variable("y", y);
                                    symbol_table.add_variable("z", z);

                                    expression_t expression;
                                    expression.register_symbol_table(
                                        symbol_table);

                                    parser_t parser;

                                    if (!parser.compile(data.expression()[0],
                                                        expression))
                                    {
                                        errorMsg("Error in the expression "
                                                 "provided for " +
                                                 this->name());
                                    }

                                    // set time value
                                    t = this->meshRef().controlsRef().time;

                                    // get fields
                                    auto& stkFieldRef = this->stkFieldRef();
                                    auto& nodeSideSTKFieldRef =
                                        this->nodeSideFieldRef().stkFieldRef();

                                    // extract geometric fields
                                    const STKScalarField* coordsSTKFieldPtr =
                                        metaData.template get_field<scalar>(
                                            stk::topology::NODE_RANK,
                                            this->getCoordinatesID_(iZone));

                                    // Find exposed node area vector (sum of
                                    // surrounding ip's)
                                    {
                                        // define some common selectors
                                        stk::mesh::Selector selAllSides =
                                            metaData.universal_part() &
                                            stk::mesh::selectUnion(
                                                boundary->parts());

                                        // Get fields
                                        STKScalarField&
                                            exposedAreaVecSTKFieldRef =
                                                *metaData.get_field<scalar>(
                                                    metaData.side_rank(),
                                                    this->getExposedAreaVectorID_(
                                                        iZone));

                                        // define vector of parent topos; should
                                        // always be UNITY in size
                                        std::vector<stk::topology> parentTopo;

                                        // setup for buckets; union parts and
                                        // ask for universal part
                                        stk::mesh::BucketVector const&
                                            sideBuckets = bulkData.get_buckets(
                                                metaData.side_rank(),
                                                selAllSides);

                                        for (stk::mesh::BucketVector::
                                                 const_iterator ib =
                                                     sideBuckets.begin();
                                             ib != sideBuckets.end();
                                             ++ib)
                                        {
                                            stk::mesh::Bucket& sideBucket =
                                                **ib;

                                            // extract master element
                                            MasterElement* meFC =
                                                MasterElementRepo::
                                                    get_surface_master_element(
                                                        sideBucket.topology());

                                            // mapping from ip to nodes for this
                                            // ordinal
                                            const label* ipNodeMap =
                                                meFC->ipNodeMap();

                                            // extract master element specifics
                                            const label numScsIp =
                                                meFC->numIntPoints_;

                                            const stk::mesh::Bucket::size_type
                                                nSidesPerBucket =
                                                    sideBucket.size();
                                            for (stk::mesh::Bucket::size_type
                                                     iSide = 0;
                                                 iSide < nSidesPerBucket;
                                                 ++iSide)
                                            {
                                                // get face
                                                stk::mesh::Entity side =
                                                    sideBucket[iSide];

                                                stk::mesh::Entity const*
                                                    sideNodeRels =
                                                        bulkData.begin_nodes(
                                                            side);

                                                // face data
                                                scalar* areaVec =
                                                    stk::mesh::field_data(
                                                        exposedAreaVecSTKFieldRef,
                                                        sideBucket,
                                                        iSide);

                                                // accumulate ip areas around
                                                // node
                                                for (label ip = 0;
                                                     ip < numScsIp;
                                                     ++ip)
                                                {
                                                    const label nearestNode =
                                                        ipNodeMap[ip];
                                                    stk::mesh::Entity node =
                                                        sideNodeRels
                                                            [nearestNode];

                                                    scalar* snvalue =
                                                        stk::mesh::field_data(
                                                            nodeSideSTKFieldRef,
                                                            node);

                                                    for (label i = 0;
                                                         i < SPATIAL_DIM;
                                                         i++)
                                                    {
                                                        snvalue[i] -= areaVec
                                                            [ip * SPATIAL_DIM +
                                                             i];
                                                    }
                                                }
                                            }
                                        }
                                    }

                                    // multiply by speed and normalize
                                    {
                                        // select all nodes relevant to the node
                                        // side field
                                        stk::mesh::Selector selAllNodes =
                                            metaData.universal_part() &
                                            stk::mesh::selectUnion(
                                                boundary->parts());

                                        stk::mesh::BucketVector const&
                                            sideNodeBuckets =
                                                bulkData.get_buckets(
                                                    stk::topology::NODE_RANK,
                                                    selAllNodes);
                                        for (stk::mesh::BucketVector::
                                                 const_iterator ib =
                                                     sideNodeBuckets.begin();
                                             ib != sideNodeBuckets.end();
                                             ++ib)
                                        {
                                            stk::mesh::Bucket& sideNodeBucket =
                                                **ib;
                                            const stk::mesh::Bucket::size_type
                                                nSideNodesPerBucket =
                                                    sideNodeBucket.size();
                                            scalar* snvalueb =
                                                stk::mesh::field_data(
                                                    nodeSideSTKFieldRef,
                                                    sideNodeBucket);
                                            scalar* valueb =
                                                stk::mesh::field_data(
                                                    stkFieldRef,
                                                    sideNodeBucket);

                                            for (stk::mesh::Bucket::size_type
                                                     iSideNode = 0;
                                                 iSideNode <
                                                 nSideNodesPerBucket;
                                                 ++iSideNode)
                                            {
                                                const scalar* coords =
                                                    stk::mesh::field_data(
                                                        *coordsSTKFieldPtr,
                                                        sideNodeBucket,
                                                        iSideNode);

#if SPATIAL_DIM == 3
                                                x = coords[0];
                                                y = coords[1];
                                                z = coords[2];
#elif SPATIAL_DIM == 2
                                                x = coords[0];
                                                y = coords[1];
#endif

                                                scalar normalSpeed =
                                                    expression.value();

                                                scalar asq = 0.0;
                                                for (label i = 0;
                                                     i < SPATIAL_DIM;
                                                     ++i)
                                                {
                                                    const scalar axi =
                                                        snvalueb[SPATIAL_DIM *
                                                                     iSideNode +
                                                                 i];
                                                    asq += axi * axi;
                                                }
                                                const scalar amag =
                                                    std::sqrt(asq);

                                                for (label i = 0;
                                                     i < SPATIAL_DIM;
                                                     i++)
                                                {
                                                    snvalueb[SPATIAL_DIM *
                                                                 iSideNode +
                                                             i] *=
                                                        normalSpeed / amag;
                                                }

                                                // override internal node values
                                                if (correctedBoundaryNodeValues_)
                                                {
                                                    for (label i = 0;
                                                         i < SPATIAL_DIM;
                                                         i++)
                                                    {
                                                        valueb[SPATIAL_DIM *
                                                                   iSideNode +
                                                               i] = snvalueb
                                                            [SPATIAL_DIM *
                                                                 iSideNode +
                                                             i];
                                                    }
                                                }
                                            }
                                        }
                                    }

                                    // Update side field on the current boundary
                                    this->sideFieldRef().interpolate(
                                        this->nodeSideFieldRef(),
                                        iZone,
                                        iBoundary,
                                        this->isShifted());
                                }
                                break;

                            case inputDataType::profileData:
                                {
                                    errorMsg("profile data not provided yet");
                                }
                                break;
                        }
                    }
                    else
                    {
                        // Reset node side field
                        this->nodeSideFieldRef().setToValue(
                            std::vector<scalar>(SPATIAL_DIM, 0.0),
                            boundary->parts());

                        // rotation data pointers
                        const scalar* p_mat = zone->transformationRef()
                                                  .rotation()
                                                  .coriolisMatrix_.data();
                        const scalar* p_ori =
                            zone->transformationRef().rotation().origin_.data();

                        switch (data.type())
                        {
                            case inputDataType::null:
                                break;

                            case inputDataType::constant:
                            case inputDataType::timeTable:
                                {
                                    scalar normalSpeed =
                                        data.type() == inputDataType::constant
                                            ? *data.value()
                                            : data.interpolate(
                                                  this->meshRef()
                                                      .controlsRef()
                                                      .time)[0];

                                    // get fields
                                    auto& stkFieldRef = this->stkFieldRef();
                                    auto& nodeSideSTKFieldRef =
                                        this->nodeSideFieldRef().stkFieldRef();

                                    // extract geometric fields
                                    const STKScalarField* coordsSTKFieldPtr =
                                        metaData.template get_field<scalar>(
                                            stk::topology::NODE_RANK,
                                            this->getCoordinatesID_(iZone));

                                    // Find exposed node area vector (sum of
                                    // surrounding ip's)
                                    {
                                        // define some common selectors
                                        stk::mesh::Selector selAllSides =
                                            metaData.universal_part() &
                                            stk::mesh::selectUnion(
                                                boundary->parts());

                                        // Get fields
                                        STKScalarField&
                                            exposedAreaVecSTKFieldRef =
                                                *metaData.get_field<scalar>(
                                                    metaData.side_rank(),
                                                    this->getExposedAreaVectorID_(
                                                        iZone));

                                        // define vector of parent topos; should
                                        // always be UNITY in size
                                        std::vector<stk::topology> parentTopo;

                                        // setup for buckets; union parts and
                                        // ask for universal part
                                        stk::mesh::BucketVector const&
                                            sideBuckets = bulkData.get_buckets(
                                                metaData.side_rank(),
                                                selAllSides);

                                        for (stk::mesh::BucketVector::
                                                 const_iterator ib =
                                                     sideBuckets.begin();
                                             ib != sideBuckets.end();
                                             ++ib)
                                        {
                                            stk::mesh::Bucket& sideBucket =
                                                **ib;

                                            // extract master element
                                            MasterElement* meFC =
                                                MasterElementRepo::
                                                    get_surface_master_element(
                                                        sideBucket.topology());

                                            // mapping from ip to nodes for this
                                            // ordinal
                                            const label* ipNodeMap =
                                                meFC->ipNodeMap();

                                            // extract master element specifics
                                            const label numScsIp =
                                                meFC->numIntPoints_;

                                            const stk::mesh::Bucket::size_type
                                                nSidesPerBucket =
                                                    sideBucket.size();
                                            for (stk::mesh::Bucket::size_type
                                                     iSide = 0;
                                                 iSide < nSidesPerBucket;
                                                 ++iSide)
                                            {
                                                // get face
                                                stk::mesh::Entity side =
                                                    sideBucket[iSide];

                                                stk::mesh::Entity const*
                                                    sideNodeRels =
                                                        bulkData.begin_nodes(
                                                            side);

                                                // face data
                                                scalar* areaVec =
                                                    stk::mesh::field_data(
                                                        exposedAreaVecSTKFieldRef,
                                                        sideBucket,
                                                        iSide);

                                                // accumulate ip areas around
                                                // node
                                                for (label ip = 0;
                                                     ip < numScsIp;
                                                     ++ip)
                                                {
                                                    const label nearestNode =
                                                        ipNodeMap[ip];
                                                    stk::mesh::Entity node =
                                                        sideNodeRels
                                                            [nearestNode];

                                                    scalar* snvalue =
                                                        stk::mesh::field_data(
                                                            nodeSideSTKFieldRef,
                                                            node);

                                                    for (label i = 0;
                                                         i < SPATIAL_DIM;
                                                         i++)
                                                    {
                                                        snvalue[i] -= areaVec
                                                            [ip * SPATIAL_DIM +
                                                             i];
                                                    }
                                                }
                                            }
                                        }
                                    }

                                    // multiply by speed and normalize
                                    {
                                        // select all nodes relevant to the node
                                        // side field
                                        stk::mesh::Selector selAllNodes =
                                            metaData.universal_part() &
                                            stk::mesh::selectUnion(
                                                boundary->parts());

                                        stk::mesh::BucketVector const&
                                            sideNodeBuckets =
                                                bulkData.get_buckets(
                                                    stk::topology::NODE_RANK,
                                                    selAllNodes);
                                        for (stk::mesh::BucketVector::
                                                 const_iterator ib =
                                                     sideNodeBuckets.begin();
                                             ib != sideNodeBuckets.end();
                                             ++ib)
                                        {
                                            stk::mesh::Bucket& sideNodeBucket =
                                                **ib;
                                            const stk::mesh::Bucket::size_type
                                                nSideNodesPerBucket =
                                                    sideNodeBucket.size();
                                            scalar* snvalue =
                                                stk::mesh::field_data(
                                                    nodeSideSTKFieldRef,
                                                    sideNodeBucket);
                                            scalar* value =
                                                stk::mesh::field_data(
                                                    stkFieldRef,
                                                    sideNodeBucket);

                                            for (stk::mesh::Bucket::size_type
                                                     iSideNode = 0;
                                                 iSideNode <
                                                 nSideNodesPerBucket;
                                                 ++iSideNode)
                                            {
                                                const scalar* coords =
                                                    stk::mesh::field_data(
                                                        *coordsSTKFieldPtr,
                                                        sideNodeBucket,
                                                        iSideNode);

                                                scalar asq = 0.0;
                                                for (label i = 0;
                                                     i < SPATIAL_DIM;
                                                     ++i)
                                                {
                                                    const scalar axi =
                                                        snvalue[SPATIAL_DIM *
                                                                    iSideNode +
                                                                i];
                                                    asq += axi * axi;
                                                }
                                                const scalar amag =
                                                    std::sqrt(asq);

                                                for (label i = 0;
                                                     i < SPATIAL_DIM;
                                                     i++)
                                                {
                                                    snvalue[SPATIAL_DIM *
                                                                iSideNode +
                                                            i] *=
                                                        normalSpeed / amag;
                                                }

                                                // add frame velocity
                                                for (label i = 0;
                                                     i < SPATIAL_DIM;
                                                     i++)
                                                {
                                                    for (label j = 0;
                                                         j < SPATIAL_DIM;
                                                         j++)
                                                    {
                                                        snvalue[SPATIAL_DIM *
                                                                    iSideNode +
                                                                i] +=
                                                            p_mat
                                                                [i * SPATIAL_DIM +
                                                                 j] *
                                                            (coords[j] -
                                                             p_ori[j]);
                                                    }
                                                }

                                                // override internal node values
                                                if (correctedBoundaryNodeValues_)
                                                {
                                                    for (label i = 0;
                                                         i < SPATIAL_DIM;
                                                         i++)
                                                    {
                                                        value[SPATIAL_DIM *
                                                                  iSideNode +
                                                              i] = snvalue
                                                            [SPATIAL_DIM *
                                                                 iSideNode +
                                                             i];
                                                    }
                                                }
                                            }
                                        }
                                    }

                                    // Update side field on the current boundary
                                    this->sideFieldRef().interpolate(
                                        this->nodeSideFieldRef(),
                                        iZone,
                                        iBoundary,
                                        this->isShifted());
                                }
                                break;

                            case inputDataType::expression:
                                {
                                    typedef exprtk::symbol_table<scalar>
                                        symbol_table_t;
                                    typedef exprtk::expression<scalar>
                                        expression_t;
                                    typedef exprtk::parser<scalar> parser_t;

                                    symbol_table_t symbol_table;
                                    symbol_table.add_constants();

                                    // declare function vars
                                    scalar t, x, y, z;
                                    symbol_table.add_variable("t", t);
                                    symbol_table.add_variable("x", x);
                                    symbol_table.add_variable("y", y);
                                    symbol_table.add_variable("z", z);

                                    expression_t expression;
                                    expression.register_symbol_table(
                                        symbol_table);

                                    parser_t parser;

                                    if (!parser.compile(data.expression()[0],
                                                        expression))
                                    {
                                        errorMsg("Error in the expression "
                                                 "provided for " +
                                                 this->name());
                                    }

                                    // set time value
                                    t = this->meshRef().controlsRef().time;

                                    // get fields
                                    auto& stkFieldRef = this->stkFieldRef();
                                    auto& nodeSideSTKFieldRef =
                                        this->nodeSideFieldRef().stkFieldRef();

                                    // extract geometric fields
                                    const STKScalarField* coordsSTKFieldPtr =
                                        metaData.template get_field<scalar>(
                                            stk::topology::NODE_RANK,
                                            this->getCoordinatesID_(iZone));

                                    // Find exposed node area vector (sum of
                                    // surrounding ip's)
                                    {
                                        // define some common selectors
                                        stk::mesh::Selector selAllSides =
                                            metaData.universal_part() &
                                            stk::mesh::selectUnion(
                                                boundary->parts());

                                        // Get fields
                                        STKScalarField&
                                            exposedAreaVecSTKFieldRef =
                                                *metaData.get_field<scalar>(
                                                    metaData.side_rank(),
                                                    this->getExposedAreaVectorID_(
                                                        iZone));

                                        // define vector of parent topos; should
                                        // always be UNITY in size
                                        std::vector<stk::topology> parentTopo;

                                        // setup for buckets; union parts and
                                        // ask for universal part
                                        stk::mesh::BucketVector const&
                                            sideBuckets = bulkData.get_buckets(
                                                metaData.side_rank(),
                                                selAllSides);

                                        for (stk::mesh::BucketVector::
                                                 const_iterator ib =
                                                     sideBuckets.begin();
                                             ib != sideBuckets.end();
                                             ++ib)
                                        {
                                            stk::mesh::Bucket& sideBucket =
                                                **ib;

                                            // extract master element
                                            MasterElement* meFC =
                                                MasterElementRepo::
                                                    get_surface_master_element(
                                                        sideBucket.topology());

                                            // mapping from ip to nodes for this
                                            // ordinal
                                            const label* ipNodeMap =
                                                meFC->ipNodeMap();

                                            // extract master element specifics
                                            const label numScsIp =
                                                meFC->numIntPoints_;

                                            const stk::mesh::Bucket::size_type
                                                nSidesPerBucket =
                                                    sideBucket.size();
                                            for (stk::mesh::Bucket::size_type
                                                     iSide = 0;
                                                 iSide < nSidesPerBucket;
                                                 ++iSide)
                                            {
                                                // get face
                                                stk::mesh::Entity side =
                                                    sideBucket[iSide];

                                                stk::mesh::Entity const*
                                                    sideNodeRels =
                                                        bulkData.begin_nodes(
                                                            side);

                                                // face data
                                                scalar* areaVec =
                                                    stk::mesh::field_data(
                                                        exposedAreaVecSTKFieldRef,
                                                        sideBucket,
                                                        iSide);

                                                // accumulate ip areas around
                                                // node
                                                for (label ip = 0;
                                                     ip < numScsIp;
                                                     ++ip)
                                                {
                                                    const label nearestNode =
                                                        ipNodeMap[ip];
                                                    stk::mesh::Entity node =
                                                        sideNodeRels
                                                            [nearestNode];

                                                    scalar* snvalue =
                                                        stk::mesh::field_data(
                                                            nodeSideSTKFieldRef,
                                                            node);

                                                    for (label i = 0;
                                                         i < SPATIAL_DIM;
                                                         i++)
                                                    {
                                                        snvalue[i] -= areaVec
                                                            [ip * SPATIAL_DIM +
                                                             i];
                                                    }
                                                }
                                            }
                                        }
                                    }

                                    // multiply by speed and normalize
                                    {
                                        // select all nodes relevant to the node
                                        // side field
                                        stk::mesh::Selector selAllNodes =
                                            metaData.universal_part() &
                                            stk::mesh::selectUnion(
                                                boundary->parts());

                                        stk::mesh::BucketVector const&
                                            sideNodeBuckets =
                                                bulkData.get_buckets(
                                                    stk::topology::NODE_RANK,
                                                    selAllNodes);
                                        for (stk::mesh::BucketVector::
                                                 const_iterator ib =
                                                     sideNodeBuckets.begin();
                                             ib != sideNodeBuckets.end();
                                             ++ib)
                                        {
                                            stk::mesh::Bucket& sideNodeBucket =
                                                **ib;
                                            const stk::mesh::Bucket::size_type
                                                nSideNodesPerBucket =
                                                    sideNodeBucket.size();
                                            scalar* snvalueb =
                                                stk::mesh::field_data(
                                                    nodeSideSTKFieldRef,
                                                    sideNodeBucket);
                                            scalar* valueb =
                                                stk::mesh::field_data(
                                                    stkFieldRef,
                                                    sideNodeBucket);

                                            for (stk::mesh::Bucket::size_type
                                                     iSideNode = 0;
                                                 iSideNode <
                                                 nSideNodesPerBucket;
                                                 ++iSideNode)
                                            {
                                                const scalar* coords =
                                                    stk::mesh::field_data(
                                                        *coordsSTKFieldPtr,
                                                        sideNodeBucket,
                                                        iSideNode);

#if SPATIAL_DIM == 3
                                                x = coords[0];
                                                y = coords[1];
                                                z = coords[2];
#elif SPATIAL_DIM == 2
                                                x = coords[0];
                                                y = coords[1];
#endif

                                                scalar normalSpeed =
                                                    expression.value();

                                                scalar asq = 0.0;
                                                for (label i = 0;
                                                     i < SPATIAL_DIM;
                                                     ++i)
                                                {
                                                    const scalar axi =
                                                        snvalueb[SPATIAL_DIM *
                                                                     iSideNode +
                                                                 i];
                                                    asq += axi * axi;
                                                }
                                                const scalar amag =
                                                    std::sqrt(asq);

                                                for (label i = 0;
                                                     i < SPATIAL_DIM;
                                                     i++)
                                                {
                                                    snvalueb[SPATIAL_DIM *
                                                                 iSideNode +
                                                             i] *=
                                                        normalSpeed / amag;
                                                }

                                                // add frame velocity
                                                for (label i = 0;
                                                     i < SPATIAL_DIM;
                                                     i++)
                                                {
                                                    for (label j = 0;
                                                         j < SPATIAL_DIM;
                                                         j++)
                                                    {
                                                        snvalueb[SPATIAL_DIM *
                                                                     iSideNode +
                                                                 i] +=
                                                            p_mat
                                                                [i * SPATIAL_DIM +
                                                                 j] *
                                                            (coords[j] -
                                                             p_ori[j]);
                                                    }
                                                }

                                                // override internal node values
                                                if (correctedBoundaryNodeValues_)
                                                {
                                                    for (label i = 0;
                                                         i < SPATIAL_DIM;
                                                         i++)
                                                    {
                                                        valueb[SPATIAL_DIM *
                                                                   iSideNode +
                                                               i] = snvalueb
                                                            [SPATIAL_DIM *
                                                                 iSideNode +
                                                             i];
                                                    }
                                                }
                                            }
                                        }
                                    }

                                    // Update side field on the current boundary
                                    this->sideFieldRef().interpolate(
                                        this->nodeSideFieldRef(),
                                        iZone,
                                        iBoundary,
                                        this->isShifted());
                                }
                                break;

                            case inputDataType::profileData:
                                {
                                    errorMsg("profile data not provided yet");
                                }
                                break;
                        }
                    }
                }
                break;
        }
    }
    else if (zone->meshTransforming())
    {
        switch (zone->transformationRef().type())
        {
            case meshMotionType::translating:
            case meshMotionType::rotating:
                {
                    if (boundary->frameType() ==
                        boundaryRelativeFrameType::absolute)
                    {
                        // Reset node side field
                        this->nodeSideFieldRef().setToValue(
                            std::vector<scalar>(SPATIAL_DIM, 0.0),
                            boundary->parts());

                        switch (data.type())
                        {
                            case inputDataType::null:
                                break;

                            case inputDataType::constant:
                            case inputDataType::timeTable:
                                {
                                    scalar normalSpeed =
                                        data.type() == inputDataType::constant
                                            ? *data.value()
                                            : data.interpolate(
                                                  this->meshRef()
                                                      .controlsRef()
                                                      .time)[0];
                                    // get fields
                                    auto& stkFieldRef = this->stkFieldRef();
                                    auto& nodeSideSTKFieldRef =
                                        this->nodeSideFieldRef().stkFieldRef();

                                    // Find exposed node area vector (sum of
                                    // surrounding ip's)
                                    {
                                        // define some common selectors
                                        stk::mesh::Selector selAllSides =
                                            metaData.universal_part() &
                                            stk::mesh::selectUnion(
                                                boundary->parts());

                                        // Get fields
                                        STKScalarField&
                                            exposedAreaVecSTKFieldRef =
                                                *metaData.get_field<scalar>(
                                                    metaData.side_rank(),
                                                    this->getExposedAreaVectorID_(
                                                        iZone));

                                        // define vector of parent topos; should
                                        // always be UNITY in size
                                        std::vector<stk::topology> parentTopo;

                                        // setup for buckets; union parts and
                                        // ask for universal part
                                        stk::mesh::BucketVector const&
                                            sideBuckets = bulkData.get_buckets(
                                                metaData.side_rank(),
                                                selAllSides);

                                        for (stk::mesh::BucketVector::
                                                 const_iterator ib =
                                                     sideBuckets.begin();
                                             ib != sideBuckets.end();
                                             ++ib)
                                        {
                                            stk::mesh::Bucket& sideBucket =
                                                **ib;

                                            // extract master element
                                            MasterElement* meFC =
                                                MasterElementRepo::
                                                    get_surface_master_element(
                                                        sideBucket.topology());

                                            // mapping from ip to nodes for this
                                            // ordinal
                                            const label* ipNodeMap =
                                                meFC->ipNodeMap();

                                            // extract master element specifics
                                            const label numScsIp =
                                                meFC->numIntPoints_;

                                            const stk::mesh::Bucket::size_type
                                                nSidesPerBucket =
                                                    sideBucket.size();
                                            for (stk::mesh::Bucket::size_type
                                                     iSide = 0;
                                                 iSide < nSidesPerBucket;
                                                 ++iSide)
                                            {
                                                // get face
                                                stk::mesh::Entity side =
                                                    sideBucket[iSide];

                                                stk::mesh::Entity const*
                                                    sideNodeRels =
                                                        bulkData.begin_nodes(
                                                            side);

                                                // face data
                                                scalar* areaVec =
                                                    stk::mesh::field_data(
                                                        exposedAreaVecSTKFieldRef,
                                                        sideBucket,
                                                        iSide);

                                                // accumulate ip areas around
                                                // node
                                                for (label ip = 0;
                                                     ip < numScsIp;
                                                     ++ip)
                                                {
                                                    const label nearestNode =
                                                        ipNodeMap[ip];
                                                    stk::mesh::Entity node =
                                                        sideNodeRels
                                                            [nearestNode];

                                                    scalar* snvalue =
                                                        stk::mesh::field_data(
                                                            nodeSideSTKFieldRef,
                                                            node);

                                                    for (label i = 0;
                                                         i < SPATIAL_DIM;
                                                         i++)
                                                    {
                                                        snvalue[i] -= areaVec
                                                            [ip * SPATIAL_DIM +
                                                             i];
                                                    }
                                                }
                                            }
                                        }
                                    }

                                    // multiply by speed and normalize
                                    {
                                        // select all nodes relevant to the node
                                        // side field
                                        stk::mesh::Selector selAllNodes =
                                            metaData.universal_part() &
                                            stk::mesh::selectUnion(
                                                boundary->parts());

                                        stk::mesh::BucketVector const&
                                            sideNodeBuckets =
                                                bulkData.get_buckets(
                                                    stk::topology::NODE_RANK,
                                                    selAllNodes);
                                        for (stk::mesh::BucketVector::
                                                 const_iterator ib =
                                                     sideNodeBuckets.begin();
                                             ib != sideNodeBuckets.end();
                                             ++ib)
                                        {
                                            stk::mesh::Bucket& sideNodeBucket =
                                                **ib;
                                            const stk::mesh::Bucket::size_type
                                                nSideNodesPerBucket =
                                                    sideNodeBucket.size();
                                            scalar* snvalue =
                                                stk::mesh::field_data(
                                                    nodeSideSTKFieldRef,
                                                    sideNodeBucket);
                                            scalar* value =
                                                stk::mesh::field_data(
                                                    stkFieldRef,
                                                    sideNodeBucket);

                                            for (stk::mesh::Bucket::size_type
                                                     iSideNode = 0;
                                                 iSideNode <
                                                 nSideNodesPerBucket;
                                                 ++iSideNode)
                                            {
                                                scalar asq = 0.0;
                                                for (label i = 0;
                                                     i < SPATIAL_DIM;
                                                     ++i)
                                                {
                                                    const scalar axi =
                                                        snvalue[SPATIAL_DIM *
                                                                    iSideNode +
                                                                i];
                                                    asq += axi * axi;
                                                }
                                                const scalar amag =
                                                    std::sqrt(asq);

                                                for (label i = 0;
                                                     i < SPATIAL_DIM;
                                                     i++)
                                                {
                                                    snvalue[SPATIAL_DIM *
                                                                iSideNode +
                                                            i] *=
                                                        normalSpeed / amag;
                                                }

                                                // override internal node values
                                                if (correctedBoundaryNodeValues_)
                                                {
                                                    for (label i = 0;
                                                         i < SPATIAL_DIM;
                                                         i++)
                                                    {
                                                        value[SPATIAL_DIM *
                                                                  iSideNode +
                                                              i] = snvalue
                                                            [SPATIAL_DIM *
                                                                 iSideNode +
                                                             i];
                                                    }
                                                }
                                            }
                                        }
                                    }

                                    // Update side field on the current boundary
                                    this->sideFieldRef().interpolate(
                                        this->nodeSideFieldRef(),
                                        iZone,
                                        iBoundary,
                                        this->isShifted());
                                }
                                break;

                            case inputDataType::expression:
                                {
                                    typedef exprtk::symbol_table<scalar>
                                        symbol_table_t;
                                    typedef exprtk::expression<scalar>
                                        expression_t;
                                    typedef exprtk::parser<scalar> parser_t;

                                    symbol_table_t symbol_table;
                                    symbol_table.add_constants();

                                    // declare function vars
                                    scalar t, x, y, z;
                                    symbol_table.add_variable("t", t);
                                    symbol_table.add_variable("x", x);
                                    symbol_table.add_variable("y", y);
                                    symbol_table.add_variable("z", z);

                                    expression_t expression;
                                    expression.register_symbol_table(
                                        symbol_table);

                                    parser_t parser;

                                    if (!parser.compile(data.expression()[0],
                                                        expression))
                                    {
                                        errorMsg("Error in the expression "
                                                 "provided for " +
                                                 this->name());
                                    }

                                    // set time value
                                    t = this->meshRef().controlsRef().time;

                                    // get fields
                                    auto& stkFieldRef = this->stkFieldRef();
                                    auto& nodeSideSTKFieldRef =
                                        this->nodeSideFieldRef().stkFieldRef();

                                    // extract geometric fields
                                    const STKScalarField* coordsSTKFieldPtr =
                                        metaData.template get_field<scalar>(
                                            stk::topology::NODE_RANK,
                                            this->getCoordinatesID_(iZone));

                                    // Find exposed node area vector (sum of
                                    // surrounding ip's)
                                    {
                                        // define some common selectors
                                        stk::mesh::Selector selAllSides =
                                            metaData.universal_part() &
                                            stk::mesh::selectUnion(
                                                boundary->parts());

                                        // Get fields
                                        STKScalarField&
                                            exposedAreaVecSTKFieldRef =
                                                *metaData.get_field<scalar>(
                                                    metaData.side_rank(),
                                                    this->getExposedAreaVectorID_(
                                                        iZone));

                                        // define vector of parent topos; should
                                        // always be UNITY in size
                                        std::vector<stk::topology> parentTopo;

                                        // setup for buckets; union parts and
                                        // ask for universal part
                                        stk::mesh::BucketVector const&
                                            sideBuckets = bulkData.get_buckets(
                                                metaData.side_rank(),
                                                selAllSides);

                                        for (stk::mesh::BucketVector::
                                                 const_iterator ib =
                                                     sideBuckets.begin();
                                             ib != sideBuckets.end();
                                             ++ib)
                                        {
                                            stk::mesh::Bucket& sideBucket =
                                                **ib;

                                            // extract master element
                                            MasterElement* meFC =
                                                MasterElementRepo::
                                                    get_surface_master_element(
                                                        sideBucket.topology());

                                            // mapping from ip to nodes for this
                                            // ordinal
                                            const label* ipNodeMap =
                                                meFC->ipNodeMap();

                                            // extract master element specifics
                                            const label numScsIp =
                                                meFC->numIntPoints_;

                                            const stk::mesh::Bucket::size_type
                                                nSidesPerBucket =
                                                    sideBucket.size();
                                            for (stk::mesh::Bucket::size_type
                                                     iSide = 0;
                                                 iSide < nSidesPerBucket;
                                                 ++iSide)
                                            {
                                                // get face
                                                stk::mesh::Entity side =
                                                    sideBucket[iSide];

                                                stk::mesh::Entity const*
                                                    sideNodeRels =
                                                        bulkData.begin_nodes(
                                                            side);

                                                // face data
                                                scalar* areaVec =
                                                    stk::mesh::field_data(
                                                        exposedAreaVecSTKFieldRef,
                                                        sideBucket,
                                                        iSide);

                                                // accumulate ip areas around
                                                // node
                                                for (label ip = 0;
                                                     ip < numScsIp;
                                                     ++ip)
                                                {
                                                    const label nearestNode =
                                                        ipNodeMap[ip];
                                                    stk::mesh::Entity node =
                                                        sideNodeRels
                                                            [nearestNode];

                                                    scalar* snvalue =
                                                        stk::mesh::field_data(
                                                            nodeSideSTKFieldRef,
                                                            node);

                                                    for (label i = 0;
                                                         i < SPATIAL_DIM;
                                                         i++)
                                                    {
                                                        snvalue[i] -= areaVec
                                                            [ip * SPATIAL_DIM +
                                                             i];
                                                    }
                                                }
                                            }
                                        }
                                    }

                                    // multiply by speed and normalize
                                    {
                                        // select all nodes relevant to the node
                                        // side field
                                        stk::mesh::Selector selAllNodes =
                                            metaData.universal_part() &
                                            stk::mesh::selectUnion(
                                                boundary->parts());

                                        stk::mesh::BucketVector const&
                                            sideNodeBuckets =
                                                bulkData.get_buckets(
                                                    stk::topology::NODE_RANK,
                                                    selAllNodes);
                                        for (stk::mesh::BucketVector::
                                                 const_iterator ib =
                                                     sideNodeBuckets.begin();
                                             ib != sideNodeBuckets.end();
                                             ++ib)
                                        {
                                            stk::mesh::Bucket& sideNodeBucket =
                                                **ib;
                                            const stk::mesh::Bucket::size_type
                                                nSideNodesPerBucket =
                                                    sideNodeBucket.size();
                                            scalar* snvalueb =
                                                stk::mesh::field_data(
                                                    nodeSideSTKFieldRef,
                                                    sideNodeBucket);
                                            scalar* valueb =
                                                stk::mesh::field_data(
                                                    stkFieldRef,
                                                    sideNodeBucket);

                                            for (stk::mesh::Bucket::size_type
                                                     iSideNode = 0;
                                                 iSideNode <
                                                 nSideNodesPerBucket;
                                                 ++iSideNode)
                                            {
                                                const scalar* coords =
                                                    stk::mesh::field_data(
                                                        *coordsSTKFieldPtr,
                                                        sideNodeBucket,
                                                        iSideNode);

#if SPATIAL_DIM == 3
                                                x = coords[0];
                                                y = coords[1];
                                                z = coords[2];
#elif SPATIAL_DIM == 2
                                                x = coords[0];
                                                y = coords[1];
#endif

                                                scalar normalSpeed =
                                                    expression.value();

                                                scalar asq = 0.0;
                                                for (label i = 0;
                                                     i < SPATIAL_DIM;
                                                     ++i)
                                                {
                                                    const scalar axi =
                                                        snvalueb[SPATIAL_DIM *
                                                                     iSideNode +
                                                                 i];
                                                    asq += axi * axi;
                                                }
                                                const scalar amag =
                                                    std::sqrt(asq);

                                                for (label i = 0;
                                                     i < SPATIAL_DIM;
                                                     i++)
                                                {
                                                    snvalueb[SPATIAL_DIM *
                                                                 iSideNode +
                                                             i] *=
                                                        normalSpeed / amag;
                                                }

                                                // override internal node values
                                                if (correctedBoundaryNodeValues_)
                                                {
                                                    for (label i = 0;
                                                         i < SPATIAL_DIM;
                                                         i++)
                                                    {
                                                        valueb[SPATIAL_DIM *
                                                                   iSideNode +
                                                               i] = snvalueb
                                                            [SPATIAL_DIM *
                                                                 iSideNode +
                                                             i];
                                                    }
                                                }
                                            }
                                        }
                                    }

                                    // Update side field on the current boundary
                                    this->sideFieldRef().interpolate(
                                        this->nodeSideFieldRef(),
                                        iZone,
                                        iBoundary,
                                        this->isShifted());
                                }
                                break;

                            case inputDataType::profileData:
                                {
                                    errorMsg("profile data not provided yet");
                                }
                                break;
                        }
                    }
                    else
                    {
                        // Reset node side field
                        this->nodeSideFieldRef().setToValue(
                            std::vector<scalar>(SPATIAL_DIM, 0.0),
                            boundary->parts());

                        // Get mesh velocity field
                        const auto& UmSTKFieldRef = this->simulationRef()
                                                        .meshMotionRef()
                                                        .UmRef()
                                                        .stkFieldRef();

                        switch (data.type())
                        {
                            case inputDataType::null:
                                break;

                            case inputDataType::constant:
                            case inputDataType::timeTable:
                                {
                                    scalar normalSpeed =
                                        data.type() == inputDataType::constant
                                            ? *data.value()
                                            : data.interpolate(
                                                  this->meshRef()
                                                      .controlsRef()
                                                      .time)[0];

                                    // get fields
                                    auto& nodeSideSTKFieldRef =
                                        this->nodeSideFieldRef().stkFieldRef();
                                    auto& stkFieldRef = this->stkFieldRef();

                                    // Find exposed node area vector (sum of
                                    // surrounding ip's)
                                    {
                                        // define some common selectors
                                        stk::mesh::Selector selAllSides =
                                            metaData.universal_part() &
                                            stk::mesh::selectUnion(
                                                boundary->parts());

                                        // Get fields
                                        STKScalarField&
                                            exposedAreaVecSTKFieldRef =
                                                *metaData.get_field<scalar>(
                                                    metaData.side_rank(),
                                                    this->getExposedAreaVectorID_(
                                                        iZone));

                                        // define vector of parent topos; should
                                        // always be UNITY in size
                                        std::vector<stk::topology> parentTopo;

                                        // setup for buckets; union parts and
                                        // ask for universal part
                                        stk::mesh::BucketVector const&
                                            sideBuckets = bulkData.get_buckets(
                                                metaData.side_rank(),
                                                selAllSides);

                                        for (stk::mesh::BucketVector::
                                                 const_iterator ib =
                                                     sideBuckets.begin();
                                             ib != sideBuckets.end();
                                             ++ib)
                                        {
                                            stk::mesh::Bucket& sideBucket =
                                                **ib;

                                            // extract master element
                                            MasterElement* meFC =
                                                MasterElementRepo::
                                                    get_surface_master_element(
                                                        sideBucket.topology());

                                            // mapping from ip to nodes for this
                                            // ordinal
                                            const label* ipNodeMap =
                                                meFC->ipNodeMap();

                                            // extract master element specifics
                                            const label numScsIp =
                                                meFC->numIntPoints_;

                                            const stk::mesh::Bucket::size_type
                                                nSidesPerBucket =
                                                    sideBucket.size();
                                            for (stk::mesh::Bucket::size_type
                                                     iSide = 0;
                                                 iSide < nSidesPerBucket;
                                                 ++iSide)
                                            {
                                                // get face
                                                stk::mesh::Entity side =
                                                    sideBucket[iSide];

                                                stk::mesh::Entity const*
                                                    sideNodeRels =
                                                        bulkData.begin_nodes(
                                                            side);

                                                // face data
                                                scalar* areaVec =
                                                    stk::mesh::field_data(
                                                        exposedAreaVecSTKFieldRef,
                                                        sideBucket,
                                                        iSide);

                                                // accumulate ip areas around
                                                // node
                                                for (label ip = 0;
                                                     ip < numScsIp;
                                                     ++ip)
                                                {
                                                    const label nearestNode =
                                                        ipNodeMap[ip];
                                                    stk::mesh::Entity node =
                                                        sideNodeRels
                                                            [nearestNode];

                                                    scalar* snvalue =
                                                        stk::mesh::field_data(
                                                            nodeSideSTKFieldRef,
                                                            node);

                                                    for (label i = 0;
                                                         i < SPATIAL_DIM;
                                                         i++)
                                                    {
                                                        snvalue[i] -= areaVec
                                                            [ip * SPATIAL_DIM +
                                                             i];
                                                    }
                                                }
                                            }
                                        }
                                    }

                                    // multiply by speed and normalize
                                    {
                                        // select all nodes relevant to the node
                                        // side field
                                        stk::mesh::Selector selAllNodes =
                                            metaData.universal_part() &
                                            stk::mesh::selectUnion(
                                                boundary->parts());

                                        stk::mesh::BucketVector const&
                                            sideNodeBuckets =
                                                bulkData.get_buckets(
                                                    stk::topology::NODE_RANK,
                                                    selAllNodes);
                                        for (stk::mesh::BucketVector::
                                                 const_iterator ib =
                                                     sideNodeBuckets.begin();
                                             ib != sideNodeBuckets.end();
                                             ++ib)
                                        {
                                            stk::mesh::Bucket& sideNodeBucket =
                                                **ib;
                                            const stk::mesh::Bucket::size_type
                                                nSideNodesPerBucket =
                                                    sideNodeBucket.size();
                                            scalar* snvalueb =
                                                stk::mesh::field_data(
                                                    nodeSideSTKFieldRef,
                                                    sideNodeBucket);
                                            scalar* valueb =
                                                stk::mesh::field_data(
                                                    stkFieldRef,
                                                    sideNodeBucket);

                                            for (stk::mesh::Bucket::size_type
                                                     iSideNode = 0;
                                                 iSideNode <
                                                 nSideNodesPerBucket;
                                                 ++iSideNode)
                                            {
                                                const scalar* Um =
                                                    stk::mesh::field_data(
                                                        UmSTKFieldRef,
                                                        sideNodeBucket,
                                                        iSideNode);

                                                scalar asq = 0.0;
                                                for (label i = 0;
                                                     i < SPATIAL_DIM;
                                                     ++i)
                                                {
                                                    const scalar axi =
                                                        snvalueb[SPATIAL_DIM *
                                                                     iSideNode +
                                                                 i];
                                                    asq += axi * axi;
                                                }
                                                const scalar amag =
                                                    std::sqrt(asq);

                                                for (label i = 0;
                                                     i < SPATIAL_DIM;
                                                     i++)
                                                {
                                                    snvalueb[SPATIAL_DIM *
                                                                 iSideNode +
                                                             i] *=
                                                        normalSpeed / amag;
                                                }

                                                // add mesh velocity
                                                for (label i = 0;
                                                     i < SPATIAL_DIM;
                                                     i++)
                                                {
                                                    snvalueb[SPATIAL_DIM *
                                                                 iSideNode +
                                                             i] += Um[i];
                                                }

                                                // override internal node values
                                                if (correctedBoundaryNodeValues_)
                                                {
                                                    for (label i = 0;
                                                         i < SPATIAL_DIM;
                                                         i++)
                                                    {
                                                        valueb[SPATIAL_DIM *
                                                                   iSideNode +
                                                               i] = snvalueb
                                                            [SPATIAL_DIM *
                                                                 iSideNode +
                                                             i];
                                                    }
                                                }
                                            }
                                        }
                                    }

                                    // Update side field on the current boundary
                                    this->sideFieldRef().interpolate(
                                        this->nodeSideFieldRef(),
                                        iZone,
                                        iBoundary,
                                        this->isShifted());
                                }
                                break;

                            case inputDataType::expression:
                                {
                                    typedef exprtk::symbol_table<scalar>
                                        symbol_table_t;
                                    typedef exprtk::expression<scalar>
                                        expression_t;
                                    typedef exprtk::parser<scalar> parser_t;

                                    symbol_table_t symbol_table;
                                    symbol_table.add_constants();

                                    // declare function vars
                                    scalar t, x, y, z;
                                    symbol_table.add_variable("t", t);
                                    symbol_table.add_variable("x", x);
                                    symbol_table.add_variable("y", y);
                                    symbol_table.add_variable("z", z);

                                    expression_t expression;
                                    expression.register_symbol_table(
                                        symbol_table);

                                    parser_t parser;

                                    if (!parser.compile(data.expression()[0],
                                                        expression))
                                    {
                                        errorMsg("Error in the expression "
                                                 "provided for " +
                                                 this->name());
                                    }

                                    // set time value
                                    t = this->meshRef().controlsRef().time;

                                    // get fields
                                    auto& nodeSideSTKFieldRef =
                                        this->nodeSideFieldRef().stkFieldRef();
                                    auto& stkFieldRef = this->stkFieldRef();

                                    // extract geometric fields
                                    const STKScalarField* coordsSTKFieldPtr =
                                        metaData.template get_field<scalar>(
                                            stk::topology::NODE_RANK,
                                            this->getCoordinatesID_(iZone));

                                    // Find exposed node area vector (sum of
                                    // surrounding ip's)
                                    {
                                        // define some common selectors
                                        stk::mesh::Selector selAllSides =
                                            metaData.universal_part() &
                                            stk::mesh::selectUnion(
                                                boundary->parts());

                                        // Get fields
                                        STKScalarField&
                                            exposedAreaVecSTKFieldRef =
                                                *metaData.get_field<scalar>(
                                                    metaData.side_rank(),
                                                    this->getExposedAreaVectorID_(
                                                        iZone));

                                        // define vector of parent topos; should
                                        // always be UNITY in size
                                        std::vector<stk::topology> parentTopo;

                                        // setup for buckets; union parts and
                                        // ask for universal part
                                        stk::mesh::BucketVector const&
                                            sideBuckets = bulkData.get_buckets(
                                                metaData.side_rank(),
                                                selAllSides);

                                        for (stk::mesh::BucketVector::
                                                 const_iterator ib =
                                                     sideBuckets.begin();
                                             ib != sideBuckets.end();
                                             ++ib)
                                        {
                                            stk::mesh::Bucket& sideBucket =
                                                **ib;

                                            // extract master element
                                            MasterElement* meFC =
                                                MasterElementRepo::
                                                    get_surface_master_element(
                                                        sideBucket.topology());

                                            // mapping from ip to nodes for this
                                            // ordinal
                                            const label* ipNodeMap =
                                                meFC->ipNodeMap();

                                            // extract master element specifics
                                            const label numScsIp =
                                                meFC->numIntPoints_;

                                            const stk::mesh::Bucket::size_type
                                                nSidesPerBucket =
                                                    sideBucket.size();
                                            for (stk::mesh::Bucket::size_type
                                                     iSide = 0;
                                                 iSide < nSidesPerBucket;
                                                 ++iSide)
                                            {
                                                // get face
                                                stk::mesh::Entity side =
                                                    sideBucket[iSide];

                                                stk::mesh::Entity const*
                                                    sideNodeRels =
                                                        bulkData.begin_nodes(
                                                            side);

                                                // face data
                                                scalar* areaVec =
                                                    stk::mesh::field_data(
                                                        exposedAreaVecSTKFieldRef,
                                                        sideBucket,
                                                        iSide);

                                                // accumulate ip areas around
                                                // node
                                                for (label ip = 0;
                                                     ip < numScsIp;
                                                     ++ip)
                                                {
                                                    const label nearestNode =
                                                        ipNodeMap[ip];
                                                    stk::mesh::Entity node =
                                                        sideNodeRels
                                                            [nearestNode];

                                                    scalar* snvalue =
                                                        stk::mesh::field_data(
                                                            nodeSideSTKFieldRef,
                                                            node);

                                                    for (label i = 0;
                                                         i < SPATIAL_DIM;
                                                         i++)
                                                    {
                                                        snvalue[i] -= areaVec
                                                            [ip * SPATIAL_DIM +
                                                             i];
                                                    }
                                                }
                                            }
                                        }
                                    }

                                    // multiply by speed and normalize
                                    {
                                        // select all nodes relevant to the node
                                        // side field
                                        stk::mesh::Selector selAllNodes =
                                            metaData.universal_part() &
                                            stk::mesh::selectUnion(
                                                boundary->parts());

                                        stk::mesh::BucketVector const&
                                            sideNodeBuckets =
                                                bulkData.get_buckets(
                                                    stk::topology::NODE_RANK,
                                                    selAllNodes);
                                        for (stk::mesh::BucketVector::
                                                 const_iterator ib =
                                                     sideNodeBuckets.begin();
                                             ib != sideNodeBuckets.end();
                                             ++ib)
                                        {
                                            stk::mesh::Bucket& sideNodeBucket =
                                                **ib;
                                            const stk::mesh::Bucket::size_type
                                                nSideNodesPerBucket =
                                                    sideNodeBucket.size();
                                            scalar* snvalueb =
                                                stk::mesh::field_data(
                                                    nodeSideSTKFieldRef,
                                                    sideNodeBucket);
                                            scalar* valueb =
                                                stk::mesh::field_data(
                                                    stkFieldRef,
                                                    sideNodeBucket);

                                            for (stk::mesh::Bucket::size_type
                                                     iSideNode = 0;
                                                 iSideNode <
                                                 nSideNodesPerBucket;
                                                 ++iSideNode)
                                            {
                                                const scalar* coords =
                                                    stk::mesh::field_data(
                                                        *coordsSTKFieldPtr,
                                                        sideNodeBucket,
                                                        iSideNode);

#if SPATIAL_DIM == 3
                                                x = coords[0];
                                                y = coords[1];
                                                z = coords[2];
#elif SPATIAL_DIM == 2
                                                x = coords[0];
                                                y = coords[1];
#endif

                                                scalar normalSpeed =
                                                    expression.value();

                                                scalar asq = 0.0;
                                                for (label i = 0;
                                                     i < SPATIAL_DIM;
                                                     ++i)
                                                {
                                                    const scalar axi =
                                                        snvalueb[SPATIAL_DIM *
                                                                     iSideNode +
                                                                 i];
                                                    asq += axi * axi;
                                                }
                                                const scalar amag =
                                                    std::sqrt(asq);

                                                for (label i = 0;
                                                     i < SPATIAL_DIM;
                                                     i++)
                                                {
                                                    snvalueb[SPATIAL_DIM *
                                                                 iSideNode +
                                                             i] *=
                                                        normalSpeed / amag;
                                                }

                                                const scalar* Um =
                                                    stk::mesh::field_data(
                                                        UmSTKFieldRef,
                                                        sideNodeBucket,
                                                        iSideNode);

                                                // add mesh velocity
                                                for (label i = 0;
                                                     i < SPATIAL_DIM;
                                                     i++)
                                                {
                                                    snvalueb[SPATIAL_DIM *
                                                                 iSideNode +
                                                             i] += Um[i];
                                                }

                                                // override internal node values
                                                if (correctedBoundaryNodeValues_)
                                                {
                                                    for (label i = 0;
                                                         i < SPATIAL_DIM;
                                                         i++)
                                                    {
                                                        valueb[SPATIAL_DIM *
                                                                   iSideNode +
                                                               i] = snvalueb
                                                            [SPATIAL_DIM *
                                                                 iSideNode +
                                                             i];
                                                    }
                                                }
                                            }
                                        }
                                    }

                                    // Update side field on the current boundary
                                    this->sideFieldRef().interpolate(
                                        this->nodeSideFieldRef(),
                                        iZone,
                                        iBoundary,
                                        this->isShifted());
                                }
                                break;

                            case inputDataType::profileData:
                                {
                                    errorMsg("profile data not provided yet");
                                }
                                break;
                        }
                    }
                }
                break;

            default:
                errorMsg("Must not reach here");
                break;
        }
    }
    else if (zone->meshDeforming())
    {
        // FIXME: Inlet boundary may be subject to mesh deformation
        assert(boundary->frameType() == boundaryRelativeFrameType::absolute);

        // Reset node side field
        this->nodeSideFieldRef().setToValue(
            std::vector<scalar>(SPATIAL_DIM, 0.0), boundary->parts());

        switch (data.type())
        {
            case inputDataType::null:
                break;

            case inputDataType::constant:
            case inputDataType::timeTable:
                {
                    scalar normalSpeed =
                        data.type() == inputDataType::constant
                            ? *data.value()
                            : data.interpolate(
                                  this->meshRef().controlsRef().time)[0];
                    // get fields
                    auto& stkFieldRef = this->stkFieldRef();
                    auto& nodeSideSTKFieldRef =
                        this->nodeSideFieldRef().stkFieldRef();

                    // Find exposed node area vector (sum of
                    // surrounding ip's)
                    {
                        // define some common selectors
                        stk::mesh::Selector selAllSides =
                            metaData.universal_part() &
                            stk::mesh::selectUnion(boundary->parts());

                        // Get fields
                        STKScalarField& exposedAreaVecSTKFieldRef =
                            *metaData.get_field<scalar>(
                                metaData.side_rank(),
                                this->getExposedAreaVectorID_(iZone));

                        // define vector of parent topos; should
                        // always be UNITY in size
                        std::vector<stk::topology> parentTopo;

                        // setup for buckets; union parts and
                        // ask for universal part
                        stk::mesh::BucketVector const& sideBuckets =
                            bulkData.get_buckets(metaData.side_rank(),
                                                 selAllSides);

                        for (stk::mesh::BucketVector::const_iterator ib =
                                 sideBuckets.begin();
                             ib != sideBuckets.end();
                             ++ib)
                        {
                            stk::mesh::Bucket& sideBucket = **ib;

                            // extract master element
                            MasterElement* meFC =
                                MasterElementRepo::get_surface_master_element(
                                    sideBucket.topology());

                            // mapping from ip to nodes for this
                            // ordinal
                            const label* ipNodeMap = meFC->ipNodeMap();

                            // extract master element specifics
                            const label numScsIp = meFC->numIntPoints_;

                            const stk::mesh::Bucket::size_type nSidesPerBucket =
                                sideBucket.size();
                            for (stk::mesh::Bucket::size_type iSide = 0;
                                 iSide < nSidesPerBucket;
                                 ++iSide)
                            {
                                // get face
                                stk::mesh::Entity side = sideBucket[iSide];

                                stk::mesh::Entity const* sideNodeRels =
                                    bulkData.begin_nodes(side);

                                // face data
                                scalar* areaVec = stk::mesh::field_data(
                                    exposedAreaVecSTKFieldRef,
                                    sideBucket,
                                    iSide);

                                // accumulate ip areas around
                                // node
                                for (label ip = 0; ip < numScsIp; ++ip)
                                {
                                    const label nearestNode = ipNodeMap[ip];
                                    stk::mesh::Entity node =
                                        sideNodeRels[nearestNode];

                                    scalar* snvalue = stk::mesh::field_data(
                                        nodeSideSTKFieldRef, node);

                                    for (label i = 0; i < SPATIAL_DIM; i++)
                                    {
                                        snvalue[i] -=
                                            areaVec[ip * SPATIAL_DIM + i];
                                    }
                                }
                            }
                        }
                    }

                    // multiply by speed and normalize
                    {
                        // select all nodes relevant to the node
                        // side field
                        stk::mesh::Selector selAllNodes =
                            metaData.universal_part() &
                            stk::mesh::selectUnion(boundary->parts());

                        stk::mesh::BucketVector const& sideNodeBuckets =
                            bulkData.get_buckets(stk::topology::NODE_RANK,
                                                 selAllNodes);
                        for (stk::mesh::BucketVector::const_iterator ib =
                                 sideNodeBuckets.begin();
                             ib != sideNodeBuckets.end();
                             ++ib)
                        {
                            stk::mesh::Bucket& sideNodeBucket = **ib;
                            const stk::mesh::Bucket::size_type
                                nSideNodesPerBucket = sideNodeBucket.size();
                            scalar* snvalue = stk::mesh::field_data(
                                nodeSideSTKFieldRef, sideNodeBucket);
                            scalar* value = stk::mesh::field_data(
                                stkFieldRef, sideNodeBucket);

                            for (stk::mesh::Bucket::size_type iSideNode = 0;
                                 iSideNode < nSideNodesPerBucket;
                                 ++iSideNode)
                            {
                                scalar asq = 0.0;
                                for (label i = 0; i < SPATIAL_DIM; ++i)
                                {
                                    const scalar axi =
                                        snvalue[SPATIAL_DIM * iSideNode + i];
                                    asq += axi * axi;
                                }
                                const scalar amag = std::sqrt(asq);

                                for (label i = 0; i < SPATIAL_DIM; i++)
                                {
                                    snvalue[SPATIAL_DIM * iSideNode + i] *=
                                        normalSpeed / amag;
                                }

                                // override internal node values
                                if (correctedBoundaryNodeValues_)
                                {
                                    for (label i = 0; i < SPATIAL_DIM; i++)
                                    {
                                        value[SPATIAL_DIM * iSideNode + i] =
                                            snvalue[SPATIAL_DIM * iSideNode +
                                                    i];
                                    }
                                }
                            }
                        }
                    }

                    // Update side field on the current boundary
                    this->sideFieldRef().interpolate(this->nodeSideFieldRef(),
                                                     iZone,
                                                     iBoundary,
                                                     this->isShifted());
                }
                break;

            case inputDataType::expression:
                {
                    typedef exprtk::symbol_table<scalar> symbol_table_t;
                    typedef exprtk::expression<scalar> expression_t;
                    typedef exprtk::parser<scalar> parser_t;

                    symbol_table_t symbol_table;
                    symbol_table.add_constants();

                    // declare function vars
                    scalar t, x, y, z;
                    symbol_table.add_variable("t", t);
                    symbol_table.add_variable("x", x);
                    symbol_table.add_variable("y", y);
                    symbol_table.add_variable("z", z);

                    expression_t expression;
                    expression.register_symbol_table(symbol_table);

                    parser_t parser;

                    if (!parser.compile(data.expression()[0], expression))
                    {
                        errorMsg("Error in the expression "
                                 "provided for " +
                                 this->name());
                    }

                    // set time value
                    t = this->meshRef().controlsRef().time;

                    // get fields
                    auto& stkFieldRef = this->stkFieldRef();
                    auto& nodeSideSTKFieldRef =
                        this->nodeSideFieldRef().stkFieldRef();

                    // extract geometric fields
                    const STKScalarField* coordsSTKFieldPtr =
                        metaData.template get_field<scalar>(
                            stk::topology::NODE_RANK,
                            this->getCoordinatesID_(iZone));

                    // Find exposed node area vector (sum of
                    // surrounding ip's)
                    {
                        // define some common selectors
                        stk::mesh::Selector selAllSides =
                            metaData.universal_part() &
                            stk::mesh::selectUnion(boundary->parts());

                        // Get fields
                        STKScalarField& exposedAreaVecSTKFieldRef =
                            *metaData.get_field<scalar>(
                                metaData.side_rank(),
                                this->getExposedAreaVectorID_(iZone));

                        // define vector of parent topos; should
                        // always be UNITY in size
                        std::vector<stk::topology> parentTopo;

                        // setup for buckets; union parts and
                        // ask for universal part
                        stk::mesh::BucketVector const& sideBuckets =
                            bulkData.get_buckets(metaData.side_rank(),
                                                 selAllSides);

                        for (stk::mesh::BucketVector::const_iterator ib =
                                 sideBuckets.begin();
                             ib != sideBuckets.end();
                             ++ib)
                        {
                            stk::mesh::Bucket& sideBucket = **ib;

                            // extract master element
                            MasterElement* meFC =
                                MasterElementRepo::get_surface_master_element(
                                    sideBucket.topology());

                            // mapping from ip to nodes for this
                            // ordinal
                            const label* ipNodeMap = meFC->ipNodeMap();

                            // extract master element specifics
                            const label numScsIp = meFC->numIntPoints_;

                            const stk::mesh::Bucket::size_type nSidesPerBucket =
                                sideBucket.size();
                            for (stk::mesh::Bucket::size_type iSide = 0;
                                 iSide < nSidesPerBucket;
                                 ++iSide)
                            {
                                // get face
                                stk::mesh::Entity side = sideBucket[iSide];

                                stk::mesh::Entity const* sideNodeRels =
                                    bulkData.begin_nodes(side);

                                // face data
                                scalar* areaVec = stk::mesh::field_data(
                                    exposedAreaVecSTKFieldRef,
                                    sideBucket,
                                    iSide);

                                // accumulate ip areas around
                                // node
                                for (label ip = 0; ip < numScsIp; ++ip)
                                {
                                    const label nearestNode = ipNodeMap[ip];
                                    stk::mesh::Entity node =
                                        sideNodeRels[nearestNode];

                                    scalar* snvalue = stk::mesh::field_data(
                                        nodeSideSTKFieldRef, node);

                                    for (label i = 0; i < SPATIAL_DIM; i++)
                                    {
                                        snvalue[i] -=
                                            areaVec[ip * SPATIAL_DIM + i];
                                    }
                                }
                            }
                        }
                    }

                    // multiply by speed and normalize
                    {
                        // select all nodes relevant to the node
                        // side field
                        stk::mesh::Selector selAllNodes =
                            metaData.universal_part() &
                            stk::mesh::selectUnion(boundary->parts());

                        stk::mesh::BucketVector const& sideNodeBuckets =
                            bulkData.get_buckets(stk::topology::NODE_RANK,
                                                 selAllNodes);
                        for (stk::mesh::BucketVector::const_iterator ib =
                                 sideNodeBuckets.begin();
                             ib != sideNodeBuckets.end();
                             ++ib)
                        {
                            stk::mesh::Bucket& sideNodeBucket = **ib;
                            const stk::mesh::Bucket::size_type
                                nSideNodesPerBucket = sideNodeBucket.size();
                            scalar* snvalueb = stk::mesh::field_data(
                                nodeSideSTKFieldRef, sideNodeBucket);
                            scalar* valueb = stk::mesh::field_data(
                                stkFieldRef, sideNodeBucket);

                            for (stk::mesh::Bucket::size_type iSideNode = 0;
                                 iSideNode < nSideNodesPerBucket;
                                 ++iSideNode)
                            {
                                const scalar* coords =
                                    stk::mesh::field_data(*coordsSTKFieldPtr,
                                                          sideNodeBucket,
                                                          iSideNode);

#if SPATIAL_DIM == 3
                                x = coords[0];
                                y = coords[1];
                                z = coords[2];
#elif SPATIAL_DIM == 2
                                x = coords[0];
                                y = coords[1];
#endif

                                scalar normalSpeed = expression.value();

                                scalar asq = 0.0;
                                for (label i = 0; i < SPATIAL_DIM; ++i)
                                {
                                    const scalar axi =
                                        snvalueb[SPATIAL_DIM * iSideNode + i];
                                    asq += axi * axi;
                                }
                                const scalar amag = std::sqrt(asq);

                                for (label i = 0; i < SPATIAL_DIM; i++)
                                {
                                    snvalueb[SPATIAL_DIM * iSideNode + i] *=
                                        normalSpeed / amag;
                                }

                                // override internal node values
                                if (correctedBoundaryNodeValues_)
                                {
                                    for (label i = 0; i < SPATIAL_DIM; i++)
                                    {
                                        valueb[SPATIAL_DIM * iSideNode + i] =
                                            snvalueb[SPATIAL_DIM * iSideNode +
                                                     i];
                                    }
                                }
                            }
                        }
                    }

                    // Update side field on the current boundary
                    this->sideFieldRef().interpolate(this->nodeSideFieldRef(),
                                                     iZone,
                                                     iBoundary,
                                                     this->isShifted());
                }
                break;

            case inputDataType::profileData:
                {
                    errorMsg("profile data not provided yet");
                }
                break;
        }
    }
    else
    {
        // Reset node side field
        this->nodeSideFieldRef().setToValue(
            std::vector<scalar>(SPATIAL_DIM, 0.0), boundary->parts());

        switch (data.type())
        {
            case inputDataType::null:
                break;

            case inputDataType::constant:
            case inputDataType::timeTable:
                {
                    scalar normalSpeed =
                        data.type() == inputDataType::constant
                            ? *data.value()
                            : data.interpolate(
                                  this->meshRef().controlsRef().time)[0];
                    // get fields
                    auto& stkFieldRef = this->stkFieldRef();
                    auto& nodeSideSTKFieldRef =
                        this->nodeSideFieldRef().stkFieldRef();

                    // Find exposed node area vector (sum of
                    // surrounding ip's)
                    {
                        // define some common selectors
                        stk::mesh::Selector selAllSides =
                            metaData.universal_part() &
                            stk::mesh::selectUnion(boundary->parts());

                        // Get fields
                        STKScalarField& exposedAreaVecSTKFieldRef =
                            *metaData.get_field<scalar>(
                                metaData.side_rank(),
                                this->getExposedAreaVectorID_(iZone));

                        // define vector of parent topos; should
                        // always be UNITY in size
                        std::vector<stk::topology> parentTopo;

                        // setup for buckets; union parts and
                        // ask for universal part
                        stk::mesh::BucketVector const& sideBuckets =
                            bulkData.get_buckets(metaData.side_rank(),
                                                 selAllSides);

                        for (stk::mesh::BucketVector::const_iterator ib =
                                 sideBuckets.begin();
                             ib != sideBuckets.end();
                             ++ib)
                        {
                            stk::mesh::Bucket& sideBucket = **ib;

                            // extract master element
                            MasterElement* meFC =
                                MasterElementRepo::get_surface_master_element(
                                    sideBucket.topology());

                            // mapping from ip to nodes for this
                            // ordinal
                            const label* ipNodeMap = meFC->ipNodeMap();

                            // extract master element specifics
                            const label numScsIp = meFC->numIntPoints_;

                            const stk::mesh::Bucket::size_type nSidesPerBucket =
                                sideBucket.size();
                            for (stk::mesh::Bucket::size_type iSide = 0;
                                 iSide < nSidesPerBucket;
                                 ++iSide)
                            {
                                // get face
                                stk::mesh::Entity side = sideBucket[iSide];

                                stk::mesh::Entity const* sideNodeRels =
                                    bulkData.begin_nodes(side);

                                // face data
                                scalar* areaVec = stk::mesh::field_data(
                                    exposedAreaVecSTKFieldRef,
                                    sideBucket,
                                    iSide);

                                // accumulate ip areas around
                                // node
                                for (label ip = 0; ip < numScsIp; ++ip)
                                {
                                    const label nearestNode = ipNodeMap[ip];
                                    stk::mesh::Entity node =
                                        sideNodeRels[nearestNode];

                                    scalar* snvalue = stk::mesh::field_data(
                                        nodeSideSTKFieldRef, node);

                                    for (label i = 0; i < SPATIAL_DIM; i++)
                                    {
                                        snvalue[i] -=
                                            areaVec[ip * SPATIAL_DIM + i];
                                    }
                                }
                            }
                        }
                    }

                    // multiply by speed and normalize
                    {
                        // select all nodes relevant to the node
                        // side field
                        stk::mesh::Selector selAllNodes =
                            metaData.universal_part() &
                            stk::mesh::selectUnion(boundary->parts());

                        stk::mesh::BucketVector const& sideNodeBuckets =
                            bulkData.get_buckets(stk::topology::NODE_RANK,
                                                 selAllNodes);
                        for (stk::mesh::BucketVector::const_iterator ib =
                                 sideNodeBuckets.begin();
                             ib != sideNodeBuckets.end();
                             ++ib)
                        {
                            stk::mesh::Bucket& sideNodeBucket = **ib;
                            const stk::mesh::Bucket::size_type
                                nSideNodesPerBucket = sideNodeBucket.size();
                            scalar* snvalue = stk::mesh::field_data(
                                nodeSideSTKFieldRef, sideNodeBucket);
                            scalar* value = stk::mesh::field_data(
                                stkFieldRef, sideNodeBucket);

                            for (stk::mesh::Bucket::size_type iSideNode = 0;
                                 iSideNode < nSideNodesPerBucket;
                                 ++iSideNode)
                            {
                                scalar asq = 0.0;
                                for (label i = 0; i < SPATIAL_DIM; ++i)
                                {
                                    const scalar axi =
                                        snvalue[SPATIAL_DIM * iSideNode + i];
                                    asq += axi * axi;
                                }
                                const scalar amag = std::sqrt(asq);

                                for (label i = 0; i < SPATIAL_DIM; i++)
                                {
                                    snvalue[SPATIAL_DIM * iSideNode + i] *=
                                        normalSpeed / amag;
                                }

                                // override internal node values
                                if (correctedBoundaryNodeValues_)
                                {
                                    for (label i = 0; i < SPATIAL_DIM; i++)
                                    {
                                        value[SPATIAL_DIM * iSideNode + i] =
                                            snvalue[SPATIAL_DIM * iSideNode +
                                                    i];
                                    }
                                }
                            }
                        }
                    }

                    // Update side field on the current boundary
                    this->sideFieldRef().interpolate(this->nodeSideFieldRef(),
                                                     iZone,
                                                     iBoundary,
                                                     this->isShifted());
                }
                break;

            case inputDataType::expression:
                {
                    typedef exprtk::symbol_table<scalar> symbol_table_t;
                    typedef exprtk::expression<scalar> expression_t;
                    typedef exprtk::parser<scalar> parser_t;

                    symbol_table_t symbol_table;
                    symbol_table.add_constants();

                    // declare function vars
                    scalar t, x, y, z;
                    symbol_table.add_variable("t", t);
                    symbol_table.add_variable("x", x);
                    symbol_table.add_variable("y", y);
                    symbol_table.add_variable("z", z);

                    expression_t expression;
                    expression.register_symbol_table(symbol_table);

                    parser_t parser;

                    if (!parser.compile(data.expression()[0], expression))
                    {
                        errorMsg("Error in the expression "
                                 "provided for " +
                                 this->name());
                    }

                    // set time value
                    t = this->meshRef().controlsRef().time;

                    // get fields
                    auto& stkFieldRef = this->stkFieldRef();
                    auto& nodeSideSTKFieldRef =
                        this->nodeSideFieldRef().stkFieldRef();

                    // extract geometric fields
                    const STKScalarField* coordsSTKFieldPtr =
                        metaData.template get_field<scalar>(
                            stk::topology::NODE_RANK,
                            this->getCoordinatesID_(iZone));

                    // Find exposed node area vector (sum of
                    // surrounding ip's)
                    {
                        // define some common selectors
                        stk::mesh::Selector selAllSides =
                            metaData.universal_part() &
                            stk::mesh::selectUnion(boundary->parts());

                        // Get fields
                        STKScalarField& exposedAreaVecSTKFieldRef =
                            *metaData.get_field<scalar>(
                                metaData.side_rank(),
                                this->getExposedAreaVectorID_(iZone));

                        // define vector of parent topos; should
                        // always be UNITY in size
                        std::vector<stk::topology> parentTopo;

                        // setup for buckets; union parts and
                        // ask for universal part
                        stk::mesh::BucketVector const& sideBuckets =
                            bulkData.get_buckets(metaData.side_rank(),
                                                 selAllSides);

                        for (stk::mesh::BucketVector::const_iterator ib =
                                 sideBuckets.begin();
                             ib != sideBuckets.end();
                             ++ib)
                        {
                            stk::mesh::Bucket& sideBucket = **ib;

                            // extract master element
                            MasterElement* meFC =
                                MasterElementRepo::get_surface_master_element(
                                    sideBucket.topology());

                            // mapping from ip to nodes for this
                            // ordinal
                            const label* ipNodeMap = meFC->ipNodeMap();

                            // extract master element specifics
                            const label numScsIp = meFC->numIntPoints_;

                            const stk::mesh::Bucket::size_type nSidesPerBucket =
                                sideBucket.size();
                            for (stk::mesh::Bucket::size_type iSide = 0;
                                 iSide < nSidesPerBucket;
                                 ++iSide)
                            {
                                // get face
                                stk::mesh::Entity side = sideBucket[iSide];

                                stk::mesh::Entity const* sideNodeRels =
                                    bulkData.begin_nodes(side);

                                // face data
                                scalar* areaVec = stk::mesh::field_data(
                                    exposedAreaVecSTKFieldRef,
                                    sideBucket,
                                    iSide);

                                // accumulate ip areas around
                                // node
                                for (label ip = 0; ip < numScsIp; ++ip)
                                {
                                    const label nearestNode = ipNodeMap[ip];
                                    stk::mesh::Entity node =
                                        sideNodeRels[nearestNode];

                                    scalar* snvalue = stk::mesh::field_data(
                                        nodeSideSTKFieldRef, node);

                                    for (label i = 0; i < SPATIAL_DIM; i++)
                                    {
                                        snvalue[i] -=
                                            areaVec[ip * SPATIAL_DIM + i];
                                    }
                                }
                            }
                        }
                    }

                    // multiply by speed and normalize
                    {
                        // select all nodes relevant to the node
                        // side field
                        stk::mesh::Selector selAllNodes =
                            metaData.universal_part() &
                            stk::mesh::selectUnion(boundary->parts());

                        stk::mesh::BucketVector const& sideNodeBuckets =
                            bulkData.get_buckets(stk::topology::NODE_RANK,
                                                 selAllNodes);
                        for (stk::mesh::BucketVector::const_iterator ib =
                                 sideNodeBuckets.begin();
                             ib != sideNodeBuckets.end();
                             ++ib)
                        {
                            stk::mesh::Bucket& sideNodeBucket = **ib;
                            const stk::mesh::Bucket::size_type
                                nSideNodesPerBucket = sideNodeBucket.size();
                            scalar* snvalueb = stk::mesh::field_data(
                                nodeSideSTKFieldRef, sideNodeBucket);
                            scalar* valueb = stk::mesh::field_data(
                                stkFieldRef, sideNodeBucket);

                            for (stk::mesh::Bucket::size_type iSideNode = 0;
                                 iSideNode < nSideNodesPerBucket;
                                 ++iSideNode)
                            {
                                const scalar* coords =
                                    stk::mesh::field_data(*coordsSTKFieldPtr,
                                                          sideNodeBucket,
                                                          iSideNode);

#if SPATIAL_DIM == 3
                                x = coords[0];
                                y = coords[1];
                                z = coords[2];
#elif SPATIAL_DIM == 2
                                x = coords[0];
                                y = coords[1];
#endif

                                scalar normalSpeed = expression.value();

                                scalar asq = 0.0;
                                for (label i = 0; i < SPATIAL_DIM; ++i)
                                {
                                    const scalar axi =
                                        snvalueb[SPATIAL_DIM * iSideNode + i];
                                    asq += axi * axi;
                                }
                                const scalar amag = std::sqrt(asq);

                                for (label i = 0; i < SPATIAL_DIM; i++)
                                {
                                    snvalueb[SPATIAL_DIM * iSideNode + i] *=
                                        normalSpeed / amag;
                                }

                                // override internal node values
                                if (correctedBoundaryNodeValues_)
                                {
                                    for (label i = 0; i < SPATIAL_DIM; i++)
                                    {
                                        valueb[SPATIAL_DIM * iSideNode + i] =
                                            snvalueb[SPATIAL_DIM * iSideNode +
                                                     i];
                                    }
                                }
                            }
                        }
                    }

                    // Update side field on the current boundary
                    this->sideFieldRef().interpolate(this->nodeSideFieldRef(),
                                                     iZone,
                                                     iBoundary,
                                                     this->isShifted());
                }
                break;

            case inputDataType::profileData:
                {
                    errorMsg("profile data not provided yet");
                }
                break;
        }
    }
}

void velocity::updateBoundarySideDirectionFields(label iZone, label iBoundary)
{
    assert(nodeSideFlowDirectionFieldPtr_);
    assert(sideFlowDirectionFieldPtr_);

    const auto& mesh = this->meshRef();
    const stk::mesh::BulkData& bulkData = mesh.bulkDataRef();
    const stk::mesh::MetaData& metaData = mesh.metaDataRef();

    const auto* zone = mesh.zonePtr(iZone);
    const auto* boundary = zone->boundaryPtr(iBoundary);

    // Get data from boundary condition
    auto& bc = this->boundaryConditionRef(iZone, iBoundary);
    auto flowDirectionOption = bc.rawStringValue("flow_direction_option");

    if (flowDirectionOption == "normal_to_boundary_condition")
    {
        STKScalarField& sideFlowDirectionSTKFieldRef =
            this->sideFlowDirectionFieldRef().stkFieldRef();

        // side field: straight forward from the unit normal vector

        // define vector of parent topos; should always be UNITY in size
        std::vector<stk::topology> parentTopo;

        // define some common selectors
        stk::mesh::Selector selAllSides =
            metaData.universal_part() &
            stk::mesh::selectUnion(boundary->parts());

        stk::mesh::BucketVector const& sideBuckets =
            bulkData.get_buckets(metaData.side_rank(), selAllSides);

        const STKScalarField& exposedAreaVecSTKFieldRef =
            *metaData.get_field<scalar>(metaData.side_rank(),
                                        this->getExposedAreaVectorID_(iZone));

        for (const stk::mesh::Bucket* bucket : sideBuckets)
        {
            const MasterElement* meFC =
                MasterElementRepo::get_surface_master_element(
                    bucket->topology());
            const label numScsBip = meFC->numIntPoints_;
            for (const stk::mesh::Entity face : *bucket)
            {
                const scalar* surface_vec =
                    stk::mesh::field_data(exposedAreaVecSTKFieldRef, face);
                scalar* svalue =
                    stk::mesh::field_data(sideFlowDirectionSTKFieldRef, face);

                for (label ip = 0; ip < numScsBip; ++ip)
                {
                    scalar surface_area2 = 0.0;
                    for (label j = 0; j < SPATIAL_DIM; ++j)
                    {
                        const scalar c = surface_vec[SPATIAL_DIM * ip + j];
                        surface_area2 += c * c;
                    }
                    assert(surface_area2 > 0.0);
                    const scalar sarea_inv = 1.0 / std::sqrt(surface_area2);
                    for (label j = 0; j < SPATIAL_DIM; ++j)
                    {
                        // negative sign because flow direction is
                        // defined by flow entering the domain
                        // (surface vectors point outwards by
                        // definition).
                        svalue[SPATIAL_DIM * ip + j] =
                            -surface_vec[SPATIAL_DIM * ip + j] * sarea_inv;
                    }
                }
            }
        }

        // interpolate to node-side field using area-average
        this->nodeSideFlowDirectionFieldRef().interpolate(
            this->sideFlowDirectionFieldRef(), iZone, iBoundary);
    }
    else if (flowDirectionOption == "cartesian_components")
    {
        auto& data = bc.data<SPATIAL_DIM>("flow_direction");

        switch (data.type())
        {
            case inputDataType::null:
                break;

            case inputDataType::constant:
            case inputDataType::timeTable:
                {
                    std::array<scalar, SPATIAL_DIM> inputValue;
                    if (data.type() == inputDataType::constant)
                    {
                        std::copy(data.value(),
                                  data.value() + SPATIAL_DIM,
                                  inputValue.begin());
                    }
                    else
                    {
                        inputValue = data.interpolate(
                            this->meshRef().controlsRef().time);
                    }

                    // select all nodes relevant to the node
                    // side field
                    stk::mesh::Selector selAllNodes =
                        metaData.universal_part() &
                        stk::mesh::selectUnion(boundary->parts());

                    // get field
                    auto& nodeSideFlowDirectionSTKFieldRef =
                        this->nodeSideFlowDirectionFieldRef().stkFieldRef();

                    stk::mesh::BucketVector const& sideNodeBuckets =
                        bulkData.get_buckets(stk::topology::NODE_RANK,
                                             selAllNodes);
                    for (stk::mesh::BucketVector::const_iterator ib =
                             sideNodeBuckets.begin();
                         ib != sideNodeBuckets.end();
                         ++ib)
                    {
                        stk::mesh::Bucket& sideNodeBucket = **ib;
                        const stk::mesh::Bucket::size_type nSideNodesPerBucket =
                            sideNodeBucket.size();
                        scalar* snvalueb = stk::mesh::field_data(
                            nodeSideFlowDirectionSTKFieldRef, sideNodeBucket);
                        for (stk::mesh::Bucket::size_type iSideNode = 0;
                             iSideNode < nSideNodesPerBucket;
                             ++iSideNode)
                        {
                            // normalize the flow direction
                            scalar norm = 0.0;
                            for (label i = 0; i < SPATIAL_DIM; i++)
                            {
                                norm += inputValue[i] * inputValue[i];
                            }
                            norm = std::sqrt(norm);

                            // fill
                            for (label i = 0; i < SPATIAL_DIM; i++)
                            {
                                snvalueb[SPATIAL_DIM * iSideNode + i] =
                                    inputValue[i] / norm;
                            }
                        }
                    }

                    this->sideFlowDirectionFieldRef().interpolate(
                        this->nodeSideFlowDirectionFieldRef(),
                        iZone,
                        iBoundary,
                        false);
                }
                break;

            case inputDataType::expression:
                {
                    // re-usable array
                    std::array<scalar, SPATIAL_DIM> inputValue;

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

                    for (label i = 0; i < SPATIAL_DIM; i++)
                    {
                        if (parser.compile(data.expression()[i],
                                           componentExpression))
                        {
                            expression_list.push_back(componentExpression);
                        }
                        else
                        {
                            errorMsg("Error in the expression "
                                     "provided "
                                     "for field " +
                                     this->name() + ": " +
                                     data.expression()[i]);
                        }
                    }

                    // set time value
                    t = this->meshRef().controlsRef().time;

                    // select all nodes relevant to the node
                    // side field
                    stk::mesh::Selector selAllNodes =
                        metaData.universal_part() &
                        stk::mesh::selectUnion(boundary->parts());

                    // get field
                    auto& nodeSideFlowDirectionSTKFieldRef =
                        this->nodeSideFlowDirectionFieldRef().stkFieldRef();

                    // extract geometric fields
                    const STKScalarField* coordsSTKFieldPtr =
                        metaData.template get_field<scalar>(
                            stk::topology::NODE_RANK,
                            this->getCoordinatesID_(iZone));

                    stk::mesh::BucketVector const& sideNodeBuckets =
                        bulkData.get_buckets(stk::topology::NODE_RANK,
                                             selAllNodes);
                    for (stk::mesh::BucketVector::const_iterator ib =
                             sideNodeBuckets.begin();
                         ib != sideNodeBuckets.end();
                         ++ib)
                    {
                        stk::mesh::Bucket& sideNodeBucket = **ib;
                        const stk::mesh::Bucket::size_type nSideNodesPerBucket =
                            sideNodeBucket.size();
                        scalar* snvalueb = stk::mesh::field_data(
                            nodeSideFlowDirectionSTKFieldRef, sideNodeBucket);
                        for (stk::mesh::Bucket::size_type iSideNode = 0;
                             iSideNode < nSideNodesPerBucket;
                             ++iSideNode)
                        {
                            const scalar* coords = stk::mesh::field_data(
                                *coordsSTKFieldPtr, sideNodeBucket, iSideNode);

#if SPATIAL_DIM == 3
                            x = coords[0];
                            y = coords[1];
                            z = coords[2];
#elif SPATIAL_DIM == 2
                            x = coords[0];
                            y = coords[1];
#endif
                            for (label i = 0; i < SPATIAL_DIM; i++)
                            {
                                inputValue[i] = expression_list[i].value();
                            }

                            // normalize the flow direction
                            scalar norm = 0.0;
                            for (label i = 0; i < SPATIAL_DIM; i++)
                            {
                                norm += inputValue[i] * inputValue[i];
                            }
                            norm = std::sqrt(norm);

                            for (label i = 0; i < SPATIAL_DIM; i++)
                            {
                                snvalueb[SPATIAL_DIM * iSideNode + i] =
                                    inputValue[i] / norm;
                            }
                        }
                    }

                    this->sideFlowDirectionFieldRef().interpolate(
                        this->nodeSideFlowDirectionFieldRef(),
                        iZone,
                        iBoundary,
                        false);
                }
                break;

            case inputDataType::profileData:
                {
                    errorMsg("profile data not provided yet");
                }
                break;
        }
    }
    else if (flowDirectionOption == "cylindrical_components")
    {
        auto& data = bc.data<SPATIAL_DIM>("flow_direction");
        auto& rotAxisData = bc.data<SPATIAL_DIM>("rotation_axis");

        switch (data.type())
        {
            case inputDataType::null:
                break;

            case inputDataType::constant:
            case inputDataType::timeTable:
                {
                    std::array<scalar, SPATIAL_DIM> inputValue;
                    if (data.type() == inputDataType::constant)
                    {
                        std::copy(data.value(),
                                  data.value() + SPATIAL_DIM,
                                  inputValue.begin());
                    }
                    else
                    {
                        inputValue = data.interpolate(
                            this->meshRef().controlsRef().time);
                    }

                    // get rotation axis data
                    std::vector<scalar> rotAxisConstValue(SPATIAL_DIM);
                    std::copy(rotAxisData.value(),
                              rotAxisData.value() + SPATIAL_DIM,
                              rotAxisConstValue.begin());

                    // fixed size arrays
                    std::vector<scalar> xPrj(SPATIAL_DIM);
                    std::vector<scalar> dxNorm(SPATIAL_DIM);
                    std::vector<scalar> er(SPATIAL_DIM);
                    std::vector<scalar> etheta(SPATIAL_DIM);

                    // get pointers
                    scalar* p_xPrj = &xPrj[0];
                    scalar* p_dxNorm = &dxNorm[0];
                    scalar* p_er = &er[0];
                    scalar* p_etheta = &etheta[0];

                    // select all nodes relevant to the node
                    // side field
                    stk::mesh::Selector selAllNodes =
                        metaData.universal_part() &
                        stk::mesh::selectUnion(boundary->parts());

                    // get field
                    auto& nodeSideFlowDirectionSTKFieldRef =
                        this->nodeSideFlowDirectionFieldRef().stkFieldRef();

                    const auto& coordsSTKFieldRef = *metaData.get_field<scalar>(
                        stk::topology::NODE_RANK,
                        this->getCoordinatesID_(iZone));

                    stk::mesh::BucketVector const& sideNodeBuckets =
                        bulkData.get_buckets(stk::topology::NODE_RANK,
                                             selAllNodes);
                    for (stk::mesh::BucketVector::const_iterator ib =
                             sideNodeBuckets.begin();
                         ib != sideNodeBuckets.end();
                         ++ib)
                    {
                        stk::mesh::Bucket& sideNodeBucket = **ib;
                        const stk::mesh::Bucket::size_type nSideNodesPerBucket =
                            sideNodeBucket.size();
                        scalar* snvalue = stk::mesh::field_data(
                            nodeSideFlowDirectionSTKFieldRef, sideNodeBucket);
                        for (stk::mesh::Bucket::size_type iSideNode = 0;
                             iSideNode < nSideNodesPerBucket;
                             ++iSideNode)
                        {
                            const scalar* coords = stk::mesh::field_data(
                                coordsSTKFieldRef, sideNodeBucket, iSideNode);

                            // ensure flow direction is unit vector
                            scalar mag = 0.0;
                            for (label i = 0; i < SPATIAL_DIM; i++)
                            {
                                mag += inputValue[i] * inputValue[i];
                            }
                            mag = std::sqrt(mag);
                            for (label i = 0; i < SPATIAL_DIM; i++)
                            {
                                inputValue[i] /= mag;
                            }

                            // ensure rotation axis is unit vector
                            mag = 0.0;
                            for (label i = 0; i < SPATIAL_DIM; i++)
                            {
                                mag +=
                                    rotAxisConstValue[i] * rotAxisConstValue[i];
                            }
                            mag = std::sqrt(mag);
                            for (label i = 0; i < SPATIAL_DIM; i++)
                            {
                                rotAxisConstValue[i] /= mag;
                            }

                            // determine projected position vector onto the
                            // rotation axis
                            scalar dot = 0.0;
                            for (label i = 0; i < SPATIAL_DIM; i++)
                            {
                                dot += coords[i] * rotAxisConstValue[i];
                            }
                            for (label i = 0; i < SPATIAL_DIM; i++)
                            {
                                p_xPrj[i] = dot * rotAxisConstValue[i];
                            }

                            // get radial unit direction
                            for (label i = 0; i < SPATIAL_DIM; i++)
                            {
                                p_dxNorm[i] = coords[i] - p_xPrj[i];
                            }
                            scalar magdxNorm = 0.0;
                            for (label i = 0; i < SPATIAL_DIM; i++)
                            {
                                magdxNorm += p_dxNorm[i] * p_dxNorm[i];
                            }
                            magdxNorm = std::sqrt(magdxNorm);
                            for (label i = 0; i < SPATIAL_DIM; i++)
                            {
                                p_er[i] = p_dxNorm[i] / magdxNorm;
                            }

                            // get circumferential unit direction
                            p_etheta[0] = rotAxisConstValue[1] * p_er[2] -
                                          rotAxisConstValue[2] * p_er[1];
                            p_etheta[1] = rotAxisConstValue[2] * p_er[0] -
                                          rotAxisConstValue[0] * p_er[2];
                            p_etheta[2] = rotAxisConstValue[0] * p_er[1] -
                                          rotAxisConstValue[1] * p_er[0];

                            snvalue[SPATIAL_DIM * iSideNode + 0] =
                                inputValue[0] * p_er[0] +
                                inputValue[1] * p_etheta[0] +
                                inputValue[2] * rotAxisConstValue[0];
                            snvalue[SPATIAL_DIM * iSideNode + 1] =
                                inputValue[0] * p_er[1] +
                                inputValue[1] * p_etheta[1] +
                                inputValue[2] * rotAxisConstValue[1];
                            snvalue[SPATIAL_DIM * iSideNode + 2] =
                                inputValue[0] * p_er[2] +
                                inputValue[1] * p_etheta[2] +
                                inputValue[2] * rotAxisConstValue[2];
                        }
                    }

                    this->sideFlowDirectionFieldRef().interpolate(
                        this->nodeSideFlowDirectionFieldRef(),
                        iZone,
                        iBoundary,
                        false);
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

                    for (label i = 0; i < SPATIAL_DIM; i++)
                    {
                        if (parser.compile(data.expression()[i],
                                           componentExpression))
                        {
                            expression_list.push_back(componentExpression);
                        }
                        else
                        {
                            errorMsg("Error in the expression "
                                     "provided "
                                     "for field " +
                                     this->name() + ": " +
                                     data.expression()[i]);
                        }
                    }

                    // set time value
                    t = this->meshRef().controlsRef().time;

                    std::vector<scalar> inputValue(SPATIAL_DIM);

                    // get rotation axis data
                    std::vector<scalar> rotAxisConstValue(SPATIAL_DIM);
                    std::copy(rotAxisData.value(),
                              rotAxisData.value() + SPATIAL_DIM,
                              rotAxisConstValue.begin());

                    // fixed size arrays
                    std::vector<scalar> xPrj(SPATIAL_DIM);
                    std::vector<scalar> dxNorm(SPATIAL_DIM);
                    std::vector<scalar> er(SPATIAL_DIM);
                    std::vector<scalar> etheta(SPATIAL_DIM);

                    // get pointers
                    scalar* p_xPrj = &xPrj[0];
                    scalar* p_dxNorm = &dxNorm[0];
                    scalar* p_er = &er[0];
                    scalar* p_etheta = &etheta[0];

                    // select all nodes relevant to the node
                    // side field
                    stk::mesh::Selector selAllNodes =
                        metaData.universal_part() &
                        stk::mesh::selectUnion(boundary->parts());

                    // get field
                    auto& nodeSideFlowDirectionSTKFieldRef =
                        this->nodeSideFlowDirectionFieldRef().stkFieldRef();

                    // extract geometric fields
                    const STKScalarField* coordsSTKFieldPtr =
                        metaData.template get_field<scalar>(
                            stk::topology::NODE_RANK,
                            this->getCoordinatesID_(iZone));

                    stk::mesh::BucketVector const& sideNodeBuckets =
                        bulkData.get_buckets(stk::topology::NODE_RANK,
                                             selAllNodes);
                    for (stk::mesh::BucketVector::const_iterator ib =
                             sideNodeBuckets.begin();
                         ib != sideNodeBuckets.end();
                         ++ib)
                    {
                        stk::mesh::Bucket& sideNodeBucket = **ib;
                        const stk::mesh::Bucket::size_type nSideNodesPerBucket =
                            sideNodeBucket.size();
                        scalar* snvalue = stk::mesh::field_data(
                            nodeSideFlowDirectionSTKFieldRef, sideNodeBucket);
                        for (stk::mesh::Bucket::size_type iSideNode = 0;
                             iSideNode < nSideNodesPerBucket;
                             ++iSideNode)
                        {
                            const scalar* coords = stk::mesh::field_data(
                                *coordsSTKFieldPtr, sideNodeBucket, iSideNode);

#if SPATIAL_DIM == 3
                            x = coords[0];
                            y = coords[1];
                            z = coords[2];
#elif SPATIAL_DIM == 2
                            x = coords[0];
                            y = coords[1];
#endif

                            for (label i = 0; i < SPATIAL_DIM; i++)
                            {
                                inputValue[i] = expression_list[i].value();
                            }

                            // ensure flow direction is unit vector
                            scalar mag = 0.0;
                            for (label i = 0; i < SPATIAL_DIM; i++)
                            {
                                mag += inputValue[i] * inputValue[i];
                            }
                            mag = std::sqrt(mag);
                            for (label i = 0; i < SPATIAL_DIM; i++)
                            {
                                inputValue[i] /= mag;
                            }

                            // ensure rotation axis is unit vector
                            mag = 0.0;
                            for (label i = 0; i < SPATIAL_DIM; i++)
                            {
                                mag +=
                                    rotAxisConstValue[i] * rotAxisConstValue[i];
                            }
                            mag = std::sqrt(mag);
                            for (label i = 0; i < SPATIAL_DIM; i++)
                            {
                                rotAxisConstValue[i] /= mag;
                            }

                            // determine projected position vector onto the
                            // rotation axis
                            scalar dot = 0.0;
                            for (label i = 0; i < SPATIAL_DIM; i++)
                            {
                                dot += coords[i] * rotAxisConstValue[i];
                            }
                            for (label i = 0; i < SPATIAL_DIM; i++)
                            {
                                p_xPrj[i] = dot * rotAxisConstValue[i];
                            }

                            // get radial unit direction
                            for (label i = 0; i < SPATIAL_DIM; i++)
                            {
                                p_dxNorm[i] = coords[i] - p_xPrj[i];
                            }
                            scalar magdxNorm = 0.0;
                            for (label i = 0; i < SPATIAL_DIM; i++)
                            {
                                magdxNorm += p_dxNorm[i] * p_dxNorm[i];
                            }
                            magdxNorm = std::sqrt(magdxNorm);
                            for (label i = 0; i < SPATIAL_DIM; i++)
                            {
                                p_er[i] = p_dxNorm[i] / magdxNorm;
                            }

                            // get circumferential unit direction
                            p_etheta[0] = rotAxisConstValue[1] * p_er[2] -
                                          rotAxisConstValue[2] * p_er[1];
                            p_etheta[1] = rotAxisConstValue[2] * p_er[0] -
                                          rotAxisConstValue[0] * p_er[2];
                            p_etheta[2] = rotAxisConstValue[0] * p_er[1] -
                                          rotAxisConstValue[1] * p_er[0];

                            snvalue[SPATIAL_DIM * iSideNode + 0] =
                                inputValue[0] * p_er[0] +
                                inputValue[1] * p_etheta[0] +
                                inputValue[2] * rotAxisConstValue[0];
                            snvalue[SPATIAL_DIM * iSideNode + 1] =
                                inputValue[0] * p_er[1] +
                                inputValue[1] * p_etheta[1] +
                                inputValue[2] * rotAxisConstValue[1];
                            snvalue[SPATIAL_DIM * iSideNode + 2] =
                                inputValue[0] * p_er[2] +
                                inputValue[1] * p_etheta[2] +
                                inputValue[2] * rotAxisConstValue[2];
                        }
                    }

                    this->sideFlowDirectionFieldRef().interpolate(
                        this->nodeSideFlowDirectionFieldRef(),
                        iZone,
                        iBoundary,
                        false);
                }
                break;

            case inputDataType::profileData:
                {
                    errorMsg("profile data not provided yet");
                }
                break;
        }
    }
    else
    {
        errorMsg("Must not reach here");
    }
}

#ifdef HAS_INTERFACE
void velocity::updateInterfaceSideField(label iInterface, bool master)
{
    const auto& mesh = this->meshRef();
    const stk::mesh::BulkData& bulkData = mesh.bulkDataRef();
    const stk::mesh::MetaData& metaData = mesh.metaDataRef();

    auto& interf = mesh.interfaceRef(iInterface);

    const interfaceSideInfo* interfaceSideInfoPtr =
        master ? interf.masterInfoPtr() : interf.slaveInfoPtr();

    if (interf.isFluidSolidType() &&
        interfaceSideInfoPtr->zonePtr()->meshDeforming())
    {
        // Get mesh velocity field
        const auto& UmSTKFieldRef =
            this->simulationRef().meshMotionRef().UmRef().stkFieldRef();

        // select all nodes relevant to the node side field
        stk::mesh::Selector selAllNodes =
            metaData.universal_part() &
            stk::mesh::selectUnion(interfaceSideInfoPtr->currentPartVec_);

        // get fields
        auto& stkFieldRef = this->stkFieldRef();
        auto& nodeSideSTKFieldRef = this->nodeSideFieldRef().stkFieldRef();

        stk::mesh::BucketVector const& sideNodeBuckets =
            bulkData.get_buckets(stk::topology::NODE_RANK, selAllNodes);
        for (stk::mesh::BucketVector::const_iterator ib =
                 sideNodeBuckets.begin();
             ib != sideNodeBuckets.end();
             ++ib)
        {
            stk::mesh::Bucket& sideNodeBucket = **ib;
            const stk::mesh::Bucket::size_type nSideNodesPerBucket =
                sideNodeBucket.size();
            scalar* snvalue =
                stk::mesh::field_data(nodeSideSTKFieldRef, sideNodeBucket);
            scalar* value = stk::mesh::field_data(stkFieldRef, sideNodeBucket);
            for (stk::mesh::Bucket::size_type iSideNode = 0;
                 iSideNode < nSideNodesPerBucket;
                 ++iSideNode)
            {
                const scalar* Um = stk::mesh::field_data(
                    UmSTKFieldRef, sideNodeBucket, iSideNode);

                for (label i = 0; i < SPATIAL_DIM; i++)
                {
                    snvalue[SPATIAL_DIM * iSideNode + i] = Um[i];
                }

                // override internal node values
                if (correctedBoundaryNodeValues_)
                {
                    for (label i = 0; i < SPATIAL_DIM; i++)
                    {
                        value[SPATIAL_DIM * iSideNode + i] =
                            snvalue[SPATIAL_DIM * iSideNode + i];
                    }
                }
            }
        }

        // Update side field on the current boundary
        this->sideFieldRef().interpolate(
            this->nodeSideFieldRef(), iInterface, master, this->isShifted());
    }
}
#endif /* HAS_INTERFACE */

void velocity::registerSideFlowDirectionFields(label iZone, label iBoundary)
{
    if (!nodeSideFlowDirectionFieldPtr_)
    {
        nodeSideFlowDirectionFieldPtr_ =
            std::make_unique<nodeSideField<scalar, SPATIAL_DIM>>(
                this->meshPtr(), "flow_direction_node_side", 1);
    }
    if (!sideFlowDirectionFieldPtr_)
    {
        sideFlowDirectionFieldPtr_ =
            std::make_unique<sideField<scalar, SPATIAL_DIM>>(
                this->meshPtr(), "flow_direction_side", 1);
    }

    // Put the side fields on the corresponding boundary part
    for (auto* part :
         this->meshPtr()->zonePtr(iZone)->boundaryRef(iBoundary).parts())
    {
        this->nodeSideFlowDirectionFieldRef().putFieldOnPart(*part);
        this->sideFlowDirectionFieldRef().putFieldOnPart(*part);
    }
}

// Access

sideField<label, 1>& velocity::reversalFlagRef()
{
    return *reversalFlagPtr_.get();
}

const sideField<label, 1>& velocity::reversalFlagRef() const
{
    return *reversalFlagPtr_.get();
}

nodeSideField<scalar, SPATIAL_DIM>& velocity::nodeSideFlowDirectionFieldRef()
{
    return *nodeSideFlowDirectionFieldPtr_.get();
}

const nodeSideField<scalar, SPATIAL_DIM>&
velocity::nodeSideFlowDirectionFieldRef() const
{
    return *nodeSideFlowDirectionFieldPtr_.get();
}

sideField<scalar, SPATIAL_DIM>& velocity::sideFlowDirectionFieldRef()
{
    return *sideFlowDirectionFieldPtr_.get();
}

const sideField<scalar, SPATIAL_DIM>&
velocity::sideFlowDirectionFieldRef() const
{
    return *sideFlowDirectionFieldPtr_.get();
}

} // namespace accel
