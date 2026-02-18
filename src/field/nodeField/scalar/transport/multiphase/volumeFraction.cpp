// File : volumeFraction.cpp
// Created : Sun Jan 26 2025 22:09:24 (+0100)
// Author : Mhamad Mahdi Alloush
// Description:
// Copyright (c) 2025 CCFNUM, Lucerne University of Applied Sciences and Arts.
// SPDX-License-Identifier: BSD-3-Clause

#include "volumeFraction.h"
#include "controls.h"
#include "realm.h"

namespace accel
{

volumeFraction::volumeFraction(realm* realmPtr,
                               const std::string name,
                               unsigned numberOfStates,
                               bool highResolution)
    : nodeScalarField(realmPtr,
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
            .volumeFractionInterpolationType_;

    gradientInterpolationScheme_ =
        meshRef()
            .controlsRef()
            .solverRef()
            .solverControl_.basicSettings_.interpolationSchemeType_
            .volumeFractionGradientInterpolationType_;

    // compressive advection scheme: allow beta > 1 for interface sharpening
    blendingFactorMax_ =
        meshRef()
            .controlsRef()
            .solverRef()
            .solverControl_.expertParameters_.volumeFractionBlendingFactorMax_;

    // volume fraction may only apply to fluid domains
    mediumIndependent_ = false;

    // set min/max accepted values for volume fraction
    minAcceptedValue_ = 0.0;
    maxAcceptedValue_ = 1.0;
}

void volumeFraction::updateBlendingFactorField(label iZone)
{
    assert(this->isZoneSet(iZone));

    if (advectionScheme_ != advectionSchemeType::highResolution)
        return;

    stk::mesh::BulkData& bulkData = this->meshRef().bulkDataRef();
    stk::mesh::MetaData& metaData = this->meshRef().metaDataRef();

    if (!limitGradient_)
    {
        this->updateMinMaxFields(iZone);
    }

    // store current beta for the purpose of relaxing it
    this->blendingFactorRef().updatePrevIterField(iZone);

    // Set the upper limit for beta: this is decisive of how compressive the
    // field is to be
    this->blendingFactorRef().setToValue(
        std::vector<scalar>(1, blendingFactorMax_),
        this->meshPtr()->zonePtr(iZone)->interiorParts());

    // get fields
    const STKScalarField* phiSTKFieldPtr = this->stkFieldPtr();
    const STKScalarField* gradPhiSTKFieldPtr = this->gradRef().stkFieldPtr();

    const STKScalarField* minValueSTKFieldPtr =
        this->minValueRef().stkFieldPtr();
    const STKScalarField* maxValueSTKFieldPtr =
        this->maxValueRef().stkFieldPtr();

    STKScalarField* blendingFactorSTKFieldPtr =
        this->blendingFactorRef().stkFieldPtr();

    // Interior ip
    {
        std::vector<scalar> ws_phi;
        std::vector<scalar> ws_gradPhi;
        std::vector<scalar> ws_min;
        std::vector<scalar> ws_max;
        std::vector<scalar> ws_coordinates;
        std::vector<scalar> ws_scs_areav;

        std::vector<scalar> ws_shape_function;
        std::vector<scalar> ws_coordinate_shape_function;

        stk::mesh::Selector selAllElements =
            metaData.universal_part() &
            stk::mesh::selectField(this->stkFieldRef()) &
            stk::mesh::selectUnion(
                this->meshPtr()->zonePtr(iZone)->interiorParts());

        STKScalarField* coordsSTKFieldPtr = metaData.get_field<scalar>(
            stk::topology::NODE_RANK, this->getCoordinatesID_(iZone));

        std::vector<scalar> coordIp(SPATIAL_DIM);
        std::vector<scalar> dxj(SPATIAL_DIM);

        scalar* p_coordIp = &coordIp[0];
        scalar* p_dxj = &dxj[0];

        stk::mesh::BucketVector const& elementBuckets =
            bulkData.get_buckets(stk::topology::ELEMENT_RANK, selAllElements);
        for (stk::mesh::BucketVector::const_iterator ib =
                 elementBuckets.begin();
             ib != elementBuckets.end();
             ++ib)
        {
            stk::mesh::Bucket& elementBucket = **ib;
            const stk::mesh::Bucket::size_type nElementsPerBucket =
                elementBucket.size();

            MasterElement* meSCS =
                MasterElementRepo::get_surface_master_element(
                    elementBucket.topology());
            const label nodesPerElement = meSCS->nodesPerElement_;
            const label numScsIp = meSCS->numIntPoints_;
            ws_shape_function.resize(numScsIp * nodesPerElement);

            scalar* p_shape_function = &ws_shape_function[0];

            if (interpolationScheme_ == interpolationSchemeType::linearLinear)
            {
                meSCS->shifted_shape_fcn(&p_shape_function[0]);
            }
            else
            {
                meSCS->shape_fcn(&p_shape_function[0]);
            }

            ws_coordinate_shape_function.resize(numScsIp * nodesPerElement);
            ws_phi.resize(nodesPerElement);
            ws_gradPhi.resize(nodesPerElement * SPATIAL_DIM);
            ws_min.resize(nodesPerElement);
            ws_max.resize(nodesPerElement);
            ws_coordinates.resize(nodesPerElement * SPATIAL_DIM);
            ws_scs_areav.resize(numScsIp * SPATIAL_DIM);

            scalar* p_coordinate_shape_function =
                &ws_coordinate_shape_function[0];
            scalar* p_phi = &ws_phi[0];
            scalar* p_gradPhi = &ws_gradPhi[0];
            scalar* p_min = &ws_min[0];
            scalar* p_max = &ws_max[0];
            scalar* p_coordinates = &ws_coordinates[0];
            scalar* p_scs_areav = &ws_scs_areav[0];

            meSCS->shape_fcn(&p_coordinate_shape_function[0]);

            const label* lrscv = meSCS->adjacentNodes();

            for (stk::mesh::Bucket::size_type iElement = 0;
                 iElement < nElementsPerBucket;
                 ++iElement)
            {
                stk::mesh::Entity const* nodeRels =
                    elementBucket.begin_nodes(iElement);
                const label nNodesPerElement =
                    elementBucket.num_nodes(iElement);

                for (label iNode = 0; iNode < nNodesPerElement; ++iNode)
                {
                    stk::mesh::Entity node = nodeRels[iNode];

                    const scalar* phi =
                        stk::mesh::field_data(*phiSTKFieldPtr, node);
                    const scalar* min =
                        stk::mesh::field_data(*minValueSTKFieldPtr, node);
                    const scalar* max =
                        stk::mesh::field_data(*maxValueSTKFieldPtr, node);
                    const scalar* gradPhi =
                        stk::mesh::field_data(*gradPhiSTKFieldPtr, node);
                    const scalar* coords =
                        stk::mesh::field_data(*coordsSTKFieldPtr, node);

                    p_phi[iNode] = phi[0];
                    p_min[iNode] = min[0];
                    p_max[iNode] = max[0];

                    for (label j = 0; j < SPATIAL_DIM; ++j)
                    {
                        p_gradPhi[iNode * SPATIAL_DIM + j] = gradPhi[j];
                    }

                    for (label i = 0; i < SPATIAL_DIM; ++i)
                    {
                        p_coordinates[iNode * SPATIAL_DIM + i] = coords[i];
                    }
                }

                // compute geometry
                scalar scs_error = 0.0;
                meSCS->determinant(
                    1, &p_coordinates[0], &p_scs_areav[0], &scs_error);

                for (label ip = 0; ip < numScsIp; ++ip)
                {
                    const label il = lrscv[2 * ip];
                    const label ir = lrscv[2 * ip + 1];

                    stk::mesh::Entity nodeL = nodeRels[il];
                    stk::mesh::Entity nodeR = nodeRels[ir];

                    // interpolate to scs ip
                    for (label j = 0; j < SPATIAL_DIM; j++)
                    {
                        p_coordIp[j] = 0.0;
                    }
                    const label offset = ip * nodesPerElement;
                    for (label ic = 0; ic < nodesPerElement; ++ic)
                    {
                        const scalar r_coord =
                            p_coordinate_shape_function[offset + ic];
                        for (label j = 0; j < SPATIAL_DIM; ++j)
                        {
                            p_coordIp[j] +=
                                r_coord * p_coordinates[ic * SPATIAL_DIM + j];
                        }
                    }

                    // left
                    {
                        for (label j = 0; j < SPATIAL_DIM; ++j)
                        {
                            dxj[j] = p_coordIp[j] -
                                     p_coordinates[il * SPATIAL_DIM + j];
                        }

                        scalar* beta = stk::mesh::field_data(
                            *blendingFactorSTKFieldPtr, nodeL);

                        scalar rgrad = 0;
                        for (label j = 0; j < SPATIAL_DIM; ++j)
                        {
                            rgrad += p_gradPhi[il * SPATIAL_DIM + j] * dxj[j];
                        }

                        if (rgrad > 0)
                        {
                            rgrad = rgrad + SMALL;
                        }
                        else
                        {
                            rgrad = rgrad - SMALL;
                        }

                        // one-sided Barth-Jespersen bound
                        scalar y = std::max((p_max[il] - p_phi[il]) / rgrad,
                                            (p_min[il] - p_phi[il]) / rgrad);

                        if (blendingFactorMax_ < 2.0 - SMALL)
                        {
                            y = ((y * y + 2.0 * y) / (y * y + y + 2.0));
                        }

                        if ((rgrad > 100 * SMALL) || (rgrad < -100 * SMALL))
                        {
                            beta[0] = std::min(beta[0], y);
                        }
                    }

                    // right
                    {
                        for (label j = 0; j < SPATIAL_DIM; ++j)
                        {
                            p_dxj[j] = p_coordIp[j] -
                                       p_coordinates[ir * SPATIAL_DIM + j];
                        }

                        scalar* beta = stk::mesh::field_data(
                            *blendingFactorSTKFieldPtr, nodeR);

                        scalar rgrad = 0;
                        for (label j = 0; j < SPATIAL_DIM; ++j)
                        {
                            rgrad += p_gradPhi[ir * SPATIAL_DIM + j] * dxj[j];
                        }

                        if (rgrad > 0)
                        {
                            rgrad = rgrad + SMALL;
                        }
                        else
                        {
                            rgrad = rgrad - SMALL;
                        }

                        // one-sided Barth-Jespersen bound
                        scalar y = std::max((p_max[ir] - p_phi[ir]) / rgrad,
                                            (p_min[ir] - p_phi[ir]) / rgrad);

                        if (blendingFactorMax_ < 2.0 - SMALL)
                        {
                            y = ((y * y + 2.0 * y) / (y * y + y + 2.0));
                        }

                        if ((rgrad > 100 * SMALL) || (rgrad < -100 * SMALL))
                        {
                            beta[0] = std::min(beta[0], y);
                        }
                    }
                }
            }
        }
    }

    // Boundary ip
    {
        STKScalarField* coordsSTKFieldPtr = metaData.get_field<scalar>(
            stk::topology::NODE_RANK, this->getCoordinatesID_(iZone));

        std::vector<scalar> ws_phi;
        std::vector<scalar> ws_gradPhi;
        std::vector<scalar> ws_min;
        std::vector<scalar> ws_max;
        std::vector<scalar> ws_coordinates;

        std::vector<scalar> ws_shape_function;
        std::vector<scalar> ws_coordinate_shape_function;

        std::vector<scalar> coordIp(SPATIAL_DIM);
        std::vector<scalar> dxj(SPATIAL_DIM);

        scalar* p_coordIp = &coordIp[0];
        scalar* p_dxj = &dxj[0];

        for (label iBoundary = 0;
             iBoundary < this->meshPtr()->zonePtr(iZone)->nBoundaries();
             iBoundary++)
        {
            std::vector<stk::topology> parentTopo;

            stk::mesh::PartVector partVec =
                this->meshPtr()->zonePtr(iZone)->boundaryRef(iBoundary).parts();

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

                MasterElement* meFC =
                    accel::MasterElementRepo::get_surface_master_element(
                        sideBucket.topology());

                const label nodesPerSide = meFC->nodesPerElement_;
                const label numScsIp = meFC->numIntPoints_;
                const label* ipNodeMap = meFC->ipNodeMap();

                ws_phi.resize(nodesPerSide);
                ws_gradPhi.resize(nodesPerSide * SPATIAL_DIM);
                ws_min.resize(nodesPerSide);
                ws_max.resize(nodesPerSide);
                ws_shape_function.resize(numScsIp * nodesPerSide);
                ws_coordinate_shape_function.resize(numScsIp * nodesPerSide);
                ws_coordinates.resize(nodesPerSide * SPATIAL_DIM);

                scalar* p_phi = &ws_phi[0];
                scalar* p_gradPhi = &ws_gradPhi[0];
                scalar* p_min = &ws_min[0];
                scalar* p_max = &ws_max[0];
                scalar* p_coordinates = &ws_coordinates[0];
                scalar* p_shape_function = &ws_shape_function[0];
                scalar* p_coordinate_shape_function =
                    &ws_coordinate_shape_function[0];

                if (interpolationScheme_ ==
                    interpolationSchemeType::linearLinear)
                {
                    meFC->shifted_shape_fcn(&p_shape_function[0]);
                }
                else
                {
                    meFC->shape_fcn(&p_shape_function[0]);
                }

                meFC->shape_fcn(&p_coordinate_shape_function[0]);

                for (stk::mesh::Bucket::size_type iSide = 0;
                     iSide < nSidesPerBucket;
                     ++iSide)
                {
                    stk::mesh::Entity const* sideNodeRels =
                        sideBucket.begin_nodes(iSide);
                    label numNodes = sideBucket.num_nodes(iSide);

                    STK_ThrowAssert(numNodes == nodesPerSide);

                    for (label ni = 0; ni < numNodes; ++ni)
                    {
                        stk::mesh::Entity node = sideNodeRels[ni];

                        const scalar* phi =
                            stk::mesh::field_data(*phiSTKFieldPtr, node);
                        const scalar* min =
                            stk::mesh::field_data(*minValueSTKFieldPtr, node);
                        const scalar* max =
                            stk::mesh::field_data(*maxValueSTKFieldPtr, node);
                        const scalar* gradPhi =
                            stk::mesh::field_data(*gradPhiSTKFieldPtr, node);
                        const scalar* coords =
                            stk::mesh::field_data(*coordsSTKFieldPtr, node);

                        p_phi[ni] = phi[0];
                        p_min[ni] = min[0];
                        p_max[ni] = max[0];

                        for (label j = 0; j < SPATIAL_DIM; ++j)
                        {
                            p_gradPhi[ni * SPATIAL_DIM + j] = gradPhi[j];
                        }

                        for (label i = 0; i < SPATIAL_DIM; ++i)
                        {
                            p_coordinates[ni * SPATIAL_DIM + i] = coords[i];
                        }
                    }

                    for (label ip = 0; ip < numScsIp; ++ip)
                    {
                        const label nn = ipNodeMap[ip];

                        stk::mesh::Entity node = sideNodeRels[nn];

                        for (label j = 0; j < SPATIAL_DIM; j++)
                        {
                            p_coordIp[j] = 0.0;
                        }
                        const label offset = ip * nodesPerSide;
                        for (label ic = 0; ic < nodesPerSide; ++ic)
                        {
                            const scalar r_coord =
                                p_coordinate_shape_function[offset + ic];
                            for (label j = 0; j < SPATIAL_DIM; ++j)
                            {
                                p_coordIp[j] +=
                                    r_coord *
                                    p_coordinates[ic * SPATIAL_DIM + j];
                            }
                        }

                        for (label j = 0; j < SPATIAL_DIM; ++j)
                        {
                            p_dxj[j] = p_coordIp[j] -
                                       p_coordinates[nn * SPATIAL_DIM + j];
                        }

                        scalar* beta = stk::mesh::field_data(
                            *blendingFactorSTKFieldPtr, node);

                        scalar rgrad = 0;
                        for (label j = 0; j < SPATIAL_DIM; ++j)
                        {
                            rgrad += p_gradPhi[nn * SPATIAL_DIM + j] * p_dxj[j];
                        }

                        if (rgrad > 0)
                        {
                            rgrad = rgrad + SMALL;
                        }
                        else
                        {
                            rgrad = rgrad - SMALL;
                        }

                        // one-sided Barth-Jespersen bound
                        scalar y = std::max((p_max[nn] - p_phi[nn]) / rgrad,
                                            (p_min[nn] - p_phi[nn]) / rgrad);

                        if (blendingFactorMax_ < 2.0 - SMALL)
                        {
                            y = ((y * y + 2.0 * y) / (y * y + y + 2.0));
                        }

                        if ((rgrad > 100 * SMALL) || (rgrad < -100 * SMALL))
                        {
                            beta[0] = std::min(beta[0], y);
                        }
                    }
                }
            }
        }
    }

    // Relax beta field
    stk::mesh::Selector selAllNodes =
        this->metaDataRef().locally_owned_part() &
        stk::mesh::selectUnion(
            this->meshPtr()->zonePtr(iZone)->interiorParts());

    const auto& blendingFactorSTKFieldRefPrev =
        this->blendingFactorRef().prevIterRef().stkFieldRef();

    stk::mesh::BucketVector const& nodeBuckets =
        this->bulkDataRef().get_buckets(stk::topology::NODE_RANK, selAllNodes);
    for (stk::mesh::BucketVector::const_iterator ib = nodeBuckets.begin();
         ib != nodeBuckets.end();
         ++ib)
    {
        stk::mesh::Bucket& nodeBucket = **ib;
        const stk::mesh::Bucket::size_type nNodesPerBucket = nodeBucket.size();

        scalar* beta =
            stk::mesh::field_data(*blendingFactorSTKFieldPtr, nodeBucket);
        const scalar* betaPrev =
            stk::mesh::field_data(blendingFactorSTKFieldRefPrev, nodeBucket);

        for (stk::mesh::Bucket::size_type iNode = 0; iNode < nNodesPerBucket;
             ++iNode)
        {
            beta[iNode] = 0.75 * betaPrev[iNode] + 0.25 * beta[iNode];
        }
    }

    // synchronize in case of parallel
    this->blendingFactorRef().synchronizeGhostedEntities(iZone);
}

} // namespace accel
