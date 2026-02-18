// File : shearStressTransportModelWallFunctions.cpp
// Created : Fri Jun 20 2024 16:48:19 (+0100)
// Author : Mhamad Mahdi Alloush
// Description:
// Copyright (c) 2024 CCFNUM, Lucerne University of Applied Sciences and Arts.
// SPDX-License-Identifier: BSD-3-Clause

#include "messager.h"
#include "realm.h"
#include "shearStressTransportModel.h"

namespace accel
{

void shearStressTransportModel::updateTurbulentEddyFrequencyAtWalls(
    const std::shared_ptr<domain> domain)
{
    stk::mesh::BulkData& bulkData = meshRef().bulkDataRef();
    stk::mesh::MetaData& metaData = meshRef().metaDataRef();

    // Get fields
    STKScalarField* omegaSTKFieldPtr = omegaRef().stkFieldPtr();
    STKScalarField* bcOmegaSTKFieldPtr =
        omegaRef().nodeSideFieldRef().stkFieldPtr();

    // Get geometric fields
    const auto& assembledWallAreaSTKFieldRef = *metaData.get_field<scalar>(
        stk::topology::NODE_RANK, mesh::assembled_wall_area_ID);

    // zero omega field at no-slip walls
    {
        stk::mesh::Selector selAllNodes =
            metaData.universal_part() &
            stk::mesh::selectUnion(collectNoSlipWallParts_(domain));

        stk::mesh::BucketVector const& nodeBuckets =
            bulkData.get_buckets(stk::topology::NODE_RANK, selAllNodes);

        for (stk::mesh::BucketVector::const_iterator ib = nodeBuckets.begin();
             ib != nodeBuckets.end();
             ++ib)
        {
            stk::mesh::Bucket& b = **ib;
            const stk::mesh::Bucket::size_type length = b.size();
            scalar* omegab = stk::mesh::field_data(*omegaSTKFieldPtr, b);
            scalar* bcOmegab = stk::mesh::field_data(*bcOmegaSTKFieldPtr, b);
            for (stk::mesh::Bucket::size_type k = 0; k < length; ++k)
            {
                omegab[k] = 0.0;
                bcOmegab[k] = 0.0;
            }
        }
    }

    // process the algorithm
    {
        // Get fields
        const STKScalarField* exposedAreaSTKFieldPtr =
            metaData.get_field<scalar>(metaData.side_rank(),
                                       this->getExposedAreaVectorID_(domain));
        const STKScalarField* wallNormalDistanceSTKFieldPtr =
            metaData.get_field<scalar>(metaData.side_rank(),
                                       mesh::wall_normal_distance_ID);

        const STKScalarField* yPlusSTKFieldPtr = yPlusRef().stkFieldPtr();
        const STKScalarField* uStarSTKFieldPtr = uStarRef().stkFieldPtr();

        const STKScalarField* rhoSTKFieldPtr = rhoRef().stkFieldPtr();
        const STKScalarField* muSTKFieldPtr = muRef().stkFieldPtr();

        // define vector of parent topos; should always
        // be UNITY in size
        std::vector<stk::topology> parentTopo;

        // nodal fields to gather
        std::vector<scalar> ws_rho;
        std::vector<scalar> ws_mu;

        // master element
        std::vector<scalar> ws_face_shape_function;

        // others
        const scalar Cmu = this->Cmu();

        // select faces: only those sitting on no-slip
        // walls
        stk::mesh::Selector selAllSides =
            metaData.universal_part() &
            stk::mesh::selectUnion(collectNoSlipWallParts_(domain));

        // shifted ip's for fields?
        const bool isUShifted = URef().isShifted();

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

            // face master element
            MasterElement* meFC = MasterElementRepo::get_surface_master_element(
                sideBucket.topology());
            const label nodesPerSide = meFC->nodesPerElement_;
            const label numScsBip = meFC->numIntPoints_;

            // mapping from ip to nodes for this
            // ordinal; face perspective (use with
            // face_node_relations)
            const label* faceIpNodeMap = meFC->ipNodeMap();

            // algorithm related; element
            ws_rho.resize(nodesPerSide);
            ws_mu.resize(nodesPerSide);
            ws_face_shape_function.resize(numScsBip * nodesPerSide);

            // pointers
            scalar* p_rho = &ws_rho[0];
            scalar* p_mu = &ws_mu[0];
            scalar* p_face_shape_function = &ws_face_shape_function[0];

            // shape functions
            if (isUShifted)
            {
                meFC->shifted_shape_fcn(&p_face_shape_function[0]);
            }
            else
            {
                meFC->shape_fcn(&p_face_shape_function[0]);
            }

            const stk::mesh::Bucket::size_type length = sideBucket.size();

            for (stk::mesh::Bucket::size_type k = 0; k < length; ++k)
            {
                // get face
                stk::mesh::Entity side = sideBucket[k];

                //======================================
                // gather nodal data off of face; n/a
                //======================================
                label numSideNodes = bulkData.num_nodes(side);

                //======================================
                // gather nodal data off of face
                //======================================
                stk::mesh::Entity const* sideNodeRels =
                    bulkData.begin_nodes(side);

                // sanity check on num nodes
                STK_ThrowAssert(numSideNodes == nodesPerSide);
                for (label ni = 0; ni < nodesPerSide; ++ni)
                {
                    stk::mesh::Entity node = sideNodeRels[ni];

                    // gather scalars
                    p_mu[ni] = *stk::mesh::field_data(*muSTKFieldPtr, node);
                    p_rho[ni] = *stk::mesh::field_data(*rhoSTKFieldPtr, node);
                }

                // pointer to face data
                const scalar* areaVec =
                    stk::mesh::field_data(*exposedAreaSTKFieldPtr, side);

                const scalar* wallNormalDistanceBip =
                    stk::mesh::field_data(*wallNormalDistanceSTKFieldPtr, side);
                const scalar* uStarBip =
                    stk::mesh::field_data(*uStarSTKFieldPtr, side);
                const scalar* yPlusBip =
                    stk::mesh::field_data(*yPlusSTKFieldPtr, side);

                // loop over face nodes
                for (label ip = 0; ip < numSideNodes; ++ip)
                {
                    const label offSetSF_face = ip * nodesPerSide;
                    const label offSetAveraVec = ip * SPATIAL_DIM;

                    const label nearestNode = faceIpNodeMap[ip];
                    stk::mesh::Entity node = sideNodeRels[nearestNode];

                    // interpolate to bip
                    scalar muBip = 0.0;
                    scalar rhoBip = 0.0;
                    for (label ic = 0; ic < nodesPerSide; ++ic)
                    {
                        const scalar r =
                            p_face_shape_function[offSetSF_face + ic];
                        rhoBip += r * p_rho[ic];
                        muBip += r * p_mu[ic];
                    }

                    // zero out vector quantities;
                    // squeeze in aMag
                    scalar aMag = 0.0;
                    for (label j = 0; j < SPATIAL_DIM; ++j)
                    {
                        const scalar axj = areaVec[offSetAveraVec + j];
                        aMag += axj * axj;
                    }
                    aMag = std::sqrt(aMag);

                    scalar* omega =
                        stk::mesh::field_data(*omegaSTKFieldPtr, node);
                    scalar* bcOmega =
                        stk::mesh::field_data(*bcOmegaSTKFieldPtr, node);

                    scalar Cmu50 = std::pow(Cmu, 0.5);

                    scalar omegaVis = scalar(6.0) * (muBip / rhoBip) /
                                      (betaOne_ * wallNormalDistanceBip[ip] *
                                       wallNormalDistanceBip[ip]);
                    scalar omegaLog = std::pow(uStarBip[ip], 2.0) * rhoBip /
                                      (kappa() * muBip * yPlusBip[ip] * Cmu50);

                    scalar omegaWall = std::sqrt(std::pow(omegaVis, 2.0) +
                                                 std::pow(omegaLog, 2.0));

                    (*omega) += omegaWall * aMag;
                    (*bcOmega) += omegaWall * aMag;
                }
            }
        }
    }

    // normalize and set assembled tef to tef bc
    {
        stk::mesh::Selector selAllNodes =
            metaData.universal_part() &
            stk::mesh::selectUnion(collectNoSlipWallParts_(domain));

        stk::mesh::BucketVector const& nodeBuckets =
            bulkData.get_buckets(stk::topology::NODE_RANK, selAllNodes);

        for (stk::mesh::BucketVector::const_iterator ib = nodeBuckets.begin();
             ib != nodeBuckets.end();
             ++ib)
        {
            stk::mesh::Bucket& b = **ib;
            const stk::mesh::Bucket::size_type length = b.size();
            const scalar* areab =
                stk::mesh::field_data(assembledWallAreaSTKFieldRef, b);
            scalar* omegab = stk::mesh::field_data(*omegaSTKFieldPtr, b);
            scalar* bcOmegab = stk::mesh::field_data(*bcOmegaSTKFieldPtr, b);
            for (stk::mesh::Bucket::size_type k = 0; k < length; ++k)
            {
                const scalar tefBnd = omegab[k] / areab[k];

                omegab[k] = tefBnd;
                bcOmegab[k] = tefBnd;
            }
        }
    }
}

} /* namespace accel */
