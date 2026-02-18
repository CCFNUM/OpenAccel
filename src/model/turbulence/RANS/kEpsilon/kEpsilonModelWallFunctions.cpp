// File : kEpsilonModelWallFunctions.cpp
// Created : Thu Feb 22 2025 13:38:51 (+0100)
// Author : Achraf Nagihi
// Description:
// Copyright (c) 2025 CCFNUM, Lucerne University of Applied Sciences and Arts.
// SPDX-License-Identifier: BSD-3-Clause

#include "kEpsilonModel.h"
#include "messager.h"
#include "realm.h"

namespace accel
{

void kEpsilonModel::updateEpsilonAtWalls(const std::shared_ptr<domain> domain)
{
    stk::mesh::BulkData& bulkData = meshRef().bulkDataRef();
    stk::mesh::MetaData& metaData = meshRef().metaDataRef();

    // Get fields
    STKScalarField* epsilonSTKFieldPtr = epsilonRef().stkFieldPtr();
    STKScalarField* bcEpsilonSTKFieldPtr =
        epsilonRef().nodeSideFieldRef().stkFieldPtr();

    const auto& assembledWallAreaSTKFieldRef = *metaData.get_field<scalar>(
        stk::topology::NODE_RANK, mesh::assembled_wall_area_ID);

    // Zero-out epsilon values at no-slip walls
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
            scalar* epsilon = stk::mesh::field_data(*epsilonSTKFieldPtr, b);
            scalar* bcEpsilon = stk::mesh::field_data(*bcEpsilonSTKFieldPtr, b);
            for (stk::mesh::Bucket::size_type k = 0; k < length; ++k)
            {
                epsilon[k] = 0.0;
                bcEpsilon[k] = 0.0;
            }
        }
    }

    // process the algorithm
    {
        // Get fields
        STKScalarField* exposedAreaSTKFieldPtr = metaData.get_field<scalar>(
            metaData.side_rank(), this->getExposedAreaVectorID_(domain));
        STKScalarField* wallNormalDistanceSTKFieldPtr =
            metaData.get_field<scalar>(metaData.side_rank(),
                                       mesh::wall_normal_distance_ID);

        STKScalarField* uStarSTKFieldPtr = uStarRef().stkFieldPtr();
        STKScalarField* rhoSTKFieldPtr = rhoRef().stkFieldPtr();
        STKScalarField* muSTKFieldPtr = muRef().stkFieldPtr();
        STKScalarField* mutSTKFieldPtr = mutRef().stkFieldPtr();
        STKScalarField& kSTKFieldRef = kRef().stkFieldRef();

        scalar Cmu = this->Cmu();
        scalar kappa = this->kappa();

        // nodal fields to gather
        std::vector<scalar> ws_rho;
        std::vector<scalar> ws_mu;
        std::vector<scalar> ws_mut;
        std::vector<scalar> ws_tke;

        // master element
        std::vector<scalar> ws_face_shape_function;

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

            // face master element
            MasterElement* meFC = MasterElementRepo::get_surface_master_element(
                sideBucket.topology());
            const label nodesPerSide = meFC->nodesPerElement_;
            const label numScsBip = meFC->numIntPoints_;

            // mapping from ip to nodes for this ordinal; face perspective
            // (use with face_node_relations)
            const label* faceIpNodeMap = meFC->ipNodeMap();

            // algorithm related; element
            ws_rho.resize(nodesPerSide);
            ws_mu.resize(nodesPerSide);
            ws_mut.resize(nodesPerSide);
            ws_tke.resize(nodesPerSide);
            ws_face_shape_function.resize(numScsBip * nodesPerSide);

            // pointers
            scalar* p_rho = &ws_rho[0];
            scalar* p_mu = &ws_mu[0];
            scalar* p_mut = &ws_mut[0];
            scalar* p_tke = &ws_tke[0];
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
                    p_mut[ni] = *stk::mesh::field_data(*mutSTKFieldPtr, node);
                    p_tke[ni] = *stk::mesh::field_data(kSTKFieldRef, node);
                }

                // pointer to face data
                const scalar* areaVec =
                    stk::mesh::field_data(*exposedAreaSTKFieldPtr, side);

                const scalar* wallNormalDistanceBip =
                    stk::mesh::field_data(*wallNormalDistanceSTKFieldPtr, side);
                const scalar* uStarBip =
                    stk::mesh::field_data(*uStarSTKFieldPtr, side);

                // loop over face nodes
                for (label ip = 0; ip < numSideNodes; ++ip)
                {
                    const label offSetSF_face = ip * nodesPerSide;
                    const label offSetAveraVec = ip * SPATIAL_DIM;

                    const label localFaceNode = faceIpNodeMap[ip];
                    stk::mesh::Entity node = sideNodeRels[localFaceNode];

                    // interpolate to bip
                    scalar muBip = 0.0;
                    scalar rhoBip = 0.0;
                    scalar mutBip = 0.0;
                    scalar tkeBip = 0.0;

                    for (label ic = 0; ic < nodesPerSide; ++ic)
                    {
                        const scalar r =
                            p_face_shape_function[offSetSF_face + ic];
                        rhoBip += r * p_rho[ic];
                        muBip += r * p_mu[ic];
                        mutBip += r * p_mut[ic];
                        tkeBip += r * p_tke[ic];
                    }

                    // zero out vector quantities; squeeze in aMag
                    scalar aMag = 0.0;
                    for (label j = 0; j < SPATIAL_DIM; ++j)
                    {
                        const scalar axj = areaVec[offSetAveraVec + j];
                        aMag += axj * axj;
                    }
                    aMag = std::sqrt(aMag);

                    scalar* epsilon =
                        stk::mesh::field_data(*epsilonSTKFieldPtr, node);
                    scalar* bcEpsilon =
                        stk::mesh::field_data(*bcEpsilonSTKFieldPtr, node);
                    scalar yStar = (rhoBip * uStarBip[ip] *
                                    wallNormalDistanceBip[ip] / 4.0) /
                                   muBip;

                    scalar yStarTilde = std::max(yStar, 11.06);
                    scalar epsilon_log =
                        (rhoBip * uStarBip[ip] / (yStarTilde * muBip)) *
                        (std::pow(Cmu, 0.75) / kappa) * std::pow(tkeBip, 1.5);

                    (*epsilon) += epsilon_log * aMag;
                    (*bcEpsilon) += epsilon_log * aMag;
                }
            }
        }
    }

    // normalize and set assembled tdr to tdr bc
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
            const scalar* area =
                stk::mesh::field_data(assembledWallAreaSTKFieldRef, b);
            scalar* epsilon = stk::mesh::field_data(*epsilonSTKFieldPtr, b);
            scalar* bcEpsilon = stk::mesh::field_data(*bcEpsilonSTKFieldPtr, b);
            for (stk::mesh::Bucket::size_type k = 0; k < length; ++k)
            {
                const scalar tdrBnd = epsilon[k] / area[k];

                epsilon[k] = tdrBnd;
                bcEpsilon[k] = tdrBnd;
            }
        }
    }
}

} /* namespace accel */
