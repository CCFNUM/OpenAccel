// File : kEpsilonModel.cpp
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

kEpsilonModel::kEpsilonModel(realm* realm) : RANSModel(realm)
{
    // Create field instances:

    // transport
    yScaleRef();
    kRef();
    epsilonRef();

    // auxiliary
    yMinRef();
    PkRef();
}

void kEpsilonModel::clipValues(const std::shared_ptr<domain> domain)
{
    const scalar clipValue = 1.0e-8;

    stk::mesh::BulkData& bulkData = this->meshRef().bulkDataRef();
    stk::mesh::MetaData& metaData = this->meshRef().metaDataRef();

    // required fields
    STKScalarField* rhoSTKFieldRef = this->rhoRef().stkFieldPtr();
    STKScalarField* muSTKFieldRef = this->muRef().stkFieldPtr();

    // required fields with state
    STKScalarField& kSTKFieldRef = this->kRef().stkFieldRef();
    STKScalarField& epsilonSTKFieldRef = this->epsilonRef().stkFieldRef();

    scalar Cmu = this->Cmu();

    // define some common selectors
    stk::mesh::Selector selAllNodes =
        metaData.universal_part() &
        stk::mesh::selectUnion(domain->zonePtr()->interiorParts());

    stk::mesh::BucketVector const& nodeBuckets =
        bulkData.get_buckets(stk::topology::NODE_RANK, selAllNodes);

    for (stk::mesh::BucketVector::const_iterator ib = nodeBuckets.begin();
         ib != nodeBuckets.end();
         ++ib)
    {
        stk::mesh::Bucket& b = **ib;
        const stk::mesh::Bucket::size_type length = b.size();

        const scalar* mu = stk::mesh::field_data(*muSTKFieldRef, b);
        const scalar* rho = stk::mesh::field_data(*rhoSTKFieldRef, b);
        scalar* tke = stk::mesh::field_data(kSTKFieldRef, b);
        scalar* tdr = stk::mesh::field_data(epsilonSTKFieldRef, b);

        for (stk::mesh::Bucket::size_type k = 0; k < length; ++k)
        {
            scalar tkeNew = tke[k];
            scalar tdrNew = tdr[k];

            if ((tkeNew >= 0.0) && (tdrNew > 0.0))
            {
                // nothing
            }
            else if ((tkeNew < 0.0) && (tdrNew < 0.0))
            {
                // both negative;
                tkeNew = clipValue;
                tdrNew = rho[k] * clipValue * clipValue * Cmu / mu[k];
            }
            else if (tkeNew < 0.0)
            {
                tkeNew = std::sqrt(mu[k] * tdrNew / (rho[k] * Cmu));
            }
            else
            {
                tdrNew = rho[k] * tke[k] * tke[k] * Cmu / mu[k];
            }

            tke[k] = tkeNew;
            tdr[k] = tdrNew;
        }
    }
}

void kEpsilonModel::updateTurbulentProduction(
    const std::shared_ptr<domain> domain)
{
    stk::mesh::BulkData& bulkData = this->meshRef().bulkDataRef();
    stk::mesh::MetaData& metaData = this->meshRef().metaDataRef();

    // Get fields
    STKScalarField& PkSTKFieldRef = PkRef().stkFieldRef();
    STKScalarField* rhoSTKFieldPtr = rhoRef().stkFieldPtr();
    STKScalarField* muSTKFieldPtr = muRef().stkFieldPtr();

    // Compute tke production
    {
        STKScalarField* mutSTKFieldPtr = mutRef().stkFieldPtr();
        STKScalarField* gradUSTKFieldPtr = URef().gradRef().stkFieldPtr();

        // select all nodes
        stk::mesh::Selector selAllNodes =
            metaData.universal_part() &
            stk::mesh::selectUnion(domain->zonePtr()->interiorParts());

        stk::mesh::BucketVector const& nodeBuckets =
            bulkData.get_buckets(stk::topology::NODE_RANK, selAllNodes);
        for (stk::mesh::BucketVector::const_iterator ib = nodeBuckets.begin();
             ib != nodeBuckets.end();
             ++ib)
        {
            stk::mesh::Bucket& b = **ib;
            const stk::mesh::Bucket::size_type length = b.size();

            // fields;
            const scalar* mut = stk::mesh::field_data(*mutSTKFieldPtr, b);
            const scalar* gradU = stk::mesh::field_data(*gradUSTKFieldPtr, b);
            scalar* Pk = stk::mesh::field_data(PkSTKFieldRef, b);

            for (stk::mesh::Bucket::size_type k = 0; k < length; ++k)
            {
                Pk[k] = 0.0;
                for (label i = 0; i < SPATIAL_DIM; ++i)
                {
                    const label offSet = SPATIAL_DIM * i;
                    for (label j = 0; j < SPATIAL_DIM; ++j)
                    {
                        Pk[k] +=
                            gradU[SPATIAL_DIM * SPATIAL_DIM * k + offSet + j] *
                            (gradU[SPATIAL_DIM * SPATIAL_DIM * k + offSet + j] +
                             gradU[SPATIAL_DIM * SPATIAL_DIM * k +
                                   SPATIAL_DIM * j + i]);
                    }
                }
                Pk[k] *= mut[k];
            }
        }
    }

    // Zero-out Pk values at no-slip wall nodes
    {
        // sel nodes: only those sitting on no-slip
        // walls
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
            for (stk::mesh::Bucket::size_type k = 0; k < length; ++k)
            {
                // get node
                stk::mesh::Entity node = b[k];

                scalar* Pk = stk::mesh::field_data(PkSTKFieldRef, node);
                (*Pk) = 0.0;
            }
        }
    }

    // Set Pk values at no-slip wall nodes
    {
        const auto& USTKFieldRef = URef().stkFieldRef();
        const auto& sideUSTKFieldRef = URef().sideFieldRef().stkFieldRef();

        STKScalarField* uPlusSTKFieldPtr = uPlusRef().stkFieldPtr();
        STKScalarField* uStarSTKFieldPtr = uStarRef().stkFieldPtr();
        STKScalarField* duPlusdyPlusSTKFieldPtr =
            duPlusdyPlusRef().stkFieldPtr();

        const auto& coordsSTKFieldRef = *metaData.get_field<scalar>(
            stk::topology::NODE_RANK, this->getCoordinatesID_(domain));
        const auto& exposedAreaVecSTKFieldRef = *metaData.get_field<scalar>(
            metaData.side_rank(), this->getExposedAreaVectorID_(domain));

        // define vector of parent topos; should always
        // be UNITY in size
        std::vector<stk::topology> parentTopo;

        // bip values
        std::vector<scalar> uBip(SPATIAL_DIM);
        std::vector<scalar> nx(SPATIAL_DIM);

        // pointers to fixed values
        scalar* p_uBip = &uBip[0];
        scalar* p_nx = &nx[0];

        // nodal fields to gather
        std::vector<scalar> ws_U;
        std::vector<scalar> ws_rho;
        std::vector<scalar> ws_mu;

        // master element
        std::vector<scalar> ws_face_shape_function;

        // select all sides: only those sitting on
        // no-slip walls
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

            // extract master element
            MasterElement* meSCS =
                MasterElementRepo::get_surface_master_element(theElemTopo);

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
            ws_U.resize(nodesPerSide * SPATIAL_DIM);
            ws_rho.resize(nodesPerSide);
            ws_mu.resize(nodesPerSide);
            ws_face_shape_function.resize(numScsBip * nodesPerSide);

            // pointers
            scalar* p_U = &ws_U[0];
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

            const stk::mesh::Bucket::size_type nSidesPerBucket =
                sideBucket.size();

            for (stk::mesh::Bucket::size_type iSide = 0;
                 iSide < nSidesPerBucket;
                 ++iSide)
            {
                // get face
                stk::mesh::Entity face = sideBucket[iSide];

                //======================================
                // gather nodal data off of face
                //======================================
                stk::mesh::Entity const* sideNodeRels =
                    bulkData.begin_nodes(face);
                label numSideNodes = bulkData.num_nodes(face);

                // sanity check on num nodes
                STK_ThrowAssert(numSideNodes == nodesPerSide);
                for (label ni = 0; ni < nodesPerSide; ++ni)
                {
                    stk::mesh::Entity node = sideNodeRels[ni];

                    // gather scalars
                    p_rho[ni] = *stk::mesh::field_data(*rhoSTKFieldPtr, node);
                    p_mu[ni] = *stk::mesh::field_data(*muSTKFieldPtr, node);

                    // gather vectors
                    scalar* U = stk::mesh::field_data(USTKFieldRef, node);
                    const label offSet = ni * SPATIAL_DIM;
                    for (label j = 0; j < SPATIAL_DIM; ++j)
                    {
                        p_U[offSet + j] = U[j];
                    }
                }

                // pointer to face data
                const scalar* areaVec =
                    stk::mesh::field_data(exposedAreaVecSTKFieldRef, face);
                const scalar* uStarBip =
                    stk::mesh::field_data(*uStarSTKFieldPtr, face);
                const scalar* uPlusBip =
                    stk::mesh::field_data(*uPlusSTKFieldPtr, face);
                const scalar* UbcVec =
                    stk::mesh::field_data(sideUSTKFieldRef, face);
                const scalar* duPlusdyPlusBip =
                    stk::mesh::field_data(*duPlusdyPlusSTKFieldPtr, face);

                // extract the connected element to this
                // exposed face; should be single in
                // size!
                const stk::mesh::Entity* face_elem_rels =
                    bulkData.begin_elements(face);
                STK_ThrowAssert(bulkData.num_elements(face) == 1);

                // get element; its face ordinal number
                stk::mesh::Entity element = face_elem_rels[0];
                const label face_ordinal =
                    bulkData.begin_element_ordinals(face)[0];

                // get the relations off of element
                stk::mesh::Entity const* elemNodeRels =
                    bulkData.begin_nodes(element);

                // loop over face nodes
                for (label ip = 0; ip < numSideNodes; ++ip)
                {
                    const label offSetSF_face = ip * nodesPerSide;
                    const label offSetAreaVec = ip * SPATIAL_DIM;

                    const label opposingNode =
                        meSCS->opposingNodes(face_ordinal, ip);
                    const label localFaceNode = faceIpNodeMap[ip];

                    // left and right nodes; right is on
                    // the face; left is the opposing
                    // node
                    stk::mesh::Entity nodeL = elemNodeRels[opposingNode];
                    stk::mesh::Entity nodeR = sideNodeRels[localFaceNode];

                    // zero out vector quantities;
                    // squeeze in aMag
                    scalar aMag = 0.0;
                    for (label j = 0; j < SPATIAL_DIM; ++j)
                    {
                        p_uBip[j] = 0.0;
                        const scalar axj = areaVec[offSetAreaVec + j];
                        aMag += axj * axj;
                    }
                    aMag = std::sqrt(aMag);

                    // interpolate to bip
                    scalar muBip = 0.0;
                    scalar rhoBip = 0.0;
                    for (label ic = 0; ic < nodesPerSide; ++ic)
                    {
                        const scalar r =
                            p_face_shape_function[offSetSF_face + ic];
                        rhoBip += r * p_rho[ic];
                        muBip += r * p_mu[ic];

                        const label offSetFN = ic * SPATIAL_DIM;
                        for (label j = 0; j < SPATIAL_DIM; ++j)
                        {
                            p_uBip[j] += r * p_U[offSetFN + j];
                        }
                    }

                    // form unit normal
                    for (label j = 0; j < SPATIAL_DIM; ++j)
                    {
                        const scalar nj = areaVec[offSetAreaVec + j] / aMag;
                        p_nx[j] = nj;
                    }

                    // determine tangential velocity
                    scalar uTangential = 0.0;
                    for (label i = 0; i < SPATIAL_DIM; ++i)
                    {
                        scalar uiTan = 0.0;
                        scalar uiBcTan = 0.0;
                        for (label j = 0; j < SPATIAL_DIM; ++j)
                        {
                            const scalar ninj = p_nx[i] * p_nx[j];
                            if (i == j)
                            {
                                const scalar om_nini = 1.0 - ninj;
                                uiTan += om_nini * p_uBip[j];
                                uiBcTan +=
                                    om_nini * UbcVec[ip * SPATIAL_DIM + j];
                            }
                            else
                            {
                                uiTan -= ninj * p_uBip[j];
                                uiBcTan -= ninj * UbcVec[ip * SPATIAL_DIM + j];
                            }
                        }
                        uTangential += (uiTan - uiBcTan) * (uiTan - uiBcTan);
                    }
                    uTangential = std::sqrt(uTangential);

                    scalar Pk_val =
                        std::pow(rhoBip, scalar(2.0)) * duPlusdyPlusBip[ip] *
                        std::pow(uStarBip[ip] / uPlusBip[ip], scalar(2.0)) *
                        pow(uTangential, scalar(2.0)) / (muBip + SMALL);

                    scalar* Pk = stk::mesh::field_data(PkSTKFieldRef, nodeR);
                    (*Pk) += Pk_val * aMag;
                }
            }
        }
    }

    // Normalize Pk values at no-slip wall nodes
    {
        const auto& assembledWallAreaSTKFieldRef = *metaData.get_field<scalar>(
            stk::topology::NODE_RANK, mesh::assembled_wall_area_ID);

        // sel nodes: only those sitting on no-slip
        // walls
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
            for (stk::mesh::Bucket::size_type k = 0; k < length; ++k)
            {
                // get node
                stk::mesh::Entity node = b[k];

                const scalar* area =
                    stk::mesh::field_data(assembledWallAreaSTKFieldRef, node);
                scalar* Pk = stk::mesh::field_data(PkSTKFieldRef, node);

                (*Pk) /= (*area);
            }
        }
    }
}

void kEpsilonModel::updateTurbulentDynamicViscosity(
    const std::shared_ptr<domain> domain)
{
    // get common data
    stk::mesh::MetaData& metaData = this->meshRef().metaDataRef();
    stk::mesh::BulkData& bulkData = this->meshRef().bulkDataRef();

    // get fields
    STKScalarField* mutSTKFieldPtr = this->mutRef().stkFieldPtr();
    const STKScalarField* rhoSTKFieldPtr = this->rhoRef().stkFieldPtr();
    const STKScalarField& kSTKFieldRef = this->kRef().stkFieldRef();
    const STKScalarField& epsilonSTKFieldRef = this->epsilonRef().stkFieldRef();

    // other
    scalar Cmu = this->Cmu();

    // define some common selectors
    stk::mesh::Selector selAllNodes =
        metaData.universal_part() &
        stk::mesh::selectUnion(domain->zonePtr()->interiorParts());

    stk::mesh::BucketVector const& nodeBuckets =
        bulkData.get_buckets(stk::topology::NODE_RANK, selAllNodes);
    for (stk::mesh::BucketVector::const_iterator ib = nodeBuckets.begin();
         ib != nodeBuckets.end();
         ++ib)
    {
        stk::mesh::Bucket& b = **ib;
        const stk::mesh::Bucket::size_type length = b.size();

        const scalar* rho = stk::mesh::field_data(*rhoSTKFieldPtr, b);
        const scalar* tke = stk::mesh::field_data(kSTKFieldRef, b);
        const scalar* tdr = stk::mesh::field_data(epsilonSTKFieldRef, b);
        scalar* mut = stk::mesh::field_data(*mutSTKFieldPtr, b);

        for (stk::mesh::Bucket::size_type k = 0; k < length; ++k)
        {
            scalar mut_val = Cmu * rho[k] * tke[k] * tke[k] / (tdr[k] + SMALL);

            mut[k] = 0.75 * mut_val + 0.25 * mut[k];
        }
    }
}

void kEpsilonModel::clipMinDistToWall(const std::shared_ptr<domain> domain)
{
    STKScalarField* minDistanceToWallSTKFieldPtr = yMinRef().stkFieldPtr();

    stk::mesh::BulkData& bulkData = meshRef().bulkDataRef();
    stk::mesh::MetaData& metaData = meshRef().metaDataRef();

    // extract fields required
    const STKScalarField* wallNormalDistanceSTKFieldPtr =
        metaData.get_field<scalar>(metaData.side_rank(),
                                   mesh::wall_normal_distance_ID);

    // Process for no-slip wall parts
    stk::mesh::Selector selAllSides =
        metaData.universal_part() &
        stk::mesh::selectUnion(collectNoSlipWallParts_(domain));

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

        // mapping from ip to nodes for this
        // ordinal; face perspective (use with
        // face_node_relations)
        const label* faceIpNodeMap = meFC->ipNodeMap();

        const stk::mesh::Bucket::size_type length = sideBucket.size();

        for (stk::mesh::Bucket::size_type k = 0; k < length; ++k)
        {
            // get face
            stk::mesh::Entity side = sideBucket[k];
            stk::mesh::Entity const* sideNodeRels = bulkData.begin_nodes(side);

            label numSideNodes = bulkData.num_nodes(side);

            // pointer to face data
            const scalar* wallNormalDistanceBip =
                stk::mesh::field_data(*wallNormalDistanceSTKFieldPtr, side);

            // loop over face nodes
            for (label ip = 0; ip < numSideNodes; ++ip)
            {
                const label localFaceNode = faceIpNodeMap[ip];

                // right is on the face
                stk::mesh::Entity node = sideNodeRels[localFaceNode];

                // assemble to nodal quantities
                scalar* minD =
                    stk::mesh::field_data(*minDistanceToWallSTKFieldPtr, node);
                (*minD) = std::max(*minD, wallNormalDistanceBip[ip] / 4.0);
            }
        }
    }
    // parallel reduce
    std::vector<const stk::mesh::FieldBase*> fieldVec;
    fieldVec.push_back(minDistanceToWallSTKFieldPtr);
    stk::mesh::parallel_max(bulkData, fieldVec);
}

} /* namespace accel */
