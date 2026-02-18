// File : shearStressTransportModel.cpp
// Created : Mon Mar 25 2024 16:48:19 (+0100)
// Author : Mhamad Mahdi Alloush
// Description:
// Copyright (c) 2024 CCFNUM, Lucerne University of Applied Sciences and Arts.
// SPDX-License-Identifier: BSD-3-Clause

#include "shearStressTransportModel.h"
#include "messager.h"
#include "realm.h"

namespace accel
{

shearStressTransportModel::shearStressTransportModel(realm* realm)
    : RANSModel(realm)
{
    // Create field instances:

    // transport
    yScaleRef();
    kRef();
    omegaRef();

    // auxiliary
    yMinRef();
    F1Ref();
    PkRef();
}

void shearStressTransportModel::clipValues(const std::shared_ptr<domain> domain)
{
    const scalar clipValue = 1.0e-8;

    stk::mesh::BulkData& bulkData = this->meshRef().bulkDataRef();
    stk::mesh::MetaData& metaData = this->meshRef().metaDataRef();

    // required fields
    const STKScalarField* rhoSTKFieldPtr = this->rhoRef().stkFieldPtr();
    const STKScalarField* muSTKFieldPtr = this->muRef().stkFieldPtr();

    // required fields with state
    STKScalarField& kSTKFieldRef = this->kRef().stkFieldRef();
    STKScalarField& omegaSTKFieldRef = this->omegaRef().stkFieldRef();

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

        const scalar* mub = stk::mesh::field_data(*muSTKFieldPtr, b);
        const scalar* rhob = stk::mesh::field_data(*rhoSTKFieldPtr, b);
        scalar* tkeb = stk::mesh::field_data(kSTKFieldRef, b);
        scalar* tefb = stk::mesh::field_data(omegaSTKFieldRef, b);

        for (stk::mesh::Bucket::size_type k = 0; k < length; ++k)
        {
            scalar tkeNew = tkeb[k];
            scalar tefNew = tefb[k];

            if ((tkeNew >= 0.0) && (tefNew > 0.0))
            {
                // nothing
            }
            else if ((tkeNew < 0.0) && (tefNew < 0.0))
            {
                // both negative;
                tkeNew = clipValue;
                tefNew = rhob[k] * clipValue / mub[k];
            }
            else if (tkeNew < 0.0)
            {
                tkeNew = mub[k] * tefNew / rhob[k];
            }
            else
            {
                tefNew = rhob[k] * tkeb[k] / mub[k];
            }

            tkeb[k] = tkeNew;
            tefb[k] = tefNew;
        }
    }
}

void shearStressTransportModel::updateFOneBlending(
    const std::shared_ptr<domain> domain)
{
    // compute fone with parameters appropriate for 2003 SST implementation
    stk::mesh::BulkData& bulkData = this->meshRef().bulkDataRef();
    stk::mesh::MetaData& metaData = this->meshRef().metaDataRef();

    const scalar CDkwClip = 1.0e-10; // 2003 SST
    const scalar sigmaWTwo = this->sigmaWTwo();
    const scalar betaStar = this->betaStar();

    // required fields with state; min_distance is fine
    STKScalarField& kSTKFieldRef = kRef().stkFieldRef();
    STKScalarField& omegaSTKFieldRef = omegaRef().stkFieldRef();

    // fields not saved off
    const STKScalarField* rhoSTKFieldPtr = rhoRef().stkFieldPtr();
    const STKScalarField* muSTKFieldPtr = muRef().stkFieldPtr();
    const STKScalarField* dkdxSTKFieldPtr = kRef().gradRef().stkFieldPtr();
    const STKScalarField* dwdxSTKFieldPtr = omegaRef().gradRef().stkFieldPtr();

    const STKScalarField* minDistanceToWallSTKFieldPtr =
        yMinRef().stkFieldPtr();
    STKScalarField* F1STKFieldPtr = F1Ref().stkFieldPtr();

    stk::mesh::Selector selAllNodes =
        metaData.universal_part() &
        stk::mesh::selectUnion(domain->zonePtr()->interiorParts());

    stk::mesh::BucketVector const& node_buckets =
        bulkData.get_buckets(stk::topology::NODE_RANK, selAllNodes);
    for (stk::mesh::BucketVector::const_iterator ib = node_buckets.begin();
         ib != node_buckets.end();
         ++ib)
    {
        stk::mesh::Bucket& b = **ib;
        const stk::mesh::Bucket::size_type length = b.size();

        // fields; supplemental and non-const fOne and ftwo
        const scalar* tefb = stk::mesh::field_data(omegaSTKFieldRef, b);
        const scalar* tkeb = stk::mesh::field_data(kSTKFieldRef, b);
        const scalar* yb =
            stk::mesh::field_data(*minDistanceToWallSTKFieldPtr, b);
        const scalar* rhob = stk::mesh::field_data(*rhoSTKFieldPtr, b);
        const scalar* mub = stk::mesh::field_data(*muSTKFieldPtr, b);
        const scalar* dkdxb = stk::mesh::field_data(*dkdxSTKFieldPtr, b);
        const scalar* dwdxb = stk::mesh::field_data(*dwdxSTKFieldPtr, b);
        scalar* fOneb = stk::mesh::field_data(*F1STKFieldPtr, b);

        for (stk::mesh::Bucket::size_type k = 0; k < length; ++k)
        {
            // compute cross diff
            scalar crossDiff = 0.0;
            for (label j = 0; j < SPATIAL_DIM; ++j)
            {
                crossDiff +=
                    dkdxb[k * SPATIAL_DIM + j] * dwdxb[k * SPATIAL_DIM + j];
            }

            // some temps
            const scalar minDSq = yb[k] * yb[k];
            const scalar trbDiss = std::sqrt(tkeb[k]) / betaStar /
                                   (tefb[k] + SMALL) / (yb[k] + SMALL);
            const scalar lamDiss =
                500.0 * mub[k] / rhob[k] / (tefb[k] + SMALL) / (minDSq + SMALL);
            const scalar CDkw = std::max(2.0 * rhob[k] * sigmaWTwo * crossDiff /
                                             (tefb[k] + SMALL),
                                         CDkwClip);

            // arguments
            const scalar fArgOne =
                std::min(std::max(trbDiss, lamDiss),
                         4.0 * rhob[k] * sigmaWTwo * tkeb[k] / (CDkw + SMALL) /
                             (minDSq + SMALL));

            // real deal
            fOneb[k] = std::tanh(fArgOne * fArgOne * fArgOne * fArgOne);
        }
    }
}

void shearStressTransportModel::updateTurbulentProduction(
    const std::shared_ptr<domain> domain)
{
    auto& mesh = this->meshRef();
    stk::mesh::BulkData& bulkData = mesh.bulkDataRef();
    stk::mesh::MetaData& metaData = mesh.metaDataRef();

    // Get fields
    STKScalarField* PkSTKFieldPtr = PkRef().stkFieldPtr();
    const STKScalarField* rhoSTKFieldPtr = rhoRef().stkFieldPtr();
    const STKScalarField* muSTKFieldPtr = muRef().stkFieldPtr();

    // Compute tke production
    {
        const STKScalarField* kSTKFieldPtr = kRef().stkFieldPtr();
        const STKScalarField* mutSTKFieldPtr = mutRef().stkFieldPtr();
        const STKScalarField* gradUSTKFieldPtr = URef().gradRef().stkFieldPtr();

        scalar comp = domain->isMaterialCompressible() ? 1.0 : 0.0;

        // MRF: Get rotation data for Coriolis production term
        const bool frameRotating = domain->zonePtr()->frameRotating();

        // Expert parameter: Check if Coriolis production is enabled
        // Use multiplier to avoid branching in hot loop
        const scalar coriolisProductionMultiplier =
            (frameRotating && controlsRef()
                                  .solverRef()
                                  .solverControl_.expertParameters_
                                  .coriolisProductionTurbulence_)
                ? 1.0
                : 0.0;

        const auto coriolisMatrix = frameRotating ? domain->zonePtr()
                                                        ->transformationRef()
                                                        .rotation()
                                                        .coriolisMatrix_
                                                  : utils::matrix::Zero();
        const scalar* p_mat = coriolisMatrix.data();

        const auto origin =
            frameRotating
                ? domain->zonePtr()->transformationRef().rotation().origin_
                : utils::vector::Zero();
        const scalar* p_ori = origin.data();

        const STKScalarField* coordinatesPtr = metaData.get_field<scalar>(
            stk::topology::NODE_RANK, this->getCoordinatesID_(domain));

        // MRF: Declare work array for Coriolis computation outside hot loop
        std::vector<scalar> omega_cross_r(SPATIAL_DIM);

        // pointers
        scalar* p_omega_cross_r = &omega_cross_r[0];

        stk::mesh::Selector selAllNodes =
            metaData.universal_part() &
            stk::mesh::selectUnion(domain->zonePtr()->interiorParts());

        stk::mesh::BucketVector const& nodeBuckets =
            bulkData.get_buckets(stk::topology::NODE_RANK, selAllNodes);
        for (stk::mesh::BucketVector::const_iterator ib = nodeBuckets.begin();
             ib != nodeBuckets.end();
             ++ib)
        {
            stk::mesh::Bucket& nodeBucket = **ib;
            const stk::mesh::Bucket::size_type nNodesPerBucket =
                nodeBucket.size();

            // get field chunks
            const scalar* kb = stk::mesh::field_data(*kSTKFieldPtr, nodeBucket);
            const scalar* rhob =
                stk::mesh::field_data(*rhoSTKFieldPtr, nodeBucket);
            const scalar* mutb =
                stk::mesh::field_data(*mutSTKFieldPtr, nodeBucket);
            const scalar* gradUb =
                stk::mesh::field_data(*gradUSTKFieldPtr, nodeBucket);
            scalar* Pkb = stk::mesh::field_data(*PkSTKFieldPtr, nodeBucket);

            for (stk::mesh::Bucket::size_type iNode = 0;
                 iNode < nNodesPerBucket;
                 ++iNode)
            {
                scalar PkByMut = 0;
                for (label i = 0; i < SPATIAL_DIM; ++i)
                {
                    const label offSet = SPATIAL_DIM * i;
                    for (label j = 0; j < SPATIAL_DIM; ++j)
                    {
                        PkByMut += gradUb[SPATIAL_DIM * SPATIAL_DIM * iNode +
                                          offSet + j] *
                                   (gradUb[SPATIAL_DIM * SPATIAL_DIM * iNode +
                                           offSet + j] +
                                    gradUb[SPATIAL_DIM * SPATIAL_DIM * iNode +
                                           SPATIAL_DIM * j + i]);
                    }
                }

                // ensure non-negative production
                PkByMut = std::max(PkByMut, 0.0);

                // Add dilation production 1
                scalar divergence = 0.0;
                for (label i = 0; i < SPATIAL_DIM; ++i)
                {
                    const label offSet = SPATIAL_DIM * i;
                    divergence +=
                        gradUb[SPATIAL_DIM * SPATIAL_DIM * iNode + offSet + i];
                }
                PkByMut = std::max(
                    PkByMut - 2.0 / 3.0 * divergence * divergence * comp, 0.0);

                // calculate Pk
                Pkb[iNode] = mutb[iNode] * PkByMut;

                // Add dilation production 2
                Pkb[iNode] = std::max(std::max(Pkb[iNode], 0.0) -
                                          2.0 / 3.0 * divergence * rhob[iNode] *
                                              kb[iNode] * comp,
                                      0.0);

                // MRF: Add Coriolis-induced turbulence production
                // P_k^rotation = 2*mu_t*|Omega x r|*|S|
                // Reference: Bardina et al., "Turbulence Modeling in Rotating
                // Channels"
                // Activated by expert parameter: coriolis_production_turbulence
                // Using multiplier to avoid branching in hot loop
                const scalar* coords =
                    stk::mesh::field_data(*coordinatesPtr, nodeBucket, iNode);

                // Compute Omega x (r - r_origin)
                // Work array omega_cross_r declared outside loop
                for (label i = 0; i < SPATIAL_DIM; i++)
                {
                    p_omega_cross_r[i] = 0.0;
                    for (label j = 0; j < SPATIAL_DIM; j++)
                    {
                        p_omega_cross_r[i] +=
                            p_mat[i * SPATIAL_DIM + j] * (coords[j] - p_ori[j]);
                    }
                }

                // Compute rotation magnitude |Omega x r|
                scalar rotation_mag = 0.0;
                for (label i = 0; i < SPATIAL_DIM; i++)
                {
                    rotation_mag += p_omega_cross_r[i] * p_omega_cross_r[i];
                }
                rotation_mag = std::sqrt(rotation_mag);

                // Strain magnitude already available from PkByMut
                // |S| = sqrt(P_k/mu_t) since P_k = mu_t * S:S
                scalar strain_mag = std::sqrt(std::max(PkByMut, 0.0) + SMALL);

                // Add rotation production with multiplier (no branching)
                // Multiplier is 0.0 if disabled, 1.0 if enabled
                scalar rotation_production = 2.0 * mutb[iNode] * rotation_mag *
                                             strain_mag *
                                             coriolisProductionMultiplier;

                Pkb[iNode] += rotation_production;
            }
        }
    }

    // Zero-out Pk values at no-slip wall nodes
    {
        const auto& assembledWallAreaSTKFieldRef = *metaData.get_field<scalar>(
            stk::topology::NODE_RANK, mesh::assembled_wall_area_ID);

        // sel nodes: only those sitting on no-slip walls
        stk::mesh::Selector selAllNodes =
            metaData.universal_part() &
            stk::mesh::selectUnion(collectNoSlipWallParts_(domain));

        stk::mesh::BucketVector const& nodeBuckets =
            bulkData.get_buckets(stk::topology::NODE_RANK, selAllNodes);
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
                // get node
                stk::mesh::Entity node = nodeBucket[iNode];

                scalar* Pk = stk::mesh::field_data(*PkSTKFieldPtr, node);
                (*Pk) = 0.0;
            }
        }
    }

    // Set Pk values at no-slip wall nodes
    {
        const STKScalarField* uStarSTKFieldPtr = uStarRef().stkFieldPtr();
        const STKScalarField* uTauSTKFieldPtr = uTauRef().stkFieldPtr();
        const STKScalarField* duPlusdyPlusSTKFieldPtr =
            duPlusdyPlusRef().stkFieldPtr();

        const auto& coordsSTKFieldRef = *metaData.get_field<scalar>(
            stk::topology::NODE_RANK, this->getCoordinatesID_(domain));
        const auto& exposedAreaVecSTKFieldRef = *metaData.get_field<scalar>(
            metaData.side_rank(), this->getExposedAreaVectorID_(domain));

        // define vector of parent topos; should always
        // be UNITY in size
        std::vector<stk::topology> parentTopo;

        // nodal fields to gather
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
            MasterElement* meFC =
                accel::MasterElementRepo::get_surface_master_element(
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

            const stk::mesh::Bucket::size_type nSidesPerBucket =
                sideBucket.size();

            for (stk::mesh::Bucket::size_type iSide = 0;
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
                for (label ni = 0; ni < nodesPerSide; ++ni)
                {
                    stk::mesh::Entity node = sideNodeRels[ni];

                    // gather scalars
                    p_rho[ni] = *stk::mesh::field_data(*rhoSTKFieldPtr, node);
                    p_mu[ni] = *stk::mesh::field_data(*muSTKFieldPtr, node);
                }

                // pointer to face data
                const scalar* areaVec =
                    stk::mesh::field_data(exposedAreaVecSTKFieldRef, side);
                const scalar* uStarBip =
                    stk::mesh::field_data(*uStarSTKFieldPtr, side);
                const scalar* uTauBip =
                    stk::mesh::field_data(*uTauSTKFieldPtr, side);
                const scalar* duPlusdyPlusBip =
                    stk::mesh::field_data(*duPlusdyPlusSTKFieldPtr, side);

                // extract the connected element to this
                // exposed face; should be single in
                // size!
                const stk::mesh::Entity* face_elem_rels =
                    bulkData.begin_elements(side);
                STK_ThrowAssert(bulkData.num_elements(side) == 1);

                // get element; its face ordinal number
                stk::mesh::Entity element = face_elem_rels[0];
                const label face_ordinal =
                    bulkData.begin_element_ordinals(side)[0];

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
                    }

                    // definition which avoids u+
                    scalar Pk_val =
                        std::pow(rhoBip, scalar(2.0)) * duPlusdyPlusBip[ip] *
                        std::pow(uStarBip[ip] * uTauBip[ip], scalar(2.0)) /
                        (muBip + SMALL);

                    scalar* Pk = stk::mesh::field_data(*PkSTKFieldPtr, nodeR);
                    (*Pk) += Pk_val * aMag;
                }
            }
        }

        // sync in case of parallel
        if (messager::parallel())
        {
            stk::mesh::communicate_field_data(bulkData, {PkSTKFieldPtr});
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
                scalar* Pk = stk::mesh::field_data(*PkSTKFieldPtr, node);

                (*Pk) /= (*area);
            }
        }
    }
}

void shearStressTransportModel::updateTurbulentDynamicViscosity(
    const std::shared_ptr<domain> domain)
{
    // get common data
    stk::mesh::MetaData& metaData = this->meshRef().metaDataRef();
    stk::mesh::BulkData& bulkData = this->meshRef().bulkDataRef();

    // get fields
    STKScalarField* mutSTKFieldPtr = this->mutRef().stkFieldPtr();
    const STKScalarField* rhoSTKFieldPtr = this->rhoRef().stkFieldPtr();
    const STKScalarField* muSTKFieldPtr = this->muRef().stkFieldPtr();
    const STKScalarField& kSTKFieldRef = this->kRef().stkFieldRef();
    const STKScalarField& omegaSTKFieldRef = this->omegaRef().stkFieldRef();
    const STKScalarField& gradUSTKFieldRef =
        this->URef().gradRef().stkFieldRef();
    const STKScalarField* minDistanceToWallSTKFieldPtr =
        yMinRef().stkFieldPtr();

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

        const scalar* rhob = stk::mesh::field_data(*rhoSTKFieldPtr, b);
        const scalar* mub = stk::mesh::field_data(*muSTKFieldPtr, b);
        const scalar* tkeb = stk::mesh::field_data(kSTKFieldRef, b);
        const scalar* tefb = stk::mesh::field_data(omegaSTKFieldRef, b);
        const scalar* minDb =
            stk::mesh::field_data(*minDistanceToWallSTKFieldPtr, b);
        scalar* mutb = stk::mesh::field_data(*mutSTKFieldPtr, b);

        for (stk::mesh::Bucket::size_type k = 0; k < length; ++k)
        {
            // compute strain rate magnitude; pull pointer within the loop to
            // make it managable
            const scalar* dudx = stk::mesh::field_data(gradUSTKFieldRef, b[k]);
            scalar sijMag = 0.0;
            for (label i = 0; i < SPATIAL_DIM; ++i)
            {
                for (label j = 0; j < SPATIAL_DIM; ++j)
                {
                    const scalar rateOfStrain =
                        0.5 *
                        (dudx[SPATIAL_DIM * i + j] + dudx[SPATIAL_DIM * j + i]);
                    sijMag += rateOfStrain * rateOfStrain;
                }
            }
            sijMag = std::sqrt(2.0 * sijMag);

            // some temps
            const scalar minDSq = minDb[k] * minDb[k];
            const scalar trbDiss = std::sqrt(tkeb[k]) / betaStar() /
                                   (tefb[k] + SMALL) / (minDb[k] + SMALL);
            const scalar lamDiss =
                500.0 * mub[k] / rhob[k] / (tefb[k] + SMALL) / (minDSq + SMALL);
            const scalar fArgTwo = std::max(2.0 * trbDiss, lamDiss);
            const scalar fTwo = std::tanh(fArgTwo * fArgTwo);

            scalar mut_val =
                aOne() * rhob[k] * tkeb[k] /
                (std::max(aOne() * tefb[k], sijMag * fTwo) + SMALL);

            mutb[k] = 0.75 * mut_val + 0.25 * mutb[k];
        }
    }
}

} /* namespace accel */
