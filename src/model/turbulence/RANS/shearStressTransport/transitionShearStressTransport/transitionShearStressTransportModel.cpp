// File       : transitionShearStressTransportModel.cpp
// Created    : Mon Jan 14 2025
// Author     : Adam Fares
// Description:
// Copyright (c) 2025 CCFNUM, Lucerne University of Applied Sciences and Arts.
// SPDX-License-Identifier: BSD-3-Clause

// code
#include "transitionShearStressTransportModel.h"
#include "messager.h"
#include "realm.h"

namespace accel
{

transitionShearStressTransportModel::transitionShearStressTransportModel(
    realm* realm)
    : shearStressTransportModel(realm)
{
    // Create field instances:

    // transport
    gammaRef();
    ReThetaRef();
}

void transitionShearStressTransportModel::
    initializeTransitionOnsetReynoldsNumber(
        const std::shared_ptr<domain> domain)
{
    updateTransitionOnsetReynoldsNumber(domain);
}

void transitionShearStressTransportModel::updateTransitionOnsetReynoldsNumber(
    const std::shared_ptr<domain> domain)
{
    fieldBroker::updateTransitionOnsetReynoldsNumber(domain);

    updateTransitionOnsetReynoldsNumberSideFields_(domain);
}

void transitionShearStressTransportModel::
    updateTransitionOnsetReynoldsNumberSideFields_(
        const std::shared_ptr<domain> domain)
{
    for (label iBoundary = 0;
         iBoundary < this->meshRef().zonePtr(domain->index())->nBoundaries();
         iBoundary++)
    {
        const auto* boundary =
            this->meshRef().zonePtr(domain->index())->boundaryPtr(iBoundary);

        boundaryPhysicalType physicalType = boundary->type();

        switch (physicalType)
        {
            case boundaryPhysicalType::inlet:
            case boundaryPhysicalType::opening:
                {
                    // boolean to translate boundary values to field
                    const bool correctedBoundaryNodeValues =
                        ReThetaRef().correctedBoundaryNodeValues();

                    // common
                    stk::mesh::BulkData& bulkData =
                        this->meshRef().bulkDataRef();
                    stk::mesh::MetaData& metaData =
                        this->meshRef().metaDataRef();

                    // Get fields
                    const STKScalarField* USTKFieldPtr =
                        this->URef().stkFieldPtr();
                    const STKScalarField* kSTKFieldPtr =
                        this->kRef().stkFieldPtr();
                    STKScalarField* nodeSideReThetaSTKFieldPtr =
                        this->ReThetaRef().nodeSideFieldRef().stkFieldPtr();

                    // get interior parts the domain is defined on
                    const stk::mesh::PartVector& partVec =
                        domain->zonePtr()->interiorParts();

                    // define some common selectors; select owned nodes
                    stk::mesh::Selector selUniversalNodes =
                        metaData.universal_part() &
                        stk::mesh::selectUnion(boundary->parts());

                    stk::mesh::BucketVector const& nodeBuckets =
                        bulkData.get_buckets(stk::topology::NODE_RANK,
                                             selUniversalNodes);

                    for (stk::mesh::BucketVector::const_iterator ib =
                             nodeBuckets.begin();
                         ib != nodeBuckets.end();
                         ++ib)
                    {
                        stk::mesh::Bucket& nodeBucket = **ib;

                        const stk::mesh::Bucket::size_type nNodesPerBucket =
                            nodeBucket.size();

                        // field chunks in bucket
                        const scalar* Ub =
                            stk::mesh::field_data(*USTKFieldPtr, nodeBucket);
                        const scalar* kb =
                            stk::mesh::field_data(*kSTKFieldPtr, nodeBucket);
                        scalar* ReThetab = stk::mesh::field_data(
                            *nodeSideReThetaSTKFieldPtr, nodeBucket);

                        for (stk::mesh::Bucket::size_type iNode = 0;
                             iNode < nNodesPerBucket;
                             ++iNode)
                        {
                            scalar Umag = 0.0;
                            for (label i = 0; i < SPATIAL_DIM; i++)
                            {
                                Umag += Ub[SPATIAL_DIM * iNode + i] *
                                        Ub[SPATIAL_DIM * iNode + i];
                            }
                            Umag = std::sqrt(Umag);

                            scalar TuInfty = 100.0 *
                                             std::sqrt(2 * kb[iNode] / 3.0) /
                                             (Umag + SMALL);
                            TuInfty = std::max(TuInfty, 0.027);

                            if (TuInfty > 1.3)
                            {
                                ReThetab[iNode] =
                                    331.50 *
                                    std::pow((TuInfty - 0.5658), -0.671);
                            }
                            else
                            {
                                ReThetab[iNode] =
                                    1173.51 - 589.428 * TuInfty +
                                    0.2196 * 1.0 / (TuInfty * TuInfty + SMALL);
                            }
                        }
                    }

                    // interpolate to node-side field
                    ReThetaRef().sideFieldRef().interpolate(
                        ReThetaRef().nodeSideFieldRef(),
                        domain->index(),
                        boundary->index(),
                        this->ReThetaRef().isShifted());

                    // translate to field
                    if (correctedBoundaryNodeValues)
                    {
                        ReThetaRef().correctBoundaryNodes(domain->index(),
                                                          boundary->index());
                    }
                }
                break;

            default:
                break;
        }
    }
}

void transitionShearStressTransportModel::updateFOneBlending(
    const std::shared_ptr<domain> domain)
{
    // compute fone with parameters appropriate for 2003 SST implementation
    stk::mesh::BulkData& bulkData = this->meshRef().bulkDataRef();
    stk::mesh::MetaData& metaData = this->meshRef().metaDataRef();

    const scalar CDkwClip = 1.0e-10; // 2003 SST
    const scalar sigmaWTwo = this->sigmaWTwo();

    // required fields with state; min_distance is fine
    const STKScalarField& kSTKFieldRef = kRef().stkFieldRef();
    const STKScalarField& omegaSTKFieldRef = omegaRef().stkFieldRef();

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
        const scalar* minDb =
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
            const scalar minDSq = minDb[k] * minDb[k];
            const scalar trbDiss = std::sqrt(tkeb[k]) / betaStar() /
                                   (tefb[k] + SMALL) / (minDb[k] + SMALL);
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

            // Transition SST Modification
            scalar Ry = rhob[k] * minDb[k] * std::sqrt(tkeb[k]) / mub[k];
            scalar F3 = std::exp(-1 * std::pow((Ry / 120.0), 8.0));

            fOneb[k] =
                std::max(std::tanh(fArgOne * fArgOne * fArgOne * fArgOne), F3);
        }
    }
}

void transitionShearStressTransportModel::updateTurbulentProduction(
    const std::shared_ptr<domain> domain)
{
    stk::mesh::BulkData& bulkData = this->meshRef().bulkDataRef();
    stk::mesh::MetaData& metaData = this->meshRef().metaDataRef();

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
            }
        }
    }

    // Zero-out Pk values at no-slip wall nodes
    {
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

                // loop over face nodes
                for (label ip = 0; ip < numSideNodes; ++ip)
                {
                    const label offSetSF_face = ip * nodesPerSide;
                    const label offSetAreaVec = ip * SPATIAL_DIM;

                    const label localFaceNode = faceIpNodeMap[ip];

                    // left and right nodes; right is on
                    // the face; left is the opposing
                    // node
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

} /* namespace accel */
