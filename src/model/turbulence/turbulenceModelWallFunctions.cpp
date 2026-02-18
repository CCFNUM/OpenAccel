// File : turbulenceModelWallFunctions.cpp
// Created : Fri Jun 20 2024 16:48:19 (+0100)
// Author : Mhamad Mahdi Alloush
// Description:
// Copyright (c) 2024 CCFNUM, Lucerne University of Applied Sciences and Arts.
// SPDX-License-Identifier: BSD-3-Clause

#include "initialConditions.h"
#include "messager.h"
#include "realm.h"
#include "turbulenceModel.h"

namespace accel
{

void turbulenceModel::updateYStar(const std::shared_ptr<domain> domain)
{
    switch (domain->turbulence_.wallFunctionType_)
    {
        case wallFunctionType::standard:
            errorMsg("Not implemented yet");
            break;

        case wallFunctionType::scalable:
            updateYStarScalable_(domain);
            break;

        case wallFunctionType::automatic:
            updateYStarAutomatic_(domain);
            break;
    }
}

void turbulenceModel::updateUStar(const std::shared_ptr<domain> domain)
{
    switch (domain->turbulence_.wallFunctionType_)
    {
        case wallFunctionType::standard:
            errorMsg("Not implemented yet");
            break;

        case wallFunctionType::scalable:
            updateUStarScalable_(domain);
            break;

        case wallFunctionType::automatic:
            updateUStarAutomatic_(domain);
            break;
    }
}

void turbulenceModel::updateYPlus(const std::shared_ptr<domain> domain)
{
    switch (domain->turbulence_.wallFunctionType_)
    {
        case wallFunctionType::standard:
            errorMsg("Not implemented yet");
            break;

        case wallFunctionType::scalable:
            updateYPlusScalable_(domain);
            break;

        case wallFunctionType::automatic:
            updateYPlusAutomatic_(domain);
            break;
    }
}

void turbulenceModel::updateUTau(const std::shared_ptr<domain> domain)
{
    switch (domain->turbulence_.wallFunctionType_)
    {
        case wallFunctionType::standard:
            errorMsg("Not implemented yet");
            break;

        case wallFunctionType::scalable:
            updateUTauScalable_(domain);
            break;

        case wallFunctionType::automatic:
            updateUTauAutomatic_(domain);
            break;
    }
}

void turbulenceModel::updateUPlus(const std::shared_ptr<domain> domain)
{
    stk::mesh::BulkData& bulkData = this->meshRef().bulkDataRef();
    stk::mesh::MetaData& metaData = this->meshRef().metaDataRef();

    // Get fields
    const auto& USTKFieldRef = URef().stkFieldRef();
    const auto& sideUSTKFieldRef = URef().sideFieldRef().stkFieldRef();

    const STKScalarField* uTauSTKFieldPtr = uTauRef().stkFieldPtr();
    STKScalarField* uPlusSTKFieldPtr = uPlusRef().stkFieldPtr();

    const auto& coordsSTKFieldRef = *metaData.get_field<scalar>(
        stk::topology::NODE_RANK, this->getCoordinatesID_(domain));
    const auto& exposedAreaVecSTKFieldRef = *metaData.get_field<scalar>(
        metaData.side_rank(), this->getExposedAreaVectorID_(domain));

    // define vector of parent topos; should always
    // be UNITY in size
    std::vector<stk::topology> parentTopo;

    // bip values
    std::vector<scalar> uBip(SPATIAL_DIM);
    std::vector<scalar> unitNormal(SPATIAL_DIM);

    // pointers to fixed values
    scalar* p_uBip = &uBip[0];
    scalar* p_unitNormal = &unitNormal[0];

    // nodal fields to gather
    std::vector<scalar> ws_U;

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

        // face master element
        MasterElement* meFC = MasterElementRepo::get_surface_master_element(
            sideBucket.topology());
        const label nodesPerSide = meFC->nodesPerElement_;
        const label numScsBip = meFC->numIntPoints_;

        ws_U.resize(nodesPerSide * SPATIAL_DIM);
        ws_face_shape_function.resize(numScsBip * nodesPerSide);

        // pointers
        scalar* p_U = &ws_U[0];
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

        const stk::mesh::Bucket::size_type nSidesPerBucket = sideBucket.size();

        for (stk::mesh::Bucket::size_type iSide = 0; iSide < nSidesPerBucket;
             ++iSide)
        {
            // get face
            stk::mesh::Entity side = sideBucket[iSide];

            //======================================
            // gather nodal data off of face
            //======================================
            stk::mesh::Entity const* sideNodeRels = bulkData.begin_nodes(side);
            label numSideNodes = bulkData.num_nodes(side);

            // sanity check on num nodes
            STK_ThrowAssert(numSideNodes == nodesPerSide);
            for (label ni = 0; ni < nodesPerSide; ++ni)
            {
                stk::mesh::Entity node = sideNodeRels[ni];

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
                stk::mesh::field_data(exposedAreaVecSTKFieldRef, side);
            const scalar* uTauBip =
                stk::mesh::field_data(*uTauSTKFieldPtr, side);
            const scalar* UbcVec =
                stk::mesh::field_data(sideUSTKFieldRef, side);
            scalar* uPlusBip = stk::mesh::field_data(*uPlusSTKFieldPtr, side);

            // extract the connected element to this
            // exposed face; should be single in
            // size!
            const stk::mesh::Entity* face_elem_rels =
                bulkData.begin_elements(side);
            STK_ThrowAssert(bulkData.num_elements(side) == 1);

            // get element; its face ordinal number
            stk::mesh::Entity element = face_elem_rels[0];

            // loop over face nodes
            for (label ip = 0; ip < numSideNodes; ++ip)
            {
                const label offSetSF_face = ip * nodesPerSide;
                const label offSetAveraVec = ip * SPATIAL_DIM;
                // zero out vector quantities;
                // squeeze in aMag
                scalar aMag = 0.0;
                for (label j = 0; j < SPATIAL_DIM; ++j)
                {
                    p_uBip[j] = 0.0;
                    const scalar axj = areaVec[offSetAveraVec + j];
                    aMag += axj * axj;
                }
                aMag = std::sqrt(aMag);

                // interpolate to bip
                for (label ic = 0; ic < nodesPerSide; ++ic)
                {
                    const scalar r = p_face_shape_function[offSetSF_face + ic];

                    const label offSetFN = ic * SPATIAL_DIM;
                    for (label j = 0; j < SPATIAL_DIM; ++j)
                    {
                        p_uBip[j] += r * p_U[offSetFN + j];
                    }
                }

                // form unit normal
                for (label j = 0; j < SPATIAL_DIM; ++j)
                {
                    const scalar nj = areaVec[offSetAveraVec + j] / aMag;
                    p_unitNormal[j] = nj;
                }

                // determine tangential velocity
                scalar uTangential = 0.0;
                for (label i = 0; i < SPATIAL_DIM; ++i)
                {
                    scalar uiTan = 0.0;
                    scalar uiBcTan = 0.0;
                    for (label j = 0; j < SPATIAL_DIM; ++j)
                    {
                        const scalar ninj = p_unitNormal[i] * p_unitNormal[j];
                        if (i == j)
                        {
                            const scalar om_nini = 1.0 - ninj;
                            uiTan += om_nini * p_uBip[j];
                            uiBcTan += om_nini * UbcVec[ip * SPATIAL_DIM + j];
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

                // For bad initialization utau could be zero because uTangential
                // could be zero
                if (uTauBip[ip] < SMALL)
                {
                    // So that uWallCoeffs field is zero
                    uPlusBip[ip] = BIG;
                }
                else
                {
                    uPlusBip[ip] = uTangential / (uTauBip[ip] + SMALL);
                }
            }
        }
    }
}

void turbulenceModel::updateTPlus(const std::shared_ptr<domain> domain)
{
    if (domain->heatTransfer_.option_ == heatTransferOption::none)
    {
        return;
    }

    switch (domain->turbulence_.wallFunctionType_)
    {
        case wallFunctionType::standard:
            errorMsg("Not implemented yet");
            break;

        case wallFunctionType::scalable:
            updateTPlusScalable_(domain);
            break;

        case wallFunctionType::automatic:
            updateTPlusAutomatic_(domain);
            break;
    }
}

void turbulenceModel::updateDuPlusDyPlus(const std::shared_ptr<domain> domain)
{
    stk::mesh::BulkData& bulkData = this->meshRef().bulkDataRef();
    stk::mesh::MetaData& metaData = this->meshRef().metaDataRef();

    // Get fields
    const STKScalarField* yPlusSTKFieldPtr = yPlusRef().stkFieldPtr();
    STKScalarField* duPlusdyPlusSTKFieldPtr = duPlusdyPlusRef().stkFieldPtr();

    // select all sides: only those sitting on
    // no-slip walls
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

        const stk::mesh::Bucket::size_type nSidesPerBucket = sideBucket.size();

        for (stk::mesh::Bucket::size_type iSide = 0; iSide < nSidesPerBucket;
             ++iSide)
        {
            // get face
            stk::mesh::Entity side = sideBucket[iSide];
            label numSideNodes = bulkData.num_nodes(side);

            // pointer to face data
            const scalar* yPlusBip =
                stk::mesh::field_data(*yPlusSTKFieldPtr, side);
            scalar* duPlusdyPlusBip =
                stk::mesh::field_data(*duPlusdyPlusSTKFieldPtr, side);

            // loop over face nodes
            for (label ip = 0; ip < numSideNodes; ++ip)
            {
                scalar blendFactor = scalar(0.01) *
                                     std::pow(yPlusBip[ip], scalar(4)) /
                                     (scalar(1.0) + scalar(5.0) * yPlusBip[ip]);
                duPlusdyPlusBip[ip] =
                    scalar(1.0) * exp(-blendFactor + SMALL) +
                    scalar(1.0) / (kappa() * yPlusBip[ip]) *
                        std::exp(-scalar(1) / (blendFactor + SMALL));
            }
        }
    }
}

void turbulenceModel::updateUWallCoeffs(const std::shared_ptr<domain> domain)
{
    stk::mesh::BulkData& bulkData = this->meshRef().bulkDataRef();
    stk::mesh::MetaData& metaData = this->meshRef().metaDataRef();

    // Get fields
    const STKScalarField* uPlusSTKFieldPtr = uPlusRef().stkFieldPtr();
    const STKScalarField* uStarSTKFieldPtr = uStarRef().stkFieldPtr();
    STKScalarField* uWallCoeffsSTKFieldPtr = uWallCoeffsRef().stkFieldPtr();

    auto& rhoSTKFieldRef = rhoRef().stkFieldRef();
    auto& muSTKFieldRef = muRef().stkFieldRef();

    const auto& wallNormalDistanceSTKFieldRef = *metaData.get_field<scalar>(
        metaData.side_rank(), mesh::wall_normal_distance_ID);
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

        const stk::mesh::Bucket::size_type nSidesPerBucket = sideBucket.size();

        for (stk::mesh::Bucket::size_type iSide = 0; iSide < nSidesPerBucket;
             ++iSide)
        {
            // get face
            stk::mesh::Entity side = sideBucket[iSide];
            label numSideNodes = bulkData.num_nodes(side);

            //======================================
            // gather nodal data off of face
            //======================================
            stk::mesh::Entity const* sideNodeRels = bulkData.begin_nodes(side);
            for (label ni = 0; ni < nodesPerSide; ++ni)
            {
                stk::mesh::Entity node = sideNodeRels[ni];

                // gather scalars
                p_rho[ni] = *stk::mesh::field_data(rhoSTKFieldRef, node);
                p_mu[ni] = *stk::mesh::field_data(muSTKFieldRef, node);
            }

            // pointer to face data
            const scalar* areaVec =
                stk::mesh::field_data(exposedAreaVecSTKFieldRef, side);
            const scalar* uPlusBip =
                stk::mesh::field_data(*uPlusSTKFieldPtr, side);
            const scalar* uStarBip =
                stk::mesh::field_data(*uStarSTKFieldPtr, side);

            scalar* uWallCoeffsBip =
                stk::mesh::field_data(*uWallCoeffsSTKFieldPtr, side);

            // loop over face nodes
            for (label ip = 0; ip < numSideNodes; ++ip)
            {
                const label offSetAreaVec = ip * SPATIAL_DIM;
                const label offSetSF_face = ip * nodesPerSide;

                // interpolate to bip
                scalar rhoBip = 0.0;
                scalar muBip = 0.0;
                for (label ic = 0; ic < nodesPerSide; ++ic)
                {
                    const scalar r = p_face_shape_function[offSetSF_face + ic];
                    rhoBip += r * p_rho[ic];
                    muBip += r * p_mu[ic];
                }

                // squeeze in aMag
                scalar aMag = 0.0;
                for (label j = 0; j < SPATIAL_DIM; ++j)
                {
                    const scalar axj = areaVec[offSetAreaVec + j];
                    aMag += axj * axj;
                }
                aMag = std::sqrt(aMag);

                uWallCoeffsBip[ip] =
                    rhoBip * uStarBip[ip] / uPlusBip[ip] * aMag;
            }
        }
    }
}

void turbulenceModel::updateTWallCoeffs(const std::shared_ptr<domain> domain)
{
    if (domain->heatTransfer_.option_ == heatTransferOption::none)
    {
        return;
    }

    // update TWallCoeffs field at no slip walls and fluid-solid interfaces, the
    // fluid side
    {
        stk::mesh::BulkData& bulkData = this->meshRef().bulkDataRef();
        stk::mesh::MetaData& metaData = this->meshRef().metaDataRef();

        // Get fields
        const STKScalarField* TPlusSTKFieldPtr = TPlusRef().stkFieldPtr();
        const STKScalarField* uStarSTKFieldPtr = uStarRef().stkFieldPtr();
        STKScalarField* TWallCoeffsSTKFieldPtr = TWallCoeffsRef().stkFieldPtr();

        const auto& rhoSTKFieldRef = rhoRef().stkFieldRef();
        const auto& cpSTKFieldRef = cpRef().stkFieldRef();

        // nodal fields to gather
        std::vector<scalar> ws_rho;
        std::vector<scalar> ws_cp;

        // master element
        std::vector<scalar> ws_face_shape_function;

        // select all sides: only those sitting on
        // no-slip walls
        stk::mesh::Selector selAllSides =
            metaData.universal_part() &
            stk::mesh::selectUnion(collectNoSlipWallParts_(domain));

        // shifted ip's for fields?
        const bool isTShifted = TRef().isShifted();

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

            // algorithm related; element
            ws_rho.resize(nodesPerSide);
            ws_cp.resize(nodesPerSide);
            ws_face_shape_function.resize(numScsBip * nodesPerSide);

            // pointers
            scalar* p_rho = &ws_rho[0];
            scalar* p_cp = &ws_cp[0];
            scalar* p_face_shape_function = &ws_face_shape_function[0];

            // shape functions
            if (isTShifted)
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
                label numSideNodes = bulkData.num_nodes(side);

                //======================================
                // gather nodal data off of face
                //======================================
                stk::mesh::Entity const* sideNodeRels =
                    bulkData.begin_nodes(side);
                for (label ni = 0; ni < nodesPerSide; ++ni)
                {
                    stk::mesh::Entity node = sideNodeRels[ni];

                    // gather scalars
                    p_rho[ni] = *stk::mesh::field_data(rhoSTKFieldRef, node);
                    p_cp[ni] = *stk::mesh::field_data(cpSTKFieldRef, node);
                }

                // pointer to face data
                const scalar* TPlusBip =
                    stk::mesh::field_data(*TPlusSTKFieldPtr, side);
                const scalar* uStarBip =
                    stk::mesh::field_data(*uStarSTKFieldPtr, side);

                scalar* TWallCoeffsBip =
                    stk::mesh::field_data(*TWallCoeffsSTKFieldPtr, side);

                // loop over face nodes
                for (label ip = 0; ip < numSideNodes; ++ip)
                {
                    const label offSetSF_face = ip * nodesPerSide;

                    // interpolate to bip
                    scalar rhoBip = 0.0;
                    scalar cpBip = 0.0;
                    for (label ic = 0; ic < nodesPerSide; ++ic)
                    {
                        const scalar r =
                            p_face_shape_function[offSetSF_face + ic];
                        rhoBip += r * p_rho[ic];
                        cpBip += r * p_cp[ic];
                    }

                    TWallCoeffsBip[ip] =
                        rhoBip * cpBip * uStarBip[ip] / TPlusBip[ip];
                }
            }
        }
    }

#ifdef HAS_INTERFACE
    // interpolate T-wall coefficient field to the solid side: use temporary
    // sideField container to make use of the transfer functionality
    sideField<scalar, 1> t_TWallField(&this->meshRef(),
                                      this->TWallCoeffsRef().stkFieldPtr());

    for (const interface* interf : domain->interfacesRef())
    {
        if (interf->isFluidSolidType())
        {
            t_TWallField.transfer(interf->index(),
                                  !interf->isMasterZone(domain->index()),
                                  TRef().isShifted());
        }
    }
#endif /* HAS_INTERFACE */
}

void turbulenceModel::updateWallShearStress(
    const std::shared_ptr<domain> domain)
{
    stk::mesh::BulkData& bulkData = this->meshRef().bulkDataRef();
    stk::mesh::MetaData& metaData = this->meshRef().metaDataRef();

    // Get fields
    STKScalarField* wallShearStressSTKFieldPtr =
        wallShearStressRef().stkFieldPtr();
    const STKScalarField* uWallCoeffsSTKFieldPtr =
        uWallCoeffsRef().stkFieldPtr();
    const auto& USTKFieldRef = URef().stkFieldRef();
    const auto& sideUSTKFieldRef = URef().sideFieldRef().stkFieldRef();

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

        // algorithm related; element
        ws_U.resize(nodesPerSide * SPATIAL_DIM);
        ws_face_shape_function.resize(numScsBip * nodesPerSide);

        // pointers
        scalar* p_U = &ws_U[0];
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

        const stk::mesh::Bucket::size_type nSidesPerBucket = sideBucket.size();

        for (stk::mesh::Bucket::size_type iSide = 0; iSide < nSidesPerBucket;
             ++iSide)
        {
            // get face
            stk::mesh::Entity side = sideBucket[iSide];

            //======================================
            // gather nodal data off of face
            //======================================
            stk::mesh::Entity const* sideNodeRels = bulkData.begin_nodes(side);
            label numSideNodes = bulkData.num_nodes(side);

            // sanity check on num nodes
            STK_ThrowAssert(numSideNodes == nodesPerSide);
            for (label ni = 0; ni < nodesPerSide; ++ni)
            {
                stk::mesh::Entity node = sideNodeRels[ni];

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
                stk::mesh::field_data(exposedAreaVecSTKFieldRef, side);
            const scalar* UbcVec =
                stk::mesh::field_data(sideUSTKFieldRef, side);
            const scalar* uWallCoeffsBip =
                stk::mesh::field_data(*uWallCoeffsSTKFieldPtr, side);
            scalar* wallShearStressBip =
                stk::mesh::field_data(*wallShearStressSTKFieldPtr, side);

            // loop over face nodes
            for (label ip = 0; ip < numSideNodes; ++ip)
            {
                const label offSetSF_face = ip * nodesPerSide;
                const label offSetAveraVec = ip * SPATIAL_DIM;

                // zero out vector quantities;
                // squeeze in aMag
                scalar aMag = 0.0;
                for (label j = 0; j < SPATIAL_DIM; ++j)
                {
                    p_uBip[j] = 0.0;
                    const scalar axj = areaVec[offSetAveraVec + j];
                    aMag += axj * axj;
                }
                aMag = std::sqrt(aMag);

                // interpolate to bip
                for (label ic = 0; ic < nodesPerSide; ++ic)
                {
                    const scalar r = p_face_shape_function[offSetSF_face + ic];

                    const label offSetFN = ic * SPATIAL_DIM;
                    for (label j = 0; j < SPATIAL_DIM; ++j)
                    {
                        p_uBip[j] += r * p_U[offSetFN + j];
                    }
                }

                // form unit normal
                for (label j = 0; j < SPATIAL_DIM; ++j)
                {
                    const scalar nj = areaVec[offSetAveraVec + j] / aMag;
                    p_nx[j] = nj;
                }

                // determine tangential velocity
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
                            uiBcTan += om_nini * UbcVec[ip * SPATIAL_DIM + j];
                        }
                        else
                        {
                            uiTan -= ninj * p_uBip[j];
                            uiBcTan -= ninj * UbcVec[ip * SPATIAL_DIM + j];
                        }
                    }

                    // calculate wall shear stress
                    wallShearStressBip[ip * SPATIAL_DIM + i] =
                        uWallCoeffsBip[ip] * (uiTan - uiBcTan) / aMag;
                }
            }
        }
    }
}

void turbulenceModel::updateYStarScalable_(const std::shared_ptr<domain> domain)
{
    stk::mesh::BulkData& bulkData = this->meshRef().bulkDataRef();
    stk::mesh::MetaData& metaData = this->meshRef().metaDataRef();

    // Get fields
    STKScalarField* yStarSTKFieldPtr = yStarRef().stkFieldPtr();
    const STKScalarField* uStarSTKFieldPtr = uStarRef().stkFieldPtr();
    const STKScalarField* wallNormalDistanceSTKFieldPtr =
        metaData.get_field<scalar>(metaData.side_rank(),
                                   mesh::wall_normal_distance_ID);

    const auto& rhoSTKFieldRef = rhoRef().stkFieldRef();
    const auto& muSTKFieldRef = muRef().stkFieldRef();
    const auto& sideUSTKFieldRef = URef().sideFieldRef().stkFieldRef();

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

        const stk::mesh::Bucket::size_type nSidesPerBucket = sideBucket.size();

        for (stk::mesh::Bucket::size_type iSide = 0; iSide < nSidesPerBucket;
             ++iSide)
        {
            // get face
            stk::mesh::Entity side = sideBucket[iSide];
            label numSideNodes = bulkData.num_nodes(side);

            //======================================
            // gather nodal data off of face
            //======================================
            stk::mesh::Entity const* sideNodeRels = bulkData.begin_nodes(side);
            for (label ni = 0; ni < nodesPerSide; ++ni)
            {
                stk::mesh::Entity node = sideNodeRels[ni];

                // gather scalars
                p_rho[ni] = *stk::mesh::field_data(rhoSTKFieldRef, node);
                p_mu[ni] = *stk::mesh::field_data(muSTKFieldRef, node);
            }

            // pointer to face data
            const scalar* wallNormalDistanceBip =
                stk::mesh::field_data(*wallNormalDistanceSTKFieldPtr, side);
            const scalar* uStarBip =
                stk::mesh::field_data(*uStarSTKFieldPtr, side);
            scalar* yStarBip = stk::mesh::field_data(*yStarSTKFieldPtr, side);

            // loop over face nodes
            for (label ip = 0; ip < numSideNodes; ++ip)
            {
                const label offSetSF_face = ip * nodesPerSide;

                // interpolate to bip
                scalar muBip = 0.0;
                scalar rhoBip = 0.0;
                for (label ic = 0; ic < nodesPerSide; ++ic)
                {
                    const scalar r = p_face_shape_function[offSetSF_face + ic];
                    rhoBip += r * p_rho[ic];
                    muBip += r * p_mu[ic];
                }

                // extract bip data
                const scalar yp = wallNormalDistanceBip[ip] / 4.0;
                const scalar uStar = uStarBip[ip];

                // determine yplus
                yStarBip[ip] = rhoBip * yp * uStar / muBip;
            }
        }
    }
}

void turbulenceModel::updateUStarScalable_(const std::shared_ptr<domain> domain)
{
    stk::mesh::BulkData& bulkData = this->meshRef().bulkDataRef();
    stk::mesh::MetaData& metaData = this->meshRef().metaDataRef();

    // Get fields
    STKScalarField* uStarSTKFieldPtr = uStarRef().stkFieldPtr();

    const auto& kSTKFieldRef = kRef().stkFieldRef();

    // nodal fields to gather
    std::vector<scalar> ws_k;

    // master element
    std::vector<scalar> ws_face_shape_function;

    // select all sides: only those sitting on
    // no-slip walls
    stk::mesh::Selector selAllSides =
        metaData.universal_part() &
        stk::mesh::selectUnion(collectNoSlipWallParts_(domain));

    // shifted ip's for fields?
    const bool isKShifted = kRef().isShifted();

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

        // algorithm related; element
        ws_k.resize(nodesPerSide);
        ws_face_shape_function.resize(numScsBip * nodesPerSide);

        // pointers
        scalar* p_k = &ws_k[0];
        scalar* p_face_shape_function = &ws_face_shape_function[0];

        // shape functions
        if (isKShifted)
        {
            meFC->shifted_shape_fcn(&p_face_shape_function[0]);
        }
        else
        {
            meFC->shape_fcn(&p_face_shape_function[0]);
        }

        const stk::mesh::Bucket::size_type nSidesPerBucket = sideBucket.size();

        for (stk::mesh::Bucket::size_type iSide = 0; iSide < nSidesPerBucket;
             ++iSide)
        {
            // get face
            stk::mesh::Entity side = sideBucket[iSide];

            //======================================
            // gather nodal data off of face
            //======================================
            stk::mesh::Entity const* sideNodeRels = bulkData.begin_nodes(side);
            label numSideNodes = bulkData.num_nodes(side);

            // sanity check on num nodes
            STK_ThrowAssert(numSideNodes == nodesPerSide);
            for (label ni = 0; ni < nodesPerSide; ++ni)
            {
                stk::mesh::Entity node = sideNodeRels[ni];

                // gather scalars
                p_k[ni] = *stk::mesh::field_data(kSTKFieldRef, node);
            }

            // pointer to face data
            scalar* uStarBip = stk::mesh::field_data(*uStarSTKFieldPtr, side);

            // extract the connected element to this
            // exposed face; should be single in
            // size!
            const stk::mesh::Entity* face_elem_rels =
                bulkData.begin_elements(side);
            STK_ThrowAssert(bulkData.num_elements(side) == 1);

            // get element; its face ordinal number
            stk::mesh::Entity element = face_elem_rels[0];

            // loop over face nodes
            for (label ip = 0; ip < numSideNodes; ++ip)
            {
                const label offSetSF_face = ip * nodesPerSide;

                // interpolate to bip
                scalar kBip = 0.0;
                for (label ic = 0; ic < nodesPerSide; ++ic)
                {
                    const scalar r = p_face_shape_function[offSetSF_face + ic];
                    kBip += r * p_k[ic];
                }

                uStarBip[ip] = std::pow(Cmu_, 0.25) * sqrt(kBip);
            }
        }
    }
}

void turbulenceModel::updateYPlusScalable_(const std::shared_ptr<domain> domain)
{
    stk::mesh::BulkData& bulkData = this->meshRef().bulkDataRef();
    stk::mesh::MetaData& metaData = this->meshRef().metaDataRef();

    // Get fields
    STKScalarField* yPlusSTKFieldPtr = yPlusRef().stkFieldPtr();
    const STKScalarField* uStarSTKFieldPtr = uStarRef().stkFieldPtr();
    const STKScalarField* wallNormalDistanceSTKFieldPtr =
        metaData.get_field<scalar>(metaData.side_rank(),
                                   mesh::wall_normal_distance_ID);

    const auto& rhoSTKFieldRef = rhoRef().stkFieldRef();
    const auto& muSTKFieldRef = muRef().stkFieldRef();
    const auto& sideUSTKFieldRef = URef().sideFieldRef().stkFieldRef();

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

        const stk::mesh::Bucket::size_type nSidesPerBucket = sideBucket.size();

        for (stk::mesh::Bucket::size_type iSide = 0; iSide < nSidesPerBucket;
             ++iSide)
        {
            // get face
            stk::mesh::Entity side = sideBucket[iSide];
            label numSideNodes = bulkData.num_nodes(side);

            //======================================
            // gather nodal data off of face
            //======================================
            stk::mesh::Entity const* sideNodeRels = bulkData.begin_nodes(side);
            for (label ni = 0; ni < nodesPerSide; ++ni)
            {
                stk::mesh::Entity node = sideNodeRels[ni];

                // gather scalars
                p_rho[ni] = *stk::mesh::field_data(rhoSTKFieldRef, node);
                p_mu[ni] = *stk::mesh::field_data(muSTKFieldRef, node);
            }

            // pointer to face data
            const scalar* wallNormalDistanceBip =
                stk::mesh::field_data(*wallNormalDistanceSTKFieldPtr, side);
            const scalar* uStarBip =
                stk::mesh::field_data(*uStarSTKFieldPtr, side);
            scalar* yPlusBip = stk::mesh::field_data(*yPlusSTKFieldPtr, side);

            // loop over face nodes
            for (label ip = 0; ip < numSideNodes; ++ip)
            {
                const label offSetSF_face = ip * nodesPerSide;

                // interpolate to bip
                scalar muBip = 0.0;
                scalar rhoBip = 0.0;
                for (label ic = 0; ic < nodesPerSide; ++ic)
                {
                    const scalar r = p_face_shape_function[offSetSF_face + ic];
                    rhoBip += r * p_rho[ic];
                    muBip += r * p_mu[ic];
                }

                // extract bip data
                const scalar yp = wallNormalDistanceBip[ip] / 4.0;
                const scalar uStar = uStarBip[ip];

                // determine yplus
                yPlusBip[ip] = std::max(rhoBip * yp * uStar / muBip, 11.06);
            }
        }
    }
}

void turbulenceModel::updateUTauScalable_(const std::shared_ptr<domain> domain)
{
    stk::mesh::BulkData& bulkData = meshRef().bulkDataRef();
    stk::mesh::MetaData& metaData = meshRef().metaDataRef();

    // bip values
    std::vector<scalar> uBip(SPATIAL_DIM);
    std::vector<scalar> unitNormal(SPATIAL_DIM);

    // pointers to fixed values
    scalar* p_uBip = &uBip[0];
    scalar* p_unitNormal = &unitNormal[0];

    // nodal fields to gather
    std::vector<scalar> ws_U;
    std::vector<scalar> ws_rho;
    std::vector<scalar> ws_mu;

    // master element
    std::vector<scalar> ws_face_shape_function;

    // deal with state
    auto& USTKFieldRef = URef().stkFieldRef();
    auto& sideUSTKFieldRef = URef().sideFieldRef().stkFieldRef();
    auto& rhoSTKFieldRef = rhoRef().stkFieldRef();
    auto& muSTKFieldRef = muRef().stkFieldRef();
    auto& uTauSTKFieldRef = uTauRef().stkFieldRef();
    auto& yPlusSTKFieldRef = yPlusRef().stkFieldRef();

    const auto& coordsSTKFieldRef = *metaData.get_field<scalar>(
        stk::topology::NODE_RANK, this->getCoordinatesID_(domain));
    const auto& wallNormalDistanceSTKFieldRef = *metaData.get_field<scalar>(
        metaData.side_rank(), mesh::wall_normal_distance_ID);
    const auto& exposedAreaVecSTKFieldRef = *metaData.get_field<scalar>(
        metaData.side_rank(), this->getExposedAreaVectorID_(domain));

    // define vector of parent topos; should always
    // be UNITY in size
    std::vector<stk::topology> parentTopo;

    // define some common selectors: only those
    // sitting on no-slip walls
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
        const label nodesPerSide = sideBucket.topology().num_nodes();
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

        const stk::mesh::Bucket::size_type length = sideBucket.size();

        for (stk::mesh::Bucket::size_type k = 0; k < length; ++k)
        {
            // get face
            stk::mesh::Entity side = sideBucket[k];

            //======================================
            // gather nodal data off of face
            //======================================
            stk::mesh::Entity const* face_node_rels =
                bulkData.begin_nodes(side);
            label num_face_nodes = bulkData.num_nodes(side);

            // sanity check on num nodes
            STK_ThrowAssert(num_face_nodes == nodesPerSide);
            for (label ni = 0; ni < num_face_nodes; ++ni)
            {
                stk::mesh::Entity node = face_node_rels[ni];

                // gather scalars
                p_rho[ni] = *stk::mesh::field_data(rhoSTKFieldRef, node);
                p_mu[ni] = *stk::mesh::field_data(muSTKFieldRef, node);

                // gather vectors
                scalar* uNp1 = stk::mesh::field_data(USTKFieldRef, node);
                const label offSet = ni * SPATIAL_DIM;
                for (label j = 0; j < SPATIAL_DIM; ++j)
                {
                    p_U[offSet + j] = uNp1[j];
                }
            }

            // pointer to face data
            const scalar* areaVec =
                stk::mesh::field_data(exposedAreaVecSTKFieldRef, side);
            const scalar* yPlusBip =
                stk::mesh::field_data(yPlusSTKFieldRef, side);
            const scalar* UbcVec =
                stk::mesh::field_data(sideUSTKFieldRef, side);
            scalar* wallFrictionVelocityBip =
                stk::mesh::field_data(uTauSTKFieldRef, side);

            // extract the connected element to this
            // exposed face; should be single in
            // size!
            const stk::mesh::Entity* face_elem_rels =
                bulkData.begin_elements(side);
            STK_ThrowAssert(bulkData.num_elements(side) == 1);

            // get element; its face ordinal number
            stk::mesh::Entity element = face_elem_rels[0];
            const label face_ordinal = bulkData.begin_element_ordinals(side)[0];

            // get the relations off of element
            stk::mesh::Entity const* elem_node_rels =
                bulkData.begin_nodes(element);

            // loop over face nodes
            for (label ip = 0; ip < numScsBip; ++ip)
            {
                const label offSetAveraVec = ip * SPATIAL_DIM;

                const label opposingNode =
                    meSCS->opposingNodes(face_ordinal, ip);
                const label localFaceNode = faceIpNodeMap[ip];

                // left and right nodes; right is on
                // the face; left is the opposing
                // node
                stk::mesh::Entity nodeL = elem_node_rels[opposingNode];
                stk::mesh::Entity nodeR = face_node_rels[localFaceNode];

                // zero out vector quantities;
                // squeeze in aMag
                scalar aMag = 0.0;
                for (label j = 0; j < SPATIAL_DIM; ++j)
                {
                    p_uBip[j] = 0.0;
                    const scalar axj = areaVec[offSetAveraVec + j];
                    aMag += axj * axj;
                }
                aMag = std::sqrt(aMag);

                // interpolate to bip
                scalar rhoBip = 0.0;
                scalar muBip = 0.0;
                const label offSetSF_face = ip * nodesPerSide;
                for (label ic = 0; ic < nodesPerSide; ++ic)
                {
                    const scalar r = p_face_shape_function[offSetSF_face + ic];
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
                    const scalar nj = areaVec[offSetAveraVec + j] / aMag;
                    p_unitNormal[j] = nj;
                }

                // determine tangential velocity
                scalar uTangential = 0.0;
                for (label i = 0; i < SPATIAL_DIM; ++i)
                {
                    scalar uiTan = 0.0;
                    scalar uiBcTan = 0.0;
                    for (label j = 0; j < SPATIAL_DIM; ++j)
                    {
                        const scalar ninj = p_unitNormal[i] * p_unitNormal[j];
                        if (i == j)
                        {
                            const scalar om_nini = 1.0 - ninj;
                            uiTan += om_nini * p_uBip[j];
                            uiBcTan += om_nini * UbcVec[ip * SPATIAL_DIM + j];
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

                scalar uPlusLog = std::log(yPlusBip[ip]) / kappa_ + B_;
                scalar uTauLog = uTangential / uPlusLog;

                wallFrictionVelocityBip[ip] = uTauLog;
            }
        }
    }
}

void turbulenceModel::updateYStarAutomatic_(
    const std::shared_ptr<domain> domain)
{
    stk::mesh::BulkData& bulkData = this->meshRef().bulkDataRef();
    stk::mesh::MetaData& metaData = this->meshRef().metaDataRef();

    // Get fields
    STKScalarField* yStarSTKFieldPtr = yStarRef().stkFieldPtr();
    const STKScalarField* uStarSTKFieldPtr = uStarRef().stkFieldPtr();
    const STKScalarField* wallNormalDistanceSTKFieldPtr =
        metaData.get_field<scalar>(metaData.side_rank(),
                                   mesh::wall_normal_distance_ID);

    const auto& rhoSTKFieldRef = rhoRef().stkFieldRef();
    const auto& muSTKFieldRef = muRef().stkFieldRef();
    const auto& sideUSTKFieldRef = URef().sideFieldRef().stkFieldRef();

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

        const stk::mesh::Bucket::size_type nSidesPerBucket = sideBucket.size();

        for (stk::mesh::Bucket::size_type iSide = 0; iSide < nSidesPerBucket;
             ++iSide)
        {
            // get face
            stk::mesh::Entity side = sideBucket[iSide];
            label numSideNodes = bulkData.num_nodes(side);

            //======================================
            // gather nodal data off of face
            //======================================
            stk::mesh::Entity const* sideNodeRels = bulkData.begin_nodes(side);
            for (label ni = 0; ni < nodesPerSide; ++ni)
            {
                stk::mesh::Entity node = sideNodeRels[ni];

                // gather scalars
                p_rho[ni] = *stk::mesh::field_data(rhoSTKFieldRef, node);
                p_mu[ni] = *stk::mesh::field_data(muSTKFieldRef, node);
            }

            // pointer to face data
            const scalar* wallNormalDistanceBip =
                stk::mesh::field_data(*wallNormalDistanceSTKFieldPtr, side);
            const scalar* uStarBip =
                stk::mesh::field_data(*uStarSTKFieldPtr, side);
            scalar* yStarBip = stk::mesh::field_data(*yStarSTKFieldPtr, side);

            // loop over face nodes
            for (label ip = 0; ip < numSideNodes; ++ip)
            {
                const label offSetSF_face = ip * nodesPerSide;

                // interpolate to bip
                scalar muBip = 0.0;
                scalar rhoBip = 0.0;
                for (label ic = 0; ic < nodesPerSide; ++ic)
                {
                    const scalar r = p_face_shape_function[offSetSF_face + ic];
                    rhoBip += r * p_rho[ic];
                    muBip += r * p_mu[ic];
                }

                // extract bip data
                const scalar yp = wallNormalDistanceBip[ip];
                const scalar uStar = uStarBip[ip];

                // determine yplus
                yStarBip[ip] = rhoBip * yp * uStar / muBip;
            }
        }
    }
}

void turbulenceModel::updateUStarAutomatic_(
    const std::shared_ptr<domain> domain)
{
    stk::mesh::BulkData& bulkData = this->meshRef().bulkDataRef();
    stk::mesh::MetaData& metaData = this->meshRef().metaDataRef();

    // Get fields
    const auto& USTKFieldRef = URef().stkFieldRef();
    const auto& sideUSTKFieldRef = URef().sideFieldRef().stkFieldRef();
    STKScalarField* uStarSTKFieldPtr = uStarRef().stkFieldPtr();

    const auto& coordsSTKFieldRef = *metaData.get_field<scalar>(
        stk::topology::NODE_RANK, this->getCoordinatesID_(domain));
    const auto& exposedAreaVecSTKFieldRef = *metaData.get_field<scalar>(
        metaData.side_rank(), this->getExposedAreaVectorID_(domain));
    STKScalarField* wallNormalDistanceSTKFieldPtr = metaData.get_field<scalar>(
        metaData.side_rank(), mesh::wall_normal_distance_ID);

    const auto& kSTKFieldRef = kRef().stkFieldRef();
    const auto& rhoSTKFieldRef = rhoRef().stkFieldRef();
    const auto& muSTKFieldRef = muRef().stkFieldRef();

    // define vector of parent topos; should always
    // be UNITY in size
    std::vector<stk::topology> parentTopo;

    // bip values
    std::vector<scalar> uBip(SPATIAL_DIM);
    std::vector<scalar> unitNormal(SPATIAL_DIM);

    // pointers to fixed values
    scalar* p_uBip = &uBip[0];
    scalar* p_unitNormal = &unitNormal[0];

    // nodal fields to gather
    std::vector<scalar> ws_U;
    std::vector<scalar> ws_k;
    std::vector<scalar> ws_rho;
    std::vector<scalar> ws_mu;

    // master element
    std::vector<scalar> ws_face_shape_function;

    const scalar Cmu = this->Cmu();

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

        // face master element
        MasterElement* meFC = MasterElementRepo::get_surface_master_element(
            sideBucket.topology());
        const label nodesPerSide = meFC->nodesPerElement_;
        const label numScsBip = meFC->numIntPoints_;

        // algorithm related; element
        ws_U.resize(nodesPerSide * SPATIAL_DIM);
        ws_k.resize(nodesPerSide);
        ws_rho.resize(nodesPerSide);
        ws_mu.resize(nodesPerSide);
        ws_face_shape_function.resize(numScsBip * nodesPerSide);

        // pointers
        scalar* p_U = &ws_U[0];
        scalar* p_k = &ws_k[0];
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

        const stk::mesh::Bucket::size_type nSidesPerBucket = sideBucket.size();

        for (stk::mesh::Bucket::size_type iSide = 0; iSide < nSidesPerBucket;
             ++iSide)
        {
            // get face
            stk::mesh::Entity side = sideBucket[iSide];

            //======================================
            // gather nodal data off of face
            //======================================
            stk::mesh::Entity const* sideNodeRels = bulkData.begin_nodes(side);
            label numSideNodes = bulkData.num_nodes(side);

            // sanity check on num nodes
            STK_ThrowAssert(numSideNodes == nodesPerSide);
            for (label ni = 0; ni < nodesPerSide; ++ni)
            {
                stk::mesh::Entity node = sideNodeRels[ni];

                // gather scalars
                p_k[ni] = *stk::mesh::field_data(kSTKFieldRef, node);
                p_rho[ni] = *stk::mesh::field_data(rhoSTKFieldRef, node);
                p_mu[ni] = *stk::mesh::field_data(muSTKFieldRef, node);

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
                stk::mesh::field_data(exposedAreaVecSTKFieldRef, side);
            const scalar* wallNormalDistanceBip =
                stk::mesh::field_data(*wallNormalDistanceSTKFieldPtr, side);
            const scalar* UbcVec =
                stk::mesh::field_data(sideUSTKFieldRef, side);
            scalar* uStarBip = stk::mesh::field_data(*uStarSTKFieldPtr, side);

            // extract the connected element to this
            // exposed face; should be single in
            // size!
            const stk::mesh::Entity* face_elem_rels =
                bulkData.begin_elements(side);
            STK_ThrowAssert(bulkData.num_elements(side) == 1);

            // get element; its face ordinal number
            stk::mesh::Entity element = face_elem_rels[0];

            // loop over face nodes
            for (label ip = 0; ip < numSideNodes; ++ip)
            {
                const label offSetSF_face = ip * nodesPerSide;
                const label offSetAveraVec = ip * SPATIAL_DIM;

                // zero out vector quantities;
                // squeeze in aMag
                scalar aMag = 0.0;
                for (label j = 0; j < SPATIAL_DIM; ++j)
                {
                    p_uBip[j] = 0.0;
                    const scalar axj = areaVec[offSetAveraVec + j];
                    aMag += axj * axj;
                }
                aMag = std::sqrt(aMag);

                // interpolate to bip
                scalar kBip = 0.0;
                scalar muBip = 0.0;
                scalar rhoBip = 0.0;
                for (label ic = 0; ic < nodesPerSide; ++ic)
                {
                    const scalar r = p_face_shape_function[offSetSF_face + ic];
                    kBip += r * p_k[ic];
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
                    const scalar nj = areaVec[offSetAveraVec + j] / aMag;
                    p_unitNormal[j] = nj;
                }

                // determine tangential velocity
                scalar uTangential = 0.0;
                for (label i = 0; i < SPATIAL_DIM; ++i)
                {
                    scalar uiTan = 0.0;
                    scalar uiBcTan = 0.0;
                    for (label j = 0; j < SPATIAL_DIM; ++j)
                    {
                        const scalar ninj = p_unitNormal[i] * p_unitNormal[j];
                        if (i == j)
                        {
                            const scalar om_nini = 1.0 - ninj;
                            uiTan += om_nini * p_uBip[j];
                            uiBcTan += om_nini * UbcVec[ip * SPATIAL_DIM + j];
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

                // extract bip data
                const scalar yp = wallNormalDistanceBip[ip];

                // determine ustar
                scalar uStar_vis =
                    std::sqrt(muBip * uTangential / (rhoBip * yp));
                scalar uStar_log = std::pow(Cmu, 0.25) * sqrt(kBip);

                uStarBip[ip] = std::pow(
                    std::pow(uStar_vis, 4.0) + std::pow(uStar_log, 4.0), 0.25);
            }
        }
    }
}

void turbulenceModel::updateYPlusAutomatic_(
    const std::shared_ptr<domain> domain)
{
    stk::mesh::BulkData& bulkData = this->meshRef().bulkDataRef();
    stk::mesh::MetaData& metaData = this->meshRef().metaDataRef();

    // Get fields
    STKScalarField* yPlusSTKFieldPtr = yPlusRef().stkFieldPtr();
    STKScalarField* uStarSTKFieldPtr = uStarRef().stkFieldPtr();
    STKScalarField* wallNormalDistanceSTKFieldPtr = metaData.get_field<scalar>(
        metaData.side_rank(), mesh::wall_normal_distance_ID);

    const auto& rhoSTKFieldRef = rhoRef().stkFieldRef();
    const auto& muSTKFieldRef = muRef().stkFieldRef();
    const auto& sideUSTKFieldRef = URef().sideFieldRef().stkFieldRef();

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

        const stk::mesh::Bucket::size_type nSidesPerBucket = sideBucket.size();

        for (stk::mesh::Bucket::size_type iSide = 0; iSide < nSidesPerBucket;
             ++iSide)
        {
            // get face
            stk::mesh::Entity side = sideBucket[iSide];
            label numSideNodes = bulkData.num_nodes(side);

            //======================================
            // gather nodal data off of face
            //======================================
            stk::mesh::Entity const* sideNodeRels = bulkData.begin_nodes(side);
            for (label ni = 0; ni < nodesPerSide; ++ni)
            {
                stk::mesh::Entity node = sideNodeRels[ni];

                // gather scalars
                p_rho[ni] = *stk::mesh::field_data(rhoSTKFieldRef, node);
                p_mu[ni] = *stk::mesh::field_data(muSTKFieldRef, node);
            }

            // pointer to face data
            const scalar* wallNormalDistanceBip =
                stk::mesh::field_data(*wallNormalDistanceSTKFieldPtr, side);
            const scalar* uStarBip =
                stk::mesh::field_data(*uStarSTKFieldPtr, side);
            scalar* yPlusBip = stk::mesh::field_data(*yPlusSTKFieldPtr, side);

            // loop over face nodes
            for (label ip = 0; ip < numSideNodes; ++ip)
            {
                const label offSetSF_face = ip * nodesPerSide;

                // interpolate to bip
                scalar muBip = 0.0;
                scalar rhoBip = 0.0;
                for (label ic = 0; ic < nodesPerSide; ++ic)
                {
                    const scalar r = p_face_shape_function[offSetSF_face + ic];
                    rhoBip += r * p_rho[ic];
                    muBip += r * p_mu[ic];
                }

                // extract bip data
                const scalar yp = wallNormalDistanceBip[ip];
                const scalar uStar = uStarBip[ip];

                // determine yplus
                yPlusBip[ip] = rhoBip * yp * uStar / muBip;
            }
        }
    }
}

void turbulenceModel::updateUTauAutomatic_(const std::shared_ptr<domain> domain)
{
    stk::mesh::BulkData& bulkData = meshRef().bulkDataRef();
    stk::mesh::MetaData& metaData = meshRef().metaDataRef();

    // bip values
    std::vector<scalar> uBip(SPATIAL_DIM);
    std::vector<scalar> unitNormal(SPATIAL_DIM);

    // pointers to fixed values
    scalar* p_uBip = &uBip[0];
    scalar* p_unitNormal = &unitNormal[0];

    // nodal fields to gather
    std::vector<scalar> ws_U;
    std::vector<scalar> ws_rho;
    std::vector<scalar> ws_mu;

    // master element
    std::vector<scalar> ws_face_shape_function;

    // deal with state
    auto& USTKFieldRef = URef().stkFieldRef();
    auto& sideUSTKFieldRef = URef().sideFieldRef().stkFieldRef();
    auto& rhoSTKFieldRef = rhoRef().stkFieldRef();
    auto& muSTKFieldRef = muRef().stkFieldRef();
    auto& uTauSTKFieldRef = uTauRef().stkFieldRef();
    auto& yPlusSTKFieldRef = yPlusRef().stkFieldRef();

    const auto& coordsSTKFieldRef = *metaData.get_field<scalar>(
        stk::topology::NODE_RANK, this->getCoordinatesID_(domain));
    const auto& wallNormalDistanceSTKFieldRef = *metaData.get_field<scalar>(
        metaData.side_rank(), mesh::wall_normal_distance_ID);
    const auto& exposedAreaVecSTKFieldRef = *metaData.get_field<scalar>(
        metaData.side_rank(), this->getExposedAreaVectorID_(domain));

    // define vector of parent topos; should always
    // be UNITY in size
    std::vector<stk::topology> parentTopo;

    // define some common selectors: only those
    // sitting on no-slip walls
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
        const label nodesPerSide = sideBucket.topology().num_nodes();
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

        const stk::mesh::Bucket::size_type length = sideBucket.size();

        for (stk::mesh::Bucket::size_type k = 0; k < length; ++k)
        {
            // get face
            stk::mesh::Entity side = sideBucket[k];

            //======================================
            // gather nodal data off of face
            //======================================
            stk::mesh::Entity const* face_node_rels =
                bulkData.begin_nodes(side);
            label num_face_nodes = bulkData.num_nodes(side);

            // sanity check on num nodes
            STK_ThrowAssert(num_face_nodes == nodesPerSide);
            for (label ni = 0; ni < num_face_nodes; ++ni)
            {
                stk::mesh::Entity node = face_node_rels[ni];

                // gather scalars
                p_rho[ni] = *stk::mesh::field_data(rhoSTKFieldRef, node);
                p_mu[ni] = *stk::mesh::field_data(muSTKFieldRef, node);

                // gather vectors
                scalar* uNp1 = stk::mesh::field_data(USTKFieldRef, node);
                const label offSet = ni * SPATIAL_DIM;
                for (label j = 0; j < SPATIAL_DIM; ++j)
                {
                    p_U[offSet + j] = uNp1[j];
                }
            }

            // pointer to face data
            const scalar* areaVec =
                stk::mesh::field_data(exposedAreaVecSTKFieldRef, side);
            const scalar* wallNormalDistanceBip =
                stk::mesh::field_data(wallNormalDistanceSTKFieldRef, side);
            const scalar* yPlusBip =
                stk::mesh::field_data(yPlusSTKFieldRef, side);
            const scalar* UbcVec =
                stk::mesh::field_data(sideUSTKFieldRef, side);
            scalar* wallFrictionVelocityBip =
                stk::mesh::field_data(uTauSTKFieldRef, side);

            // extract the connected element to this
            // exposed face; should be single in
            // size!
            const stk::mesh::Entity* face_elem_rels =
                bulkData.begin_elements(side);
            STK_ThrowAssert(bulkData.num_elements(side) == 1);

            // get element; its face ordinal number
            stk::mesh::Entity element = face_elem_rels[0];
            const label face_ordinal = bulkData.begin_element_ordinals(side)[0];

            // get the relations off of element
            stk::mesh::Entity const* elem_node_rels =
                bulkData.begin_nodes(element);

            // loop over face nodes
            for (label ip = 0; ip < numScsBip; ++ip)
            {
                const label offSetAveraVec = ip * SPATIAL_DIM;

                const label opposingNode =
                    meSCS->opposingNodes(face_ordinal, ip);
                const label localFaceNode = faceIpNodeMap[ip];

                // left and right nodes; right is on
                // the face; left is the opposing
                // node
                stk::mesh::Entity nodeL = elem_node_rels[opposingNode];
                stk::mesh::Entity nodeR = face_node_rels[localFaceNode];

                // zero out vector quantities;
                // squeeze in aMag
                scalar aMag = 0.0;
                for (label j = 0; j < SPATIAL_DIM; ++j)
                {
                    p_uBip[j] = 0.0;
                    const scalar axj = areaVec[offSetAveraVec + j];
                    aMag += axj * axj;
                }
                aMag = std::sqrt(aMag);

                // interpolate to bip
                scalar rhoBip = 0.0;
                scalar muBip = 0.0;
                const label offSetSF_face = ip * nodesPerSide;
                for (label ic = 0; ic < nodesPerSide; ++ic)
                {
                    const scalar r = p_face_shape_function[offSetSF_face + ic];
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
                    const scalar nj = areaVec[offSetAveraVec + j] / aMag;
                    p_unitNormal[j] = nj;
                }

                // determine tangential velocity
                scalar uTangential = 0.0;
                for (label i = 0; i < SPATIAL_DIM; ++i)
                {
                    scalar uiTan = 0.0;
                    scalar uiBcTan = 0.0;
                    for (label j = 0; j < SPATIAL_DIM; ++j)
                    {
                        const scalar ninj = p_unitNormal[i] * p_unitNormal[j];
                        if (i == j)
                        {
                            const scalar om_nini = 1.0 - ninj;
                            uiTan += om_nini * p_uBip[j];
                            uiBcTan += om_nini * UbcVec[ip * SPATIAL_DIM + j];
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

                scalar uTauVis = sqrt(muBip * uTangential /
                                      (rhoBip * wallNormalDistanceBip[ip]));

                scalar yPlusMin = 0.17871;
                scalar uPlusLog = std::log(yPlusBip[ip]) / kappa_ + B_;
                scalar uPlusLogMin =
                    std::log(yPlusMin) / kappa_ + B_; // slightly more than one
                scalar uTauLog = uTangential / std::max(uPlusLog, uPlusLogMin);

                scalar uTau = std::pow(
                    std::pow(uTauVis, 4.0) + std::pow(uTauLog, 4.0), 0.25);

                wallFrictionVelocityBip[ip] = uTau;
            }
        }
    }
}

void turbulenceModel::updateTPlusScalable_(const std::shared_ptr<domain> domain)
{
    stk::mesh::BulkData& bulkData = this->meshRef().bulkDataRef();
    stk::mesh::MetaData& metaData = this->meshRef().metaDataRef();

    // Get fields
    const auto& cpSTKFieldRef = cpRef().stkFieldRef();
    const auto& muSTKFieldRef = muRef().stkFieldRef();
    const auto& lambdaSTKFieldRef = lambdaRef().stkFieldRef();

    const STKScalarField* yStarSTKFieldPtr = yStarRef().stkFieldPtr();
    STKScalarField* TPlusSTKFieldPtr = TPlusRef().stkFieldPtr();

    // define vector of parent topos; should always
    // be UNITY in size
    std::vector<stk::topology> parentTopo;

    // nodal fields to gather
    std::vector<scalar> ws_cp;
    std::vector<scalar> ws_mu;
    std::vector<scalar> ws_lambda;

    // master element
    std::vector<scalar> ws_face_shape_function;

    // select all sides: only those sitting on
    // no-slip walls
    stk::mesh::Selector selAllSides =
        metaData.universal_part() &
        stk::mesh::selectUnion(collectNoSlipWallParts_(domain));

    // shifted ip's for fields?
    const bool isTShifted = TRef().isShifted();

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

        ws_cp.resize(nodesPerSide);
        ws_mu.resize(nodesPerSide);
        ws_face_shape_function.resize(numScsBip * nodesPerSide);

        // pointers
        scalar* p_cp = &ws_cp[0];
        scalar* p_mu = &ws_mu[0];
        scalar* p_lambda = &ws_lambda[0];
        scalar* p_face_shape_function = &ws_face_shape_function[0];

        // shape functions
        if (isTShifted)
        {
            meFC->shifted_shape_fcn(&p_face_shape_function[0]);
        }
        else
        {
            meFC->shape_fcn(&p_face_shape_function[0]);
        }

        const stk::mesh::Bucket::size_type nSidesPerBucket = sideBucket.size();

        for (stk::mesh::Bucket::size_type iSide = 0; iSide < nSidesPerBucket;
             ++iSide)
        {
            // get face
            stk::mesh::Entity side = sideBucket[iSide];

            //======================================
            // gather nodal data off of face
            //======================================
            stk::mesh::Entity const* sideNodeRels = bulkData.begin_nodes(side);
            label numSideNodes = bulkData.num_nodes(side);

            // sanity check on num nodes
            STK_ThrowAssert(numSideNodes == nodesPerSide);
            for (label ni = 0; ni < nodesPerSide; ++ni)
            {
                stk::mesh::Entity node = sideNodeRels[ni];

                // gather scalars
                p_cp[ni] = *stk::mesh::field_data(cpSTKFieldRef, node);
                p_mu[ni] = *stk::mesh::field_data(muSTKFieldRef, node);
                p_lambda[ni] = *stk::mesh::field_data(lambdaSTKFieldRef, node);
            }

            // pointer to face data
            const scalar* yStarBip =
                stk::mesh::field_data(*yStarSTKFieldPtr, side);
            scalar* TPlusBip = stk::mesh::field_data(*TPlusSTKFieldPtr, side);

            // extract the connected element to this
            // exposed face; should be single in
            // size!
            const stk::mesh::Entity* face_elem_rels =
                bulkData.begin_elements(side);
            STK_ThrowAssert(bulkData.num_elements(side) == 1);

            // get element; its face ordinal number
            stk::mesh::Entity element = face_elem_rels[0];

            // loop over face nodes
            for (label ip = 0; ip < numSideNodes; ++ip)
            {
                const label offSetSF_face = ip * nodesPerSide;

                // interpolate to bip
                scalar cpBip = 0;
                scalar muBip = 0;
                scalar lambdaBip = 0;
                for (label ic = 0; ic < nodesPerSide; ++ic)
                {
                    const scalar r = p_face_shape_function[offSetSF_face + ic];

                    cpBip += r * p_cp[ic];
                    muBip += r * p_mu[ic];
                    lambdaBip += r * p_lambda[ic];
                }

                scalar Pr = cpBip * muBip / lambdaBip;
                scalar beta =
                    std::pow(3.85 * std::pow(Pr, 1.0 / 3.0) - 1.3, 2.0) +
                    2.12 * std::log(Pr);

                TPlusBip[ip] = 2.12 * std::log(yStarBip[ip]) + beta;
            }
        }
    }
}

void turbulenceModel::updateTPlusAutomatic_(
    const std::shared_ptr<domain> domain)
{
    stk::mesh::BulkData& bulkData = this->meshRef().bulkDataRef();
    stk::mesh::MetaData& metaData = this->meshRef().metaDataRef();

    // Get fields
    const auto& cpSTKFieldRef = cpRef().stkFieldRef();
    const auto& muSTKFieldRef = muRef().stkFieldRef();
    const auto& lambdaSTKFieldRef = lambdaRef().stkFieldRef();

    const STKScalarField* yStarSTKFieldPtr = yStarRef().stkFieldPtr();
    STKScalarField* TPlusSTKFieldPtr = TPlusRef().stkFieldPtr();

    // define vector of parent topos; should always
    // be UNITY in size
    std::vector<stk::topology> parentTopo;

    // nodal fields to gather
    std::vector<scalar> ws_cp;
    std::vector<scalar> ws_mu;
    std::vector<scalar> ws_lambda;

    // master element
    std::vector<scalar> ws_face_shape_function;

    // select all sides: only those sitting on
    // no-slip walls
    stk::mesh::Selector selAllSides =
        metaData.universal_part() &
        stk::mesh::selectUnion(collectNoSlipWallParts_(domain));

    // shifted ip's for fields?
    const bool isTShifted = TRef().isShifted();

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

        ws_cp.resize(nodesPerSide);
        ws_mu.resize(nodesPerSide);
        ws_lambda.resize(nodesPerSide);
        ws_face_shape_function.resize(numScsBip * nodesPerSide);

        // pointers
        scalar* p_cp = &ws_cp[0];
        scalar* p_mu = &ws_mu[0];
        scalar* p_lambda = &ws_lambda[0];
        scalar* p_face_shape_function = &ws_face_shape_function[0];

        // shape functions
        if (isTShifted)
        {
            meFC->shifted_shape_fcn(&p_face_shape_function[0]);
        }
        else
        {
            meFC->shape_fcn(&p_face_shape_function[0]);
        }

        const stk::mesh::Bucket::size_type nSidesPerBucket = sideBucket.size();

        for (stk::mesh::Bucket::size_type iSide = 0; iSide < nSidesPerBucket;
             ++iSide)
        {
            // get face
            stk::mesh::Entity side = sideBucket[iSide];

            //======================================
            // gather nodal data off of face
            //======================================
            stk::mesh::Entity const* sideNodeRels = bulkData.begin_nodes(side);
            label numSideNodes = bulkData.num_nodes(side);

            // sanity check on num nodes
            STK_ThrowAssert(numSideNodes == nodesPerSide);
            for (label ni = 0; ni < nodesPerSide; ++ni)
            {
                stk::mesh::Entity node = sideNodeRels[ni];

                // gather scalars
                p_cp[ni] = *stk::mesh::field_data(cpSTKFieldRef, node);
                p_mu[ni] = *stk::mesh::field_data(muSTKFieldRef, node);
                p_lambda[ni] = *stk::mesh::field_data(lambdaSTKFieldRef, node);
            }

            // pointer to face data
            const scalar* yStarBip =
                stk::mesh::field_data(*yStarSTKFieldPtr, side);
            scalar* TPlusBip = stk::mesh::field_data(*TPlusSTKFieldPtr, side);

            // extract the connected element to this
            // exposed face; should be single in
            // size!
            const stk::mesh::Entity* face_elem_rels =
                bulkData.begin_elements(side);
            STK_ThrowAssert(bulkData.num_elements(side) == 1);

            // get element; its face ordinal number
            stk::mesh::Entity element = face_elem_rels[0];

            // loop over face nodes
            for (label ip = 0; ip < numSideNodes; ++ip)
            {
                const label offSetSF_face = ip * nodesPerSide;

                // interpolate to bip
                scalar cpBip = 0;
                scalar muBip = 0;
                scalar lambdaBip = 0;
                for (label ic = 0; ic < nodesPerSide; ++ic)
                {
                    const scalar r = p_face_shape_function[offSetSF_face + ic];

                    cpBip += r * p_cp[ic];
                    muBip += r * p_mu[ic];
                    lambdaBip += r * p_lambda[ic];
                }

                scalar Pr = cpBip * muBip / lambdaBip;
                scalar beta =
                    std::pow(3.85 * std::pow(Pr, 1.0 / 3.0) - 1.3, 2.0) +
                    2.12 * std::log(Pr);
                scalar Gamma = 0.01 * std::pow(Pr * yStarBip[ip], 4.0) /
                               (1.5 * pow(Pr, 3.0) * yStarBip[ip]);

                TPlusBip[ip] = Pr * yStarBip[ip] * exp(-Gamma) +
                               (2.12 * std::log(yStarBip[ip]) + beta) *
                                   std::exp(-1.0 / Gamma);
            }
        }
    }
}

} /* namespace accel */
