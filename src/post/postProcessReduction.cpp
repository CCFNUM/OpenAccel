// File : postProcessReduction.cpp
// Created : Tue Aug 05 2025 19:49:24 (+0100)
// Author : Mhamad Mahdi Alloush
// Description:
// Copyright (c) 2025 CCFNUM, Lucerne University of Applied Sciences and Arts.
// SPDX-License-Identifier: BSD-3-Clause

#include "postProcess.h"

namespace accel
{

reductionObject::reductionObject(postProcess* postProcessPtr,
                                 std::string name,
                                 postProcessType type,
                                 std::vector<std::string> location,
                                 label frequency,
                                 bool writeToFile,
                                 reductionType rType,
                                 std::string field)
    : postProcessObject(postProcessPtr,
                        name,
                        type,
                        location,
                        frequency,
                        writeToFile),
      reductionType_(rType), field_(field)
{
    if (messager::master())
    {
        const stk::mesh::MetaData& metaData =
            postProcessPtr_->meshRef().metaDataRef();
        const stk::mesh::BulkData& bulkData =
            postProcessPtr_->meshRef().bulkDataRef();

        if (writeToFile_)
        {
            if (fs::is_regular_file(postProcessPtr->directory() / name_))
            {
                fs::remove_all(postProcessPtr->directory() / name_);
            }
            std::string fileName(postProcessPtr->directory() / name_);
            std::ofstream file(fileName);

            // Check if a node field
            const STKScalarField* STKFieldPtr =
                metaData.get_field<scalar>(stk::topology::NODE_RANK, field_);

            if (STKFieldPtr)
            {
                file << postProcessPtr_->instanceHeader() + "\t";

                if (STKFieldPtr->max_size() == 1)
                {
                    file << "value";
                }
                else
                {
                    for (label i = 0; i < STKFieldPtr->max_size(); i++)
                    {
                        file << "value_" << std::to_string(i) << "\t";
                    }
                }

                file << "\n";
                file.close();
            }
            else
            {
                // Check if element field
                const STKScalarField* elementSTKFieldPtr =
                    metaData.get_field<scalar>(stk::topology::ELEMENT_RANK,
                                               field_);

                if (elementSTKFieldPtr)
                {
                    // a header will be made at execution: for
                    // now just close the file
                    file.close();
                }
                else
                {
                    std::string msg = field_ + " not valid for post-process. "
                                               "Available fields:\n";
                    const auto& fields = metaData.get_fields();
                    for (auto field : fields)
                    {
                        msg += field->name() + "\n";
                    }

                    errorMsg(msg);
                }
            }
        }
        else
        {
            // Check if a node field
            const STKScalarField* STKFieldPtr =
                metaData.get_field<scalar>(stk::topology::NODE_RANK, field_);

            if (STKFieldPtr)
            {
                // do nothing ...
            }
            else
            {
                // Check if element field
                const STKScalarField* elementSTKFieldPtr =
                    metaData.get_field<scalar>(stk::topology::ELEMENT_RANK,
                                               field_);

                if (elementSTKFieldPtr)
                {
                    // do nothing ...
                }
                else
                {
                    std::string msg = field_ + " not valid for post-process. "
                                               "Available fields:\n";
                    const auto& fields = metaData.get_fields();
                    for (auto field : fields)
                    {
                        msg += field->name() + "\n";
                    }

                    errorMsg(msg);
                }
            }
        }
    }
}

void reductionObject::update()
{
    const stk::mesh::MetaData& metaData =
        postProcessPtr_->meshRef().metaDataRef();
    const stk::mesh::BulkData& bulkData =
        postProcessPtr_->meshRef().bulkDataRef();

    switch (reductionType_)
    {
        case reductionType::sum:
            {
                // Check if a node field
                const STKScalarField* STKFieldPtr = metaData.get_field<scalar>(
                    stk::topology::NODE_RANK, field_);

                if (STKFieldPtr)
                {
                    // Collect parts
                    stk::mesh::PartVector partVec;
                    for (const auto& locationName : location_)
                    {
                        stk::mesh::Part* part = metaData.get_part(locationName);

                        partVec.push_back(part);
                    }

                    label fieldDim = STKFieldPtr->max_size();

                    std::vector<scalar> total(fieldDim, 0.0);

                    const stk::mesh::Selector selOwnedNodes =
                        metaData.locally_owned_part() &
                        stk::mesh::selectUnion(partVec);
                    const stk::mesh::BucketVector& nodeBuckets =
                        bulkData.get_buckets(stk::topology::NODE_RANK,
                                             selOwnedNodes);

                    for (stk::mesh::Bucket::size_type ib = 0;
                         ib < nodeBuckets.size();
                         ++ib)
                    {
                        const stk::mesh::Bucket& bucket = *nodeBuckets[ib];
                        const stk::mesh::Bucket::size_type n_entities =
                            bucket.size();

                        for (stk::mesh::Bucket::size_type i = 0; i < n_entities;
                             ++i)
                        {
                            const stk::mesh::Entity node = bucket[i];
                            const scalar* val =
                                stk::mesh::field_data(*STKFieldPtr, node);

                            for (label i = 0; i < fieldDim; i++)
                            {
                                total[i] += val[i];
                            }
                        }
                    }

                    if (messager::parallel())
                    {
                        messager::sumReduce(total);
                    }

                    // print
                    if (messager::master())
                    {
                        std::cout << "Object name: " << name_ << ", value: (";
                        std::cout << total[0];
                        for (label i = 1; i < fieldDim; i++)
                        {
                            std::cout << "\t" << total[i];
                        }
                        std::cout << ")" << std::endl;

                        if (writeToFile_)
                        {
                            std::string fileName(postProcessPtr_->directory() /
                                                 name_);
                            std::ofstream file(fileName, std::ios_base::app);
                            file << postProcessPtr_->instance() << "\t";
                            for (label i = 0; i < fieldDim; i++)
                            {
                                file << total[i] << "\t";
                            }
                            file << "\n";
                            file.close();
                        }
                    }
                }
                else
                {
                    // Check if a side field
                    const STKScalarField* sideSTKFieldPtr =
                        metaData.get_field<scalar>(metaData.side_rank(),
                                                   field_ + "_side");

                    if (sideSTKFieldPtr)
                    {
                        // Collect parts
                        stk::mesh::PartVector partVec;
                        for (const auto& locationName : location_)
                        {
                            stk::mesh::Part* part =
                                metaData.get_part(locationName);

                            for (stk::mesh::Part* subPart : part->subsets())
                            {
                                // check if field registered
                                // for all required parts
                                if (!sideSTKFieldPtr->defined_on(*subPart))
                                {
                                    errorMsg(field_ +
                                             "_side is not "
                                             "defined "
                                             "for part " +
                                             subPart->name());
                                }
                                else
                                {
                                    partVec.push_back(subPart);
                                }
                            }
                        }

                        // max dim
                        label fieldDim = sideSTKFieldPtr->max_size();

                        const stk::mesh::Selector selOwnedSides =
                            metaData.locally_owned_part() &
                            stk::mesh::selectUnion(partVec);
                        const stk::mesh::BucketVector& sideBuckets =
                            bulkData.get_buckets(metaData.side_rank(),
                                                 selOwnedSides);

                        // get real field dimensions
                        {
                            // first find maximum number of
                            // intergration points per side
                            label maxNbip = 1;
                            for (stk::mesh::Bucket::size_type ib = 0;
                                 ib < sideBuckets.size();
                                 ++ib)
                            {
                                const stk::mesh::Bucket& bucket =
                                    *sideBuckets[ib];
                                const MasterElement* meFC = MasterElementRepo::
                                    get_surface_master_element(
                                        bucket.topology());
                                const label n_bip = meFC->numIntPoints_;

                                maxNbip = std::max(maxNbip, n_bip);
                            }

                            // max across partitions
                            messager::maxReduce(maxNbip);

                            // get real field dim
                            fieldDim /= maxNbip;
                        }

                        std::vector<scalar> total(fieldDim, 0.0);

                        for (stk::mesh::Bucket::size_type ib = 0;
                             ib < sideBuckets.size();
                             ++ib)
                        {
                            const stk::mesh::Bucket& bucket = *sideBuckets[ib];
                            const stk::mesh::Bucket::size_type n_entities =
                                bucket.size();
                            const MasterElement* meFC =
                                MasterElementRepo::get_surface_master_element(
                                    bucket.topology());
                            const label n_bip = meFC->numIntPoints_;
                            for (stk::mesh::Bucket::size_type i = 0;
                                 i < n_entities;
                                 ++i)
                            {
                                const stk::mesh::Entity side = bucket[i];
                                const scalar* val = stk::mesh::field_data(
                                    *sideSTKFieldPtr, side);

                                for (label bip = 0; bip < n_bip; bip++)
                                {
                                    for (label i = 0; i < fieldDim; i++)
                                    {
                                        total[i] += val[bip * fieldDim + i];
                                    }
                                }
                            }
                        }

                        if (messager::parallel())
                        {
                            messager::sumReduce(total);
                        }

                        // print
                        if (messager::master())
                        {
                            std::cout << "Object name: " << name_
                                      << ", value: (";
                            std::cout << total[0];
                            for (label i = 1; i < fieldDim; i++)
                            {
                                std::cout << "\t" << total[i];
                            }
                            std::cout << ")" << std::endl;

                            if (writeToFile_)
                            {
                                std::string fileName(
                                    postProcessPtr_->directory() / name_);
                                std::ofstream file(fileName,
                                                   std::ios_base::app);

                                // write header
                                if (postProcessPtr_->iter() == 1)
                                {
                                    file << postProcessPtr_->instanceHeader() +
                                                "\t";

                                    if (fieldDim == 1)
                                    {
                                        file << "value";
                                    }
                                    else
                                    {
                                        for (label i = 0; i < fieldDim; i++)
                                        {
                                            file << "value_"
                                                 << std::to_string(i) << "\t";
                                        }
                                    }

                                    file << "\n";
                                }

                                file << postProcessPtr_->instance() << "\t";
                                for (label i = 0; i < fieldDim; i++)
                                {
                                    file << total[i] << "\t";
                                }
                                file << "\n";
                                file.close();
                            }
                        }
                    }
                    else
                    {
                        errorMsg("A side field with the name " +
                                 sideSTKFieldPtr->name() + " does not exist");
                    }
                }
            }
            break;

        case reductionType::areaAverage:
            {
                // Check if a node field
                const STKScalarField* STKFieldPtr = metaData.get_field<scalar>(
                    stk::topology::NODE_RANK, field_);

                if (STKFieldPtr)
                {
                    STKScalarField* coordinates = metaData.get_field<scalar>(
                        stk::topology::NODE_RANK,
                        postProcessPtr_->meshRef().getCoordinateFieldName());

                    // Collect parts
                    stk::mesh::PartVector partVec;
                    for (const auto& locationName : location_)
                    {
                        stk::mesh::Part* part = metaData.get_part(locationName);

                        partVec.push_back(part);
                    }

                    label fieldDim = STKFieldPtr->max_size();

                    std::vector<scalar> average(fieldDim, 0.0);

                    // ip values; both boundary and opposing
                    // surface
                    std::vector<scalar> phiBip(fieldDim);

                    // pointers to fixed values
                    scalar* p_phiBip = &phiBip[0];

                    // nodal fields to gather
                    std::vector<scalar> ws_phi;

                    // master element
                    std::vector<scalar> ws_shape_function;

                    const stk::mesh::Selector selOwnedSides =
                        metaData.locally_owned_part() &
                        stk::mesh::selectUnion(partVec);
                    const stk::mesh::BucketVector& sideBuckets =
                        bulkData.get_buckets(metaData.side_rank(),
                                             selOwnedSides);

                    scalar totalBoundaryArea = 0.0;

                    for (stk::mesh::Bucket::size_type ib = 0;
                         ib < sideBuckets.size();
                         ++ib)
                    {
                        const stk::mesh::Bucket& bucket = *sideBuckets[ib];

                        // face master element
                        MasterElement* meFC =
                            MasterElementRepo::get_surface_master_element(
                                bucket.topology());
                        const label nodesPerSide =
                            bucket.topology().num_nodes();
                        const label numScsBip = meFC->numIntPoints_;

                        // define scratch field
                        std::vector<scalar> ws_coordinates(nodesPerSide *
                                                           SPATIAL_DIM);
                        std::vector<scalar> ws_scs_areav(numScsBip *
                                                         SPATIAL_DIM);

                        // algorithm related; element
                        // (exposed face and element)
                        ws_phi.resize(nodesPerSide * fieldDim);
                        ws_shape_function.resize(numScsBip * nodesPerSide);

                        // pointers
                        scalar* p_phi = &ws_phi[0];
                        scalar* p_face_shape_function = &ws_shape_function[0];

                        // shape functions; boundary
                        meFC->shape_fcn(&p_face_shape_function[0]);

                        const stk::mesh::Bucket::size_type nSidesPerBucket =
                            bucket.size();

                        for (stk::mesh::Bucket::size_type iSide = 0;
                             iSide < nSidesPerBucket;
                             ++iSide)
                        {
                            // get face
                            stk::mesh::Entity side = bucket[iSide];

                            //======================================
                            // gather nodal data off of face
                            //======================================
                            stk::mesh::Entity const* sideNodeRels =
                                bulkData.begin_nodes(side);
                            label numSideNodes = bulkData.num_nodes(side);

                            // sanity check on num nodes
                            STK_ThrowAssert(numSideNodes == nodesPerSide);
                            for (label ni = 0; ni < numSideNodes; ++ni)
                            {
                                stk::mesh::Entity node = sideNodeRels[ni];

                                // gather vectors
                                scalar* phi =
                                    stk::mesh::field_data(*STKFieldPtr, node);

                                const label offSet = ni * fieldDim;
                                for (label j = 0; j < fieldDim; ++j)
                                {
                                    p_phi[offSet + j] = phi[j];
                                }

                                scalar* coords =
                                    stk::mesh::field_data(*coordinates, node);
                                for (label j = 0; j < SPATIAL_DIM; ++j)
                                {
                                    ws_coordinates[ni * SPATIAL_DIM + j] =
                                        coords[j];
                                }
                            }

                            // compute scs integration point
                            // areavec
                            scalar scs_error = 0.0;
                            meFC->determinant(1,
                                              &ws_coordinates[0],
                                              &ws_scs_areav[0],
                                              &scs_error);

                            // loop over boundary ips
                            for (label ip = 0; ip < numScsBip; ++ip)
                            {
                                // calc vector quantities
                                scalar asq = 0.0;
                                for (label j = 0; j < SPATIAL_DIM; ++j)
                                {
                                    const scalar axj =
                                        ws_scs_areav[ip * SPATIAL_DIM + j];
                                    asq += axj * axj;
                                }
                                const scalar amag = std::sqrt(asq);

                                // interpolate to scs point;
                                // operate on saved off
                                // ws_field
                                for (label j = 0; j < fieldDim; ++j)
                                {
                                    p_phiBip[j] = 0.0;
                                }

                                // interpolate to bip
                                const label offSetSF_face = ip * nodesPerSide;
                                for (label ic = 0; ic < nodesPerSide; ++ic)
                                {
                                    const scalar r =
                                        p_face_shape_function[offSetSF_face +
                                                              ic];

                                    const label icNdim = ic * fieldDim;
                                    for (label j = 0; j < fieldDim; ++j)
                                    {
                                        p_phiBip[j] += r * p_phi[icNdim + j];
                                    }
                                }

                                // calculate average
                                for (label i = 0; i < fieldDim; i++)
                                {
                                    average[i] += amag * p_phiBip[i];
                                }

                                totalBoundaryArea += amag;
                            }
                        }
                    }

                    if (messager::parallel())
                    {
                        messager::sumReduce(totalBoundaryArea);
                        messager::sumReduce(average);
                    }

                    // calculate area average
                    for (label i = 0; i < fieldDim; i++)
                    {
                        average[i] /= totalBoundaryArea;
                    }

                    // print
                    if (messager::master())
                    {
                        std::cout << "Object name: " << name_ << ", value: (";
                        std::cout << average[0];
                        for (label i = 1; i < fieldDim; i++)
                        {
                            std::cout << "\t" << average[i];
                        }
                        std::cout << ")" << std::endl;

                        if (writeToFile_)
                        {
                            std::string fileName(postProcessPtr_->directory() /
                                                 name_);
                            std::ofstream file(fileName, std::ios_base::app);
                            file << postProcessPtr_->instance() << "\t";
                            for (label i = 0; i < fieldDim; i++)
                            {
                                file << average[i] << "\t";
                            }
                            file << "\n";
                            file.close();
                        }
                    }
                }
                else
                {
                    // Check if a side field
                    const STKScalarField* sideSTKFieldPtr =
                        metaData.get_field<scalar>(metaData.side_rank(),
                                                   field_ + "_side");

                    if (sideSTKFieldPtr)
                    {
                        STKScalarField* coordinates =
                            metaData.get_field<scalar>(
                                stk::topology::NODE_RANK,
                                postProcessPtr_->meshRef()
                                    .getCoordinateFieldName());

                        // Collect parts
                        stk::mesh::PartVector partVec;
                        for (const auto& locationName : location_)
                        {
                            stk::mesh::Part* part =
                                metaData.get_part(locationName);

                            for (stk::mesh::Part* subPart : part->subsets())
                            {
                                // check if field registered
                                // for all required parts
                                if (!sideSTKFieldPtr->defined_on(*subPart))
                                {
                                    errorMsg(field_ +
                                             "_side is not "
                                             "defined "
                                             "for part " +
                                             subPart->name());
                                }
                                else
                                {
                                    partVec.push_back(subPart);
                                }
                            }
                        }

                        label fieldDim = sideSTKFieldPtr->max_size();

                        const stk::mesh::Selector selOwnedSides =
                            metaData.locally_owned_part() &
                            stk::mesh::selectUnion(partVec);
                        const stk::mesh::BucketVector& sideBuckets =
                            bulkData.get_buckets(metaData.side_rank(),
                                                 selOwnedSides);

                        // get real field dimensions
                        {
                            // first find maximum number of
                            // intergration points per side
                            label maxNbip = 1;
                            for (stk::mesh::Bucket::size_type ib = 0;
                                 ib < sideBuckets.size();
                                 ++ib)
                            {
                                const stk::mesh::Bucket& bucket =
                                    *sideBuckets[ib];
                                const MasterElement* meFC = MasterElementRepo::
                                    get_surface_master_element(
                                        bucket.topology());
                                const label n_bip = meFC->numIntPoints_;

                                maxNbip = std::max(maxNbip, n_bip);
                            }

                            // max across partitions
                            messager::maxReduce(maxNbip);

                            // get real field dim
                            fieldDim /= maxNbip;
                        }

                        std::vector<scalar> average(fieldDim, 0.0);

                        // master element
                        std::vector<scalar> ws_shape_function;

                        scalar totalBoundaryArea = 0.0;

                        for (stk::mesh::Bucket::size_type ib = 0;
                             ib < sideBuckets.size();
                             ++ib)
                        {
                            const stk::mesh::Bucket& bucket = *sideBuckets[ib];

                            // face master element
                            MasterElement* meFC =
                                MasterElementRepo::get_surface_master_element(
                                    bucket.topology());
                            const label nodesPerSide =
                                bucket.topology().num_nodes();
                            const label numScsBip = meFC->numIntPoints_;

                            // define scratch field
                            std::vector<scalar> ws_coordinates(nodesPerSide *
                                                               SPATIAL_DIM);
                            std::vector<scalar> ws_scs_areav(numScsBip *
                                                             SPATIAL_DIM);

                            for (stk::mesh::Bucket::size_type i = 0;
                                 i < sideBuckets.size();
                                 ++i)
                            {
                                const stk::mesh::Entity side = bucket[i];

                                //======================================
                                // gather nodal data off of
                                // face
                                //======================================
                                stk::mesh::Entity const* sideNodeRels =
                                    bulkData.begin_nodes(side);
                                label numSideNodes = bulkData.num_nodes(side);

                                // sanity check on num nodes
                                STK_ThrowAssert(numSideNodes == nodesPerSide);
                                for (label ni = 0; ni < numSideNodes; ++ni)
                                {
                                    stk::mesh::Entity node = sideNodeRels[ni];

                                    scalar* coords = stk::mesh::field_data(
                                        *coordinates, node);
                                    for (label j = 0; j < SPATIAL_DIM; ++j)
                                    {
                                        ws_coordinates[ni * SPATIAL_DIM + j] =
                                            coords[j];
                                    }
                                }

                                // compute scs integration
                                // point areavec
                                scalar scs_error = 0.0;
                                meFC->determinant(1,
                                                  &ws_coordinates[0],
                                                  &ws_scs_areav[0],
                                                  &scs_error);

                                const scalar* val = stk::mesh::field_data(
                                    *sideSTKFieldPtr, side);

                                for (label ip = 0; ip < numScsBip; ip++)
                                {
                                    // calc vector
                                    // quantities
                                    scalar asq = 0.0;
                                    for (label j = 0; j < SPATIAL_DIM; ++j)
                                    {
                                        const scalar axj =
                                            ws_scs_areav[ip * SPATIAL_DIM + j];
                                        asq += axj * axj;
                                    }
                                    const scalar amag = std::sqrt(asq);

                                    for (label i = 0; i < fieldDim; i++)
                                    {
                                        average[i] +=
                                            amag * val[ip * fieldDim + i];
                                    }

                                    totalBoundaryArea += amag;
                                }
                            }
                        }

                        if (messager::parallel())
                        {
                            messager::sumReduce(totalBoundaryArea);
                            messager::sumReduce(average);
                        }

                        // calculate area average
                        for (label i = 0; i < fieldDim; i++)
                        {
                            average[i] /= totalBoundaryArea;
                        }

                        // print
                        if (messager::master())
                        {
                            std::cout << "Object name: " << name_
                                      << ", value: (";
                            std::cout << average[0];
                            for (label i = 1; i < fieldDim; i++)
                            {
                                std::cout << "\t" << average[i];
                            }
                            std::cout << ")" << std::endl;

                            if (writeToFile_)
                            {
                                std::string fileName(
                                    postProcessPtr_->directory() / name_);
                                std::ofstream file(fileName,
                                                   std::ios_base::app);

                                // write header
                                if (postProcessPtr_->iter() == 1)
                                {
                                    file << postProcessPtr_->instanceHeader() +
                                                "\t";

                                    if (fieldDim == 1)
                                    {
                                        file << "value";
                                    }
                                    else
                                    {
                                        for (label i = 0; i < fieldDim; i++)
                                        {
                                            file << "value_"
                                                 << std::to_string(i) << "\t";
                                        }
                                    }

                                    file << "\n";
                                }

                                file << postProcessPtr_->instance() << "\t";
                                for (label i = 0; i < fieldDim; i++)
                                {
                                    file << average[i] << "\t";
                                }
                                file << "\n";
                                file.close();
                            }
                        }
                    }
                    else
                    {
                        errorMsg("A side field with the name " +
                                 sideSTKFieldPtr->name() + " does not exist");
                    }
                }
            }
            break;

        case reductionType::average:
            {
                // Check if a node field
                const STKScalarField* STKFieldPtr = metaData.get_field<scalar>(
                    stk::topology::NODE_RANK, field_);

                if (STKFieldPtr)
                {
                    // Collect parts
                    stk::mesh::PartVector partVec;
                    for (const auto& locationName : location_)
                    {
                        stk::mesh::Part* part = metaData.get_part(locationName);

                        partVec.push_back(part);
                    }

                    label fieldDim = STKFieldPtr->max_size();

                    std::vector<scalar> average(fieldDim, 0.0);

                    // ip values; both boundary and opposing
                    // surface
                    std::vector<scalar> phiBip(fieldDim);

                    // pointers to fixed values
                    scalar* p_phiBip = &phiBip[0];

                    // nodal fields to gather
                    std::vector<scalar> ws_phi;

                    // master element
                    std::vector<scalar> ws_shape_function;

                    const stk::mesh::Selector selOwnedSides =
                        metaData.locally_owned_part() &
                        stk::mesh::selectUnion(partVec);
                    const stk::mesh::BucketVector& sideBuckets =
                        bulkData.get_buckets(metaData.side_rank(),
                                             selOwnedSides);

                    label totalIntPoints = 0.0;

                    for (stk::mesh::Bucket::size_type ib = 0;
                         ib < sideBuckets.size();
                         ++ib)
                    {
                        const stk::mesh::Bucket& bucket = *sideBuckets[ib];

                        // face master element
                        MasterElement* meFC =
                            MasterElementRepo::get_surface_master_element(
                                bucket.topology());
                        const label nodesPerSide =
                            bucket.topology().num_nodes();
                        const label numScsBip = meFC->numIntPoints_;

                        // algorithm related; element
                        // (exposed face and element)
                        ws_phi.resize(nodesPerSide * fieldDim);
                        ws_shape_function.resize(numScsBip * nodesPerSide);

                        // pointers
                        scalar* p_phi = &ws_phi[0];
                        scalar* p_face_shape_function = &ws_shape_function[0];

                        // shape functions; boundary
                        meFC->shape_fcn(&p_face_shape_function[0]);

                        const stk::mesh::Bucket::size_type nSidesPerBucket =
                            bucket.size();

                        for (stk::mesh::Bucket::size_type iSide = 0;
                             iSide < nSidesPerBucket;
                             ++iSide)
                        {
                            // get face
                            stk::mesh::Entity side = bucket[iSide];

                            //======================================
                            // gather nodal data off of face
                            //======================================
                            stk::mesh::Entity const* sideNodeRels =
                                bulkData.begin_nodes(side);
                            label numSideNodes = bulkData.num_nodes(side);

                            // sanity check on num nodes
                            STK_ThrowAssert(numSideNodes == nodesPerSide);
                            for (label ni = 0; ni < numSideNodes; ++ni)
                            {
                                stk::mesh::Entity node = sideNodeRels[ni];

                                // gather vectors
                                scalar* phi =
                                    stk::mesh::field_data(*STKFieldPtr, node);

                                const label offSet = ni * fieldDim;
                                for (label j = 0; j < fieldDim; ++j)
                                {
                                    p_phi[offSet + j] = phi[j];
                                }
                            }

                            // loop over boundary ips
                            for (label ip = 0; ip < numScsBip; ++ip)
                            {
                                // interpolate to scs point;
                                // operate on saved off
                                // ws_field
                                for (label j = 0; j < fieldDim; ++j)
                                {
                                    p_phiBip[j] = 0.0;
                                }

                                // interpolate to bip
                                const label offSetSF_face = ip * nodesPerSide;
                                for (label ic = 0; ic < nodesPerSide; ++ic)
                                {
                                    const scalar r =
                                        p_face_shape_function[offSetSF_face +
                                                              ic];

                                    const label icNdim = ic * fieldDim;
                                    for (label j = 0; j < fieldDim; ++j)
                                    {
                                        p_phiBip[j] += r * p_phi[icNdim + j];
                                    }
                                }

                                // calculate average
                                for (label i = 0; i < fieldDim; i++)
                                {
                                    average[i] += p_phiBip[i];
                                }

                                // increment
                                totalIntPoints++;
                            }
                        }
                    }

                    if (messager::parallel())
                    {
                        messager::sumReduce(totalIntPoints);
                        messager::sumReduce(average);
                    }

                    // calculate average
                    for (label i = 0; i < fieldDim; i++)
                    {
                        average[i] /= static_cast<scalar>(totalIntPoints);
                    }

                    // print
                    if (messager::master())
                    {
                        std::cout << "Object name: " << name_ << ", value: (";
                        std::cout << average[0];
                        for (label i = 1; i < fieldDim; i++)
                        {
                            std::cout << "\t" << average[i];
                        }
                        std::cout << ")" << std::endl;

                        if (writeToFile_)
                        {
                            std::string fileName(postProcessPtr_->directory() /
                                                 name_);
                            std::ofstream file(fileName, std::ios_base::app);
                            file << postProcessPtr_->instance() << "\t";
                            for (label i = 0; i < fieldDim; i++)
                            {
                                file << average[i] << "\t";
                            }
                            file << "\n";
                            file.close();
                        }
                    }
                }
                else
                {
                    // Check if a side field
                    const STKScalarField* sideSTKFieldPtr =
                        metaData.get_field<scalar>(metaData.side_rank(),
                                                   field_ + "_side");

                    if (sideSTKFieldPtr)
                    {
                        // Collect parts
                        stk::mesh::PartVector partVec;
                        for (const auto& locationName : location_)
                        {
                            stk::mesh::Part* part =
                                metaData.get_part(locationName);

                            for (stk::mesh::Part* subPart : part->subsets())
                            {
                                // check if field registered
                                // for all required parts
                                if (!sideSTKFieldPtr->defined_on(*subPart))
                                {
                                    errorMsg(field_ +
                                             "_side is not "
                                             "defined "
                                             "for part " +
                                             subPart->name());
                                }
                                else
                                {
                                    partVec.push_back(subPart);
                                }
                            }
                        }

                        label fieldDim = sideSTKFieldPtr->max_size();

                        const stk::mesh::Selector selOwnedSides =
                            metaData.locally_owned_part() &
                            stk::mesh::selectUnion(partVec);
                        const stk::mesh::BucketVector& sideBuckets =
                            bulkData.get_buckets(metaData.side_rank(),
                                                 selOwnedSides);

                        // get real field dimensions
                        {
                            // first find maximum number of
                            // intergration points per side
                            label maxNbip = 1;
                            for (stk::mesh::Bucket::size_type ib = 0;
                                 ib < sideBuckets.size();
                                 ++ib)
                            {
                                const stk::mesh::Bucket& bucket =
                                    *sideBuckets[ib];
                                const MasterElement* meFC = MasterElementRepo::
                                    get_surface_master_element(
                                        bucket.topology());
                                const label n_bip = meFC->numIntPoints_;

                                maxNbip = std::max(maxNbip, n_bip);
                            }

                            // max across partitions
                            messager::maxReduce(maxNbip);

                            // get real field dim
                            fieldDim /= maxNbip;
                        }

                        std::vector<scalar> average(fieldDim, 0.0);

                        label totalIntPoints = 0.0;

                        for (stk::mesh::Bucket::size_type ib = 0;
                             ib < sideBuckets.size();
                             ++ib)
                        {
                            const stk::mesh::Bucket& bucket = *sideBuckets[ib];
                            const stk::mesh::Bucket::size_type n_entities =
                                bucket.size();
                            const MasterElement* meFC =
                                MasterElementRepo::get_surface_master_element(
                                    bucket.topology());
                            const label n_bip = meFC->numIntPoints_;
                            for (stk::mesh::Bucket::size_type i = 0;
                                 i < n_entities;
                                 ++i)
                            {
                                const stk::mesh::Entity side = bucket[i];
                                const scalar* val = stk::mesh::field_data(
                                    *sideSTKFieldPtr, side);

                                for (label bip = 0; bip < n_bip; bip++)
                                {
                                    for (label i = 0; i < fieldDim; i++)
                                    {
                                        average[i] += val[bip * fieldDim + i];
                                    }

                                    // increment
                                    totalIntPoints++;
                                }
                            }
                        }

                        if (messager::parallel())
                        {
                            messager::sumReduce(totalIntPoints);
                            messager::sumReduce(average);
                        }

                        // calculate average
                        for (label i = 0; i < fieldDim; i++)
                        {
                            average[i] /= static_cast<scalar>(totalIntPoints);
                        }

                        // print
                        if (messager::master())
                        {
                            std::cout << "Object name: " << name_
                                      << ", value: (";
                            std::cout << average[0];
                            for (label i = 1; i < fieldDim; i++)
                            {
                                std::cout << "\t" << average[i];
                            }
                            std::cout << ")" << std::endl;

                            if (writeToFile_)
                            {
                                std::string fileName(
                                    postProcessPtr_->directory() / name_);
                                std::ofstream file(fileName,
                                                   std::ios_base::app);

                                // write header
                                if (postProcessPtr_->iter() == 1)
                                {
                                    file << postProcessPtr_->instanceHeader() +
                                                "\t";

                                    if (fieldDim == 1)
                                    {
                                        file << "value";
                                    }
                                    else
                                    {
                                        for (label i = 0; i < fieldDim; i++)
                                        {
                                            file << "value_"
                                                 << std::to_string(i) << "\t";
                                        }
                                    }

                                    file << "\n";
                                }

                                file << postProcessPtr_->instance() << "\t";
                                for (label i = 0; i < fieldDim; i++)
                                {
                                    file << average[i] << "\t";
                                }
                                file << "\n";
                                file.close();
                            }
                        }
                    }
                    else
                    {
                        errorMsg("A side field with the name " +
                                 sideSTKFieldPtr->name() + " does not exist");
                    }
                }
            }
            break;
    }
}

} // namespace accel
