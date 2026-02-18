// File : postProcessForce.cpp
// Created : Tue Aug 05 2025 19:49:24 (+0100)
// Author : Mhamad Mahdi Alloush
// Description:
// Copyright (c) 2025 CCFNUM, Lucerne University of Applied Sciences and Arts.
// SPDX-License-Identifier: BSD-3-Clause

#include "postProcess.h"

namespace accel
{

forceObject::forceObject(postProcess* postProcessPtr,
                         std::string name,
                         postProcessType type,
                         std::vector<std::string> location,
                         label frequency,
                         bool writeToFile,
                         bool calculateMoment,
                         std::array<scalar, SPATIAL_DIM> momentCenter,
                         bool totalPrint)
    : postProcessObject(postProcessPtr,
                        name,
                        type,
                        location,
                        frequency,
                        writeToFile),
      calculateMoment_(calculateMoment), momentCenter_(momentCenter),
      totalPrint_(totalPrint)
{
    if (writeToFile_)
    {
        if (fs::is_regular_file(postProcessPtr->directory() / name_))
        {
            fs::remove_all(postProcessPtr->directory() / name_);
        }
        std::string fileName(postProcessPtr->directory() / name_);
        std::ofstream file(fileName);

        file << postProcessPtr_->instanceHeader() + "\t";
        if (totalPrint_)
        {
            file << "force_x\t";
            file << "force_y\t";
#if SPATIAL_DIM == 3
            file << "force_z\t";
#endif
            if (calculateMoment_)
            {
#if SPATIAL_DIM == 3
                file << "moment_x\t";
                file << "moment_y\t";
                file << "moment_z\t";
#elif SPATIAL_DIM == 2
                file << "moment_z\t";
#endif
            }
        }
        else
        {
            // pressure force
            file << "p_force_x\t";
            file << "p_force_y\t";
#if SPATIAL_DIM == 3
            file << "p_force_z\t";
#endif

            // viscous force
            file << "v_force_x\t";
            file << "v_force_y\t";
#if SPATIAL_DIM == 3
            file << "v_force_z\t";
#endif

            if (calculateMoment_)
            {
                // pressure moment
#if SPATIAL_DIM == 3
                file << "p_moment_x\t";
                file << "p_moment_y\t";
                file << "p_moment_z\t";
#elif SPATIAL_DIM == 2
                file << "p_moment_z\t";
#endif

                // viscous moment
#if SPATIAL_DIM == 3
                file << "v_moment_x\t";
                file << "v_moment_y\t";
                file << "v_moment_z\t";
#elif SPATIAL_DIM == 2
                file << "v_moment_z\t";
#endif
            }
        }
        file << "\n";
        file.close();
    }
}

void forceObject::update()
{
    const stk::mesh::MetaData& metaData =
        postProcessPtr_->meshRef().metaDataRef();
    const stk::mesh::BulkData& bulkData =
        postProcessPtr_->meshRef().bulkDataRef();

    // nodal fields to gather
    std::vector<scalar> ws_p;

    // master element
    std::vector<scalar> ws_face_shape_function;

    // get fields
    const auto& pSTKFieldRef = postProcessPtr_->pRef().stkFieldRef();
    const auto& wallShearStressSTKFieldRef =
        postProcessPtr_->wallShearStressRef().stkFieldRef();

    // Get geometric fields
    const auto& exposedAreaVecSTKFieldRef = *metaData.get_field<scalar>(
        metaData.side_rank(), mesh::exposed_area_vector_ID);
    const auto& coordsSTKFieldRef = *metaData.get_field<scalar>(
        stk::topology::NODE_RANK,
        postProcessPtr_->meshRef().getCoordinateFieldName());

    // local force
    scalar l_force_moment[12] = {0.0};

    // individual forces
    scalar ws_p_force[3] = {};
    scalar ws_v_force[3] = {};
    scalar ws_t_force[3] = {};
    scalar ws_p_moment[3] = {};
    scalar ws_v_moment[3] = {};
    scalar ws_t_moment[3] = {};
    scalar ws_radius[3] = {};

    // Collect parts
    stk::mesh::PartVector partVec;
    for (const auto& locationName : location_)
    {
        stk::mesh::Part* part = metaData.get_part(locationName);
        partVec.push_back(part);
    }

    // define some common selectors
    stk::mesh::Selector selOwnedSides =
        metaData.locally_owned_part() & stk::mesh::selectUnion(partVec);

    // shifted ip's for field?
    const bool isPShifted = postProcessPtr_->pRef().isShifted();

    stk::mesh::BucketVector const& sideBuckets =
        bulkData.get_buckets(metaData.side_rank(), selOwnedSides);
    for (stk::mesh::BucketVector::const_iterator ib = sideBuckets.begin();
         ib != sideBuckets.end();
         ++ib)
    {
        stk::mesh::Bucket& sideBucket = **ib;

        // face master element
        MasterElement* meFC =
            accel::MasterElementRepo::get_surface_master_element(
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

                // zero out vector quantities; squeeze in
                // aMag
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

                // extract nodal fields
                stk::mesh::Entity node = sideNodeRels[localFaceNode];
                const scalar* coord =
                    stk::mesh::field_data(coordsSTKFieldRef, node);

                // load radius; assemble force
                // -sigma_ij*njdS
                for (label i = 0; i < SPATIAL_DIM; ++i)
                {
                    const scalar ai = areaVec[offSetAveraVec + i];
                    ws_radius[i] = coord[i] - momentCenter_[i];
                    ws_p_force[i] = pBip * ai;
                    ws_v_force[i] =
                        wallShearStressBip[SPATIAL_DIM * ip + i] * aMag;
                    ws_t_force[i] = ws_p_force[i] + ws_v_force[i];
                }

                ws_p_moment[0] =
                    ws_radius[1] * ws_p_force[2] - ws_radius[2] * ws_p_force[1];
                ws_p_moment[1] = -(ws_radius[0] * ws_p_force[2] -
                                   ws_radius[2] * ws_p_force[0]);
                ws_p_moment[2] =
                    ws_radius[0] * ws_p_force[1] - ws_radius[1] * ws_p_force[0];

                ws_v_moment[0] =
                    ws_radius[1] * ws_v_force[2] - ws_radius[2] * ws_v_force[1];
                ws_v_moment[1] = -(ws_radius[0] * ws_v_force[2] -
                                   ws_radius[2] * ws_v_force[0]);
                ws_v_moment[2] =
                    ws_radius[0] * ws_v_force[1] - ws_radius[1] * ws_v_force[0];

                ws_t_moment[0] = ws_p_moment[0] + ws_v_moment[0];
                ws_t_moment[1] = ws_p_moment[1] + ws_v_moment[1];
                ws_t_moment[2] = ws_p_moment[2] + ws_v_moment[2];

                // assemble force and moment
                for (label j = 0; j < 3; ++j)
                {
                    l_force_moment[j + 0] += ws_p_force[j];
                    l_force_moment[j + 3] += ws_v_force[j];
                    l_force_moment[j + 6] += ws_p_moment[j];
                    l_force_moment[j + 9] += ws_v_moment[j];
                }
            }
        }
    }

    // parallel assemble and output
    scalar g_force_moment[12] = {};

    // Parallel assembly of L2
    stk::all_reduce_sum(
        bulkData.parallel(), &l_force_moment[0], &g_force_moment[0], 12);

    // print
    if (messager::master())
    {
        std::cout << "Object name: " << name_ << std::endl;
        if (totalPrint_)
        {
            // Copy to arg force and moment
            scalar force[3] = {0.0};
            scalar moment[3] = {0.0};
            for (label j = 0; j < 3; ++j)
            {
                force[j] = g_force_moment[j + 0] + g_force_moment[j + 3];
                moment[j] = g_force_moment[j + 6] + g_force_moment[j + 9];
            }

            std::cout << "\tforce_x: " << force[0] << std::endl;
            std::cout << "\tforce_y: " << force[1] << std::endl;
#if SPATIAL_DIM == 3
            std::cout << "\tforce_z: " << force[2] << std::endl;
#endif
            if (calculateMoment_)
            {
#if SPATIAL_DIM == 3
                std::cout << "\tmoment_x: " << moment[0] << std::endl;
                std::cout << "\tmoment_y: " << moment[1] << std::endl;
                std::cout << "\tmoment_z: " << moment[2] << std::endl;
#elif SPATIAL_DIM == 2
                std::cout << "\tmoment_z: " << moment[2] << std::endl;
#endif
            }

            if (writeToFile_)
            {
                std::string fileName(postProcessPtr_->directory() / name_);
                std::ofstream file(fileName, std::ios_base::app);
                file << postProcessPtr_->instance() << "\t";
                file << force[0] << "\t";
                file << force[1] << "\t";
#if SPATIAL_DIM == 3
                file << force[2] << "\t";
#endif
                if (calculateMoment_)
                {
#if SPATIAL_DIM == 3
                    file << moment[0] << "\t";
                    file << moment[1] << "\t";
                    file << moment[2] << "\t";
#elif SPATIAL_DIM == 2
                    file << moment[2] << "\t";
#endif
                }
                file << "\n";
                file.close();
            }
        }
        else
        {
            // Copy to arg force and moment
            scalar p_force[3] = {0.0};
            scalar v_force[3] = {0.0};
            scalar p_moment[3] = {0.0};
            scalar v_moment[3] = {0.0};
            for (label j = 0; j < 3; ++j)
            {
                p_force[j] = g_force_moment[j + 0];
                v_force[j] = g_force_moment[j + 3];
                p_moment[j] = g_force_moment[j + 6];
                v_moment[j] = g_force_moment[j + 9];
            }

            // pressure force
            std::cout << "\tp_force_x: " << p_force[0] << std::endl;
            std::cout << "\tp_force_y: " << p_force[1] << std::endl;
#if SPATIAL_DIM == 3
            std::cout << "\tp_force_z: " << p_force[2] << std::endl;
#endif

            // viscous force
            std::cout << "\tv_force_x: " << v_force[0] << std::endl;
            std::cout << "\tv_force_y: " << v_force[1] << std::endl;
#if SPATIAL_DIM == 3
            std::cout << "\tv_force_z: " << v_force[2] << std::endl;
#endif

            if (calculateMoment_)
            {
                // pressure moment
#if SPATIAL_DIM == 3
                std::cout << "\tp_moment_x: " << p_moment[0] << std::endl;
                std::cout << "\tp_moment_y: " << p_moment[1] << std::endl;
                std::cout << "\tp_moment_z: " << p_moment[2] << std::endl;
#elif SPATIAL_DIM == 2
                std::cout << "\tv_moment_z: " << p_moment[2] << std::endl;
#endif

                // viscous moment
#if SPATIAL_DIM == 3
                std::cout << "\tv_moment_x: " << v_moment[0] << std::endl;
                std::cout << "\tv_moment_y: " << v_moment[1] << std::endl;
                std::cout << "\tv_moment_z: " << v_moment[2] << std::endl;
#elif SPATIAL_DIM == 2
                std::cout << "\tv_moment_z: " << v_moment[2] << std::endl;
#endif
            }

            if (writeToFile_)
            {
                std::string fileName(postProcessPtr_->directory() / name_);
                std::ofstream file(fileName, std::ios_base::app);
                file << postProcessPtr_->instance() << "\t";

                // pressure force
                file << p_force[0] << "\t";
                file << p_force[1] << "\t";
#if SPATIAL_DIM == 3
                file << p_force[2] << "\t";
#endif

                // viscous force
                file << v_force[0] << "\t";
                file << v_force[1] << "\t";
#if SPATIAL_DIM == 3
                file << v_force[2] << "\t";
#endif

                if (calculateMoment_)
                {
                    // pressure moment
#if SPATIAL_DIM == 3
                    file << p_moment[0] << "\t";
                    file << p_moment[1] << "\t";
                    file << p_moment[2] << "\t";
#elif SPATIAL_DIM == 2
                    file << p_moment[2] << "\t";
#endif

                    // viscous moment
#if SPATIAL_DIM == 3
                    file << v_moment[0] << "\t";
                    file << v_moment[1] << "\t";
                    file << v_moment[2] << "\t";
#elif SPATIAL_DIM == 2
                    file << v_moment[2] << "\t";
#endif
                }
                file << "\n";
                file.close();
            }
        }
    }
}

} // namespace accel
