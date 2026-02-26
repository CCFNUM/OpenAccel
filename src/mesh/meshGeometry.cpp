// File       : meshGeometry.cpp
// Created    : Fri Aug 25 2023 12:55:24 (+0100)
// Author     : Mhamad Mahdi Alloush
// Description:
// Copyright (c) 2023 CCFNUM, Lucerne University of Applied Sciences and Arts.
// SPDX-License-Identifier: BSD-3-Clause

// code
#include "boundary.h"
#ifdef HAS_INTERFACE
#include "dgInfo.h"
#include "interface.h"
#include "interfaceSideInfo.h"
#endif /* HAS_INTERFACE */
#include "mesh.h"
#include "messager.h"
#include "zone.h"

namespace accel
{

void mesh::setupCoordinateField_()
{
    // put field on moving zones
    for (label iZone = 0; iZone < this->nZones(); iZone++)
    {
        if (this->zoneRef(iZone).meshMoving())
        {
#ifndef NDEBUG
            if (messager::master())
            {
                std::cout
                    << "\tsetting up original cooridinates field for zone: "
                    << this->zoneRef(iZone).name() << std::endl;
            }
#endif
            auto& originalCoordinatesSTKFieldRef =
                metaDataRef().declare_field<scalar>(stk::topology::NODE_RANK,
                                                    original_coordinates_ID);

            stk::mesh::put_field_on_mesh(
                originalCoordinatesSTKFieldRef,
                stk::mesh::selectUnion(this->zoneRef(iZone).interiorParts()),
                SPATIAL_DIM,
                nullptr);
        }
    }
}

void mesh::setupGeometricFields_()
{
#ifndef NDEBUG
    {
        auto& rankIDSTK = metaDataRef().declare_field<int>(
            stk::topology::ELEMENT_RANK, rank_ID);
        stk::mesh::put_field_on_mesh(
            rankIDSTK,
            stk::mesh::selectUnion(interiorActiveParts_),
            1,
            nullptr);

        stk::io::set_field_output_type(rankIDSTK,
                                       stk::io::FieldOutputType::SCALAR);
    }

    // SCL (Space Conservation Law) check field for mesh-moving zones
    for (label iZone = 0; iZone < this->nZones(); iZone++)
    {
        if (this->zoneRef(iZone).meshMoving())
        {
            auto& sclCheckSTK = metaDataRef().declare_field<scalar>(
                stk::topology::NODE_RANK, scl_check_ID);
            stk::mesh::put_field_on_mesh(
                sclCheckSTK,
                stk::mesh::selectUnion(this->zoneRef(iZone).interiorParts()),
                1,
                nullptr);

            stk::io::set_field_output_type(sclCheckSTK,
                                           stk::io::FieldOutputType::SCALAR);
            break; // field declared once, put on all moving zones
        }
    }
    // Put field on remaining moving zones
    for (label iZone = 0; iZone < this->nZones(); iZone++)
    {
        if (this->zoneRef(iZone).meshMoving())
        {
            auto* sclCheckSTK = metaDataRef().get_field<scalar>(
                stk::topology::NODE_RANK, scl_check_ID);
            if (sclCheckSTK)
            {
                stk::mesh::put_field_on_mesh(
                    *sclCheckSTK,
                    stk::mesh::selectUnion(
                        this->zoneRef(iZone).interiorParts()),
                    1,
                    nullptr);
            }
        }
    }
#endif /* NDEBUG */

    // Dual nodal volume
    {
        auto& dualNodalVolumeSTKFieldRef = metaDataRef().declare_field<scalar>(
            stk::topology::NODE_RANK, dual_nodal_volume_ID);
        stk::mesh::put_field_on_mesh(
            dualNodalVolumeSTKFieldRef,
            stk::mesh::selectUnion(interiorActiveParts_),
            1,
            nullptr);

        // put original field on deforming zones
        for (label iZone = 0; iZone < this->nZones(); iZone++)
        {
            if (this->zoneRef(iZone).meshDeforming())
            {
#ifndef NDEBUG
                if (messager::master())
                {
                    std::cout << "\tsetting up original vol field for zone: "
                              << this->zoneRef(iZone).name() << std::endl;
                }
#endif
                auto& originalDualNodalVolumeSTKFieldRef =
                    metaDataRef().declare_field<scalar>(
                        stk::topology::NODE_RANK,
                        original_dual_nodal_volume_ID);

                stk::mesh::put_field_on_mesh(
                    originalDualNodalVolumeSTKFieldRef,
                    stk::mesh::selectUnion(
                        this->zoneRef(iZone).interiorParts()),
                    1,
                    nullptr);
            }
        }
    }

    // Exposed area vector
    {
        auto& exposedAreaVecSTKFieldRef = metaDataRef().declare_field<scalar>(
            metaDataRef().side_rank(), exposed_area_vector_ID);

        // Put on all boundaries
        for (const stk::mesh::Part* part : boundaryActiveParts_)
        {
            // A part could have different topologies, by getting the subparts
            // we make sure we take into account all topologies.
            for (const stk::mesh::Part* subPart : part->subsets())
            {
                // Determine number of integration points in the face
                MasterElement* meFC =
                    MasterElementRepo::get_surface_master_element(
                        subPart->topology());
                const label numScsIp = meFC->numIntPoints_;

                // Put the field on mesh
                stk::mesh::put_field_on_mesh(exposedAreaVecSTKFieldRef,
                                             *subPart,
                                             SPATIAL_DIM * numScsIp,
                                             nullptr);
            }
        }

        // put original field on moving zones
        for (label iZone = 0; iZone < this->nZones(); iZone++)
        {
            if (this->zoneRef(iZone).meshMoving())
            {
#ifndef NDEBUG
                if (messager::master())
                {
                    std::cout
                        << "\tsetting up original exposed area vector field "
                           "for bondaries touching zone: "
                        << this->zoneRef(iZone).name() << std::endl;
                }
#endif
                auto& originalExposedAreaVecSTKFieldRef =
                    metaDataRef().declare_field<scalar>(
                        metaDataRef().side_rank(),
                        original_exposed_area_vector_ID);

#ifdef HAS_INTERFACE
                // Interface
                for (auto interf : this->zonePtr(iZone)->interfacesRef())
                {
                    if (interf->isInternal())
                    {
                        for (const stk::mesh::Part* part :
                             interf->masterInfoPtr()->currentPartVec_)
                        {
                            // A part could have different topologies, by
                            // getting the subparts we make sure we take into
                            // account all topologies.
                            for (const stk::mesh::Part* subPart :
                                 part->subsets())
                            {
                                // Determine number of integration points in the
                                // face
                                MasterElement* meFC = MasterElementRepo::
                                    get_surface_master_element(
                                        subPart->topology());
                                const label numScsIp = meFC->numIntPoints_;

                                // Put the field on mesh
                                stk::mesh::put_field_on_mesh(
                                    originalExposedAreaVecSTKFieldRef,
                                    *subPart,
                                    SPATIAL_DIM * numScsIp,
                                    nullptr);
                            }
                        }

                        for (const stk::mesh::Part* part :
                             interf->slaveInfoPtr()->currentPartVec_)
                        {
                            // A part could have different topologies, by
                            // getting the subparts we make sure we take into
                            // account all topologies.
                            for (const stk::mesh::Part* subPart :
                                 part->subsets())
                            {
                                // Determine number of integration points in the
                                // face
                                MasterElement* meFC = MasterElementRepo::
                                    get_surface_master_element(
                                        subPart->topology());
                                const label numScsIp = meFC->numIntPoints_;

                                // Put the field on mesh
                                stk::mesh::put_field_on_mesh(
                                    originalExposedAreaVecSTKFieldRef,
                                    *subPart,
                                    SPATIAL_DIM * numScsIp,
                                    nullptr);
                            }
                        }
                    }
                    else
                    {
                        for (const stk::mesh::Part* part :
                             interf->interfaceSideInfoPtr(iZone)
                                 ->currentPartVec_)
                        {
                            // A part could have different topologies, by
                            // getting the subparts we make sure we take into
                            // account all topologies.
                            for (const stk::mesh::Part* subPart :
                                 part->subsets())
                            {
                                // Determine number of integration points in the
                                // face
                                MasterElement* meFC = MasterElementRepo::
                                    get_surface_master_element(
                                        subPart->topology());
                                const label numScsIp = meFC->numIntPoints_;

                                // Put the field on mesh
                                stk::mesh::put_field_on_mesh(
                                    originalExposedAreaVecSTKFieldRef,
                                    *subPart,
                                    SPATIAL_DIM * numScsIp,
                                    nullptr);
                            }
                        }
                    }
                }
#endif /* HAS_INTERFACE */

                // Boundary
                for (label iBoundary = 0;
                     iBoundary < this->zoneRef(iZone).nBoundaries();
                     iBoundary++)
                {
                    for (const stk::mesh::Part* part :
                         this->zonePtr(iZone)->boundaryPtr(iBoundary)->parts())
                    {
                        // A part could have different topologies, by getting
                        // the subparts we make sure we take into account all
                        // topologies.
                        for (const stk::mesh::Part* subPart : part->subsets())
                        {
                            // Determine number of integration points in the
                            // face
                            MasterElement* meFC =
                                MasterElementRepo::get_surface_master_element(
                                    subPart->topology());
                            const label numScsIp = meFC->numIntPoints_;

                            // Put the field on mesh
                            stk::mesh::put_field_on_mesh(
                                originalExposedAreaVecSTKFieldRef,
                                *subPart,
                                SPATIAL_DIM * numScsIp,
                                nullptr);
                        }
                    }
                }
            }
        }
    }

    // Wall normal distance
    {
        auto& wallNormalDistanceSTKFieldRef =
            (metaDataRef().declare_field<scalar>(metaDataRef().side_rank(),
                                                 wall_normal_distance_ID));

        // Put on all wall boundaries
        for (const stk::mesh::Part* part : wallBoundaryActiveParts_)
        {
            // A part could have different topologies, by getting the subparts
            // we make sure we take into account all topologies.
            for (const stk::mesh::Part* subPart : part->subsets())
            {
                // Determine number of integration points in
                // the face
                MasterElement* meFC =
                    MasterElementRepo::get_surface_master_element(
                        subPart->topology());
                const label numScsIp = meFC->numIntPoints_;

                // Put the field on mesh
                stk::mesh::put_field_on_mesh(
                    wallNormalDistanceSTKFieldRef, *subPart, numScsIp, nullptr);
            }
        }
    }

    // Assembled wall area
    {
        auto& assembledWallAreaSTKFieldRef =
            (metaDataRef().declare_field<scalar>(stk::topology::NODE_RANK,
                                                 assembled_wall_area_ID));

        // initialize to 0
        scalar initialValue = 0.0;

        // Put on all wall boundaries
        for (const stk::mesh::Part* part : wallBoundaryActiveParts_)
        {
            stk::mesh::put_field_on_mesh(
                assembledWallAreaSTKFieldRef, *part, 1, &initialValue);
        }
    }

    // Assembled symm area
    {
        auto& assembledSymmAreaSTKFieldRef =
            (metaDataRef().declare_field<scalar>(stk::topology::NODE_RANK,
                                                 assembled_symm_area_ID));

        // initialize to 0
        scalar initialValue = 0.0;

        // Put on all symmetry boundaries
        for (const stk::mesh::Part* part : symmetryBoundaryActiveParts_)
        {
            stk::mesh::put_field_on_mesh(assembledSymmAreaSTKFieldRef,
                                         *part,
                                         SPATIAL_DIM,
                                         &initialValue);
        }
    }

    // Assembled wall normal distance
    {
        auto& assembledWallNormalDistanceSTKFieldRef =
            (metaDataRef().declare_field<scalar>(
                stk::topology::NODE_RANK, assembled_wall_normal_distance_ID));

        // initialize to 0
        scalar initialValue = 0.0;

        for (const stk::mesh::Part* part : wallBoundaryActiveParts_)
        {
            stk::mesh::put_field_on_mesh(assembledWallNormalDistanceSTKFieldRef,
                                         *part,
                                         1,
                                         &initialValue);
        }
    }
}

void mesh::setupZones_()
{
    for (label iZone = 0; iZone < this->nZones(); iZone++)
    {
        zonePtr(iZone)->setup();
    }
}

#ifdef HAS_INTERFACE
void mesh::setupInterfaces_()
{
    for (label iInterface = 0; iInterface < nInterfaces(); iInterface++)
    {
        interfaceRef(iInterface).setup();
    }
}
#endif /* HAS_INTERFACE */

void mesh::initializeCoordinateField_()
{
    // set custom name for coordinates field
    this->metaDataRef().set_coordinate_field_name(coordinates_ID);

    // synchronize coordinates field for ghosted nodes if parallel run
    if (messager::parallel())
    {
        STKScalarField* coordinates = this->metaDataRef().get_field<scalar>(
            stk::topology::NODE_RANK, coordinates_ID);
        stk::mesh::communicate_field_data(this->bulkDataRef(), {coordinates});
    }

    // scaling if required
    if (scaleMesh_)
    {
        scale(scaleVector_);
    }

    // copy coordinates to the original coordinates in case of a
    // deforming/moving mesh. Only done at initialization
    for (label iZone = 0; iZone < this->nZones(); iZone++)
    {
        if (this->zoneRef(iZone).meshMoving())
        {
            // get required components
            stk::mesh::BulkData& bulkData = this->bulkDataRef();
            stk::mesh::MetaData& metaData = this->metaDataRef();

            // get fields
            const STKScalarField* coordinates = metaData.get_field<scalar>(
                stk::topology::NODE_RANK, coordinates_ID);
            STKScalarField* orgCoordinates = metaData.get_field<scalar>(
                stk::topology::NODE_RANK, original_coordinates_ID);

            // if org field is not yet created, do it now. This might happen in
            // cases were inherent mesh deformation occurs
            if (!orgCoordinates)
            {
                orgCoordinates = &metaDataRef().declare_field<scalar>(
                    stk::topology::NODE_RANK, original_coordinates_ID);

                stk::mesh::put_field_on_mesh(
                    *orgCoordinates,
                    stk::mesh::selectUnion(
                        this->zoneRef(iZone).interiorParts()),
                    SPATIAL_DIM,
                    nullptr);
            }
            else
            {
                // check if defined over the current zone parts
                if (!orgCoordinates->defined_on_any(
                        this->zoneRef(iZone).interiorParts()))
                {
                    stk::mesh::put_field_on_mesh(
                        *orgCoordinates,
                        stk::mesh::selectUnion(
                            this->zoneRef(iZone).interiorParts()),
                        SPATIAL_DIM,
                        nullptr);
                }
            }

#ifndef NDEBUG
            if (messager::master())
            {
                std::cout << "copying coordinates field to original in zone: "
                          << this->zoneRef(iZone).name() << std::endl;
            }
#endif
            // copy
            ops::copy<scalar>(coordinates,
                              orgCoordinates,
                              this->zoneRef(iZone).interiorParts());
        }
    }
}

void mesh::initializeZones_()
{
    for (label iZone = 0; iZone < this->nZones(); iZone++)
    {
        zoneRef(iZone).initialize();
    }
}

#ifdef HAS_INTERFACE
void mesh::initializeInterfaces_()
{
    if (hasInterfaces())
    {
        // update interface info for master and slave of every pair
        for (label iInterface = 0; iInterface < nInterfaces(); iInterface++)
        {
            interfaceRef(iInterface).initialize();
        }
    }
}
#endif /* HAS_INTERFACE */

void mesh::initializeGeometricFields_()
{
    // Now compute geometric quantities (volume, exposed area, etc) on the
    // entire mesh

    // Dual nodal volume
    updateDualNodalVolumeField_();

    // Original dual nodal volume
    for (label iZone = 0; iZone < this->nZones(); iZone++)
    {
        if (this->zoneRef(iZone).meshDeforming())
        {
            // get fields
            const auto* volSTKFieldPtr = this->metaDataPtr()->get_field<scalar>(
                stk::topology::NODE_RANK, mesh::dual_nodal_volume_ID);
            auto* orgVolSTKFieldPtr = this->metaDataPtr()->get_field<scalar>(
                stk::topology::NODE_RANK, mesh::original_dual_nodal_volume_ID);

            // if org field is not yet created, do it now. This might happen in
            // cases were inherent mesh deformation occurs
            if (!orgVolSTKFieldPtr)
            {
                orgVolSTKFieldPtr = &metaDataRef().declare_field<scalar>(
                    stk::topology::NODE_RANK, original_dual_nodal_volume_ID);

                stk::mesh::put_field_on_mesh(
                    *orgVolSTKFieldPtr,
                    stk::mesh::selectUnion(
                        this->zoneRef(iZone).interiorParts()),
                    1,
                    nullptr);
            }
            else
            {
                // check if defined over the current zone parts
                if (!orgVolSTKFieldPtr->defined_on_any(
                        this->zoneRef(iZone).interiorParts()))
                {
                    stk::mesh::put_field_on_mesh(
                        *orgVolSTKFieldPtr,
                        stk::mesh::selectUnion(
                            this->zoneRef(iZone).interiorParts()),
                        1,
                        nullptr);
                }
            }

#ifndef NDEBUG
            if (messager::master())
            {
                std::cout << "copying vol field to original in zone: "
                          << this->zoneRef(iZone).name() << std::endl;
            }
#endif
            // copy
            ops::copy<scalar>(volSTKFieldPtr,
                              orgVolSTKFieldPtr,
                              this->zonePtr(iZone)->interiorParts());
        }
    }

    // exposed area vector
    updateExposedAreaVectorField_();

    // Original exposed area vector
    for (label iZone = 0; iZone < this->nZones(); iZone++)
    {
        if (this->zoneRef(iZone).meshMoving())
        {
            const auto* exposedAreaVectorKFieldPtr =
                this->metaDataPtr()->get_field<scalar>(
                    this->metaDataPtr()->side_rank(),
                    mesh::exposed_area_vector_ID);
            auto* orgExposedAreaVectorKFieldPtr =
                this->metaDataPtr()->get_field<scalar>(
                    this->metaDataPtr()->side_rank(),
                    mesh::original_exposed_area_vector_ID);

            // if org field is not yet created, do it now. This might happen in
            // cases were inherent mesh deformation occurs
            if (!orgExposedAreaVectorKFieldPtr)
            {
                orgExposedAreaVectorKFieldPtr =
                    &metaDataRef().declare_field<scalar>(
                        metaDataRef().side_rank(),
                        original_exposed_area_vector_ID);

#ifdef HAS_INTERFACE
                // Interface
                for (auto interf : this->zonePtr(iZone)->interfacesRef())
                {
                    if (interf->isInternal())
                    {
                        for (const stk::mesh::Part* part :
                             interf->masterInfoPtr()->currentPartVec_)
                        {
                            // A part could have different topologies, by
                            // getting the subparts we make sure we take into
                            // account all topologies.
                            for (const stk::mesh::Part* subPart :
                                 part->subsets())
                            {
                                // Determine number of integration points in the
                                // face
                                MasterElement* meFC = MasterElementRepo::
                                    get_surface_master_element(
                                        subPart->topology());
                                const label numScsIp = meFC->numIntPoints_;

                                // Put the field on mesh
                                stk::mesh::put_field_on_mesh(
                                    *orgExposedAreaVectorKFieldPtr,
                                    *subPart,
                                    SPATIAL_DIM * numScsIp,
                                    nullptr);
                            }
                        }

                        for (const stk::mesh::Part* part :
                             interf->slaveInfoPtr()->currentPartVec_)
                        {
                            // A part could have different topologies, by
                            // getting the subparts we make sure we take into
                            // account all topologies.
                            for (const stk::mesh::Part* subPart :
                                 part->subsets())
                            {
                                // Determine number of integration points in the
                                // face
                                MasterElement* meFC = MasterElementRepo::
                                    get_surface_master_element(
                                        subPart->topology());
                                const label numScsIp = meFC->numIntPoints_;

                                // Put the field on mesh
                                stk::mesh::put_field_on_mesh(
                                    *orgExposedAreaVectorKFieldPtr,
                                    *subPart,
                                    SPATIAL_DIM * numScsIp,
                                    nullptr);
                            }
                        }
                    }
                    else
                    {
                        for (const stk::mesh::Part* part :
                             interf->interfaceSideInfoPtr(iZone)
                                 ->currentPartVec_)
                        {
                            // A part could have different topologies, by
                            // getting the subparts we make sure we take into
                            // account all topologies.
                            for (const stk::mesh::Part* subPart :
                                 part->subsets())
                            {
                                // Determine number of integration points in the
                                // face
                                MasterElement* meFC = MasterElementRepo::
                                    get_surface_master_element(
                                        subPart->topology());
                                const label numScsIp = meFC->numIntPoints_;

                                // Put the field on mesh
                                stk::mesh::put_field_on_mesh(
                                    *orgExposedAreaVectorKFieldPtr,
                                    *subPart,
                                    SPATIAL_DIM * numScsIp,
                                    nullptr);
                            }
                        }
                    }
                }
#endif /* HAS_INTERFACE */

                // Boundary
                for (label iBoundary = 0;
                     iBoundary < this->zoneRef(iZone).nBoundaries();
                     iBoundary++)
                {
                    for (const stk::mesh::Part* part :
                         this->zonePtr(iZone)->boundaryPtr(iBoundary)->parts())
                    {
                        // A part could have different topologies, by getting
                        // the subparts we make sure we take into account all
                        // topologies.
                        for (const stk::mesh::Part* subPart : part->subsets())
                        {
                            // Determine number of integration points in the
                            // face
                            MasterElement* meFC =
                                MasterElementRepo::get_surface_master_element(
                                    subPart->topology());
                            const label numScsIp = meFC->numIntPoints_;

                            // Put the field on mesh
                            stk::mesh::put_field_on_mesh(
                                *orgExposedAreaVectorKFieldPtr,
                                *subPart,
                                SPATIAL_DIM * numScsIp,
                                nullptr);
                        }
                    }
                }
            }
            else
            {
#ifdef HAS_INTERFACE
                // Interface
                for (auto interf : this->zonePtr(iZone)->interfacesRef())
                {
                    if (interf->isInternal())
                    {
                        // check if defined over the current interface parts
                        if (!orgExposedAreaVectorKFieldPtr->defined_on_any(
                                interf->masterInfoPtr()->currentPartVec_))
                        {
                            for (const stk::mesh::Part* part :
                                 interf->masterInfoPtr()->currentPartVec_)
                            {
                                // A part could have different topologies, by
                                // getting the subparts we make sure we take
                                // into account all topologies.
                                for (const stk::mesh::Part* subPart :
                                     part->subsets())
                                {
                                    // Determine number of integration points in
                                    // the face
                                    MasterElement* meFC = MasterElementRepo::
                                        get_surface_master_element(
                                            subPart->topology());
                                    const label numScsIp = meFC->numIntPoints_;

                                    // Put the field on mesh
                                    stk::mesh::put_field_on_mesh(
                                        *orgExposedAreaVectorKFieldPtr,
                                        *subPart,
                                        SPATIAL_DIM * numScsIp,
                                        nullptr);
                                }
                            }
                        }

                        // check if defined over the opposing interface parts
                        if (!orgExposedAreaVectorKFieldPtr->defined_on_any(
                                interf->slaveInfoPtr()->currentPartVec_))
                        {
                            for (const stk::mesh::Part* part :
                                 interf->slaveInfoPtr()->currentPartVec_)
                            {
                                // A part could have different topologies, by
                                // getting the subparts we make sure we take
                                // into account all topologies.
                                for (const stk::mesh::Part* subPart :
                                     part->subsets())
                                {
                                    // Determine number of integration points in
                                    // the face
                                    MasterElement* meFC = MasterElementRepo::
                                        get_surface_master_element(
                                            subPart->topology());
                                    const label numScsIp = meFC->numIntPoints_;

                                    // Put the field on mesh
                                    stk::mesh::put_field_on_mesh(
                                        *orgExposedAreaVectorKFieldPtr,
                                        *subPart,
                                        SPATIAL_DIM * numScsIp,
                                        nullptr);
                                }
                            }
                        }
                    }
                    else
                    {
                        // check if defined over the current interface parts
                        if (!orgExposedAreaVectorKFieldPtr->defined_on_any(
                                interf->interfaceSideInfoPtr(iZone)
                                    ->currentPartVec_))
                        {
                            for (const stk::mesh::Part* part :
                                 interf->interfaceSideInfoPtr(iZone)
                                     ->currentPartVec_)
                            {
                                // A part could have different topologies, by
                                // getting the subparts we make sure we take
                                // into account all topologies.
                                for (const stk::mesh::Part* subPart :
                                     part->subsets())
                                {
                                    // Determine number of integration points in
                                    // the face
                                    MasterElement* meFC = MasterElementRepo::
                                        get_surface_master_element(
                                            subPart->topology());
                                    const label numScsIp = meFC->numIntPoints_;

                                    // Put the field on mesh
                                    stk::mesh::put_field_on_mesh(
                                        *orgExposedAreaVectorKFieldPtr,
                                        *subPart,
                                        SPATIAL_DIM * numScsIp,
                                        nullptr);
                                }
                            }
                        }
                    }
                }
#endif /* HAS_INTERFACE */

                // Boundary
                for (label iBoundary = 0;
                     iBoundary < this->zoneRef(iZone).nBoundaries();
                     iBoundary++)
                {
                    // check if defined over the current boundary parts
                    if (!orgExposedAreaVectorKFieldPtr->defined_on_any(
                            this->zonePtr(iZone)
                                ->boundaryPtr(iBoundary)
                                ->parts()))
                    {
                        for (const stk::mesh::Part* part :
                             this->zonePtr(iZone)
                                 ->boundaryPtr(iBoundary)
                                 ->parts())
                        {
                            // A part could have different topologies, by
                            // getting the subparts we make sure we take into
                            // account all topologies.
                            for (const stk::mesh::Part* subPart :
                                 part->subsets())
                            {
                                // Determine number of integration points in the
                                // face
                                MasterElement* meFC = MasterElementRepo::
                                    get_surface_master_element(
                                        subPart->topology());
                                const label numScsIp = meFC->numIntPoints_;

                                // Put the field on mesh
                                stk::mesh::put_field_on_mesh(
                                    *orgExposedAreaVectorKFieldPtr,
                                    *subPart,
                                    SPATIAL_DIM * numScsIp,
                                    nullptr);
                            }
                        }
                    }
                }
            }

#ifndef NDEBUG
            if (messager::master())
            {
                std::cout << "copying exposed area vector field to original on "
                             "boundaries touching zone: "
                          << this->zoneRef(iZone).name() << std::endl;
            }
#endif

#ifdef HAS_INTERFACE
            // Interface
            for (auto interf : this->zonePtr(iZone)->interfacesRef())
            {
                if (interf->isInternal())
                {
                    ops::copy<scalar>(exposedAreaVectorKFieldPtr,
                                      orgExposedAreaVectorKFieldPtr,
                                      interf->masterInfoPtr()->currentPartVec_);
                    ops::copy<scalar>(exposedAreaVectorKFieldPtr,
                                      orgExposedAreaVectorKFieldPtr,
                                      interf->slaveInfoPtr()->currentPartVec_);
                }
                else
                {
                    ops::copy<scalar>(
                        exposedAreaVectorKFieldPtr,
                        orgExposedAreaVectorKFieldPtr,
                        interf->interfaceSideInfoPtr(iZone)->currentPartVec_);
                }
            }
#endif /* HAS_INTERFACE */

            // Boundary
            for (label iBoundary = 0;
                 iBoundary < this->zoneRef(iZone).nBoundaries();
                 iBoundary++)
            {
                ops::copy<scalar>(
                    exposedAreaVectorKFieldPtr,
                    orgExposedAreaVectorKFieldPtr,
                    this->zonePtr(iZone)->boundaryPtr(iBoundary)->parts());
            }
        }
    }

    // Wall normal distance
    updateWallNormalDistanceField_();

    // Assembled wall area
    updateAssembledWallAreaField_();

    // Assembled symmetry area
    updateAssembledSymmetryAreaField_();

    // Assembled wall normal distance
    updateAssembledWallNormalDistanceField_();

#ifndef NDEBUG
    {
        // get required components
        stk::mesh::BulkData& bulkData = this->bulkDataRef();
        stk::mesh::MetaData& metaData = this->metaDataRef();

        const int myrank = messager::myProcNo();
        stk::mesh::Field<int>& rankIDSTK =
            *metaData.get_field<int>(stk::topology::ELEMENT_RANK, rank_ID);
        for (stk::mesh::Bucket* bucket : bulkData.get_buckets(
                 stk::topology::ELEMENT_RANK,
                 metaData.universal_part() &
                     stk::mesh::selectUnion(interiorActiveParts_)))
        {
            for (const stk::mesh::Entity e : *bucket)
            {
                int* r = stk::mesh::field_data(rankIDSTK, e);
                *r = myrank;
            }
        }
    }
#endif /* NDEBUG */
}

void mesh::updateZones_(bool force)
{
    for (label iZone = 0; iZone < nZones(); iZone++)
    {
        if (this->zonePtr(iZone)->meshMoving() || force)
        {
#ifndef NDEBUG
            if (messager::master())
            {
                std::cout << "updating zone: " << zoneRef(iZone).name()
                          << std::endl;
            }
#endif
            zoneRef(iZone).update();
        }
    }
}

#ifdef HAS_INTERFACE
void mesh::updateInterfaces_(bool force)
{
    // only update if any of the two connected zones is moving (deforming or
    // transforming)
    for (label iInterface = 0; iInterface < nInterfaces(); iInterface++)
    {
        if (this->interfaceRef(iInterface)
                .masterInfoPtr()
                ->zonePtr()
                ->meshMoving() ||
            this->interfaceRef(iInterface)
                .slaveInfoPtr()
                ->zonePtr()
                ->meshMoving() ||
            force)
        {
#ifndef NDEBUG
            if (messager::master())
            {
                std::cout << "updating interface: "
                          << this->interfaceRef(iInterface).name() << std::endl;
            }
#endif
            interfaceRef(iInterface).update();
        }
    }
}
#endif /* HAS_INTERFACE */

void mesh::updateGeometricFields_(bool force)
{
    // get required components
    stk::mesh::BulkData& bulkData = this->bulkDataRef();
    stk::mesh::MetaData& metaData = this->metaDataRef();

    // update dual nodal volume
    {
        // collect all interior parts of deforming zones: volume will not change
        // otherwise
        stk::mesh::ConstPartVector parts;
        for (label iZone = 0; iZone < nZones(); iZone++)
        {
            if (this->zonePtr(iZone)->meshDeforming() || force)
            {
#ifndef NDEBUG
                if (messager::master())
                {
                    std::cout << "updating vol field for zone: "
                              << zoneRef(iZone).name() << std::endl;
                }
#endif
                for (auto part : this->zonePtr(iZone)->interiorParts())
                {
                    parts.push_back(part);
                }
            }
        }

        // update
        updateDualNodalVolumeField_(parts);
    }

    // update exposed area vector
    {
        // collect all boundary/interface parts of moving zones (deforming or
        // transforming): exposed area vector will not change otherwise
        stk::mesh::ConstPartVector parts;
        for (label iZone = 0; iZone < nZones(); iZone++)
        {
            if (this->zonePtr(iZone)->meshMoving() || force)
            {
#ifndef NDEBUG
                if (messager::master())
                {
                    std::cout << "updating exposed area vector field for "
                                 "boundaries touching zone: "
                              << zoneRef(iZone).name() << std::endl;
                }
#endif

#ifdef HAS_INTERFACE
                for (auto interf : this->zonePtr(iZone)->interfacesRef())
                {
                    if (interf->isInternal())
                    {
                        for (auto part :
                             interf->masterInfoPtr()->currentPartVec_)
                        {
                            parts.push_back(part);
                        }
                        for (auto part :
                             interf->slaveInfoPtr()->currentPartVec_)
                        {
                            parts.push_back(part);
                        }
                    }
                    else
                    {
                        for (auto part : interf->interfaceSideInfoPtr(iZone)
                                             ->currentPartVec_)
                        {
                            parts.push_back(part);
                        }
                    }
                }
#endif /* HAS_INTERFACE */

                for (label iBoundary = 0;
                     iBoundary < this->zonePtr(iZone)->nBoundaries();
                     iBoundary++)
                {
                    for (auto part :
                         this->zonePtr(iZone)->boundaryPtr(iBoundary)->parts())
                    {
                        parts.push_back(part);
                    }
                }
            }
        }

        // update
        updateExposedAreaVectorField_(parts);
    }

    // update wall normal distance
    {
        // collect all wall/fluid-solid interface parts of deforming zones: wall
        // normal distance will not change otherwise
        stk::mesh::ConstPartVector parts;
        for (label iZone = 0; iZone < nZones(); iZone++)
        {
            if (this->zonePtr(iZone)->meshDeforming() || force)
            {
                for (label iBoundary = 0;
                     iBoundary < this->zonePtr(iZone)->nBoundaries();
                     iBoundary++)
                {
                    auto type = zoneRef(iZone).boundaryRef(iBoundary).type();
                    switch (type)
                    {
                        case boundaryPhysicalType::wall:
                            {
                                for (auto part : this->zonePtr(iZone)
                                                     ->boundaryPtr(iBoundary)
                                                     ->parts())
                                {
                                    parts.push_back(part);
                                }
                            }
                            break;

                        default:
                            break;
                    }
                }

#ifdef HAS_INTERFACE
                for (auto interf : this->zonePtr(iZone)->interfacesRef())
                {
                    if (interf->isFluidSolidType())
                    {
                        for (auto part : interf->interfaceSideInfoPtr(iZone)
                                             ->currentPartVec_)
                        {
                            parts.push_back(part);
                        }
                    }
                }
#endif /* HAS_INTERFACE */
            }
        }

        // update
        updateWallNormalDistanceField_(parts);
    }

    // update assembled wall area
    {
        // collect all wall/fluid-solid interface parts of deforming zones:
        // assembled wall area will not change otherwise
        stk::mesh::ConstPartVector parts;
        for (label iZone = 0; iZone < nZones(); iZone++)
        {
            if (this->zonePtr(iZone)->meshDeforming() || force)
            {
                for (label iBoundary = 0;
                     iBoundary < this->zonePtr(iZone)->nBoundaries();
                     iBoundary++)
                {
                    auto type = zoneRef(iZone).boundaryRef(iBoundary).type();
                    switch (type)
                    {
                        case boundaryPhysicalType::wall:
                            {
                                for (auto part : this->zonePtr(iZone)
                                                     ->boundaryPtr(iBoundary)
                                                     ->parts())
                                {
                                    parts.push_back(part);
                                }
                            }
                            break;

                        default:
                            break;
                    }
                }

#ifdef HAS_INTERFACE
                for (auto interf : this->zonePtr(iZone)->interfacesRef())
                {
                    if (interf->isFluidSolidType())
                    {
                        for (auto part : interf->interfaceSideInfoPtr(iZone)
                                             ->currentPartVec_)
                        {
                            parts.push_back(part);
                        }
                    }
                }
#endif /* HAS_INTERFACE */
            }
        }

        // update
        updateAssembledWallAreaField_(parts);
    }

    // update assembled symmetry area
    {
        // collect all symmetry parts of deforming zones: assembled symmetry
        // area will not change otherwise
        stk::mesh::ConstPartVector parts;
        for (label iZone = 0; iZone < nZones(); iZone++)
        {
            if (this->zonePtr(iZone)->meshDeforming() || force)
            {
                for (label iBoundary = 0;
                     iBoundary < this->zonePtr(iZone)->nBoundaries();
                     iBoundary++)
                {
                    auto type = zoneRef(iZone).boundaryRef(iBoundary).type();
                    switch (type)
                    {
                        case boundaryPhysicalType::symmetry:
                            {
                                for (auto part : this->zonePtr(iZone)
                                                     ->boundaryPtr(iBoundary)
                                                     ->parts())
                                {
                                    parts.push_back(part);
                                }
                            }
                            break;

                        default:
                            break;
                    }
                }
            }
        }

        // update
        updateAssembledSymmetryAreaField_(parts);
    }

    // update assembled wall normal distance
    {
        // collect all wall/fluid-solid interface parts of deforming zones:
        // assembled wall normal distance will not change otherwise
        stk::mesh::ConstPartVector parts;
        for (label iZone = 0; iZone < nZones(); iZone++)
        {
            if (this->zonePtr(iZone)->meshDeforming() || force)
            {
                for (label iBoundary = 0;
                     iBoundary < this->zonePtr(iZone)->nBoundaries();
                     iBoundary++)
                {
                    auto type = zoneRef(iZone).boundaryRef(iBoundary).type();
                    switch (type)
                    {
                        case boundaryPhysicalType::wall:
                            {
                                for (auto part : this->zonePtr(iZone)
                                                     ->boundaryPtr(iBoundary)
                                                     ->parts())
                                {
                                    parts.push_back(part);
                                }
                            }
                            break;

                        default:
                            break;
                    }
                }

#ifdef HAS_INTERFACE
                for (auto interf : this->zonePtr(iZone)->interfacesRef())
                {
                    if (interf->isFluidSolidType())
                    {
                        for (auto part : interf->interfaceSideInfoPtr(iZone)
                                             ->currentPartVec_)
                        {
                            parts.push_back(part);
                        }
                    }
                }
#endif /* HAS_INTERFACE */
            }
        }

        // update
        updateAssembledWallNormalDistanceField_(parts);
    }
}

void mesh::updateDualNodalVolumeField_(stk::mesh::ConstPartVector interiorParts)
{
    // if no interior parts assigned, consider all active ones
    auto parts = interiorParts.empty() ? interiorActiveParts_ : interiorParts;

    // get required components
    stk::mesh::BulkData& bulkData = this->bulkDataRef();
    stk::mesh::MetaData& metaData = this->metaDataRef();

    STKScalarField& dualNodalVolumeSTKFieldRef =
        *metaDataRef().get_field<scalar>(stk::topology::NODE_RANK,
                                         dual_nodal_volume_ID);
    // zero volume field first
    ops::zero<scalar>(&dualNodalVolumeSTKFieldRef);

    // calculate volume of control volume
    {
        // extract coordinates field
        STKScalarField* coordinates = metaData.get_field<scalar>(
            stk::topology::NODE_RANK, this->getCoordinateFieldName());

        // setup for buckets; union parts and ask for all element
        // buckets
        stk::mesh::BucketVector const& elementBuckets = bulkData.get_buckets(
            stk::topology::ELEMENT_RANK,
            metaData.universal_part() & stk::mesh::selectUnion(parts) &
                stk::mesh::selectField(dualNodalVolumeSTKFieldRef));

        for (stk::mesh::BucketVector::const_iterator ib =
                 elementBuckets.begin();
             ib != elementBuckets.end();
             ++ib)
        {
            stk::mesh::Bucket& elementBucket = **ib;

            // extract master element
            MasterElement* meSCV = MasterElementRepo::get_volume_master_element(
                elementBucket.topology());

            // extract master element specifics
            const label nodesPerElement = meSCV->nodesPerElement_;
            const label numScvIp = meSCV->numIntPoints_;
            const label* ipNodeMap = meSCV->ipNodeMap();

            // define scratch field
            std::vector<scalar> ws_coordinates(nodesPerElement * SPATIAL_DIM);
            std::vector<scalar> ws_scv_volume(numScvIp);

            const stk::mesh::Bucket::size_type nElementsPerBucket =
                elementBucket.size();
            for (stk::mesh::Bucket::size_type iElement = 0;
                 iElement < nElementsPerBucket;
                 ++iElement)
            {
                //===============================================
                // gather nodal data; this is how we do it now..
                //===============================================
                stk::mesh::Entity const* nodeRels =
                    elementBucket.begin_nodes(iElement);
                label num_nodes = elementBucket.num_nodes(iElement);

                // sanity check on num nodes
                STK_ThrowAssert(num_nodes == nodesPerElement);

                for (label ni = 0; ni < num_nodes; ++ni)
                {
                    stk::mesh::Entity node = nodeRels[ni];
                    scalar* coords = stk::mesh::field_data(*coordinates, node);
                    const label offSet = ni * SPATIAL_DIM;
                    for (label j = 0; j < SPATIAL_DIM; ++j)
                    {
                        ws_coordinates[offSet + j] = coords[j];
                    }
                }

                // compute integration point volume
                scalar scv_error = 0.0;
                meSCV->determinant(
                    1, &ws_coordinates[0], &ws_scv_volume[0], &scv_error);

                // assemble dual volume
                for (label ip = 0; ip < numScvIp; ++ip)
                {
                    // nearest node for this ip
                    const label nn = ipNodeMap[ip];
                    stk::mesh::Entity node = nodeRels[nn];

                    scalar* dualcv =
                        stk::mesh::field_data(dualNodalVolumeSTKFieldRef, node);

                    // augment nodal dual volume
                    *dualcv += ws_scv_volume[ip];
                }
            }
        }
    }

#ifndef NDEBUG
    // Check for negative volumes in debug mode
    {
        stk::mesh::Selector selAllNodes =
            metaData.universal_part() & stk::mesh::selectUnion(parts) &
            stk::mesh::selectField(dualNodalVolumeSTKFieldRef);

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
                stk::mesh::Entity node = nodeBucket[iNode];

                scalar* dualcv =
                    stk::mesh::field_data(dualNodalVolumeSTKFieldRef, node);

                if (*dualcv < 0.0)
                {
                    stk::mesh::EntityId nodeId = bulkData.identifier(node);

                    // Collect surrounding element IDs
                    std::vector<stk::mesh::EntityId> elementIds;
                    stk::mesh::Entity const* elemRels =
                        bulkData.begin_elements(node);
                    const unsigned numElements = bulkData.num_elements(node);

                    for (unsigned i = 0; i < numElements; ++i)
                    {
                        elementIds.push_back(bulkData.identifier(elemRels[i]));
                    }

                    // Print error message with rank information for
                    // parallel runs
                    std::ostringstream msg;
                    msg << "[Rank " << messager::myProcNo()
                        << "] Negative dual volume detected!\n"
                        << "  Node ID: " << nodeId << "\n"
                        << "  Volume: " << *dualcv << "\n"
                        << "  Surrounding elements: ";

                    for (size_t i = 0; i < elementIds.size(); ++i)
                    {
                        if (i > 0)
                            msg << ", ";
                        msg << elementIds[i];
                    }
                    msg << std::endl;

                    std::cout << msg.str();
                }
            }
        }
    }
#endif /* NDEBUG */
}

void mesh::updateExposedAreaVectorField_(
    stk::mesh::ConstPartVector boundaryParts)
{
    // if no boundary parts assigned, consider all active ones
    auto parts = boundaryParts.empty() ? boundaryActiveParts_ : boundaryParts;

    // get required components
    stk::mesh::BulkData& bulkData = this->bulkDataRef();
    stk::mesh::MetaData& metaData = this->metaDataRef();

    // Get field
    STKScalarField& exposedAreaVecSTKFieldRef = *metaData.get_field<scalar>(
        metaData.side_rank(), exposed_area_vector_ID);

    // extract coordinates field
    STKScalarField* coordinates = metaData.get_field<scalar>(
        stk::topology::NODE_RANK, this->getCoordinateFieldName());

    // setup for buckets; union parts and ask for universal part
    stk::mesh::Selector selAllSides =
        metaData.universal_part() & stk::mesh::selectUnion(parts) &
        stk::mesh::selectField(exposedAreaVecSTKFieldRef);
    stk::mesh::BucketVector const& sideBuckets =
        bulkData.get_buckets(metaData.side_rank(), selAllSides);

    for (stk::mesh::BucketVector::const_iterator ib = sideBuckets.begin();
         ib != sideBuckets.end();
         ++ib)
    {
        stk::mesh::Bucket& sideBucket = **ib;

        // extract master element
        MasterElement* meFC = MasterElementRepo::get_surface_master_element(
            sideBucket.topology());

        // extract master element specifics
        const label nodesPerSide = meFC->nodesPerElement_;
        const label numScsBip = meFC->numIntPoints_;

        // define scratch field
        std::vector<scalar> ws_coordinates(nodesPerSide * SPATIAL_DIM);
        std::vector<scalar> ws_scs_areav(numScsBip * SPATIAL_DIM);

        const stk::mesh::Bucket::size_type nSidesPerBucket = sideBucket.size();
        for (stk::mesh::Bucket::size_type iSide = 0; iSide < nSidesPerBucket;
             ++iSide)
        {
            // face data
            scalar* areaVec = stk::mesh::field_data(
                exposedAreaVecSTKFieldRef, sideBucket, iSide);

            // face node relations for nodal gather
            stk::mesh::Entity const* sideNodeRels =
                sideBucket.begin_nodes(iSide);

            //===============================================
            // gather nodal data; this is how we do it now..
            //===============================================
            label num_nodes = sideBucket.num_nodes(iSide);
            for (label ni = 0; ni < num_nodes; ++ni)
            {
                stk::mesh::Entity node = sideNodeRels[ni];
                scalar* coords = stk::mesh::field_data(*coordinates, node);
                const label offSet = ni * SPATIAL_DIM;
                for (label j = 0; j < SPATIAL_DIM; ++j)
                {
                    ws_coordinates[offSet + j] = coords[j];
                }
            }

            // compute scs integration point areavec
            scalar scs_error = 0.0;
            meFC->determinant(
                1, &ws_coordinates[0], &ws_scs_areav[0], &scs_error);

            // scarrer to area vector
            for (label ip = 0; ip < numScsBip; ++ip)
            {
                const label offSet = ip * SPATIAL_DIM;
                for (label j = 0; j < SPATIAL_DIM; ++j)
                {
                    areaVec[offSet + j] = ws_scs_areav[offSet + j];
                }
            }
        }
    }
}

void mesh::updateWallNormalDistanceField_(stk::mesh::ConstPartVector wallParts)
{
    // if no wall parts assigned, consider all active ones
    auto parts = wallParts.empty() ? wallBoundaryActiveParts_ : wallParts;

    // get required components
    stk::mesh::BulkData& bulkData = this->bulkDataRef();
    stk::mesh::MetaData& metaData = this->metaDataRef();

    // Get fields
    STKScalarField& wallNormalDistanceSTKFieldRef = *metaData.get_field<scalar>(
        metaData.side_rank(), wall_normal_distance_ID);
    STKScalarField& exposedAreaVecSTKFieldRef = *metaData.get_field<scalar>(
        metaData.side_rank(), exposed_area_vector_ID);
    STKScalarField& coordsSTKFieldRef = *metaData.get_field<scalar>(
        stk::topology::NODE_RANK, this->getCoordinateFieldName());

    // define vector of parent topos; should always be UNITY in size
    std::vector<stk::topology> parentTopo;

    // setup for buckets; union parts and ask for universal part
    stk::mesh::Selector selAllSides =
        metaData.universal_part() & stk::mesh::selectUnion(parts) &
        stk::mesh::selectField(wallNormalDistanceSTKFieldRef);
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

        // mapping from ip to nodes for this ordinal; face
        // perspective (use with face_node_relations)
        const label* faceIpNodeMap = meFC->ipNodeMap();

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

            // pointer to face data
            const scalar* areaVec =
                stk::mesh::field_data(exposedAreaVecSTKFieldRef, side);

            scalar* wallNormalDistanceBip =
                stk::mesh::field_data(wallNormalDistanceSTKFieldRef, side);

            // extract the connected element to this exposed face;
            // should be single in size!
            const stk::mesh::Entity* face_elem_rels =
                bulkData.begin_elements(side);
            STK_ThrowAssert(bulkData.num_elements(side) == 1);

            // get element; its face ordinal number
            stk::mesh::Entity element = face_elem_rels[0];
            const label face_ordinal = bulkData.begin_element_ordinals(side)[0];

            // get the relations off of element
            stk::mesh::Entity const* elemNodeRels =
                bulkData.begin_nodes(element);

            // loop over face nodes
            for (label ip = 0; ip < numScsBip; ++ip)
            {
                const label offSetAveraVec = ip * SPATIAL_DIM;

                const label opposingNode =
                    meSCS->opposingNodes(face_ordinal, ip);
                const label localFaceNode = faceIpNodeMap[ip];

                // left and right nodes; right is on the face; left
                // is the opposing node
                stk::mesh::Entity nodeL = elemNodeRels[opposingNode];
                stk::mesh::Entity nodeR = sideNodeRels[localFaceNode];

                // extract nodal fields
                const scalar* coordL =
                    stk::mesh::field_data(coordsSTKFieldRef, nodeL);
                const scalar* coordR =
                    stk::mesh::field_data(coordsSTKFieldRef, nodeR);

                // zero out vector quantities; squeeze in aMag
                scalar aMag = 0.0;
                for (label j = 0; j < SPATIAL_DIM; ++j)
                {
                    const scalar axj = areaVec[offSetAveraVec + j];
                    aMag += axj * axj;
                }
                aMag = std::sqrt(aMag);

                // determine yp: we adopt full edge
                // important: An alternative approach takes the average of the
                // components and adopts 1/4 edge
                // NOTE: this algorithm gives inconsistent results,
                // affecting the turbulence solutions.
                // scalar ypBip = 0.0;
                // for (label j = 0; j < SPATIAL_DIM; ++j)
                // {
                // const scalar nj = areaVec[offSetAveraVec + j]
                // / aMag; const scalar ej = coordR[j] -
                // coordL[j]; ypBip += nj * ej * nj * ej;
                // }
                // wallNormalDistanceBip[ip] = std::sqrt(ypBip);

                // this is a simple projection of the distance
                // vector from the boundary node to the first
                // internal node, projected onto the reverse face
                // normal. The lower bound should be small enough to
                // allow for micron-scale distances.
                scalar ypBip = 0.0;
                for (label j = 0; j < SPATIAL_DIM; ++j)
                {
                    const scalar nj = -areaVec[offSetAveraVec + j] / aMag;
                    const scalar ej = coordL[j] - coordR[j];
                    ypBip += nj * ej;
                }

                wallNormalDistanceBip[ip] = std::max(ypBip, SMALL);
            }
        }
    }
}

void mesh::updateAssembledWallAreaField_(stk::mesh::ConstPartVector wallParts)
{
    // if no wall parts assigned, consider all active ones
    auto parts = wallParts.empty() ? wallBoundaryActiveParts_ : wallParts;

    // get required components
    stk::mesh::BulkData& bulkData = this->bulkDataRef();
    stk::mesh::MetaData& metaData = this->metaDataRef();

    // Get fields
    STKScalarField& assembledWallAreaSTKFieldRef = *metaData.get_field<scalar>(
        stk::topology::NODE_RANK, assembled_wall_area_ID);

    // zero assembled fields first
    ops::zero<scalar>(&assembledWallAreaSTKFieldRef);

    // calculate assembled quantities
    {
        const STKScalarField& exposedAreaVecSTKFieldRef =
            *metaData.get_field<scalar>(metaData.side_rank(),
                                        exposed_area_vector_ID);

        // define vector of parent topos; should always be UNITY in
        // size
        std::vector<stk::topology> parentTopo;

        // define some common selectors
        stk::mesh::Selector selAllSides =
            metaData.universal_part() & stk::mesh::selectUnion(parts) &
            stk::mesh::selectField(assembledWallAreaSTKFieldRef);

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

            // mapping from ip to nodes for this ordinal; face
            // perspective (use with face_node_relations)
            const label* faceIpNodeMap = meFC->ipNodeMap();

            const stk::mesh::Bucket::size_type length = sideBucket.size();

            for (stk::mesh::Bucket::size_type k = 0; k < length; ++k)
            {
                // get face
                stk::mesh::Entity side = sideBucket[k];

                //======================================
                // gather nodal data off of face
                //======================================
                stk::mesh::Entity const* sideNodeRels =
                    bulkData.begin_nodes(side);
                label numSideNodes = bulkData.num_nodes(side);

                // pointer to face data
                const scalar* areaVec =
                    stk::mesh::field_data(exposedAreaVecSTKFieldRef, side);

                // loop over face nodes
                for (label ip = 0; ip < numSideNodes; ++ip)
                {
                    const label offSetAveraVec = ip * SPATIAL_DIM;

                    const label localFaceNode = faceIpNodeMap[ip];

                    // left and right nodes; right is on the face;
                    // left is the opposing node
                    stk::mesh::Entity node = sideNodeRels[localFaceNode];

                    // zero out vector quantities; squeeze in aMag
                    scalar aMag = 0.0;
                    for (label j = 0; j < SPATIAL_DIM; ++j)
                    {
                        const scalar axj = areaVec[offSetAveraVec + j];
                        aMag += axj * axj;
                    }
                    aMag = std::sqrt(aMag);

                    // assemble to nodal quantities
                    scalar* assembledArea = stk::mesh::field_data(
                        assembledWallAreaSTKFieldRef, node);

                    *assembledArea += aMag;
                }
            }
        }

        if (messager::parallel())
        {
            stk::mesh::communicate_field_data(bulkData,
                                              {&assembledWallAreaSTKFieldRef});
        }
    }
}

void mesh::updateAssembledSymmetryAreaField_(
    stk::mesh::ConstPartVector symmParts)
{
    // if no symmetry parts assigned, consider all active ones
    auto parts = symmParts.empty() ? symmetryBoundaryActiveParts_ : symmParts;

    // get required components
    stk::mesh::BulkData& bulkData = this->bulkDataRef();
    stk::mesh::MetaData& metaData = this->metaDataRef();

    // Get fields
    STKScalarField& assembledSymmAreaSTKFieldRef = *metaData.get_field<scalar>(
        stk::topology::NODE_RANK, assembled_symm_area_ID);

    // zero assembled fields first
    ops::zero<scalar>(&assembledSymmAreaSTKFieldRef);

    // calculate assembled quantities
    {
        const STKScalarField& exposedAreaVecSTKFieldRef =
            *metaData.get_field<scalar>(metaData.side_rank(),
                                        exposed_area_vector_ID);

        // define vector of parent topos; should always be UNITY in
        // size
        std::vector<stk::topology> parentTopo;

        // define some common selectors
        stk::mesh::Selector selAllSides =
            metaData.universal_part() & stk::mesh::selectUnion(parts) &
            stk::mesh::selectField(assembledSymmAreaSTKFieldRef);
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

            // mapping from ip to nodes for this ordinal; face
            // perspective (use with face_node_relations)
            const label* faceIpNodeMap = meFC->ipNodeMap();

            const stk::mesh::Bucket::size_type length = sideBucket.size();

            for (stk::mesh::Bucket::size_type k = 0; k < length; ++k)
            {
                // get face
                stk::mesh::Entity side = sideBucket[k];

                //======================================
                // gather nodal data off of face
                //======================================
                stk::mesh::Entity const* sideNodeRels =
                    bulkData.begin_nodes(side);
                label numSideNodes = bulkData.num_nodes(side);

                // pointer to face data
                const scalar* areaVec =
                    stk::mesh::field_data(exposedAreaVecSTKFieldRef, side);

                // loop over face nodes
                for (label ip = 0; ip < numSideNodes; ++ip)
                {
                    const label offSetAveraVec = ip * SPATIAL_DIM;

                    const label localFaceNode = faceIpNodeMap[ip];

                    // left and right nodes; right is on the face;
                    // left is the opposing node
                    stk::mesh::Entity node = sideNodeRels[localFaceNode];

                    scalar* aarea = stk::mesh::field_data(
                        assembledSymmAreaSTKFieldRef, node);

                    for (label i = 0; i < SPATIAL_DIM; i++)
                    {
                        aarea[i] += areaVec[offSetAveraVec + i];
                    }
                }
            }
        }

        if (messager::parallel())
        {
            stk::mesh::communicate_field_data(bulkData,
                                              {&assembledSymmAreaSTKFieldRef});
        }
    }
}

void mesh::updateAssembledWallNormalDistanceField_(
    stk::mesh::ConstPartVector wallParts)
{
    // if no wall parts assigned, consider all active ones
    auto parts = wallParts.empty() ? wallBoundaryActiveParts_ : wallParts;

    // get required components
    stk::mesh::BulkData& bulkData = this->bulkDataRef();
    stk::mesh::MetaData& metaData = this->metaDataRef();

    // Get fields
    STKScalarField& assembledWallNormalDistanceSTKFieldRef =
        *metaData.get_field<scalar>(stk::topology::NODE_RANK,
                                    assembled_wall_normal_distance_ID);

    STKScalarField* assembledWallAreaSTKFieldPtr = metaData.get_field<scalar>(
        stk::topology::NODE_RANK, assembled_wall_area_ID);

    // zero assembled fields first
    ops::zero<scalar>(&assembledWallNormalDistanceSTKFieldRef);

    // calculate assembled quantities
    {
        const STKScalarField& wallNormalDistanceSTKFieldRef =
            *metaData.get_field<scalar>(metaData.side_rank(),
                                        wall_normal_distance_ID);
        const STKScalarField& exposedAreaVecSTKFieldRef =
            *metaData.get_field<scalar>(metaData.side_rank(),
                                        exposed_area_vector_ID);

        // define some common selectors
        stk::mesh::Selector selAllSides =
            metaData.universal_part() & stk::mesh::selectUnion(parts) &
            stk::mesh::selectField(wallNormalDistanceSTKFieldRef);
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

            // mapping from ip to nodes for this ordinal; face
            // perspective (use with face_node_relations)
            const label* faceIpNodeMap = meFC->ipNodeMap();

            const stk::mesh::Bucket::size_type length = sideBucket.size();

            for (stk::mesh::Bucket::size_type k = 0; k < length; ++k)
            {
                // get face
                stk::mesh::Entity side = sideBucket[k];

                //======================================
                // gather nodal data off of face
                //======================================
                stk::mesh::Entity const* sideNodeRels =
                    bulkData.begin_nodes(side);
                label numSideNodes = bulkData.num_nodes(side);

                // pointer to face data
                const scalar* areaVec =
                    stk::mesh::field_data(exposedAreaVecSTKFieldRef, side);
                const scalar* wallNormalDistanceBip =
                    stk::mesh::field_data(wallNormalDistanceSTKFieldRef, side);

                // loop over face nodes
                for (label ip = 0; ip < numSideNodes; ++ip)
                {
                    const label offSetAveraVec = ip * SPATIAL_DIM;

                    const label localFaceNode = faceIpNodeMap[ip];

                    // left and right nodes; right is on the face;
                    // left is the opposing node
                    stk::mesh::Entity node = sideNodeRels[localFaceNode];

                    // zero out vector quantities; squeeze in aMag
                    scalar aMag = 0.0;
                    for (label j = 0; j < SPATIAL_DIM; ++j)
                    {
                        const scalar axj = areaVec[offSetAveraVec + j];
                        aMag += axj * axj;
                    }
                    aMag = std::sqrt(aMag);

                    // assemble to nodal quantities
                    scalar* assembledWallNormalDistance = stk::mesh::field_data(
                        assembledWallNormalDistanceSTKFieldRef, node);

                    scalar* assembledArea = stk::mesh::field_data(
                        *assembledWallAreaSTKFieldPtr, node);

                    *assembledWallNormalDistance +=
                        wallNormalDistanceBip[ip] * aMag / (*assembledArea);
                }
            }
        }

        if (messager::parallel())
        {
            stk::mesh::communicate_field_data(
                bulkData, {&assembledWallNormalDistanceSTKFieldRef});
        }
    }
}

void mesh::scale(std::array<scalar, SPATIAL_DIM> scaleVector)
{
    // get required components
    stk::mesh::BulkData& bulkData = this->bulkDataRef();
    stk::mesh::MetaData& metaData = this->metaDataRef();

    STKScalarField* coordinates =
        metaData.get_field<scalar>(stk::topology::NODE_RANK, coordinates_ID);

    stk::mesh::Selector selAllNodes = metaData.universal_part();

    stk::mesh::BucketVector const& nodeBuckets =
        bulkData.get_buckets(stk::topology::NODE_RANK, selAllNodes);
    for (stk::mesh::BucketVector::const_iterator ib = nodeBuckets.begin();
         ib != nodeBuckets.end();
         ++ib)
    {
        stk::mesh::Bucket& nodeBucket = **ib;
        const stk::mesh::Bucket::size_type nNodesPerBucket = nodeBucket.size();

        for (stk::mesh::Bucket::size_type iNode = 0; iNode < nNodesPerBucket;
             ++iNode)
        {
            stk::mesh::Entity node = nodeBucket[iNode];

            scalar* coords = stk::mesh::field_data(*coordinates, node);

            // translate to scale origin first
            for (label j = 0; j < SPATIAL_DIM; ++j)
            {
                coords[j] -= scaleOrigin_[j];
            }

            // scale
            for (label j = 0; j < SPATIAL_DIM; ++j)
            {
                coords[j] *= scaleVector[j];
            }

            // translate to original coordinates
            for (label j = 0; j < SPATIAL_DIM; ++j)
            {
                coords[j] += scaleOrigin_[j];
            }
        }
    }
}

void mesh::scale(label iZone, std::array<scalar, SPATIAL_DIM> scaleVector)
{
    // get required components
    stk::mesh::BulkData& bulkData = this->bulkDataRef();
    stk::mesh::MetaData& metaData = this->metaDataRef();

    STKScalarField* coordinates =
        metaData.get_field<scalar>(stk::topology::NODE_RANK, coordinates_ID);

    stk::mesh::Selector selAllNodes =
        metaData.universal_part() &
        stk::mesh::selectUnion(this->zoneRef(iZone).interiorParts());

    stk::mesh::BucketVector const& nodeBuckets =
        bulkData.get_buckets(stk::topology::NODE_RANK, selAllNodes);
    for (stk::mesh::BucketVector::const_iterator ib = nodeBuckets.begin();
         ib != nodeBuckets.end();
         ++ib)
    {
        stk::mesh::Bucket& nodeBucket = **ib;
        const stk::mesh::Bucket::size_type nNodesPerBucket = nodeBucket.size();

        for (stk::mesh::Bucket::size_type iNode = 0; iNode < nNodesPerBucket;
             ++iNode)
        {
            stk::mesh::Entity node = nodeBucket[iNode];

            scalar* coords = stk::mesh::field_data(*coordinates, node);

            // translate to scale origin first
            for (label j = 0; j < SPATIAL_DIM; ++j)
            {
                coords[j] -= scaleOrigin_[j];
            }

            // scale
            for (label j = 0; j < SPATIAL_DIM; ++j)
            {
                coords[j] *= scaleVector[j];
            }

            // translate to original coordinates
            for (label j = 0; j < SPATIAL_DIM; ++j)
            {
                coords[j] += scaleOrigin_[j];
            }
        }
    }
}

void mesh::scale(label iZone,
                 label iBoundary,
                 std::array<scalar, SPATIAL_DIM> scaleVector)
{
    // get required components
    stk::mesh::BulkData& bulkData = this->bulkDataRef();
    stk::mesh::MetaData& metaData = this->metaDataRef();

    STKScalarField* coordinates =
        metaData.get_field<scalar>(stk::topology::NODE_RANK, coordinates_ID);

    stk::mesh::Selector selAllNodes =
        metaData.universal_part() &
        stk::mesh::selectUnion(
            this->zoneRef(iZone).boundaryRef(iBoundary).parts());

    stk::mesh::BucketVector const& nodeBuckets =
        bulkData.get_buckets(stk::topology::NODE_RANK, selAllNodes);
    for (stk::mesh::BucketVector::const_iterator ib = nodeBuckets.begin();
         ib != nodeBuckets.end();
         ++ib)
    {
        stk::mesh::Bucket& nodeBucket = **ib;
        const stk::mesh::Bucket::size_type nNodesPerBucket = nodeBucket.size();

        for (stk::mesh::Bucket::size_type iNode = 0; iNode < nNodesPerBucket;
             ++iNode)
        {
            stk::mesh::Entity node = nodeBucket[iNode];

            scalar* coords = stk::mesh::field_data(*coordinates, node);

            // translate to scale origin first
            for (label j = 0; j < SPATIAL_DIM; ++j)
            {
                coords[j] -= scaleOrigin_[j];
            }

            // scale
            for (label j = 0; j < SPATIAL_DIM; ++j)
            {
                coords[j] *= scaleVector[j];
            }

            // translate to original coordinates
            for (label j = 0; j < SPATIAL_DIM; ++j)
            {
                coords[j] += scaleOrigin_[j];
            }
        }
    }
}

} // namespace accel
