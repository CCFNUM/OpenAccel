// File : mesh.cpp
// Created : Fri Aug 25 2023 12:55:24 (+0100)
// Author : Mhamad Mahdi Alloush
// Description:
// Copyright (c) 2023 CCFNUM, Lucerne University of Applied Sciences and Arts.
// SPDX-License-Identifier: BSD-3-Clause

// code
#include "mesh.h"
#include "controls.h"
#ifdef HAS_INTERFACE
#include "interface.h"
#include "interfaceSideInfo.h"
#endif /* HAS_INTERFACE */
#include "messager.h"
#include "realm.h"
#include "zone.h"

#if SPATIAL_DIM == 3
#include "elementValidator.h"
#endif // SPATIAL_DIM == 3

namespace accel
{

#if SPATIAL_DIM == 3
// Custom deleter implementation
void ElementValidatorDeleter::operator()(elementValidator* ptr)
{
    delete ptr;
}
#endif // SPATIAL_DIM == 3

mesh::mesh(controls* controlsPtr) : controlsPtr_(controlsPtr)
{
    stencil_ = controlsPtr->isReducedStencil()
                   ? ::linearSolver::GraphLayout::Stencil__Reduced
                   : ::linearSolver::GraphLayout::Stencil__Full;
}

mesh::~mesh()
{
    // Explicit destruction order for proper cleanup
#if SPATIAL_DIM == 3
    elementValidatorPtr_.reset();
#endif // SPATIAL_DIM == 3
    delete ioBrokerPtr_;
}

// Operations

void mesh::reset()
{
    // set current coordinates to original ones
    for (label iZone = 0; iZone < this->nZones(); iZone++)
    {
        if (this->zoneRef(iZone).meshMoving())
        {
            stk::mesh::BulkData& bulkData = this->bulkDataRef();
            stk::mesh::MetaData& metaData = this->metaDataRef();

            // Get fields
            const auto& originalCoordsSTKFieldRef = *metaData.get_field<scalar>(
                stk::topology::NODE_RANK, mesh::original_coordinates_ID);
            auto& coordsSTKFieldRef = *metaData.get_field<scalar>(
                stk::topology::NODE_RANK, mesh::coordinates_ID);

            // get interior parts the domain is defined on
            const stk::mesh::PartVector& partVec =
                this->zoneRef(iZone).interiorParts();

            // define some common selectors; select all nodes
            stk::mesh::Selector selUniversalNodes =
                metaData.universal_part() & stk::mesh::selectUnion(partVec);

            stk::mesh::BucketVector const& nodeBuckets = bulkData.get_buckets(
                stk::topology::NODE_RANK, selUniversalNodes);
            for (stk::mesh::BucketVector::const_iterator ib =
                     nodeBuckets.begin();
                 ib != nodeBuckets.end();
                 ++ib)
            {
                stk::mesh::Bucket& nodeBucket = **ib;

                const stk::mesh::Bucket::size_type nNodesPerBucket =
                    nodeBucket.size();

                // field chunks in bucket
                const scalar* orgCoordsb = stk::mesh::field_data(
                    originalCoordsSTKFieldRef, nodeBucket);
                scalar* coordsb =
                    stk::mesh::field_data(coordsSTKFieldRef, nodeBucket);

                for (stk::mesh::Bucket::size_type iNode = 0;
                     iNode < nNodesPerBucket;
                     ++iNode)
                {
                    for (label i = 0; i < SPATIAL_DIM; i++)
                    {
                        // update coords
                        coordsb[SPATIAL_DIM * iNode + i] =
                            orgCoordsb[SPATIAL_DIM * iNode + i];
                    }
                }
            }
        }
    }

    // update geometry
    update();
}

void mesh::setup()
{
    // An original coordinate field for cases of deforming/moving mesh has to be
    // stored on all mesh interior parts
    setupCoordinateField_();

    // Create fundamental STK geometric fields: dual volumes, exposed area
    // vector, wall distance, etc., and put on all the mesh and/or boundaries
    setupGeometricFields_();

    // Create and allocate containers
    setupZones_();

#ifdef HAS_INTERFACE
    // Create and allocate structures and containers
    setupInterfaces_();
#endif /* HAS_INTERFACE */

#if SPATIAL_DIM == 3
    // Create element validator for quality checking and correction (if enabled)
    setupElementValidation_();
#endif // SPATIAL_DIM == 3
}

void mesh::initialize()
{
    // Populate_mesh fills in the entities (nodes/elements/etc) and
    // connectivities, but no field-data. Field-data is not allocated yet.
    ioBrokerPtr_->populate_mesh();

    // Now the mesh is fully populated, so we're ready to populate
    // field-data including coordinates, and attributes and/or distribution
    // factors if those exist on the input mesh file.
    ioBrokerPtr_->populate_field_data();

    // After population, it is essential to make sure that the coordinate
    // field is filled and ready to use
    initializeCoordinateField_();

    // Initialize zones: boundaries, etc.
    initializeZones_();

#ifdef HAS_INTERFACE
    // Initialize interfaces: make required searches
    initializeInterfaces_();
#endif /* HAS_INTERFACE */

    // Calculate geometric quantities (volume, exposed area, etc.)
    initializeGeometricFields_();

    // Build crs node graph
    lazyInitializeNodeGraph_();

#if SPATIAL_DIM == 3
    // Initialize element validator (if enabled)
    initializeElementValidation_();
#endif // SPATIAL_DIM == 3

    // load imbalance
    reportLoadImbalance_();
}

void mesh::update()
{
    // Update mesh structures in case of possible mesh motion/deformation
    updateZones_();
#ifdef HAS_INTERFACE
    updateInterfaces_();
#endif /* HAS_INTERFACE */
    updateGeometricFields_();
    updateNodeGraph_();
}

// Access

controls& mesh::controlsRef()
{
    return *controlsPtr_;
}

const controls& mesh::controlsRef() const
{
    return *controlsPtr_;
}

stk::io::StkMeshIoBroker* mesh::ioBrokerPtr()
{
    return ioBrokerPtr_;
}

const stk::io::StkMeshIoBroker* mesh::ioBrokerPtr() const
{
    return ioBrokerPtr_;
}

stk::io::StkMeshIoBroker& mesh::ioBrokerRef()
{
    return *ioBrokerPtr_;
}

const stk::io::StkMeshIoBroker& mesh::ioBrokerRef() const
{
    return *ioBrokerPtr_;
}

stk::mesh::ConstPartVector mesh::interiorActiveParts()
{
    return interiorActiveParts_;
}

const stk::mesh::ConstPartVector mesh::interiorActiveParts() const
{
    return interiorActiveParts_;
}

stk::mesh::ConstPartVector mesh::boundaryActiveParts()
{
    return boundaryActiveParts_;
}

const stk::mesh::ConstPartVector mesh::boundaryActiveParts() const
{
    return boundaryActiveParts_;
}

stk::mesh::ConstPartVector mesh::wallBoundaryActiveParts()
{
    return wallBoundaryActiveParts_;
}

const stk::mesh::ConstPartVector mesh::wallBoundaryActiveParts() const
{
    return wallBoundaryActiveParts_;
}

stk::mesh::ConstPartVector mesh::symmetryBoundaryActiveParts()
{
    return symmetryBoundaryActiveParts_;
}

const stk::mesh::ConstPartVector mesh::symmetryBoundaryActiveParts() const
{
    return symmetryBoundaryActiveParts_;
}

#ifdef HAS_INTERFACE
// Interfaces

label mesh::nInterfaces() const
{
    return interfaceVector_.size();
}

std::vector<interface*> mesh::interfaceVector() const
{
    std::vector<interface*> interfaceVec(nInterfaces(), nullptr);
    for (label iInterface = 0; iInterface < nInterfaces(); iInterface++)
    {
        interfaceVec[iInterface] = interfaceVector_[iInterface].get();
    }
    return interfaceVec;
}

interface& mesh::interfaceRef(label iInterface)
{
    return *interfaceVector_[iInterface].get();
}

const interface& mesh::interfaceRef(label iInterface) const
{
    return *interfaceVector_[iInterface].get();
}
#endif /* HAS_INTERFACE */

// Zones

label mesh::nZones() const
{
    return zoneVector_.size();
}

std::vector<zone*> mesh::zoneVector() const
{
    std::vector<zone*> zoneVec(nZones(), nullptr);
    for (label iZone = 0; iZone < nZones(); iZone++)
    {
        zoneVec[iZone] = zoneVector_[iZone].get();
    }
    return zoneVec;
}

zone& mesh::zoneRef(label iZone)
{
    return *zoneVector_[iZone].get();
}

const zone& mesh::zoneRef(label iZone) const
{
    return *zoneVector_[iZone].get();
}

zone* mesh::zonePtr(label iZone)
{
    return zoneVector_[iZone].get();
}

const zone* mesh::zonePtr(label iZone) const
{
    return zoneVector_[iZone].get();
}

// private helper

void mesh::reportLoadImbalance_() const
{
    const stk::mesh::MetaData& meta_data = this->metaDataRef();
    const stk::mesh::BulkData& bulk_data = this->bulkDataRef();

    stk::mesh::EntityVector owned_nodes;
    stk::mesh::Selector local_selector = meta_data.locally_owned_part();
    bulk_data.get_entities(
        stk::topology::NODE_RANK, local_selector, owned_nodes);

    const label size = messager::nProcs();
    const label myrank = messager::myProcNo();
    const label n_locally_owned = owned_nodes.size();

    struct LocalNodes
    {
        label nodes;
        label rank;
    };

    LocalNodes max_locally_owned = {n_locally_owned, myrank};
    LocalNodes min_locally_owned = {n_locally_owned, myrank};
    std::vector<label> all_locally_owned(size, n_locally_owned);

    // clang-format off
    if (messager::parallel())
    {
        LocalNodes gmax, gmin;
        MPI_Reduce(&max_locally_owned, &gmax, 1, MPI_2INT, MPI_MAXLOC, 0, messager::comm());
        MPI_Reduce(&min_locally_owned, &gmin, 1, MPI_2INT, MPI_MINLOC, 0, messager::comm());
        MPI_Gather(&n_locally_owned, 1, MPI_INT, all_locally_owned.data(), 1, MPI_INT, 0, messager::comm());
        max_locally_owned.nodes = gmax.nodes;
        max_locally_owned.rank = gmax.rank;
        min_locally_owned.nodes = gmin.nodes;
        min_locally_owned.rank = gmin.rank;
    }
    // clang-format on
    const double global_load_imbalance =
        static_cast<double>(max_locally_owned.nodes - min_locally_owned.nodes) /
        min_locally_owned.nodes;

    if (messager::master())
    {
        label total_nodes = 0;
        for (const label n : all_locally_owned)
        {
            total_nodes += n;
        }
        std::printf("\nLOAD SUMMARY:\n");
        std::printf("Total number of nodes: %d\n", total_nodes);
    }

    if (messager::parallel())
    {
        if (messager::master())
        {
            for (label rank = 0; rank < size; rank++)
            {
                std::printf("\tNumber of nodes on rank %d: %d\n",
                            rank,
                            all_locally_owned[rank]);
            }
        }
    }

    if (messager::master())
    {
        if (size > 1)
        {
            std::printf("Maximum locally owned nodes: %d (rank %d)\n",
                        max_locally_owned.nodes,
                        max_locally_owned.rank);
            std::printf("Minimum locally owned nodes: %d (rank %d)\n",
                        min_locally_owned.nodes,
                        min_locally_owned.rank);
        }
        std::printf("Load imbalance: %.2f%%\n", global_load_imbalance * 100);
    }
}

#if SPATIAL_DIM == 3
// Element validation and correction
elementValidator& mesh::elementValidatorRef()
{
    if (!elementValidatorPtr_)
    {
        errorMsg("Element validator is not initialized. Check if mesh "
                 "validation is enabled (check_mesh: true).");
    }
    return *elementValidatorPtr_.get();
}

const elementValidator& mesh::elementValidatorRef() const
{
    if (!elementValidatorPtr_)
    {
        errorMsg("Element validator is not initialized. Check if mesh "
                 "validation is enabled (check_mesh: true).");
    }
    return *elementValidatorPtr_.get();
}

elementValidator* mesh::elementValidatorPtr()
{
    return elementValidatorPtr_.get();
}

const elementValidator* mesh::elementValidatorPtr() const
{
    return elementValidatorPtr_.get();
}

// Private element validation methods
void mesh::setupElementValidation_()
{
    if (checkMesh_)
    {
        elementValidatorPtr_ =
            std::unique_ptr<elementValidator, ElementValidatorDeleter>(
                new elementValidator(this));
        elementValidatorPtr_->setup();

        // Configure correction based on YAML setting
        elementValidatorPtr_->enableElementCorrection(enableCorrection_);
    }
}

void mesh::initializeElementValidation_()
{
    if (checkMesh_ && elementValidatorPtr_)
    {
        elementValidatorPtr_->initialize();

        // Perform initial element validation
        elementValidatorPtr_->validateAllElements();

        // Print initial mesh quality report
        if (messager::master())
        {
            elementValidatorPtr_->printValidationReport();
        }

        // Apply element corrections if needed (flash CVFEM pattern)
        // In flash, this is called after mesh quality assessment
        auto correctionResults = elementValidatorPtr_->correctInvalidElements();

        if (!correctionResults.empty())
        {
            // Count actual geometric corrections vs. just marking elements
            label geometricCorrections = 0;
            for (const auto& result : correctionResults)
            {
                if (result.correctionApplied)
                    geometricCorrections++;
            }

            if (geometricCorrections > 0)
            {
                infoMsg("Element correction in progress...");

                // Re-validate after corrections to update statistics
                elementValidatorPtr_->validateAllElements();

                // Print final mesh quality report after corrections
                if (messager::master())
                {
                    infoMsg("Mesh quality after corrections:");
                    elementValidatorPtr_->printValidationReport();
                }
            }
        }
    }
}
#endif // SPATIAL_DIM == 3

} // namespace accel
