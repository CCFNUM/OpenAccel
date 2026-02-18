// File : mesh.h
// Created : Fri Aug 25 2023 12:55:24 (+0100)
// Author : Mhamad Mahdi Alloush
// Description: Finite volume mesh container with topology, zones, and geometric
// fields
// Copyright (c) 2023 CCFNUM, Lucerne University of Applied Sciences and
// Arts.
// SPDX-License-Identifier: BSD-3-Clause

#ifndef MESH_H
#define MESH_H

// code
#include "nodeGraph.h"
#include "types.h"

namespace accel
{

class controls;
#ifdef HAS_INTERFACE
class interface;
#endif /* HAS_INTERFACE */
class zone;

#if SPATIAL_DIM == 3
class elementValidator;

// Custom deleter for elementValidator (forward declared)
struct ElementValidatorDeleter
{
    void operator()(elementValidator* ptr);
};
#endif // SPATIAL_DIM == 3

class mesh
{
private:
    // Access to realm: required for more advanced mesh activities

    controls* controlsPtr_;

    // stk containers

    std::shared_ptr<stk::mesh::BulkData> bulkDataPtr_ = nullptr;

    stk::io::StkMeshIoBroker* ioBrokerPtr_ = nullptr;

    stk::mesh::Field<label>* globalIdentityFieldPtr_;

    // Attributes

    label nNodes_ = 0;

    label nShadowNodes_ = 0;

    label nUselessNodes_ = 0;

    label nActiveNodes_ = 0;

    label nAllNodes_ = 0;

    stk::mesh::ConstPartVector interiorActiveParts_;

    stk::mesh::ConstPartVector boundaryActiveParts_;

    stk::mesh::ConstPartVector wallBoundaryActiveParts_;

    stk::mesh::ConstPartVector symmetryBoundaryActiveParts_;

    ::linearSolver::GraphLayout stencil_;

    std::unique_ptr<nodeGraph> globalOrderGraphPtr_ = nullptr;

    std::unique_ptr<nodeGraph> localOrderGraphPtr_ = nullptr;

    bool anyZoneFrameRotating_ = false;

    bool anyZoneMeshTransforming_ = false;

    bool anyZoneMeshDeforming_ = false;

    // a container which maps the local node id (from STK) to the entity. This
    // is because we re-define the local id's of the nodes
    std::vector<stk::mesh::Entity> localNodeIDToEntity_;

#ifdef HAS_INTERFACE
    // Interfaces

    bool hasInterfaces_ = false;

    // all interfaces in the simulation (translational periodic, rotational
    // periodic and general connection)
    std::vector<std::unique_ptr<interface>> interfaceVector_;
#endif /* HAS_INTERFACE */

    // Zones

    // all zones in the mesh
    std::vector<std::unique_ptr<zone>> zoneVector_;

#if SPATIAL_DIM == 3
    // Element validator for quality checking and correction
    std::unique_ptr<elementValidator, ElementValidatorDeleter>
        elementValidatorPtr_ = nullptr;
#endif // SPATIAL_DIM == 3

    // Private methods

    void validateYamlAgainstExodus_(const YAML::Node& inputNode);

    void setupCoordinateField_();

    void setupGeometricFields_();

    void setupZones_();

#ifdef HAS_INTERFACE
    void setupInterfaces_();
#endif /* HAS_INTERFACE */

    void initializeCoordinateField_();

    void initializeZones_();

#ifdef HAS_INTERFACE
    void initializeInterfaces_();
#endif /* HAS_INTERFACE */

    void initializeGeometricFields_();

    void lazyInitializeNodeGraph_();

    void updateZones_(bool force = false);

#ifdef HAS_INTERFACE
    void updateInterfaces_(bool force = false);
#endif /* HAS_INTERFACE */

    void updateGeometricFields_(bool force = false);

    void
    updateDualNodalVolumeField_(stk::mesh::ConstPartVector interiorParts = {});

    void updateExposedAreaVectorField_(
        stk::mesh::ConstPartVector boundaryParts = {});

    void
    updateWallNormalDistanceField_(stk::mesh::ConstPartVector wallParts = {});

    void
    updateAssembledWallAreaField_(stk::mesh::ConstPartVector wallParts = {});

    void updateAssembledSymmetryAreaField_(
        stk::mesh::ConstPartVector symmParts = {});

    void updateAssembledWallNormalDistanceField_(
        stk::mesh::ConstPartVector wallParts = {});

    std::unique_ptr<nodeGraph>
    createNodeGraph_(const ::linearSolver::GraphLayout layout);

    void updateNodeGraph_();

    void initializeLocalNodeIDs_();

    void updateLocalNodeIDs_();

    void reportLoadImbalance_() const;

#if SPATIAL_DIM == 3
    // Element validation private methods
    void setupElementValidation_();
    void initializeElementValidation_();
#endif // SPATIAL_DIM == 3

    // transformation

    bool scaleMesh_ = false;

    std::array<scalar, SPATIAL_DIM> scaleVector_{1.0};

    std::array<scalar, SPATIAL_DIM> scaleOrigin_{0.0};

    bool exportOriginal_ = false;

#if SPATIAL_DIM == 3
    // element validation
    bool checkMesh_ = false;
    bool enableCorrection_ = false;
#endif // SPATIAL_DIM == 3

public:
    static constexpr char coordinates_ID[] = "coordinates";
    static constexpr char original_coordinates_ID[] = "original_coordinates";
    static constexpr char dual_nodal_volume_ID[] = "dual_nodal_volume";
    static constexpr char original_dual_nodal_volume_ID[] =
        "original_dual_nodal_volume";
    static constexpr char exposed_area_vector_ID[] = "exposed_area_vector";
    static constexpr char original_exposed_area_vector_ID[] =
        "original_exposed_area_vector";
    static constexpr char wall_normal_distance_ID[] = "wall_normal_distance";
    static constexpr char assembled_wall_normal_distance_ID[] =
        "assembled_wall_normal_distance";
    static constexpr char assembled_wall_area_ID[] = "assembled_wall_area";
    static constexpr char assembled_symm_area_ID[] = "assembled_symm_area";
    static constexpr char aux[] = "aux";
#ifndef NDEBUG
    static constexpr char rank_ID[] = "rank_id";
    static constexpr char scl_check_ID[] = "scl_check";
#endif /* NDEBUG */

    // Constructors

    mesh(controls* controlsPtr);

    // Destructor (explicit definition needed for unique_ptr with incomplete
    // type)

    ~mesh();

    // Operations

    void read(const YAML::Node& inputNode);

    void reset();

    void setup();

    void initialize();

    void update();

    void scale(std::array<scalar, SPATIAL_DIM> scaleVector);

    void scale(label iZone, std::array<scalar, SPATIAL_DIM> scaleVector);

    void scale(label iZone,
               label iBoundary,
               std::array<scalar, SPATIAL_DIM> scaleVector);

    void write(size_t resultsFileIndex, scalar writeTime);

    // Access

    void setAnyZoneFrameRotating(bool state)
    {
        anyZoneFrameRotating_ = state;
    }

    const bool anyZoneFrameRotating() const
    {
        return anyZoneFrameRotating_;
    }

    void setAnyZoneMeshTransforming(bool state)
    {
        anyZoneMeshTransforming_ = state;
    }

    const bool anyZoneMeshTransforming() const
    {
        return anyZoneMeshTransforming_;
    }

    void setAnyZoneMeshDeforming(bool state)
    {
        anyZoneMeshDeforming_ = state;
    }

    const bool anyZoneMeshDeforming() const
    {
        return anyZoneMeshDeforming_;
    }

    const std::string getCoordinateFieldName() const
    {
        return coordinates_ID;
    }

    controls& controlsRef();

    const controls& controlsRef() const;

    inline stk::mesh::BulkData* bulkDataPtr()
    {
        return bulkDataPtr_.get();
    };

    inline const stk::mesh::BulkData* bulkDataPtr() const
    {
        return bulkDataPtr_.get();
    };

    inline stk::mesh::BulkData& bulkDataRef()
    {
        return *bulkDataPtr_.get();
    };

    inline const stk::mesh::BulkData& bulkDataRef() const
    {
        return *bulkDataPtr_.get();
    };

    inline stk::mesh::MetaData* metaDataPtr()
    {
        return bulkDataRef().mesh_meta_data_ptr().get();
    }

    inline const stk::mesh::MetaData* metaDataPtr() const
    {
        return bulkDataRef().mesh_meta_data_ptr().get();
    }

    inline stk::mesh::MetaData& metaDataRef()
    {
        return bulkDataRef().mesh_meta_data();
    }

    inline const stk::mesh::MetaData& metaDataRef() const
    {
        return bulkDataRef().mesh_meta_data();
    }

    stk::io::StkMeshIoBroker* ioBrokerPtr();

    const stk::io::StkMeshIoBroker* ioBrokerPtr() const;

    stk::io::StkMeshIoBroker& ioBrokerRef();

    const stk::io::StkMeshIoBroker& ioBrokerRef() const;

    stk::mesh::ConstPartVector interiorActiveParts();

    const stk::mesh::ConstPartVector interiorActiveParts() const;

    stk::mesh::ConstPartVector boundaryActiveParts();

    const stk::mesh::ConstPartVector boundaryActiveParts() const;

    stk::mesh::ConstPartVector wallBoundaryActiveParts();

    const stk::mesh::ConstPartVector wallBoundaryActiveParts() const;

    stk::mesh::ConstPartVector symmetryBoundaryActiveParts();

    const stk::mesh::ConstPartVector symmetryBoundaryActiveParts() const;

    inline label nNodes() const
    {
        return nNodes_;
    };

    inline label nShadowNodes() const
    {
        return nShadowNodes_;
    };

    inline label nUselessNodes() const
    {
        return nUselessNodes_;
    };

    inline label nActiveNodes() const
    {
        return nActiveNodes_;
    };

    inline label nAllNodes() const
    {
        return nAllNodes_;
    };

    // Selectors

    stk::mesh::Selector locallyOwnedInteriorPartsSelector() const
    {
        return this->metaDataRef().locally_owned_part() &
               stk::mesh::selectUnion(this->interiorActiveParts());
    };

    stk::mesh::Selector globallySharedInteriorPartsSelector() const
    {
        return this->metaDataRef().globally_shared_part() &
               stk::mesh::selectUnion(this->interiorActiveParts());
    };

    stk::mesh::Selector auraInteriorPartsSelector() const
    {
        return this->metaDataRef().aura_part() &
               stk::mesh::selectUnion(this->interiorActiveParts());
    };

    stk::mesh::Selector universalInteriorPartsSelector() const
    {
        return this->metaDataRef().universal_part() &
               stk::mesh::selectUnion(this->interiorActiveParts());
    };

    stk::mesh::Selector locallyOwnedBoundaryPartsSelector() const
    {
        return this->metaDataRef().locally_owned_part() &
               stk::mesh::selectUnion(this->boundaryActiveParts());
    };

    stk::mesh::Selector globallySharedBoundaryPartsSelector() const
    {
        return this->metaDataRef().globally_shared_part() &
               stk::mesh::selectUnion(this->boundaryActiveParts());
    };

    stk::mesh::Selector auraBoundaryPartsSelector() const
    {
        return this->metaDataRef().aura_part() &
               stk::mesh::selectUnion(this->boundaryActiveParts());
    };

    stk::mesh::Selector universalBoundaryPartsSelector() const
    {
        return this->metaDataRef().universal_part() &
               stk::mesh::selectUnion(this->boundaryActiveParts());
    };

#ifdef HAS_INTERFACE
    // Interfaces

    label nInterfaces() const;

    std::vector<interface*> interfaceVector() const;

    interface& interfaceRef(label iInterface);

    const interface& interfaceRef(label iInterface) const;

    bool hasInterfaces() const
    {
        return hasInterfaces_;
    }
#endif /* HAS_INTERFACE */

    // Zones

    label nZones() const;

    std::vector<zone*> zoneVector() const;

    zone& zoneRef(label iZone);

    const zone& zoneRef(label iZone) const;

    zone* zonePtr(label iZone);

    const zone* zonePtr(label iZone) const;

#if SPATIAL_DIM == 3
    // Element validation and correction
    elementValidator& elementValidatorRef();

    const elementValidator& elementValidatorRef() const;

    elementValidator* elementValidatorPtr();

    const elementValidator* elementValidatorPtr() const;

    // Element validation control
    bool isMeshValidationEnabled() const
    {
        return checkMesh_;
    }

    void enableMeshValidation(bool enable)
    {
        checkMesh_ = enable;
    }

    bool isCorrectionEnabled() const
    {
        return enableCorrection_;
    }

    void enableCorrection(bool enable)
    {
        enableCorrection_ = enable;
    }
#endif // SPATIAL_DIM == 3

    // CRS node graph

    nodeGraph* getGlobalOrderGraphPtr()
    {
        // build graph on demand (lazy)
        if (!globalOrderGraphPtr_)
        {
            globalOrderGraphPtr_ = createNodeGraph_(
                ::linearSolver::GraphLayout::ColumnIndexOrder__Global);

            assert(globalOrderGraphPtr_);
            globalOrderGraphPtr_->buildGraph();
        }
        return globalOrderGraphPtr_.get();
    }

    nodeGraph* getLocalOrderGraphPtr()
    {
        // build graph on demand (lazy)
        if (!localOrderGraphPtr_)
        {
            localOrderGraphPtr_ = createNodeGraph_(
                ::linearSolver::GraphLayout::ColumnIndexOrder__Local);

            assert(localOrderGraphPtr_);
            localOrderGraphPtr_->buildGraph();
        }
        return localOrderGraphPtr_.get();
    }

    std::vector<stk::mesh::Entity>& localNodeIDToEntity()
    {
        return localNodeIDToEntity_;
    }

    const std::vector<stk::mesh::Entity>& localNodeIDToEntity() const
    {
        return localNodeIDToEntity_;
    }
};

} // namespace accel

#endif // MESH_H
