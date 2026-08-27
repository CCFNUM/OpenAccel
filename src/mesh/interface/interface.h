// File       : interface.h
// Created    : Tue Apr 20 2024 12:55:24 (+0100)
// Author     : Mhamad Mahdi Alloush
// Description: Manipulation of an interface
// Copyright 2024 CCFNUM HSLU T&A. All Rights Reserved.

#ifndef INTERFACE_H
#define INTERFACE_H

// code
#include "interfaceSideInfo.h"
#include "mesh.h"

namespace accel
{

class interface
{
public:
    // row-merge position maps, one per graph (column orders differ)
    struct conformalGraphMap
    {
        const ::linearSolver::CRSNodeGraph* graph = nullptr;

        // per pair: stencil-2 column -> stencil-1 column
        std::vector<std::vector<label>> map;

        // per pair: does that column carry an interface-slave unknown that was
        // remapped onto its master partner?  This cannot be inferred from
        // `map[i] == i`: a slave column and its master partner may legitimately
        // occupy the same stencil position, and such a column still needs the
        // right (column) transform applied.
        std::vector<std::vector<uint8_t>> remapped;
    };

protected:
    mesh* meshPtr_;

    std::string name_;

    interfaceModelOption option_;

    nonconformalMethod nonconformalMethod_ =
        nonconformalMethod::discontinuousGalerkin; // default

    interfaceType type_{interfaceType::fluid_fluid}; // default interface type

    label index_ = -1;

    std::pair<std::string, std::string> pairZoneNames_;

    std::pair<label, label> pairZoneIndices_;

    std::unique_ptr<interfaceSideInfo> masterInfoPtr_ = nullptr;

    std::unique_ptr<interfaceSideInfo> slaveInfoPtr_ = nullptr;

    // General aspects of interface

    // is the interface connecting parts of the same domain?
    bool isInternal_ = true;

    bool isOverlap_ = true;

    bool isConformal_ = false;

    bool isForceNonconformalTreatment_ = false;

    bool isSlipNonOverlap_ = true;

    scalar overlapTolerance_ = 1e-3;

    scalar conformalityTolerance_ = 1e-6;

    // Penalty factor multiplier

    scalar penaltyFactor_ = 1.0;

    // Geometric transformations measured from salve to master

    utils::vector axisLocation_;

    utils::vector rotationAxis_;

    // Connectivities - conformal

    std::vector<std::pair<stk::mesh::Entity, stk::mesh::Entity>>
        matchingNodePairVector_;

    // stable id+owner record of matched pairs; survives ghosting churn, used to
    // owner-route ghosting and re-resolve the entity pairs
    struct conformalPairId
    {
        stk::mesh::EntityId id1, id2;
        int owner1, owner2;
    };

    std::vector<conformalPairId> conformalPairIds_;

    mutable std::vector<conformalGraphMap> conformalRowToRowMaps_;

    void initializeGhostings_();

    void updateGhostings_();

    void determineGeometricRelations_();

    const conformalGraphMap&
    populateConformalRowToRowMapping_(const ::linearSolver::CRSNodeGraph* g);

    // ghost each conformal pair's opposing node + its supporting elements to
    // the other side's rank, so split pairs are visible for graph/mapping
    void populateConformalElemsToGhost_();

    // re-resolve matchingNodePairVector_ entity handles from conformalPairIds_
    // (call after the persistent ghosting makes both sides visible)
    void rebuildMatchingPairsFromIds_();

public:
    // Constructors

    interface(mesh* meshPtr, label index, std::string name);

    // Methods

    void read(const YAML::Node& inputNode);

    void setup();

    void initialize();

    void update();

    // Operations

    void provideDiagnosis();

    void errorCheck();

    // Setters and Getters

    void setInterfaceType(interfaceType type)
    {
        type_ = type;
    }

    interfaceType type() const
    {
        return type_;
    }

    void setInternal(bool state)
    {
        isInternal_ = state;
    };

    bool isInternal() const
    {
        return isInternal_;
    };

    scalar penaltyFactor() const
    {
        return penaltyFactor_;
    }

    // Access

    label index() const
    {
        return index_;
    }

    std::string name() const
    {
        return name_;
    };

    interfaceModelOption option() const
    {
        return option_;
    };

    nonconformalMethod ncMethod() const
    {
        return nonconformalMethod_;
    };

    bool isFluidSolidType() const
    {
        return type_ == interfaceType::fluid_solid;
    };

    bool isSlipNonOverlap() const
    {
        return isSlipNonOverlap_;
    };

    bool isOverlap() const
    {
        return isOverlap_;
    };

    bool isConformal() const
    {
        return isConformal_;
    };

    bool isConformalTreatment() const
    {
        return isConformal_ && !isForceNonconformalTreatment_;
    };

    std::vector<std::pair<stk::mesh::Entity, stk::mesh::Entity>>&
    matchingNodePairVector()
    {
        return matchingNodePairVector_;
    }

    const std::vector<std::pair<stk::mesh::Entity, stk::mesh::Entity>>&
    matchingNodePairVector() const
    {
        return matchingNodePairVector_;
    }

    // graphs rebuild in place (mesh motion): drop the stale position maps
    void clearConformalRowToRowMaps() const
    {
        conformalRowToRowMaps_.clear();
    }

    // position map for the given node graph (built lazily, cached per graph)
    const conformalGraphMap&
    conformalRowToRowMap(const ::linearSolver::CRSNodeGraph* graph) const
    {
        for (const auto& entry : conformalRowToRowMaps_)
        {
            if (entry.graph == graph)
            {
                return entry;
            }
        }
        return const_cast<interface*>(this)->populateConformalRowToRowMapping_(
            graph);
    }

    interfaceSideInfo* interfaceSideInfoPtr(label zoneIndex)
    {
        return pairZoneIndices_.first == zoneIndex ? masterInfoPtr_.get()
                                                   : slaveInfoPtr_.get();
    };

    const interfaceSideInfo* interfaceSideInfoPtr(label zoneIndex) const
    {
        return pairZoneIndices_.first == zoneIndex ? masterInfoPtr_.get()
                                                   : slaveInfoPtr_.get();
    };

    interfaceSideInfo* masterInfoPtr()
    {
        return masterInfoPtr_.get();
    };

    const interfaceSideInfo* masterInfoPtr() const
    {
        return masterInfoPtr_.get();
    };

    interfaceSideInfo& masterInfoRef()
    {
        return *masterInfoPtr_.get();
    };

    const interfaceSideInfo& masterInfoRef() const
    {
        return *masterInfoPtr_.get();
    };

    interfaceSideInfo* slaveInfoPtr()
    {
        return slaveInfoPtr_.get();
    };

    const interfaceSideInfo* slaveInfoPtr() const
    {
        return slaveInfoPtr_.get();
    };

    interfaceSideInfo& slaveInfoRef()
    {
        return *slaveInfoPtr_.get();
    };

    const interfaceSideInfo& slaveInfoRef() const
    {
        return *slaveInfoPtr_.get();
    };

    std::pair<std::string, std::string> pairZoneNames() const
    {
        return pairZoneNames_;
    };

    std::string masterZoneName() const
    {
        return pairZoneNames_.first;
    };

    std::string slaveZoneName() const
    {
        return pairZoneNames_.second;
    };

    std::pair<label, label> pairZoneIndices() const
    {
        return pairZoneIndices_;
    };

    label masterZoneIndex() const
    {
        return pairZoneIndices_.first;
    };

    label slaveZoneIndex() const
    {
        return pairZoneIndices_.second;
    };

    bool isMasterZone(label zoneIndex) const
    {
        return pairZoneIndices_.first == zoneIndex;
    }

    mesh* meshPtr();

    const mesh* meshPtr() const;

    mesh& meshRef();

    const mesh& meshRef() const;

    // Public members

    // Ghostings

    stk::mesh::EntityProcVec elemsToGhost_;

    stk::mesh::Ghosting* interfaceGhosting_;

    std::vector<label> ghostCommProcs_;
};

} // namespace accel

#endif // INTERFACE_H
