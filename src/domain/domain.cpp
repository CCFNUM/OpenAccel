// File : domain.cpp
// Created : Wed Jan 03 2024 13:38:51 (+0100)
// Author : Mhamad Mahdi Alloush
// Description:
// Copyright (c) 2024 CCFNUM, Lucerne University of Applied Sciences and Arts.
// SPDX-License-Identifier: BSD-3-Clause

#include "domain.h"
#include "boundary.h"
#include "realm.h"
#include "simulation.h"
#include "zone.h"

namespace accel
{

domain::domain(simulation* simulationPtr,
               zone* zonePtr,
               const YAML::Node& domain_conf)
    : domain_conf_(domain_conf), simulationPtr_(simulationPtr),
      zonePtr_(zonePtr), energySource_(1), momentumSource_(SPATIAL_DIM)
{
    this->read_();
}

// Methods

void domain::setup()
{
    // find reference pressure attributes
    setupPressureLevelInformation_();
}

void domain::setupPressureLevelInformation_()
{
    // Check if domain requires a pressure reference: only applies if domain is
    // a fluid domain
    if (type_ != domainType::fluid)
        return;

    const auto& mesh = meshRef();
    const stk::mesh::MetaData& metaData = mesh.metaDataRef();
    const stk::mesh::BulkData& bulkData = mesh.bulkDataRef();

    // set default to true
    pressureLevelRequired_ = true;

    // general search for pressure boundaries (inlet, outlet, opening)
    for (label iBoundary = 0; iBoundary < this->zonePtr()->nBoundaries();
         iBoundary++)
    {
        boundaryPhysicalType physicalType =
            this->zonePtr()->boundaryRef(iBoundary).type();

        switch (physicalType)
        {
            case boundaryPhysicalType::inlet:
            case boundaryPhysicalType::outlet:
            case boundaryPhysicalType::opening:
                pressureLevelRequired_ = false;
                break;

            default:
                break;
        }
    }

    // return if no pressure reference is required
    if (!pressureLevelRequired_)
        return;

    // now the domain requires a pressure reference ... proceed to finding a
    // reference node

    // check if the user has specified a reference pressure point
    if (this->meshRef()
            .controlsRef()
            .solverRef()
            .solverControl_.advancedOptions_.pressureLevelInformation_
            .option_ ==
        pressureLevelInformationSpecification::cartesianCoordinates)
    {
        auto refLoc = this->meshRef()
                          .controlsRef()
                          .solverRef()
                          .solverControl_.advancedOptions_
                          .pressureLevelInformation_.cartesianCoordinates_;

        // To be added: ensure that the provided reference point is in the
        // active domains BODGEEE

        // reference node to be assigned
        stk::mesh::EntityId g_nodeID = 0;
        label g_owningProc = -1;

        // Determine the nearest node in this processor
        stk::mesh::Entity nearestNode;
        stk::mesh::Selector sel =
            metaData.locally_owned_part() &
            stk::mesh::selectUnion(zonePtr()->interiorParts());
        const auto& buckets =
            bulkData.get_buckets(stk::topology::NODE_RANK, sel);

        const auto* coordsSTKFieldPtr = metaData.get_field<scalar>(
            stk::topology::NODE_RANK, this->meshRef().getCoordinateFieldName());

        scalar distSqr = std::numeric_limits<scalar>::max();
        for (auto b : buckets)
        {
            auto length = b->size();

            for (size_t i = 0; i < length; i++)
            {
                auto node = (*b)[i];
                const scalar* coords =
                    stk::mesh::field_data(*coordsSTKFieldPtr, node);

                scalar dist = 0.0;
                for (label j = 0; j < SPATIAL_DIM; j++)
                {
                    scalar xdiff = (refLoc[j] - coords[j]);
                    dist += xdiff * xdiff;
                }
                if (dist < distSqr)
                {
                    distSqr = dist;
                    nearestNode = node;
                }
            }
        }

        // Determine the global minimum
        std::vector<scalar> minDistList(bulkData.parallel_size());
        MPI_Allgather(&distSqr,
                      1,
                      MPI_DOUBLE,
                      minDistList.data(),
                      1,
                      MPI_DOUBLE,
                      bulkData.parallel());
        scalar minDist = std::numeric_limits<scalar>::max();
        for (label i = 0; i < bulkData.parallel_size(); i++)
        {
            if (minDistList[i] < minDist)
            {
                minDist = minDistList[i];
                g_owningProc = i;
            }
        }

        // Communicate the nearest node ID to all processors.
        stk::mesh::EntityId nodeID = 0;
        if (g_owningProc == bulkData.parallel_rank())
        {
            nodeID = bulkData.identifier(nearestNode);
        }
        stk::all_reduce_max(bulkData.parallel(), &nodeID, &g_nodeID, 1);

        // Set the reference node
        if (pressureLevelRequired_)
        {
            assert(g_nodeID != 0);
            assert(g_owningProc != -1);

            pressureLevelNodeId_ = g_nodeID;
            associatedPartitionRankForPressureLevelNode_ = g_owningProc;
        }
    }
    else
    {
        // We need to select an arbitrary reference node: we search for the
        // first node in each zone as sorted by stk (node id = 1)
        stk::mesh::Selector selOwnedNodes =
            metaData.locally_owned_part() &
            stk::mesh::selectUnion(this->zonePtr()->interiorParts());

        stk::mesh::BucketVector const& nodeBuckets =
            bulkData.get_buckets(stk::topology::NODE_RANK, selOwnedNodes);

        static_assert(std::is_integral_v<stk::mesh::EntityId> &&
                      std::is_unsigned_v<stk::mesh::EntityId>);
        stk::mesh::EntityId lowest_local_node_id = ~0;
        for (stk::mesh::BucketVector::const_iterator ib = nodeBuckets.begin();
             ib != nodeBuckets.end();
             ++ib)
        {
            for (size_t iNode = 0; iNode < (**ib).size(); iNode++)
            {
                const stk::mesh::EntityId nid =
                    bulkData.identifier((**ib)[iNode]);
                lowest_local_node_id =
                    (nid < lowest_local_node_id) ? nid : lowest_local_node_id;
            }
        }

        // 2.) find global lowest node identifier on iZone
        struct STKNodeIDRankPair
        {
            scalar node_id; // MPI does not provide an unsigned type
            // wider than 32 bits, using explicit
            // scalar cast instead
            int rank_id;
        };

        STKNodeIDRankPair local = {
            .node_id = static_cast<scalar>(lowest_local_node_id),
            .rank_id = messager::myProcNo()};
        STKNodeIDRankPair global = {.node_id = local.node_id,
                                    .rank_id = local.rank_id};
        MPI_Allreduce(
            &local, &global, 1, MPI_DOUBLE_INT, MPI_MINLOC, messager::comm());

        pressureLevelNodeId_ = static_cast<stk::mesh::EntityId>(global.node_id);
        associatedPartitionRankForPressureLevelNode_ = global.rank_id;
    }
}

bool domain::isWallDistanceRequired() const
{
    bool state = false;

    if (this->zonePtr()->meshDeforming())
    {
        if (this->zonePtr()
                ->deformationRef()
                .displacementDiffusion()
                .meshStiffnessSpecification_ ==
            meshStiffnessSpecificationType::increaseNearBoundaries)
        {
            state = true;
        }
    }

    if (this->turbulence_.option_ != turbulenceOption::laminar)
    {
        state = true;
    }

    return state;
}

// Access

YAML::Node domain::getYAMLBoundaryConditions() const
{
    if (!domain_conf_["boundaries"])
    {
        errorMsg("domain: `" + this->name() +
                 "` does not define a `boundaries` sequence");
    }
    return domain_conf_["boundaries"];
}

YAML::Node domain::getYAMLInitialConditions() const
{
    if (!initialization_)
    {
        errorMsg("domain: `" + this->name() +
                 "` does not define valid initialization data");
    }
    return initialization_;
}

YAML::Node domain::getYAMLMaterial(std::string materialName) const
{
    for (const auto& materialBlock : materialBlockVector_)
    {
        if (materialBlock["name"].template as<std::string>() == materialName)
        {
            return materialBlock;
        }
    }

    errorMsg("material " + materialName + " does not exist");

    return materialBlockVector_[0];
}

YAML::Node domain::getYAMLMaterial(label iMaterial) const
{
    assert(iMaterial < materialBlockVector_.size());
    return materialBlockVector_[iMaterial == -1 ? 0 : iMaterial];
}

simulation* domain::simulationPtr()
{
    return simulationPtr_;
}

const simulation* domain::simulationPtr() const
{
    return simulationPtr_;
}

simulation& domain::simulationRef()
{
    return *simulationPtr_;
}

const simulation& domain::simulationRef() const
{
    return *simulationPtr_;
}

zone* domain::zonePtr()
{
    return zonePtr_;
}

const zone* domain::zonePtr() const
{
    return zonePtr_;
}

zone& domain::zoneRef()
{
    return *zonePtr_;
}

const zone& domain::zoneRef() const
{
    return *zonePtr_;
}

mesh* domain::meshPtr()
{
    return simulationRef().meshPtr();
}

const mesh* domain::meshPtr() const
{
    return simulationRef().meshPtr();
}

mesh& domain::meshRef()
{
    return simulationRef().meshRef();
}

const mesh& domain::meshRef() const
{
    return simulationRef().meshRef();
}

label domain::index() const
{
    return this->zonePtr()->index();
}

std::string domain::name() const
{
    return this->zonePtr()->name();
}

} // namespace accel
