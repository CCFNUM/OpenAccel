// File       : realm.cpp
// Created    :
// Author     :
// Description:
// Copyright 2023 CCFNUM HSLU T&A. All Rights Reserved.

// code
#include "realm.h"
#include "domain.h"
#include "mesh.h"
#include "simulation.h"

#include <stk_io/MeshField.hpp>
#include <stk_io/StkMeshIoBroker.hpp>
#include <stk_mesh/base/FieldParallel.hpp>

namespace accel
{

// Constructors

realm::realm(simulation* simulationPtr, std::string name)
    : simulationPtr_(simulationPtr), meshPtr_(simulationPtr_->meshPtr()),
      name_(name), tRealm_(std::make_unique<turbRealm>()),
      smRealm_(std::make_unique<smRealm>())
{
}

// Operations

void realm::initialize()
{
    if (meshPtr_->controlsRef().solverRef().restartControl_.isRestart_)
    {
        restoreWallFields_();
    }
}

void realm::restoreWallFields_()
{
    const auto& restart_ctrl =
        meshPtr_->controlsRef().solverRef().restartControl_;

    scalar restart_time = restart_ctrl.restartTime_;
    if (restart_time == 0.0)
    {
        restart_time = meshPtr_->ioBrokerRef().get_max_time();
    }

    // lagged wall state: last, as the equations rebuild it while initializing
    for (const std::string& name :
         {std::string(U_ID) + "_side", std::string(turbRealm::uWallCoeffs_ID)})
    {
        stk::mesh::FieldBase* theField =
            stk::mesh::get_field_by_name(name, meshPtr_->metaDataRef());
        if (theField == nullptr)
        {
            continue;
        }

        stk::io::MeshField mf(theField,
                              std::to_string(std::hash<std::string>{}(name)),
                              restart_ctrl.timeMatchOption_);
        mf.set_read_time(restart_time);
        mf.set_single_state(false);
        meshPtr_->ioBrokerRef().read_input_field(mf);
        if (mf.field_restored())
        {
            stk::mesh::communicate_field_data(meshPtr_->bulkDataRef(),
                                              {theField});
        }
    }
}

void realm::registerRestartField(const std::string& fieldName)
{
    simulationPtr_->registerRestartField(fieldName);
}

// Access

mesh* realm::meshPtr()
{
    return meshPtr_;
}

const mesh* realm::meshPtr() const
{
    return meshPtr_;
}

mesh& realm::meshRef()
{
    return *meshPtr_;
}

const mesh& realm::meshRef() const
{
    return *meshPtr_;
}

simulation* realm::simulationPtr()
{
    return simulationPtr_;
}

const simulation* realm::simulationPtr() const
{
    return simulationPtr_;
}

simulation& realm::simulationRef()
{
    return *simulationPtr_;
}

const simulation& realm::simulationRef() const
{
    return *simulationPtr_;
}

} // namespace accel
