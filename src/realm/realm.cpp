// code
#include "realm.h"
#include "domain.h"
#include "mesh.h"
#include "simulation.h"

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
