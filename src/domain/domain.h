// File : domain.h
// Created : Wed Jan 03 2024 13:38:51 (+0100)
// Author : Mhamad Mahdi Alloush
// Description: Physical domain definition with material properties and equation
// models
// Copyright (c) 2024 CCFNUM, Lucerne University of Applied Sciences and
// Arts.
// SPDX-License-Identifier: BSD-3-Clause

#ifndef DOMAIN_H
#define DOMAIN_H

// code
#include "types.h"

namespace accel
{

class realm;
class simulation;
#ifdef HAS_INTERFACE
class interface;
#endif /* HAS_INTERFACE */
class zone;
class mesh;

// define possible fluid model
struct heatTransfer
{
    heatTransferOption option_ = heatTransferOption::none;
    scalar fluidTemperature_ = 300;
    bool includeViscousWork_ = true;
    bool includePressureWork_ = true;
    bool includeLowSpeedCompressibility_ = true;
};

struct turbulence
{
    struct turbulentFluxClosureForHeatTransfer
    {
        turbulentFluxClosureForHeatTransferOption option_ =
            turbulentFluxClosureForHeatTransferOption::eddyDiffusivity;
        scalar turbulentPrandtlNumber_ = 0.9;
    };

    turbulenceOption option_ = turbulenceOption::laminar;
    bool transitional_ = false;
    wallFunctionType wallFunctionType_ = wallFunctionType::standard;
    turbulentFluxClosureForHeatTransfer turbulentFluxClosureForHeatTransfer_;
};

struct multiphase
{
    struct freeSurfaceModel
    {
        freeSurfaceModelOption option_ = freeSurfaceModelOption::none;
        label interfaceCompressionLevel_ = 0;
        bool fluxCorrectedTransport_ = false;
        label nAlphaCorrections_ = 1; // FCT outer corrections
    };

    bool homogeneous_ = true;
    freeSurfaceModel freeSurfaceModel_;
};

struct solidMechanics
{
    solidMechanicsOption option_ = solidMechanicsOption::none;
    kinematicFormulationType formulation_ =
        kinematicFormulationType::totalLagrangian;
    bool planeStress_ = false;
    bool lumpedMass_ =
        true; // true = lumped (diagonal), false = consistent (full mass matrix)
};

struct buoyancy
{
    buoyancyOption option_ = buoyancyOption::nonBuoyant;
    buoyancyModel model_ = buoyancyModel::full;
    std::array<scalar, SPATIAL_DIM> gravity_ = {0};
    scalar referenceDensity_ = 0.0;
    scalar referenceTemperature_ = 0.0;
    std::array<scalar, SPATIAL_DIM> referenceLocation_ = {0};
};

struct fluidPairModel
{
    std::string materialA_;
    std::string materialB_;
    label globalIndexA_ = -1;
    label globalIndexB_ = -1;

    struct surfaceTension
    {
        surfaceTensionModelOption option_ = surfaceTensionModelOption::none;
        scalar coefficient_ = 0.0;
    };

    surfaceTension surfaceTension_;
};

struct material
{
    struct thermodynamicProperties
    {
        struct equationOfState
        {
            equationOfStateOption option_ = equationOfStateOption::value;
            scalar molarMass_ = 1;
        };

        struct specificHeatCapacity
        {
            specificHeatCapacityOption option_ =
                specificHeatCapacityOption::value;

            // coefficients for zero-pressure polynomial (up to 8 coeffs)
            std::array<scalar, 8> coeffs_ =
                {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
        };

        equationOfState equationOfState_;
        specificHeatCapacity specificHeatCapacity_;
    };

    struct transportProperties
    {
        struct dynamicViscosity
        {
            dynamicViscosityOption option_ = dynamicViscosityOption::value;

            // Sutherland's parameters
            scalar referenceTemperature_ = 273.0;
            scalar referenceViscosity_ = 0.0;
            scalar sutherlandsConstant_ = 0.0;
            scalar temperatureExponent_ = 0.0;
        };

        struct thermalConductivity
        {
            thermalConductivityOption option_ =
                thermalConductivityOption::value;

            // Sutherland's parameters
            scalar referenceTemperature_ = 273.0;
            scalar referenceThermalConductivity_ = 0.0;
            scalar sutherlandsConstant_ = 0.0;
            scalar temperatureExponent_ = 0.0;

            // Kinetic theory model parameters
            scalar c1_ = 0.0;
            scalar c2_ = 0.0;
        };

        dynamicViscosity dynamicViscosity_;
        thermalConductivity thermalConductivity_;
    };

    struct mechanicalProperties
    {
        struct youngModulus
        {
            youngModulusOption option_ = youngModulusOption::value;
        };

        struct poissonRatio
        {
            poissonRatioOption option_ = poissonRatioOption::value;
        };

        youngModulus youngModulus_;
        poissonRatio poissonRatio_;
    };

    struct buoyancyProperties
    {
        struct thermalExpansivity
        {
            thermalExpansivityOption option_ = thermalExpansivityOption::value;
        };

        thermalExpansivity thermalExpansivity_;
    };

    std::string name_ = "";
    thermodynamicStateOption thermodynamicStateOption_ =
        thermodynamicStateOption::gas;
    thermodynamicProperties thermodynamicProperties_;
    transportProperties transportProperties_;
    mechanicalProperties mechanicalProperties_;
    buoyancyProperties buoyancyProperties_;
};

struct source
{
    sourceOption option_ = sourceOption::source;
    std::vector<scalar> value_;

    source() = default;

    source(label dim) : value_(dim, 0.0)
    {
    }
};

class domain
{
    const YAML::Node domain_conf_; // reference to this domain config

protected:
    // A pointer to the mesh necessary to access bulk and meta data
    simulation* simulationPtr_;

    // The mesh zone corresponding to this domain: the domain will
    // have a single zone attached to it.
    zone* zonePtr_;

    domainType type_; // physical type of domain

    std::vector<YAML::Node> materialBlockVector_; // materials defined on domain

    std::vector<material> materialVector_;

    std::map<label, label> localToGlobalMaterialIndexMap_;

    std::map<label, label> globalToLocalMaterialIndexMap_;

    YAML::Node initialization_; // initialization data for fields used by domain

    // equations declared on this domain
    std::array<bool, static_cast<size_t>(equationID::numberOfDeclaredEquations)>
        equations_;

    // domain reference pressure

    scalar referencePressure_ = 0.0;

    // domain pressure level information

    bool pressureLevelRequired_ = false;

    stk::mesh::EntityId pressureLevelNodeId_ = stk::mesh::InvalidEntityId;

    label associatedPartitionRankForPressureLevelNode_ = -1;

    void setupPressureLevelInformation_();

    // domain uniform body force
    std::array<scalar, SPATIAL_DIM> uniformBodyForce_{0};

    // equation sources

    source energySource_;

    source momentumSource_;

    void read_();

public:
    domain(simulation* simulation_ptr,
           zone* zonePtr,
           const YAML::Node& domain_conf);

    // Methods

    void setup();

    // Access

    domainType type() const
    {
        return type_;
    }

    label index() const;

    std::string name() const;

    scalar referencePressure() const
    {
        assert(type_ == domainType::fluid);
        return referencePressure_;
    }

    bool pressureLevelRequired() const
    {
        assert(type_ == domainType::fluid);
        return pressureLevelRequired_;
    }

    stk::mesh::EntityId pressureLevelNodeId() const
    {
        assert(type_ == domainType::fluid);
        return pressureLevelNodeId_;
    }

    label associatedPartitionRankForPressureLevelNode() const
    {
        assert(type_ == domainType::fluid);
        return associatedPartitionRankForPressureLevelNode_;
    }

    std::array<scalar, SPATIAL_DIM> uniformBodyForce() const
    {
        assert(type_ == domainType::fluid);
        return uniformBodyForce_;
    }

    void addEquation(const equationID eq_id)
    {
        equations_[static_cast<int>(eq_id)] = true;
    }

    // acess to the realm: contains mesh, settings, etc.

    YAML::Node getYAMLBoundaryConditions() const;

    YAML::Node getYAMLInitialConditions() const;

    YAML::Node getYAMLMaterial(std::string materialName) const;

    // iMaterial is local to the domain
    YAML::Node getYAMLMaterial(label iMaterial = 0) const;

    bool hasEquation(const equationID id) const
    {
        return equations_[static_cast<int>(id)];
    }

#ifdef HAS_INTERFACE
    bool hasInterfaces() const;

    std::vector<interface*>& interfacesRef();

    const std::vector<interface*>& interfacesRef() const;
#endif /* HAS_INTERFACE */

    simulation* simulationPtr();

    const simulation* simulationPtr() const;

    simulation& simulationRef();

    const simulation& simulationRef() const;

    // zones: associated mesh part of this domain

    zone* zonePtr();

    const zone* zonePtr() const;

    zone& zoneRef();

    const zone& zoneRef() const;

    mesh* meshPtr();

    const mesh* meshPtr() const;

    mesh& meshRef();

    const mesh& meshRef() const;

    const source& energySource() const
    {
        return energySource_;
    }

    const source& momentumSource() const
    {
        return momentumSource_;
    }

    // materials

    material& materialRef(label localMaterialIndex = 0)
    {
        return materialVector_[localMaterialIndex];
    }

    const material& materialRef(label localMaterialIndex = 0) const
    {
        return materialVector_[localMaterialIndex];
    }

    const std::vector<material>& materialVector() const
    {
        return materialVector_;
    }

    label nMaterials() const
    {
        return materialVector_.size();
    }

    label localToGlobalMaterialIndex(label localMaterialIndex) const
    {
        return const_cast<std::map<label, label>&>(
            localToGlobalMaterialIndexMap_)[localMaterialIndex];
    }

    label globalToLocalMaterialIndex(label globalMaterialIndex) const
    {
        return const_cast<std::map<label, label>&>(
            globalToLocalMaterialIndexMap_)[globalMaterialIndex];
    }

    bool hasMaterial(std::string materialName)
    {
        for (const auto& material : materialVector_)
        {
            if (material.name_ == materialName)
            {
                return true;
            }
        }

        return false;
    }

    // shortcuts to material properties

    bool isMaterialCompressible(label localMaterialIndex = 0) const
    {
        return this->materialRef(localMaterialIndex)
                   .thermodynamicProperties_.equationOfState_.option_ ==
               equationOfStateOption::idealGas;
    }

    bool isWallDistanceRequired() const;

    // Public members

    // models enabled in the domain

    heatTransfer heatTransfer_;

    turbulence turbulence_;

    multiphase multiphase_;

    solidMechanics solidMechanics_;

    buoyancy buoyancy_;

    std::vector<fluidPairModel> fluidPairModels_;

    const std::vector<fluidPairModel>& fluidPairModels() const
    {
        return fluidPairModels_;
    }
};

} // namespace accel

#endif // DOMAIN_H
