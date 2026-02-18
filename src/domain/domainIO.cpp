// File : domain.cpp
// Created : Wed Jan 03 2024 13:38:51 (+0100)
// Author : Mhamad Mahdi Alloush
// Description:
// Copyright (c) 2024 CCFNUM, Lucerne University of Applied Sciences and Arts.
// SPDX-License-Identifier: BSD-3-Clause

// code
#include "domain.h"
#include "macros.h"
#include "realm.h"
#include "simulation.h"

namespace accel
{

// IO

void domain::read_()
{
    // domain primary checks
    assert(domain_conf_["name"].template as<std::string>() == this->name());

    // type of domain
    if (domain_conf_["type"])
    {
        std::string type_name = domain_conf_["type"].template as<std::string>();
        ::accel::tolower(type_name);
        type_ = domainType::undefined;
        if (type_name == "solid")
        {
            type_ = domainType::solid;
        }
        else if (type_name == "fluid")
        {
            type_ = domainType::fluid;
        }
    }
    else
    {
        errorMsg("domain: `type` field not set for domain index " +
                 this->name());
    }

    // NOTE: [2024-06-10] Safety check: currently
    // `domainType::undefined` is not a valid state. Should that be relaxed
    // some time in the future, remove this sanity check.
    if (type_ == domainType::undefined)
    {
        errorMsg("domain: undefined type for domain index" + this->name());
    }

    // material compositions
    if (domain_conf_["materials"])
    {
        // user-defined materials in the domain
        std::vector<std::string> materialNameList =
            domain_conf_["materials"].template as<std::vector<std::string>>();

        // get material_library from simulation node
        YAML::Node materialLibraryBlock =
            simulationPtr_->getYAMLSimulationNode()["material_library"];

        // loop over all user-specified material and register them
        label localMaterialIndex = 0;
        for (std::string materialName : materialNameList)
        {
            // flag to be used
            bool foundMaterial = false;

            ::accel::tolower(materialName);

            for (const auto& materialBlock : materialLibraryBlock)
            {
                if (materialBlock["name"].template as<std::string>() ==
                    materialName)
                {
                    // material found in library: check if already registered
                    for (const auto& material : materialBlockVector_)
                    {
                        if (material["name"].template as<std::string>() ==
                            materialName)
                        {
                            errorMsg("domain: material `" + materialName +
                                     "` set for domain " + this->name() +
                                     " is already defined.");
                        }
                    }

                    // register valid material for this domain
                    materialBlockVector_.push_back(materialBlock);

                    // turn to true
                    foundMaterial = true;

                    // register material in the global map
                    simulationRef().registerMaterial(materialName);

                    // store phase index map
                    localToGlobalMaterialIndexMap_[localMaterialIndex] =
                        simulationRef().materialIndex(materialName);
                    globalToLocalMaterialIndexMap_
                        [simulationRef().materialIndex(materialName)] =
                            localMaterialIndex;

                    // must break
                    break;
                }
            }

            // raise an error
            if (!foundMaterial)
            {
                errorMsg("domain: material `" + materialName +
                         "` set for domain " + this->name() +
                         " does not exist in material library.");
            }

            // store material data

            // decalare default material
            material mat;

            // name
            mat.name_ = materialName;

            // thermodynamic properties
            if (materialBlockVector_[localMaterialIndex]
                                    ["thermodynamic_properties"])
            {
                const auto& thermodynamicPropertiesBlock =
                    materialBlockVector_[localMaterialIndex]
                                        ["thermodynamic_properties"];

                // set equation of state option
                mat.thermodynamicProperties_.equationOfState_
                    .option_ = convertEquationOfStateOptionFromString(
                    thermodynamicPropertiesBlock["equation_of_state"]["option"]
                        .template as<std::string>());

                // ensure that the density is constant value if a solid domain
                if (type_ == domainType::solid)
                {
                    assert(
                        mat.thermodynamicProperties_.equationOfState_.option_ ==
                        equationOfStateOption::value);
                }

                // set molar mass if available
                if (thermodynamicPropertiesBlock["equation_of_state"]
                                                ["molar_mass"])
                {
                    mat.thermodynamicProperties_.equationOfState_.molarMass_ =
                        thermodynamicPropertiesBlock["equation_of_state"]
                                                    ["molar_mass"]
                                                        .template as<scalar>();
                }

                // set specific heat capacity option
                if (thermodynamicPropertiesBlock["specific_heat_capacity"])
                {
                    mat.thermodynamicProperties_.specificHeatCapacity_.option_ =
                        convertSpecificHeatCapacityOptionFromString(
                            thermodynamicPropertiesBlock
                                ["specific_heat_capacity"]["option"]
                                    .template as<std::string>());

                    if (mat.thermodynamicProperties_.specificHeatCapacity_
                            .option_ ==
                        specificHeatCapacityOption::zeroPressurePolynomial)
                    {
                        mat.thermodynamicProperties_.specificHeatCapacity_
                            .coeffs_[0] = thermodynamicPropertiesBlock
                                              ["specific_heat_capacity"]["a1"]
                                                  .template as<scalar>();
                        mat.thermodynamicProperties_.specificHeatCapacity_
                            .coeffs_[1] = thermodynamicPropertiesBlock
                                              ["specific_heat_capacity"]["a2"]
                                                  .template as<scalar>();
                        mat.thermodynamicProperties_.specificHeatCapacity_
                            .coeffs_[2] = thermodynamicPropertiesBlock
                                              ["specific_heat_capacity"]["a3"]
                                                  .template as<scalar>();
                        mat.thermodynamicProperties_.specificHeatCapacity_
                            .coeffs_[3] = thermodynamicPropertiesBlock
                                              ["specific_heat_capacity"]["a4"]
                                                  .template as<scalar>();
                        mat.thermodynamicProperties_.specificHeatCapacity_
                            .coeffs_[4] = thermodynamicPropertiesBlock
                                              ["specific_heat_capacity"]["a5"]
                                                  .template as<scalar>();

                        if (thermodynamicPropertiesBlock
                                ["specific_heat_capacity"]["a6"])
                        {
                            mat.thermodynamicProperties_.specificHeatCapacity_
                                .coeffs_[5] =
                                thermodynamicPropertiesBlock
                                    ["specific_heat_capacity"]["a6"]
                                        .template as<scalar>();

                            if (thermodynamicPropertiesBlock
                                    ["specific_heat_capacity"]["a7"])
                            {
                                mat.thermodynamicProperties_
                                    .specificHeatCapacity_.coeffs_[6] =
                                    thermodynamicPropertiesBlock
                                        ["specific_heat_capacity"]["a7"]
                                            .template as<scalar>();

                                if (thermodynamicPropertiesBlock
                                        ["specific_heat_capacity"]["a8"])
                                {
                                    mat.thermodynamicProperties_
                                        .specificHeatCapacity_.coeffs_[7] =
                                        thermodynamicPropertiesBlock
                                            ["specific_heat_capacity"]["a8"]
                                                .template as<scalar>();
                                }
                            }
                        }
                    }
                }
            }

            // transport properties
            if (materialBlockVector_[localMaterialIndex]
                                    ["transport_properties"])
            {
                const auto& transportPropertiesBlock =
                    materialBlockVector_[localMaterialIndex]
                                        ["transport_properties"];

                // dynamic viscosity only for fluid domains
                if (type_ == domainType::fluid)
                {
                    // set dynamic viscosity option
                    mat.transportProperties_.dynamicViscosity_
                        .option_ = convertDynamicViscosityOptionFromString(
                        transportPropertiesBlock["dynamic_viscosity"]["option"]
                            .template as<std::string>());

                    if (mat.transportProperties_.dynamicViscosity_.option_ ==
                        dynamicViscosityOption::sutherlandsFormula)
                    {
                        mat.transportProperties_.dynamicViscosity_
                            .referenceTemperature_ =
                            transportPropertiesBlock["dynamic_viscosity"]
                                                    ["reference_temperature"]
                                                        .template as<scalar>();
                        mat.transportProperties_.dynamicViscosity_
                            .referenceViscosity_ =
                            transportPropertiesBlock["dynamic_viscosity"]
                                                    ["reference_viscosity"]
                                                        .template as<scalar>();
                        mat.transportProperties_.dynamicViscosity_
                            .sutherlandsConstant_ =
                            transportPropertiesBlock["dynamic_viscosity"]
                                                    ["sutherlands_constant"]
                                                        .template as<scalar>();
                        mat.transportProperties_.dynamicViscosity_
                            .temperatureExponent_ =
                            transportPropertiesBlock["dynamic_viscosity"]
                                                    ["temperature_exponent"]
                                                        .template as<scalar>();
                    }
                }

                // set thermal conductivity option
                if (transportPropertiesBlock["thermal_conductivity"])
                {
                    mat.transportProperties_.thermalConductivity_.option_ =
                        convertThermalConductivityOptionFromString(
                            transportPropertiesBlock
                                ["thermal_conductivity"]["option"]
                                    .template as<std::string>());

                    if (mat.transportProperties_.thermalConductivity_.option_ ==
                        thermalConductivityOption::sutherlandsFormula)
                    {
                        mat.transportProperties_.thermalConductivity_
                            .referenceTemperature_ =
                            transportPropertiesBlock["thermal_conductivity"]
                                                    ["reference_temperature"]
                                                        .template as<scalar>();
                        mat.transportProperties_.thermalConductivity_
                            .referenceThermalConductivity_ =
                            transportPropertiesBlock
                                ["thermal_conductivity"]
                                ["reference_thermal_conductivity"]
                                    .template as<scalar>();
                        mat.transportProperties_.thermalConductivity_
                            .sutherlandsConstant_ =
                            transportPropertiesBlock["thermal_conductivity"]
                                                    ["sutherlands_constant"]
                                                        .template as<scalar>();
                        mat.transportProperties_.thermalConductivity_
                            .temperatureExponent_ =
                            transportPropertiesBlock["thermal_conductivity"]
                                                    ["temperature_exponent"]
                                                        .template as<scalar>();
                    }
                    else if (mat.transportProperties_.thermalConductivity_
                                 .option_ ==
                             thermalConductivityOption::kineticTheoryModel)
                    {
                        mat.transportProperties_.thermalConductivity_.c1_ =
                            transportPropertiesBlock["thermal_conductivity"]
                                                    ["c1"]
                                                        .template as<scalar>();
                        mat.transportProperties_.thermalConductivity_.c2_ =
                            transportPropertiesBlock["thermal_conductivity"]
                                                    ["c2"]
                                                        .template as<scalar>();
                    }
                }
            }

            // buoyancy properties
            if (materialBlockVector_[localMaterialIndex]["buoyancy_properties"])
            {
                const auto& buoyancyPropertiesBlock =
                    materialBlockVector_[localMaterialIndex]
                                        ["buoyancy_properties"];

                // set thermal conductivity option
                if (buoyancyPropertiesBlock["thermal_expansivity"])
                {
                    mat.buoyancyProperties_.thermalExpansivity_
                        .option_ = convertThermalExpansivityOptionFromString(
                        buoyancyPropertiesBlock["thermal_expansivity"]["option"]
                            .template as<std::string>());
                }
            }

            // mechanical properties
            if (materialBlockVector_[localMaterialIndex]
                                    ["mechanical_properties"])
            {
                const auto& mechanicalPropertiesBlock =
                    materialBlockVector_[localMaterialIndex]
                                        ["mechanical_properties"];

                mat.mechanicalProperties_.youngModulus_.option_ =
                    convertYoungModulusOptionFromString(
                        mechanicalPropertiesBlock["young_modulus"]["option"]
                            .template as<std::string>());

                mat.mechanicalProperties_.poissonRatio_.option_ =
                    convertPoissonRatioOptionFromString(
                        mechanicalPropertiesBlock["poisson_ratio"]["option"]
                            .template as<std::string>());
            }

            // move ownership to material vector
            materialVector_.push_back(std::move(mat));

            // increment
            localMaterialIndex++;
        }
    }
    else
    {
        errorMsg("domain: `materials` field not set for domain " +
                 this->name());
    }

    // domain specific initialization
    if (domain_conf_["initialization"])
    {
        if (!initialization_)
        {
            // if no default initialization defined -> use this
            initialization_ = domain_conf_["initialization"]; // reference
        }
        else
        {
            // update customized nodes
            using Iter = YAML::iterator;
            YAML::Node custom = domain_conf_["initialization"];
            for (Iter item = custom.begin(); item != custom.end(); ++item)
            {
                const std::string key = item->first.template as<std::string>();
                YAML::Node value = item->second;
                initialization_[key] = value;
            }
        }
    }
    else if (simulationPtr_
                 ->getYAMLSimulationNode()["physical_analysis"]["domains"]
                                          ["initialization"])
    {
        // we refer to the default global initialization
        YAML::Node default_initialization =
            simulationPtr_
                ->getYAMLSimulationNode()["physical_analysis"]["domains"]
                                         ["initialization"];

        // default initialization
        initialization_.reset(YAML::Clone(default_initialization));
    }
    else
    {
        errorMsg("Initialization block from domain " + this->name() +
                 " not provided");
    }

    // set all equations to false initially
    for (size_t i = 0; i < equations_.size(); i++)
    {
        equations_[i] = false;
    }

    // check for valid initialization data is currently done in
    // domain::getYAMLInitialConditions() -- if a valid domain can be defined
    // without initialization data in a future extension, the check in
    // domain::getYAMLInitialConditions() can be removed.
    if (type_ == domainType::fluid)
    {
        if (!domain_conf_["fluid_models"])
        {
            errorMsg("fluid_models block not provided in fluid domain " +
                     name());
        }

        const auto& fluidModelsBlock = domain_conf_["fluid_models"];

        // flow model: always ON if no multiphase model
        equations_[static_cast<int>(equationID::segregatedFlow)] = true;

        // query multiphase model: must exist in case more than a material is
        // assigned
        if (this->nMaterials() > 1)
        {
            if (fluidModelsBlock["multiphase"])
            {
                const auto& multiphaseBlock = fluidModelsBlock["multiphase"];

                if (multiphaseBlock["homogeneous"])
                {
                    multiphase_.homogeneous_ =
                        multiphaseBlock["homogeneous"].template as<bool>();

                    if (!multiphase_.homogeneous_)
                    {
                        errorMsg("Inhomogeneous multiphase is not supported");
                    }
                }
                else
                {
                    errorMsg("option is not provided for multiphase block in "
                             "fluid models of domain");
                }

                if (multiphaseBlock["free_surface_model"])
                {
                    const auto& freeSurfaceModelBlock =
                        multiphaseBlock["free_surface_model"];

                    multiphase_.freeSurfaceModel_.option_ =
                        convertFreeSurfaceModelOptionFromString(
                            freeSurfaceModelBlock["option"]
                                .template as<std::string>());

                    if (multiphase_.freeSurfaceModel_.option_ !=
                        freeSurfaceModelOption::standard)
                    {
                        errorMsg("Only free surface model is supported for "
                                 "multiphase");
                    }
                    else
                    {
                        equations_[static_cast<int>(
                            equationID::segregatedFreeSurface)] = true;

                        // turn off single-phase flow solver in the domain
                        equations_[static_cast<int>(
                            equationID::segregatedFlow)] = false;

                        // assert there are more than one phase in the domain
                        assert(nMaterials() > 1);

                        if (freeSurfaceModelBlock
                                ["interface_compression_level"])
                        {
                            multiphase_.freeSurfaceModel_
                                .interfaceCompressionLevel_ =
                                freeSurfaceModelBlock
                                    ["interface_compression_level"]
                                        .template as<label>();
                        }

                        if (freeSurfaceModelBlock["flux_corrected_transport"])
                        {
                            multiphase_.freeSurfaceModel_
                                .fluxCorrectedTransport_ =
                                freeSurfaceModelBlock
                                    ["flux_corrected_transport"]
                                        .template as<bool>();
                        }

                        if (freeSurfaceModelBlock["n_alpha_corrections"])
                        {
                            multiphase_.freeSurfaceModel_.nAlphaCorrections_ =
                                freeSurfaceModelBlock["n_alpha_corrections"]
                                    .template as<label>();
                        }
                    }
                }
                else
                {
                    errorMsg("free surface model option is not provided for "
                             "multiphase block in "
                             "fluid models of domain");
                }
            }
            else
            {
                errorMsg("multiphase model must be enabled when more than a "
                         "fluid exist in the domain");
            }
        }

        // query heat transfer model
        if (fluidModelsBlock["heat_transfer"])
        {
            const auto& heatTransferBlock = fluidModelsBlock["heat_transfer"];

            if (heatTransferBlock["option"])
            {
                heatTransfer_.option_ = convertHeatTransferOptionFromString(
                    heatTransferBlock["option"].template as<std::string>());

                switch (heatTransfer_.option_)
                {
                    case heatTransferOption::none:
                        break;

                    case heatTransferOption::isothermal:
                        {
                            // query the fluid temperature
                            heatTransfer_.fluidTemperature_ =
                                heatTransferBlock["fluid_temperature"]
                                    .template as<scalar>();
                        }
                        break;

                    case heatTransferOption::thermalEnergy:
                        {
                            equations_[static_cast<int>(
                                equationID::thermalEnergy)] = true;
                        }
                        break;

                    case heatTransferOption::totalEnergy:
                        {
                            equations_[static_cast<int>(
                                equationID::totalEnergy)] = true;
                        }
                        break;
                }
            }
            else
            {
                errorMsg("option is not provided for heat transfer model in "
                         "fluid models of domain");
            }

            if (heatTransferBlock["include_viscous_work"])
            {
                heatTransfer_.includeViscousWork_ =
                    heatTransferBlock["include_viscous_work"]
                        .template as<bool>();
            }

            if (heatTransferBlock["include_pressure_work"])
            {
                heatTransfer_.includePressureWork_ =
                    heatTransferBlock["include_pressure_work"]
                        .template as<bool>();
            }

            if (heatTransferBlock["include_low_speed_compressibility"])
            {
                heatTransfer_.includeLowSpeedCompressibility_ =
                    heatTransferBlock["include_low_speed_compressibility"]
                        .template as<bool>();
            }
        }

        // ensure that a heat transfer model is enabled in case at least one
        // material in the domain exhibits a non-const density
        for (const auto& material : materialVector_)
        {
            if (type_ == domainType::fluid &&
                material.thermodynamicProperties_.equationOfState_.option_ !=
                    equationOfStateOption::value)
            {
                if (heatTransfer_.option_ == heatTransferOption::none)
                {
                    errorMsg("For non-const density material, a heat transfer "
                             "model must be enabled");
                }
            }
        }

        // query turbulence model
        if (fluidModelsBlock["turbulence"])
        {
            const auto& turbulenceBlock = fluidModelsBlock["turbulence"];

            if (turbulenceBlock["option"])
            {
                turbulence_.option_ = convertTurbulenceOptionFromString(
                    turbulenceBlock["option"].template as<std::string>());
            }
            else
            {
                errorMsg("option is not provided for turbulence model in fluid "
                         "models of domain");
            }

            switch (turbulence_.option_)
            {
                case turbulenceOption::laminar:
                    break;

                case turbulenceOption::kEpsilon:
                    {
                        equations_[static_cast<int>(
                            equationID::segregatedKEpsilon)] = true;

                        // query wall function
                        if (turbulenceBlock["wall_function"])
                        {
                            turbulence_.wallFunctionType_ =
                                convertWallFunctionTypeFromString(
                                    turbulenceBlock["wall_function"]
                                        .template as<std::string>());

                            // ensure the wall function is scalable
                            if (turbulence_.wallFunctionType_ !=
                                wallFunctionType::scalable)
                            {
                                errorMsg("For k-epsilon model, only scalable "
                                         "wall function is allowed");
                            }
                        }
                        else
                        {
                            // default to scalable wall function
                            turbulence_.wallFunctionType_ =
                                wallFunctionType::scalable;
                        }
                    }
                    break;

                case turbulenceOption::shearStressTransport:
                    {
                        // query whether transitional SST
                        if (turbulenceBlock["transitional_turbulence"])
                        {
                            if (turbulenceBlock["transitional_turbulence"]
                                    .template as<bool>() == true)
                            {
                                turbulence_.transitional_ = true;

                                // query whether correlation-based transition
                                // model (Menter 2015) or full gamma-ReTheta
                                // model (Langtry-Menter 2009)
                                if (turbulenceBlock["correlation_based"])
                                {
                                    if (turbulenceBlock["correlation_based"]
                                            .template as<bool>() == true)
                                    {
                                        // correlation-based transition model
                                        // (Menter 2015) - solves only gamma
                                        equations_[static_cast<int>(
                                            equationID::
                                                segregatedCorrelationTransitionShearStressTransport)] =
                                            true;
                                    }
                                    else
                                    {
                                        // full transition model
                                        // (Langtry-Menter 2009) - solves gamma
                                        // and ReTheta
                                        equations_[static_cast<int>(
                                            equationID::
                                                segregatedTransitionShearStressTransport)] =
                                            true;
                                    }
                                }
                                else
                                {
                                    // default to full transition model
                                    equations_[static_cast<int>(
                                        equationID::
                                            segregatedTransitionShearStressTransport)] =
                                        true;
                                }
                            }
                            else
                            {
                                equations_[static_cast<int>(
                                    equationID::
                                        segregatedShearStressTransport)] = true;
                            }
                        }
                        else
                        {
                            equations_[static_cast<int>(
                                equationID::segregatedShearStressTransport)] =
                                true;
                        }

                        // query wall function
                        if (turbulenceBlock["wall_function"])
                        {
                            turbulence_.wallFunctionType_ =
                                convertWallFunctionTypeFromString(
                                    turbulenceBlock["wall_function"]
                                        .template as<std::string>());

                            // ensure the wall function is scalable
                            if (turbulence_.wallFunctionType_ !=
                                wallFunctionType::automatic)
                            {
                                errorMsg("For SST model, only automatic wall "
                                         "function is allowed");
                            }
                        }
                        else
                        {
                            // default to automatic wall function
                            turbulence_.wallFunctionType_ =
                                wallFunctionType::automatic;
                        }
                    }
                    break;
            }

            // if heat transfer model enabled, a turbulent flux closure for heat
            // transfer is required
            if (turbulence_.option_ != turbulenceOption::laminar &&
                (heatTransfer_.option_ == heatTransferOption::thermalEnergy ||
                 heatTransfer_.option_ == heatTransferOption::totalEnergy))
            {
                if (turbulenceBlock["turbulent_flux_closure_for_heat_transfer"])
                {
                    const auto& block = turbulenceBlock
                        ["turbulent_flux_closure_for_heat_transfer"];

                    turbulence_.turbulentFluxClosureForHeatTransfer_.option_ =
                        convertTurbulentFluxClosureForHeatTransferOptionFromString(
                            block["option"].template as<std::string>());

                    switch (turbulence_.turbulentFluxClosureForHeatTransfer_
                                .option_)
                    {
                        case turbulentFluxClosureForHeatTransferOption::
                            eddyDiffusivity:
                            {
                                turbulence_.turbulentFluxClosureForHeatTransfer_
                                    .turbulentPrandtlNumber_ =
                                    block["turbulent_prandtl_number"]
                                        .template as<scalar>();
                            }
                            break;
                    }
                }
            }
        }

        // query fluid pair models (surface tension between fluid pairs)
        if (domain_conf_["fluid_pair_models"])
        {
            const auto& fluidPairModelsBlock =
                domain_conf_["fluid_pair_models"];

            for (const auto& pairBlock : fluidPairModelsBlock)
            {
                fluidPairModel fpm;

                // read pair: [materialA, materialB]
                if (!pairBlock["pair"])
                {
                    errorMsg("fluid_pair_models: `pair` key is required");
                }

                auto pairList =
                    pairBlock["pair"].template as<std::vector<std::string>>();
                if (pairList.size() != 2)
                {
                    errorMsg("fluid_pair_models: `pair` must contain exactly "
                             "two material names");
                }

                ::accel::tolower(pairList[0]);
                ::accel::tolower(pairList[1]);
                fpm.materialA_ = pairList[0];
                fpm.materialB_ = pairList[1];

                // validate both materials exist in this domain
                if (!this->hasMaterial(fpm.materialA_))
                {
                    errorMsg("fluid_pair_models: material `" + fpm.materialA_ +
                             "` not found in domain " + this->name());
                }
                if (!this->hasMaterial(fpm.materialB_))
                {
                    errorMsg("fluid_pair_models: material `" + fpm.materialB_ +
                             "` not found in domain " + this->name());
                }

                // resolve global material indices
                fpm.globalIndexA_ =
                    simulationRef().materialIndex(fpm.materialA_);
                fpm.globalIndexB_ =
                    simulationRef().materialIndex(fpm.materialB_);

                // validate no duplicate pairs (order-independent)
                for (const auto& existing : fluidPairModels_)
                {
                    if ((existing.globalIndexA_ == fpm.globalIndexA_ &&
                         existing.globalIndexB_ == fpm.globalIndexB_) ||
                        (existing.globalIndexA_ == fpm.globalIndexB_ &&
                         existing.globalIndexB_ == fpm.globalIndexA_))
                    {
                        errorMsg("fluid_pair_models: duplicate pair [" +
                                 fpm.materialA_ + ", " + fpm.materialB_ +
                                 "] in domain " + this->name());
                    }
                }

                // read surface tension model
                if (pairBlock["surface_tension"])
                {
                    const auto& stBlock = pairBlock["surface_tension"];

                    if (stBlock["option"])
                    {
                        fpm.surfaceTension_.option_ =
                            convertSurfaceTensionModelOptionFromString(
                                stBlock["option"].template as<std::string>());
                    }

                    if (stBlock["surface_tension_coefficient"])
                    {
                        fpm.surfaceTension_.coefficient_ =
                            stBlock["surface_tension_coefficient"]
                                .template as<scalar>();
                    }
                }

                fluidPairModels_.push_back(std::move(fpm));
            }
        }
    }
    else if (type_ == domainType::solid)
    {
        if (!domain_conf_["solid_models"])
        {
            errorMsg("solid_models block not provided in solid domain " +
                     name());
        }

        const auto& solidModelsBlock = domain_conf_["solid_models"];

        // query heat transfer model
        if (solidModelsBlock["heat_transfer"])
        {
            const auto& heatTransferBlock = solidModelsBlock["heat_transfer"];

            if (heatTransferBlock["option"])
            {
                heatTransfer_.option_ = convertHeatTransferOptionFromString(
                    heatTransferBlock["option"].template as<std::string>());

                switch (heatTransfer_.option_)
                {
                    case heatTransferOption::none:
                        break;

                    case heatTransferOption::isothermal:
                        {
                            // query the fluid temperature
                            heatTransfer_.fluidTemperature_ =
                                heatTransferBlock["fluid_temperature"]
                                    .template as<scalar>();
                        }
                        break;

                    case heatTransferOption::thermalEnergy:
                        {
                            equations_[static_cast<int>(
                                equationID::thermalEnergy)] = true;
                        }
                        break;

                    case heatTransferOption::totalEnergy:
                        {
                            equations_[static_cast<int>(
                                equationID::totalEnergy)] = true;
                        }
                        break;
                }
            }
            else
            {
                errorMsg("option is not provided for heat transfer model in "
                         "solid models of domain");
            }
        }

        // query solid mechanics model
        if (solidModelsBlock["solid_mechanics"])
        {
            const auto& solidMechanicsBlock =
                solidModelsBlock["solid_mechanics"];

            if (solidMechanicsBlock["option"])
            {
                solidMechanics_.option_ = convertSolidMechanicsOptionFromString(
                    solidMechanicsBlock["option"].template as<std::string>());
            }

            if (solidMechanicsBlock["formulation"])
            {
                solidMechanics_.formulation_ =
                    convertKinematicFormulationTypeFromString(
                        solidMechanicsBlock["formulation"]
                            .template as<std::string>());

                switch (solidMechanics_.formulation_)
                {
                    case kinematicFormulationType::totalLagrangian:
                        {
                            if (this->zonePtr()
                                    ->meshPtr()
                                    ->controlsRef()
                                    .isTransient())
                            {
                                // set deforming zone
                                this->zonePtr()->setMeshDeforming(true);
                                this->zonePtr()
                                    ->meshPtr()
                                    ->setAnyZoneMeshDeforming(true);

                                // call deformation
                                auto& deformation =
                                    this->zonePtr()->deformationRef();

                                // force inherent mesh deformation
                                deformation.setSpecification(
                                    meshDeformationSpecificationType::inherent);

                                // force total displacement
                                deformation
                                    .setDisplacementRelativeToPreviousMesh(
                                        false);
                            }

                            equations_[static_cast<int>(
                                equationID::solidDisplacement)] = true;
                        }
                        break;

                    case kinematicFormulationType::updatedLagrangian:
                        {
                            errorMsg("updated_lagrangian not implemented yet");

                            if (this->zonePtr()
                                    ->meshPtr()
                                    ->controlsRef()
                                    .isTransient())
                            {
                                // set deforming zone
                                this->zonePtr()->setMeshDeforming(true);
                                this->zonePtr()
                                    ->meshPtr()
                                    ->setAnyZoneMeshDeforming(true);

                                // call deformation
                                auto& deformation =
                                    this->zonePtr()->deformationRef();

                                // force inherent mesh deformation
                                deformation.setSpecification(
                                    meshDeformationSpecificationType::inherent);

                                // force relative displacement
                                deformation
                                    .setDisplacementRelativeToPreviousMesh(
                                        true);
                            }
                        }
                        break;
                }
            }

            if (solidMechanicsBlock["plane_stress"])
            {
                solidMechanics_.planeStress_ =
                    solidMechanicsBlock["plane_stress"].template as<bool>();
            }

            if (solidMechanicsBlock["lumped_mass"])
            {
                solidMechanics_.lumpedMass_ =
                    solidMechanicsBlock["lumped_mass"].template as<bool>();

                if (!solidMechanics_.lumpedMass_)
                {
                    errorMsg("false lumped mass -> consistent mass appraoch is "
                             "not implemented");
                }
            }
        }
    }
    else
    {
        errorMsg("No models are provided for domain `" + this->name() + "`");
    }

    // domain models
    if (domain_conf_["domain_models"])
    {
        const auto& domainModelsBlock = domain_conf_["domain_models"];

        if (domainModelsBlock["reference_pressure"])
        {
            // ensure this is a fluid domain
            if (type_ != domainType::fluid)
            {
                errorMsg("assigning reference_pressure to a solid domain is "
                         "invalid");
            }

            referencePressure_ =
                domainModelsBlock["reference_pressure"].template as<scalar>();
        }

        if (domainModelsBlock["uniform_body_force"])
        {
            uniformBodyForce_ =
                domainModelsBlock["uniform_body_force"]
                    .template as<std::array<scalar, SPATIAL_DIM>>();
        }

        // query buoyancy model
        if (domainModelsBlock["buoyancy_model"])
        {
            // ensure this is a fluid domain
            if (type_ != domainType::fluid)
            {
                errorMsg(
                    "assigning buoyancy_model to a solid domain is invalid");
            }

            const auto& buoyancyModelBlock =
                domainModelsBlock["buoyancy_model"];

            if (buoyancyModelBlock["option"])
            {
                buoyancy_.option_ = convertBuoyancyOptionFromString(
                    buoyancyModelBlock["option"].template as<std::string>());

                switch (buoyancy_.option_)
                {
                    case buoyancyOption::nonBuoyant:
                        break;

                    case buoyancyOption::buoyant:
                        {
                            if (buoyancyModelBlock["gravity"])
                            {
                                buoyancy_.gravity_ =
                                    buoyancyModelBlock["gravity"]
                                        .template as<
                                            std::array<scalar, SPATIAL_DIM>>();
                            }
                            else
                            {
                                errorMsg("gravity key is not provided for "
                                         "buoyancy model");
                            }

                            if (nMaterials() > 1 ||
                                this->isMaterialCompressible())
                            {
                                buoyancy_.model_ = buoyancyModel::full;
                            }
                            else
                            {
                                buoyancy_.model_ = buoyancyModel::boussinesq;

                                // ensure heat transfer is enabled
                                if (heatTransfer_.option_ ==
                                    heatTransferOption::none)
                                {
                                    errorMsg("Boussinesq buoyancy model "
                                             "requires heat transfer enabled");
                                }
                            }

                            switch (buoyancy_.model_)
                            {
                                case buoyancyModel::full:
                                    {
                                        if (buoyancyModelBlock
                                                ["buoyancy_reference_density"])
                                        {
                                            buoyancy_.referenceDensity_ =
                                                buoyancyModelBlock
                                                    ["buoyancy_reference_"
                                                     "density"]
                                                        .template as<scalar>();
                                        }
                                        else
                                        {
                                            errorMsg("buoyancy_reference_"
                                                     "density key is "
                                                     "not provided for "
                                                     "buoyancy model");
                                        }

                                        if (buoyancyModelBlock
                                                ["reference_location"])
                                        {
                                            buoyancy_.referenceLocation_ =
                                                buoyancyModelBlock
                                                    ["reference_location"]
                                                        .template as<std::array<
                                                            scalar,
                                                            SPATIAL_DIM>>();
                                        }
                                    }
                                    break;

                                case buoyancyModel::boussinesq:
                                    {
                                        // ensure that the thermal energy
                                        // equation is turned on
                                        assert(equations_[static_cast<int>(
                                            equationID::thermalEnergy)]);

                                        if (buoyancyModelBlock
                                                ["buoyancy_reference_"
                                                 "temperature"])
                                        {
                                            buoyancy_.referenceTemperature_ =
                                                buoyancyModelBlock
                                                    ["buoyancy_reference_"
                                                     "temperature"]
                                                        .template as<scalar>();
                                        }
                                        else
                                        {
                                            errorMsg("buoyancy_reference_"
                                                     "density key is "
                                                     "not provided for "
                                                     "buoyancy model");
                                        }
                                    }
                                    break;
                            }
                        }
                        break;
                };
            }
        }

        // query mesh deformation
        if (domainModelsBlock["mesh_deformation"])
        {
            // this information must have already been collected in the
            // zone-read stage. Here, we only activate the displacement
            // diffusion equation. This equation will not be solved with the
            // global set of equations, but rather, will be solved as a
            // standalone in the mesh motion
            const auto& meshDeformation = this->zonePtr()->deformationRef();

            switch (meshDeformation.specification())
            {
                case meshDeformationSpecificationType::none:
                case meshDeformationSpecificationType::inherent:
                    break;

                case meshDeformationSpecificationType::regionsOfMotionSpecified:
                    {
                        // ensure this is a fluid domain
                        if (type_ != domainType::fluid)
                        {
                            errorMsg("assigning mesh_deformation specification "
                                     "to a solid domain as a "
                                     "region_of_motion_specified is invalid");
                        }

                        equations_[static_cast<int>(
                            equationID::displacementDiffusion)] = true;
                    }
                    break;
            }
        }
    }

    // equation sources
    if (domain_conf_["sources"])
    {
        const auto& sourceBlocks = domain_conf_["sources"];
        if (sourceBlocks["energy"])
        {
            // read value
            const auto& energySourceBlock = sourceBlocks["energy"];
            if (energySourceBlock["option"])
            {
                energySource_.option_ = convertSourceOptionFromString(
                    energySourceBlock["option"].template as<std::string>());
                switch (energySource_.option_)
                {
                    case sourceOption::source:
                        {
                            scalar sourceValue = energySourceBlock["source"]
                                                     .template as<scalar>();
                            energySource_.value_[0] = sourceValue;
                        }
                        break;

                    case sourceOption::totalSource:
                        {
                            scalar totalSourceValue =
                                energySourceBlock["total_source"]
                                    .template as<scalar>();
                            energySource_.value_[0] = totalSourceValue;
                        }
                        break;
                }
            }
            else
            {
                errorMsg(
                    "option key not provided for energy source in domain " +
                    this->name());
            }
        }
        else if (sourceBlocks["momentum"])
        {
            errorMsg("momentum source for domain " + this->name() +
                     " not implemented");
        }
        else
        {
            errorMsg("Unrecognised source for domain " + this->name());
        }
    }
}

} // namespace accel
