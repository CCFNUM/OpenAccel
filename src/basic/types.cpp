// File : types.cpp
// Created : Wed Jan 03 2024 13:38:51 (+0100)
// Author : Mhamad Mahdi Alloush
// Description:
// Copyright (c) 2024 CCFNUM, Lucerne University of Applied Sciences and Arts.
// SPDX-License-Identifier: BSD-3-Clause

// code
#include "types.h"

namespace accel
{

// Field-type related

std::unordered_map<size_t, stk::io::FieldOutputType> fieldType{
    {1, stk::io::FieldOutputType::SCALAR},
    {2, stk::io::FieldOutputType::VECTOR_2D},
    {3, stk::io::FieldOutputType::VECTOR_3D},
    {4, stk::io::FieldOutputType::FULL_TENSOR_22},
    {9, stk::io::FieldOutputType::FULL_TENSOR_36}};

// Specialisations for pTraits

const std::string pTraits<label>::typeName = "label";
const std::string pTraits<scalar>::typeName = "scalar";

// equation identifiers
std::unordered_map<std::string, equationID> equationIDMap{
    {"segregated_navier_stokes", equationID::segregatedFlow},
    {"segregated_shear_stress_transport",
     equationID::segregatedShearStressTransport},
    {"segregated_transition_shear_stress_transport",
     equationID::segregatedTransitionShearStressTransport},
    {"segregated_correlation_transition_shear_stress_transport",
     equationID::segregatedCorrelationTransitionShearStressTransport},
    {"segregated_k_epsilon", equationID::segregatedKEpsilon},
    {"segregated_free_surface", equationID::segregatedFreeSurface},
    {"thermal_energy", equationID::thermalEnergy},
    {"total_energy", equationID::totalEnergy},
    {"displacement_diffusion", equationID::displacementDiffusion},
    {"volume_fraction", equationID::volumeFraction},
    {"solid_displacement", equationID::solidDisplacement},
    {"wallScaleDiffusion", equationID::wallScaleDiffusion},
};

equationID convertEquationIDFromString(std::string s)
{
    ::accel::tolower(s);
    std::replace(s.begin(), s.end(), '-', '_'); // replace hyphen with
    // underscore
    const auto it = equationIDMap.find(s); // check that `s` exists
    if (it == equationIDMap.end())
    {
        errorMsg("No equation ID found for `" + s + "`");
    }
    return equationIDMap[s];
}

// Heat transfer option
std::unordered_map<std::string, heatTransferOption> heatTransferOptionMap{
    {"none", heatTransferOption::none},
    {"isothermal", heatTransferOption::isothermal},
    {"thermal_energy", heatTransferOption::thermalEnergy},
    {"total_energy", heatTransferOption::totalEnergy}};

heatTransferOption convertHeatTransferOptionFromString(std::string s)
{
    ::accel::tolower(s);
    std::replace(s.begin(), s.end(), '-', '_'); // replace hyphen with
    // underscore
    const auto it = heatTransferOptionMap.find(s); // check that `s` exists
    if (it == heatTransferOptionMap.end())
    {
        errorMsg("No heat transfer option found for `" + s + "`");
    }
    return heatTransferOptionMap[s];
}

// Solid mechanics option
std::unordered_map<std::string, solidMechanicsOption> solidMechanicsOptionMap{
    {"none", solidMechanicsOption::none},
    {"linear_elastic", solidMechanicsOption::linearElastic},
    {"simplified_new_hookean", solidMechanicsOption::simplifiedNeoHookean}};

solidMechanicsOption convertSolidMechanicsOptionFromString(std::string s)
{
    ::accel::tolower(s);
    std::replace(s.begin(), s.end(), '-', '_'); // replace hyphen with
    // underscore
    const auto it = solidMechanicsOptionMap.find(s); // check that `s` exists
    if (it == solidMechanicsOptionMap.end())
    {
        errorMsg("No solid mechanics option found for `" + s + "`");
    }
    return solidMechanicsOptionMap[s];
}

// Turbulence option
std::unordered_map<std::string, turbulenceOption> turbulenceOptionMap{
    {"laminar", turbulenceOption::laminar},
    {"shear_stress_transport", turbulenceOption::shearStressTransport},
    {"k_epsilon", turbulenceOption::kEpsilon}};

turbulenceOption convertTurbulenceOptionFromString(std::string s)
{
    ::accel::tolower(s);
    std::replace(s.begin(), s.end(), '-', '_'); // replace hyphen with
    // underscore
    const auto it = turbulenceOptionMap.find(s); // check that `s` exists
    if (it == turbulenceOptionMap.end())
    {
        errorMsg("No turbulence option found for `" + s + "`");
    }
    return turbulenceOptionMap[s];
}

// Free surface model option
std::unordered_map<std::string, freeSurfaceModelOption>
    freeSurfaceModelOptionMap{{"none", freeSurfaceModelOption::none},
                              {"standard", freeSurfaceModelOption::standard}};

freeSurfaceModelOption convertFreeSurfaceModelOptionFromString(std::string s)
{
    ::accel::tolower(s);
    std::replace(s.begin(), s.end(), '-', '_'); // replace hyphen with
    // underscore
    const auto it = freeSurfaceModelOptionMap.find(s); // check that `s` exists
    if (it == freeSurfaceModelOptionMap.end())
    {
        errorMsg("No free surface model option found for `" + s + "`");
    }
    return freeSurfaceModelOptionMap[s];
}

// Surface tension model option
std::unordered_map<std::string, surfaceTensionModelOption>
    surfaceTensionModelOptionMap{
        {"none", surfaceTensionModelOption::none},
        {"continuum_surface_force",
         surfaceTensionModelOption::continuumSurfaceForce}};

surfaceTensionModelOption
convertSurfaceTensionModelOptionFromString(std::string s)
{
    ::accel::tolower(s);
    std::replace(s.begin(), s.end(), '-', '_');
    const auto it = surfaceTensionModelOptionMap.find(s);
    if (it == surfaceTensionModelOptionMap.end())
    {
        errorMsg("No surface tension model option found for `" + s + "`");
    }
    return surfaceTensionModelOptionMap[s];
}

// Source option
std::unordered_map<std::string, sourceOption> sourceOptionMap{
    {"source", sourceOption::source},
    {"total_source", sourceOption::totalSource}};

sourceOption convertSourceOptionFromString(std::string s)
{
    ::accel::tolower(s);
    std::replace(s.begin(), s.end(), '-', '_'); // replace hyphen with
    // underscore
    const auto it = sourceOptionMap.find(s); // check that `s` exists
    if (it == sourceOptionMap.end())
    {
        errorMsg("No source option found for `" + s + "`");
    }
    return sourceOptionMap[s];
}

// Turbulent flux closure for heat transfer option
std::unordered_map<std::string, turbulentFluxClosureForHeatTransferOption>
    turbulentFluxClosureForHeatTransferOptionMap{
        {"eddy_diffusivity",
         turbulentFluxClosureForHeatTransferOption::eddyDiffusivity}};

turbulentFluxClosureForHeatTransferOption
convertTurbulentFluxClosureForHeatTransferOptionFromString(std::string s)
{
    ::accel::tolower(s);
    std::replace(s.begin(), s.end(), '-', '_'); // replace hyphen with
    // underscore
    const auto it = turbulentFluxClosureForHeatTransferOptionMap.find(
        s); // check that `s` exists
    if (it == turbulentFluxClosureForHeatTransferOptionMap.end())
    {
        errorMsg(
            "No turbulent flux closure for heat transfer option found for `" +
            s + "`");
    }
    return turbulentFluxClosureForHeatTransferOptionMap[s];
}

// Buoyancy option
std::unordered_map<std::string, buoyancyOption> buoyancyOptionMap{
    {"non_buoyant", buoyancyOption::nonBuoyant},
    {"buoyant", buoyancyOption::buoyant}};

buoyancyOption convertBuoyancyOptionFromString(std::string s)
{
    ::accel::tolower(s);
    std::replace(s.begin(), s.end(), '-', '_'); // replace hyphen with
    // underscore
    const auto it = buoyancyOptionMap.find(s); // check that `s` exists
    if (it == buoyancyOptionMap.end())
    {
        errorMsg("No buoyancy option found for `" + s + "`");
    }
    return buoyancyOptionMap[s];
}

// Thermodynamic state option
std::unordered_map<std::string, thermodynamicStateOption>
    thermodynamicStateOptionMap{{"liquid", thermodynamicStateOption::liquid},
                                {"gas", thermodynamicStateOption::gas}};

thermodynamicStateOption
convertThermodynamicStateOptionFromString(std::string s)
{
    ::accel::tolower(s);
    std::replace(s.begin(), s.end(), '-', '_'); // replace hyphen with
    // underscore
    const auto it =
        thermodynamicStateOptionMap.find(s); // check that `s` exists
    if (it == thermodynamicStateOptionMap.end())
    {
        errorMsg("No thermodynamic state option found for `" + s + "`");
    }
    return thermodynamicStateOptionMap[s];
}

// Equation of state option
std::unordered_map<std::string, equationOfStateOption> equationOfStateOptionMap{
    {"value", equationOfStateOption::value},
    {"ideal_gas", equationOfStateOption::idealGas}};

equationOfStateOption convertEquationOfStateOptionFromString(std::string s)
{
    ::accel::tolower(s);
    std::replace(s.begin(), s.end(), '-', '_'); // replace hyphen with
    // underscore
    const auto it = equationOfStateOptionMap.find(s); // check that `s` exists
    if (it == equationOfStateOptionMap.end())
    {
        errorMsg("No equation of state option found for `" + s + "`");
    }
    return equationOfStateOptionMap[s];
}

// Specific heat option
std::unordered_map<std::string, specificHeatCapacityOption>
    specificHeatCapacityOptionMap{
        {"value", specificHeatCapacityOption::value},
        {"zero_pressure_polynomial",
         specificHeatCapacityOption::zeroPressurePolynomial}};

specificHeatCapacityOption
convertSpecificHeatCapacityOptionFromString(std::string s)
{
    ::accel::tolower(s);
    std::replace(s.begin(), s.end(), '-', '_'); // replace hyphen with
    // underscore
    const auto it =
        specificHeatCapacityOptionMap.find(s); // check that `s` exists
    if (it == specificHeatCapacityOptionMap.end())
    {
        errorMsg("No specific heat capacity option found for `" + s + "`");
    }
    return specificHeatCapacityOptionMap[s];
}

// Dynamic viscosity option
std::unordered_map<std::string, dynamicViscosityOption>
    dynamicViscosityOptionMap{
        {"value", dynamicViscosityOption::value},
        {"sutherlands_formula", dynamicViscosityOption::sutherlandsFormula}};

dynamicViscosityOption convertDynamicViscosityOptionFromString(std::string s)
{
    ::accel::tolower(s);
    std::replace(s.begin(), s.end(), '-', '_'); // replace hyphen with
    // underscore
    const auto it = dynamicViscosityOptionMap.find(s); // check that `s` exists
    if (it == dynamicViscosityOptionMap.end())
    {
        errorMsg("No dynamic viscosity option found for `" + s + "`");
    }
    return dynamicViscosityOptionMap[s];
}

// Thermal conductivity option
std::unordered_map<std::string, thermalConductivityOption>
    thermalConductivityOptionMap{
        {"value", thermalConductivityOption::value},
        {"kinetic_theory_model", thermalConductivityOption::kineticTheoryModel},
        {"sutherlands_formula", thermalConductivityOption::sutherlandsFormula}};

thermalConductivityOption
convertThermalConductivityOptionFromString(std::string s)
{
    ::accel::tolower(s);
    std::replace(s.begin(), s.end(), '-', '_'); // replace hyphen with
    // underscore
    const auto it =
        thermalConductivityOptionMap.find(s); // check that `s` exists
    if (it == thermalConductivityOptionMap.end())
    {
        errorMsg("No thermal conductivity option found for `" + s + "`");
    }
    return thermalConductivityOptionMap[s];
}

// Thermal expansivity option
std::unordered_map<std::string, thermalExpansivityOption>
    thermalExpansivityOptionMap{{"value", thermalExpansivityOption::value}};

thermalExpansivityOption
convertThermalExpansivityOptionFromString(std::string s)
{
    ::accel::tolower(s);
    std::replace(s.begin(), s.end(), '-', '_'); // replace hyphen with
    // underscore
    const auto it =
        thermalExpansivityOptionMap.find(s); // check that `s` exists
    if (it == thermalExpansivityOptionMap.end())
    {
        errorMsg("No thermal expansivity option found for `" + s + "`");
    }
    return thermalExpansivityOptionMap[s];
}

// Young modulus option
std::unordered_map<std::string, youngModulusOption> youngModulusOptionMap{
    {"value", youngModulusOption::value}};

youngModulusOption convertYoungModulusOptionFromString(std::string s)
{
    ::accel::tolower(s);
    std::replace(s.begin(), s.end(), '-', '_'); // replace hyphen with
    // underscore
    const auto it = youngModulusOptionMap.find(s); // check that `s` exists
    if (it == youngModulusOptionMap.end())
    {
        errorMsg("No Young modulus option found for `" + s + "`");
    }
    return youngModulusOptionMap[s];
}

// Poisson ratio option
std::unordered_map<std::string, poissonRatioOption> poissonRatioOptionMap{
    {"value", poissonRatioOption::value}};

poissonRatioOption convertPoissonRatioOptionFromString(std::string s)
{
    ::accel::tolower(s);
    std::replace(s.begin(), s.end(), '-', '_'); // replace hyphen with
    // underscore
    const auto it = poissonRatioOptionMap.find(s); // check that `s` exists
    if (it == poissonRatioOptionMap.end())
    {
        errorMsg("No Poisson ratio option found for `" + s + "`");
    }
    return poissonRatioOptionMap[s];
}

// Transient scheme

std::unordered_map<std::string, transientSchemeType> transientSchemeTypeMap{
    {"first_order_backward_euler",
     transientSchemeType::firstOrderBackwardEuler},
    {"second_order_backward_euler",
     transientSchemeType::secondOrderBackwardEuler}};

transientSchemeType convertTransientSchemeTypeFromString(std::string s)
{
    auto it = transientSchemeTypeMap.find(s);
    if (it != transientSchemeTypeMap.end())
    {
        return it->second;
    }
    else
    {
        errorMsg("Invalid transient scheme type: " + s);
    }

    return transientSchemeType::firstOrderBackwardEuler; // useless
}

std::array<scalar, 2> BDF1::coeff()
{
    // interpolation coefficients for BDF1
    return {1.0, -1.0};
}

std::array<scalar, 3> BDF2::coeff(const scalar dt_curr, const scalar dt_prev)
{
    // NOTE [faw 2025-10-24]: The reference implementation uses BDF1
    // coefficients for first step where dt_prev does not exist.
    //
    // generic interpolation coefficients for BDF2 (non-const timesteps)
    assert(std::abs(dt_prev) != 0);
    const scalar omega_n = dt_curr / dt_prev;
    constexpr scalar one = 1.0;
    constexpr scalar two = 2.0;
    return {(one + two * omega_n) / (one + omega_n),
            -(one + omega_n),
            omega_n * omega_n / (one + omega_n)};
}

std::unordered_map<std::string, timestepMode> timestepModeMap{
    {"constant", timestepMode::constant},
    {"periodic_interval", timestepMode::periodicInterval},
    {"specified_interval", timestepMode::specifiedInterval},
    {"adaptive", timestepMode::adaptive}};

timestepMode convertTimestepModeFromString(std::string s)
{
    auto it = timestepModeMap.find(s);
    if (it != timestepModeMap.end())
    {
        return it->second;
    }
    else
    {
        errorMsg("Invalid timestep mode: " + s);
    }

    return timestepMode::constant;
}

std::unordered_map<std::string, timestepAdaptationType>
    timestepAdaptationTypeMap{
        {"max_courant", timestepAdaptationType::maxCourant},
        {"rms_courant", timestepAdaptationType::rmsCourant}};

timestepAdaptationType convertTimestepAdaptationTypeFromString(std::string s)
{
    auto it = timestepAdaptationTypeMap.find(s);
    if (it != timestepAdaptationTypeMap.end())
    {
        return it->second;
    }
    else
    {
        errorMsg("Invalid timestep adaptation option: " + s);
    }

    return timestepAdaptationType::maxCourant;
}

std::unordered_map<std::string, outputFrequencyType> outputFrequencyTypeMap{
    {"time_interval", outputFrequencyType::timeInterval},
    {"timestep_interval", outputFrequencyType::timestepInterval}};

outputFrequencyType convertOutputFrequencyTypeFromString(std::string s)
{
    auto it = outputFrequencyTypeMap.find(s);
    if (it != outputFrequencyTypeMap.end())
    {
        return it->second;
    }
    else
    {
        errorMsg("Invalid output frequency option: " + s);
    }

    return outputFrequencyType::timestepInterval;
}

// Advection scheme

std::unordered_map<std::string, advectionSchemeType> advectionSchemeTypeMap{
    {"upwind", advectionSchemeType::upwind},
    {"high_resolution", advectionSchemeType::highResolution}};

advectionSchemeType convertAdvectionSchemeTypeFromString(std::string s)
{
    auto it = advectionSchemeTypeMap.find(s);
    if (it != advectionSchemeTypeMap.end())
    {
        return it->second;
    }
    else
    {
        errorMsg("Invalid advection scheme type: " + s);
    }

    return advectionSchemeType::upwind; // useless
}

// Physical boundary type

std::unordered_map<std::string, boundaryPhysicalType> physicalBoundaryTypeMap{
    {"inlet", boundaryPhysicalType::inlet},
    {"outlet", boundaryPhysicalType::outlet},
    {"opening", boundaryPhysicalType::opening},
    {"wall", boundaryPhysicalType::wall},
    {"symmetry", boundaryPhysicalType::symmetry}};

boundaryPhysicalType convertBoundaryPhysicalTypeFromString(std::string s)
{
    auto it = physicalBoundaryTypeMap.find(s);
    if (it != physicalBoundaryTypeMap.end())
    {
        return it->second;
    }
    else
    {
        errorMsg("Invalid physical boundary type: " + s);
    }

    return boundaryPhysicalType::wall; // useless
}

std::string toString(boundaryPhysicalType type)
{
    for (const auto& pair : physicalBoundaryTypeMap)
    {
        if (pair.second == type)
        {
            return pair.first;
        }
    }
    errorMsg("physical type not available");
    return "";
}

// Boundary condition type

std::unordered_map<std::string, boundaryConditionType> boundaryConditionTypeMap{
    {"specified_value", boundaryConditionType::specifiedValue},
    {"specified_flux", boundaryConditionType::specifiedFlux},
    {"zero_gradient", boundaryConditionType::zeroGradient},
    {"symmetry", boundaryConditionType::symmetry},
    {"mixed", boundaryConditionType::mixed},
    {"no_slip", boundaryConditionType::noSlip},
    {"slip", boundaryConditionType::slip},
    {"normal_speed", boundaryConditionType::normalSpeed},
    {"mass_flow_rate", boundaryConditionType::massFlowRate},
    {"static_pressure", boundaryConditionType::staticPressure},
    {"average_static_pressure", boundaryConditionType::averageStaticPressure},
    {"total_pressure", boundaryConditionType::totalPressure},
    {"static_temperature", boundaryConditionType::staticTemperature},
    {"total_temperature", boundaryConditionType::totalTemperature},
    {"specified_direction", boundaryConditionType::specifiedDirection},
    {"periodic_displacement", boundaryConditionType::periodicDisplacement},
    {"rigid_body_solution", boundaryConditionType::rigidBodySolution},
    {"intensity_and_length_scale",
     boundaryConditionType::intensityAndLengthScale},
    {"intensity_and_eddy_viscosity_ratio",
     boundaryConditionType::intensityAndEddyViscosityRatio}};

boundaryConditionType convertBoundaryConditionTypeFromString(std::string s)
{
    auto it = boundaryConditionTypeMap.find(s);
    if (it != boundaryConditionTypeMap.end())
    {
        return it->second;
    }
    else
    {
        errorMsg("Invalid boundary condition type: " + s);
    }

    return boundaryConditionType::zeroGradient; // useless
}

std::string toString(boundaryConditionType type)
{
    for (const auto& pair : boundaryConditionTypeMap)
    {
        if (pair.second == type)
        {
            return pair.first;
        }
    }
    errorMsg("type not available");
    return "";
}

// Data input type from yaml node

std::unordered_map<std::string, inputDataType> inputTypeMap{
    {"constant", inputDataType::constant},
    {"expression", inputDataType::expression},
    {"time_table", inputDataType::timeTable},
    {"profile_data", inputDataType::profileData}};

inputDataType convertInputTypeFromString(std::string s)
{
    auto it = inputTypeMap.find(s);
    if (it != inputTypeMap.end())
    {
        return it->second;
    }
    else
    {
        errorMsg("Invalid input data type: " + s);
    }

    return inputDataType::constant; // useless
}

std::string toString(inputDataType type)
{
    for (const auto& pair : inputTypeMap)
    {
        if (pair.second == type)
        {
            return pair.first;
        }
    }
    errorMsg("type not available");
    return "";
}

// Initial condition option

std::unordered_map<std::string, initialConditionOption>
    initialConditionOptionMap{{"null", initialConditionOption::null},
                              {"value", initialConditionOption::value}};

initialConditionOption convertInitialConditionOptionFromString(std::string s)
{
    auto it = initialConditionOptionMap.find(s);
    if (it != initialConditionOptionMap.end())
    {
        return it->second;
    }
    else
    {
        errorMsg("Invalid initial condition option: " + s);
    }

    return initialConditionOption::null; // useless
}

// Wall-function type

std::unordered_map<std::string, wallFunctionType> wallFunctionTypeMap{
    {"standard", wallFunctionType::standard},
    {"scalable", wallFunctionType::scalable},
    {"automatic", wallFunctionType::automatic}};

wallFunctionType convertWallFunctionTypeFromString(std::string s)
{
    auto it = wallFunctionTypeMap.find(s);
    if (it != wallFunctionTypeMap.end())
    {
        return it->second;
    }
    else
    {
        errorMsg("Invalid wall-function type type: " + s);
    }

    return wallFunctionType::standard; // useless
}

// Domain type

std::unordered_map<std::string, domainType> domainTypeMap{
    {"fluid", domainType::fluid},
    {"solid", domainType::solid}};

domainType convertDomainTypeFromString(std::string s)
{
    auto it = domainTypeMap.find(s);
    if (it != domainTypeMap.end())
    {
        return it->second;
    }
    else
    {
        errorMsg("Invalid domain type: " + s);
    }

    return domainType::fluid; // useless
}

// Interpolation type

std::unordered_map<std::string, interpolationSchemeType>
    interpolationSchemeTypeMap{
        {"trilinear", interpolationSchemeType::trilinear},
        {"linear_linear", interpolationSchemeType::linearLinear}};

std::unordered_map<std::string, timeInterpolationSchemeType>
    timeInterpolationSchemeTypeMap{
        {"closest", timeInterpolationSchemeType::closest},
        {"piecewise_linear", timeInterpolationSchemeType::piecewiseLinear},
        {"b_spline", timeInterpolationSchemeType::BSpline}};

interpolationSchemeType convertInterpolationSchemeTypeFromString(std::string s)
{
    auto it = interpolationSchemeTypeMap.find(s);
    if (it != interpolationSchemeTypeMap.end())
    {
        return it->second;
    }
    else
    {
        errorMsg("Invalid interpolation scheme type: " + s);
    }

    return interpolationSchemeType::trilinear; // useless
}

timeInterpolationSchemeType
convertTimeInterpolationSchemeTypeFromString(std::string s)
{
    auto it = timeInterpolationSchemeTypeMap.find(s);
    if (it != timeInterpolationSchemeTypeMap.end())
    {
        return it->second;
    }
    else
    {
        errorMsg("Invalid time interpolation scheme type: " + s);
    }

    return timeInterpolationSchemeType::closest; // useless
}

// Residual type

std::unordered_map<std::string, residualType> residualTypeMap{
    {"RMS", residualType::RMS}};

residualType convertResidualTypeFromString(std::string s)
{
    auto it = residualTypeMap.find(s);
    if (it != residualTypeMap.end())
    {
        return it->second;
    }
    else
    {
        errorMsg("Invalid residual type: " + s);
    }

    return residualType::RMS; // useless
}

// Pressure level information specification

std::unordered_map<std::string, pressureLevelInformationSpecification>
    pressureLevelInformationSpecificationMap{
        {"automatic", pressureLevelInformationSpecification::automatic},
        {"cartesian_coordinates",
         pressureLevelInformationSpecification::cartesianCoordinates}};

pressureLevelInformationSpecification
convertPressureLevelInformationSpecificationFromString(std::string s)
{
    auto it = pressureLevelInformationSpecificationMap.find(s);
    if (it != pressureLevelInformationSpecificationMap.end())
    {
        return it->second;
    }
    else
    {
        errorMsg("Invalid pressure level information specification: " + s);
    }

    return pressureLevelInformationSpecification::automatic; // useless
}

// Mesh motion type

std::unordered_map<std::string, meshMotionType> meshMotionTypeMap{
    {"stationary", meshMotionType::stationary},
    {"translating", meshMotionType::translating},
    {"rotating", meshMotionType::rotating}};

meshMotionType convertMeshMotionTypeFromString(std::string s)
{
    auto it = meshMotionTypeMap.find(s);
    if (it != meshMotionTypeMap.end())
    {
        return it->second;
    }
    else
    {
        errorMsg("Invalid mesh motion model: " + s);
    }

    return meshMotionType::stationary; // useless
}

// Mesh deformation option

std::unordered_map<std::string, meshDeformationSpecificationType>
    meshDeformationOptionMap{
        {"none", meshDeformationSpecificationType::none},
        {"regions_of_motion_specified",
         meshDeformationSpecificationType::regionsOfMotionSpecified},
        {"inherent", meshDeformationSpecificationType::inherent}};

meshDeformationSpecificationType
convertMeshDeformationSpecificationFromString(std::string s)
{
    auto it = meshDeformationOptionMap.find(s);
    if (it != meshDeformationOptionMap.end())
    {
        return it->second;
    }
    else
    {
        errorMsg("Invalid mesh deformation option: " + s);
    }

    return meshDeformationSpecificationType::none; // useless
}

// Mesh deformation model

std::unordered_map<std::string, meshDeformationModel>
    meshDeformationModelTypeMap{{"displacement_diffusion",
                                 meshDeformationModel::displacementDiffusion}};

meshDeformationModel convertMeshDeformationModelTypeFromString(std::string s)
{
    auto it = meshDeformationModelTypeMap.find(s);
    if (it != meshDeformationModelTypeMap.end())
    {
        return it->second;
    }
    else
    {
        errorMsg("Invalid mesh deformation model: " + s);
    }

    return meshDeformationModel::displacementDiffusion; // useless
}

// Mesh stiffness specification

std::unordered_map<std::string, meshStiffnessSpecificationType>
    meshStiffnessSpecificationTypeMap{
        {"value", meshStiffnessSpecificationType::value},
        {"increase_near_small_volumes",
         meshStiffnessSpecificationType::increaseNearSmallVolumes},
        {"increase_near_boundaries",
         meshStiffnessSpecificationType::increaseNearBoundaries},
        {"blended_distance_and_small_volumes",
         meshStiffnessSpecificationType::blendedDistanceAndSmallVolumes}};

meshStiffnessSpecificationType
convertMeshStiffnessSpecificationTypeFromString(std::string s)
{
    auto it = meshStiffnessSpecificationTypeMap.find(s);
    if (it != meshStiffnessSpecificationTypeMap.end())
    {
        return it->second;
    }
    else
    {
        errorMsg("Invalid mesh stiffness specification: " + s);
    }

    return meshStiffnessSpecificationType::value; // useless
}

// Frame type

std::unordered_map<std::string, boundaryRelativeFrameType> frameTypeMap{
    {"stationary", boundaryRelativeFrameType::absolute},
    {"translating", boundaryRelativeFrameType::relative},
    {"rotating", boundaryRelativeFrameType::relative}};

boundaryRelativeFrameType
convertBoundaryRelativeFrameTypeFromString(std::string s)
{
    auto it = frameTypeMap.find(s);
    if (it != frameTypeMap.end())
    {
        return it->second;
    }
    else
    {
        errorMsg("Invalid relative frame type: " + s);
    }

    return boundaryRelativeFrameType::absolute; // useless
}

// Kinematic formulation type

std::unordered_map<std::string, kinematicFormulationType>
    kinematicFormulationTypeMap{
        {"total_lagrangian", kinematicFormulationType::totalLagrangian},
        {"updated_lagrangian", kinematicFormulationType::updatedLagrangian}};

kinematicFormulationType
convertKinematicFormulationTypeFromString(std::string s)
{
    auto it = kinematicFormulationTypeMap.find(s);
    if (it != kinematicFormulationTypeMap.end())
    {
        return it->second;
    }
    else
    {
        errorMsg("Invalid kinematic formulation type: " + s);
    }

    return kinematicFormulationType::totalLagrangian; // useless
}

// Post process type

std::unordered_map<std::string, postProcessType> postProcessTypeMap{
    {"reduction", postProcessType::reduction},
    {"force", postProcessType::force},
    {"probe", postProcessType::probe}};

postProcessType convertPostProcessTypeFromString(std::string s)
{
    auto it = postProcessTypeMap.find(s);
    if (it != postProcessTypeMap.end())
    {
        return it->second;
    }
    else
    {
        errorMsg("Invalid post-process type: " + s);
    }

    return postProcessType::reduction; // useless
}

// Monitor type

std::unordered_map<std::string, reductionType> statisticsTypeMap{
    {"sum", reductionType::sum},
    {"average", reductionType::average},
    {"area_average", reductionType::areaAverage}};

reductionType convertStatisticsTypeFromString(std::string s)
{
    auto it = statisticsTypeMap.find(s);
    if (it != statisticsTypeMap.end())
    {
        return it->second;
    }
    else
    {
        errorMsg("Invalid statistics type: " + s);
    }

    return reductionType::sum; // useless
}

} // namespace accel
