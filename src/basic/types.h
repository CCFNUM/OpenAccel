// File       : types.h
// Created    : Wed Jan 03 2024 13:38:51 (+0100)
// Author     : Mhamad Mahdi Alloush
// Description: Fundamental type aliases, enumerations, and constants for the
// solver
// Copyright (c) 2024 CCFNUM, Lucerne University of Applied Sciences and
// Arts.
// SPDX-License-Identifier: BSD-3-Clause

#ifndef TYPES_H
#define TYPES_H

// stk_util
#include <stk_util/environment/WallTime.hpp>
#include <stk_util/environment/perf_util.hpp>
#include <stk_util/parallel/Parallel.hpp>

// stk_mesh/base/fem
#include "stk_search/Box.hpp" // for Box
#include "stk_search/FilterCoarseSearch.hpp"
#include "stk_search/IdentProc.hpp" // for IdentProc
#include "stk_search/Point.hpp"     // for Point
#include "stk_search/Sphere.hpp"
#include <stk_io/FillMesh.hpp>
#include <stk_mesh/base/BulkData.hpp>
#include <stk_mesh/base/Comm.hpp>
#include <stk_mesh/base/CoordinateSystems.hpp>
#include <stk_mesh/base/CreateEdges.hpp>
#include <stk_mesh/base/CreateFaces.hpp>
#include <stk_mesh/base/Entity.hpp>
#include <stk_mesh/base/Field.hpp>
#include <stk_mesh/base/FieldBLAS.hpp>
#include <stk_mesh/base/FieldParallel.hpp>
#include <stk_mesh/base/GetEntities.hpp>
#include <stk_mesh/base/Ghosting.hpp>
#include <stk_mesh/base/MeshBuilder.hpp>
#include <stk_mesh/base/MetaData.hpp>
#include <stk_mesh/base/Relation.hpp>
#include <stk_mesh/base/SkinBoundary.hpp>
#include <stk_search/BoundingBox.hpp>
#include <stk_search/CoarseSearch.hpp>
#include <stk_search/SearchMethod.hpp>
#include <stk_transfer/GeometricTransfer.hpp>

// stk_io
#include <stk_io/InputFile.hpp>
#include <stk_io/IossBridge.hpp>
#include <stk_io/StkMeshIoBroker.hpp>

// stk option parsing
#include <stk_util/environment/OptionsSpecification.hpp>
#include <stk_util/environment/ParseCommandLineArgs.hpp>
#include <stk_util/environment/ParsedOptions.hpp>
#include <stk_util/util/ParameterList.hpp>

// stk_util
#include <stk_util/parallel/ParallelReduce.hpp>

// basic c++
#include <algorithm> // std::sort, std::stable_sort
#include <array>
#include <cassert>
#include <cmath>
#include <ctime>
#include <filesystem>
#include <fstream>
#include <functional>
#include <iomanip>
#include <iostream>
#include <limits>
#include <map>
#include <memory>
#include <mpi.h>
#include <numeric> // std::iota
#include <set>
#include <span>
#include <sstream>
#include <string>
#include <unordered_map>
#include <utility>
#include <vector>

// External
#include "master_element/MasterElement.h"
#include "master_element/MasterElementFactory.h"
#include "yaml-cpp/yaml.h"
#include <Eigen/Core>
#include <Eigen/Geometry>

// code
#include "macros.h"

namespace fs = std::filesystem;

namespace accel
{

// Bring external master element types into accel namespace
using sierra::nalu::MasterElement;
using sierra::nalu::MasterElementRepo;

typedef double scalar;
typedef int label;
typedef uint64_t ulabel;
using Vector = std::vector<scalar>;

// using Eigen::RowMajor for conceptual compatibility,
// default would be Eigen::ColumMajor
using rvector2 = Eigen::Matrix<scalar, 2, 1>;
using rvector3 = Eigen::Matrix<scalar, 3, 1>;
using rvectorD = Eigen::Matrix<scalar, SPATIAL_DIM, 1>;

using ivector2 = Eigen::Matrix<label, 2, 1>;
using ivector3 = Eigen::Matrix<label, 3, 1>;
using ivectorD = Eigen::Matrix<label, SPATIAL_DIM, 1>;

using rtensor2 = Eigen::Matrix<scalar, 2, 2, Eigen::RowMajor>;
using rtensor3 = Eigen::Matrix<scalar, 3, 3, Eigen::RowMajor>;
using rtensorD =
    Eigen::Matrix<scalar, SPATIAL_DIM, SPATIAL_DIM, Eigen::RowMajor>;

using itensor2 = Eigen::Matrix<label, 2, 2, Eigen::RowMajor>;
using itensor3 = Eigen::Matrix<label, 3, 3, Eigen::RowMajor>;
using itensorD =
    Eigen::Matrix<label, SPATIAL_DIM, SPATIAL_DIM, Eigen::RowMajor>;

// View typedefs to use eigen functionality from pointer-type arrays
// avoids stack allocation and data copies
// usage:
// rvector3ViewC A(someStdVectorA.data());
// rvector3ViewC B(someStdVectorB.data());
// scalar dotProd = A.dot(B);
using rvector2View = Eigen::Map<rvector2>;
using rvector2ViewC = Eigen::Map<const rvector2>;
using rvector3View = Eigen::Map<rvector3>;
using rvector3ViewC = Eigen::Map<const rvector3>;
using rvectorDView = Eigen::Map<rvectorD>;
using rvectorDViewC = Eigen::Map<const rvectorD>;

using ivector2View = Eigen::Map<ivector2>;
using ivector2ViewC = Eigen::Map<const ivector2>;
using ivector3View = Eigen::Map<ivector3>;
using ivector3ViewC = Eigen::Map<const ivector3>;
using ivectorDView = Eigen::Map<ivectorD>;
using ivectorDViewC = Eigen::Map<const ivectorD>;

using rvector2View = Eigen::Map<rvector2>;
using rvector2ViewC = Eigen::Map<const rvector2>;
using rvector3View = Eigen::Map<rvector3>;
using rvector3ViewC = Eigen::Map<const rvector3>;
using rvectorDView = Eigen::Map<rvectorD>;
using rvectorDViewC = Eigen::Map<const rvectorD>;

using itensor2View = Eigen::Map<itensor2>;
using itensor2ViewC = Eigen::Map<const itensor2>;
using itensor3View = Eigen::Map<itensor3>;
using itensor3ViewC = Eigen::Map<const itensor3>;
using itensorDView = Eigen::Map<itensorD>;
using itensorDViewC = Eigen::Map<const itensorD>;

using rtensor2View = Eigen::Map<rtensor2>;
using rtensor2ViewC = Eigen::Map<const rtensor2>;
using rtensor3View = Eigen::Map<rtensor3>;
using rtensor3ViewC = Eigen::Map<const rtensor3>;
using rtensorDView = Eigen::Map<rtensorD>;
using rtensorDViewC = Eigen::Map<const rtensorD>;

using itensor2View = Eigen::Map<itensor2>;
using itensor2ViewC = Eigen::Map<const itensor2>;
using itensor3View = Eigen::Map<itensor3>;
using itensor3ViewC = Eigen::Map<const itensor3>;
using itensorDView = Eigen::Map<itensorD>;
using itensorDViewC = Eigen::Map<const itensorD>;

constexpr scalar BIG =
    static_cast<scalar>(1.0) / std::numeric_limits<scalar>::epsilon();
constexpr scalar VBIG = std::numeric_limits<scalar>::max();
constexpr scalar SMALL = std::numeric_limits<scalar>::epsilon();
constexpr scalar VSMALL = std::numeric_limits<scalar>::min();

// Field type map
extern std::unordered_map<size_t, stk::io::FieldOutputType> fieldType;

// typedefs for STK types
using STKScalarField = stk::mesh::Field<scalar>;

// Template class acting as a zero
template <typename T>
class zero
{
public:
    zero() = default;

    // Operator overloading to mimic zero behavior
    operator T() const
    {
        return static_cast<T>(0);
    }
};

// Templates class acting as properties (traits) for types
template <class Base>
class pTraits : public Base
{
public:
    // Copy construct from base class
    explicit pTraits(const Base& obj) : Base(obj)
    {
    }
};

// Template specialisation for pTraits<label>
template <>
class pTraits<label>
{
public:
    static constexpr int nComponents = 1;
    static const std::string typeName;
};

// Template specialisation for pTraits<scalar>
template <>
class pTraits<scalar>
{
public:
    static constexpr int nComponents = 1;
    static const std::string typeName;
};

// equation identifier used in equation factories as well as linking domains to
// equations
enum class equationID
{
    // fluid equations
    segregatedFlow = 0,
    segregatedShearStressTransport,
    segregatedTransitionShearStressTransport,
    segregatedCorrelationTransitionShearStressTransport,
    segregatedKEpsilon,
    segregatedFreeSurface,
    // solid equations
    //
    // other equations
    thermalEnergy,
    totalEnergy,
    displacementDiffusion,
    volumeFraction,
    solidDisplacement,
    wallScaleDiffusion,
    // default for equations that do not require an ID
    noID,
    // total number of equation IDs
    numberOfDeclaredEquations
};

equationID convertEquationIDFromString(std::string s);

// Heat Transfer option
enum class heatTransferOption
{
    none,
    isothermal,
    thermalEnergy,
    totalEnergy
};

heatTransferOption convertHeatTransferOptionFromString(std::string s);

// Solid mechanics option
enum class solidMechanicsOption
{
    none,
    linearElastic,
    simplifiedNeoHookean
};

solidMechanicsOption convertSolidMechanicsOptionFromString(std::string s);

// Turbulence model option
enum class turbulenceOption
{
    laminar,
    shearStressTransport,
    kEpsilon
};

turbulenceOption convertTurbulenceOptionFromString(std::string s);

// Free surface model option
enum class freeSurfaceModelOption
{
    none,
    standard
};

freeSurfaceModelOption convertFreeSurfaceModelOptionFromString(std::string s);

// Surface tension model option
enum class surfaceTensionModelOption
{
    none,
    continuumSurfaceForce
};

surfaceTensionModelOption
convertSurfaceTensionModelOptionFromString(std::string s);

// Source option
enum class sourceOption
{
    source,
    totalSource
};

sourceOption convertSourceOptionFromString(std::string s);

// Turbulent flux closure for heat transfer option
enum class turbulentFluxClosureForHeatTransferOption
{
    eddyDiffusivity
};

turbulentFluxClosureForHeatTransferOption
convertTurbulentFluxClosureForHeatTransferOptionFromString(std::string s);

// Buoyancy option
enum class buoyancyOption
{
    nonBuoyant,
    buoyant
};

buoyancyOption convertBuoyancyOptionFromString(std::string s);

// Buoyancy model
enum class buoyancyModel
{
    full,
    boussinesq
};

// Thermodynamic state option
enum class thermodynamicStateOption
{
    liquid,
    gas
};

thermodynamicStateOption
convertThermodynamicStateOptionFromString(std::string s);

// Equation of state option
enum class equationOfStateOption
{
    value,
    idealGas
};

equationOfStateOption convertEquationOfStateOptionFromString(std::string s);

// Specific heat option
enum class specificHeatCapacityOption
{
    value,
    zeroPressurePolynomial
};

specificHeatCapacityOption
convertSpecificHeatCapacityOptionFromString(std::string s);

// Dynamic viscosity option
enum class dynamicViscosityOption
{
    value,
    sutherlandsFormula
};

dynamicViscosityOption convertDynamicViscosityOptionFromString(std::string s);

// Thermal conductivity option
enum class thermalConductivityOption
{
    value,
    kineticTheoryModel,
    sutherlandsFormula
};

thermalConductivityOption
convertThermalConductivityOptionFromString(std::string s);

// Thermal expansivity option
enum class thermalExpansivityOption
{
    value
};

thermalExpansivityOption
convertThermalExpansivityOptionFromString(std::string s);

// Young modulus option
enum class youngModulusOption
{
    value
};

youngModulusOption convertYoungModulusOptionFromString(std::string s);

// Poisson ratio option
enum class poissonRatioOption
{
    value
};

poissonRatioOption convertPoissonRatioOptionFromString(std::string s);

// Term type
enum class term
{
    transient,
    falseTransient,
    advection,
    diffusion,
    source,
    stress,
    pressureGradient,
    rhieChow
};

// Transient schemes
enum class transientSchemeType
{
    firstOrderBackwardEuler, // BDF1
    secondOrderBackwardEuler // BDF2
};

transientSchemeType convertTransientSchemeTypeFromString(std::string s);

struct BDF1
{
    static std::array<scalar, 2> coeff();
};

struct BDF2
{
    static std::array<scalar, 3> coeff(const scalar dt_curr,
                                       const scalar dt_prev);
};

enum class timestepMode
{
    constant,
    periodicInterval,
    specifiedInterval,
    adaptive
};

timestepMode convertTimestepModeFromString(std::string s);

enum class timestepAdaptationType
{
    maxCourant,
    rmsCourant
};

timestepAdaptationType convertTimestepAdaptationTypeFromString(std::string s);

enum class outputFrequencyType
{
    timeInterval,
    timestepInterval
};

outputFrequencyType convertOutputFrequencyTypeFromString(std::string s);

// Advection schemes
enum class advectionSchemeType
{
    upwind,
    highResolution
};

advectionSchemeType convertAdvectionSchemeTypeFromString(std::string s);

// Physical bounadaries
enum class boundaryPhysicalType
{
    inlet,
    outlet,
    opening,
    wall,
    symmetry
};

boundaryPhysicalType convertBoundaryPhysicalTypeFromString(std::string s);

std::string toString(boundaryPhysicalType type);

// Boundary condition types
enum class boundaryConditionType
{
    specifiedValue,                // generic
    specifiedFlux,                 // generic
    zeroGradient,                  // generic
    symmetry,                      // generic
    mixed,                         // generic
    noSlip,                        // specific to velocity
    slip,                          // specific to velocity
    normalSpeed,                   // specific to velocity
    massFlowRate,                  // specific to velocity
    specifiedDirection,            // specific to velocity
    staticPressure,                // specific to pressure
    averageStaticPressure,         // specific to pressure
    totalPressure,                 // specific to pressure
    staticTemperature,             // specific to temperature
    totalTemperature,              // specific to temperature
    periodicDisplacement,          // specific to displacement
    rigidBodySolution,             // specific to displacement
    intensityAndLengthScale,       // specific to k, omega or epsilon
    intensityAndEddyViscosityRatio // specific to k, omega or epsilon
};

boundaryConditionType convertBoundaryConditionTypeFromString(std::string s);

std::string toString(boundaryConditionType type);

// Boundary condition input type
enum class inputDataType
{
    null,
    constant,
    expression,
    timeTable,
    profileData
};

inputDataType convertInputTypeFromString(std::string s);

std::string toString(inputDataType type);

// Initialization type
enum class initialConditionOption
{
    null,
    value
};

initialConditionOption convertInitialConditionOptionFromString(std::string s);

#ifdef HAS_INTERFACE
// Interface model option
enum class interfaceModelOption
{
    translationalPeriodicity,
    rotationalPeriodicity,
    generalConnection
};

interfaceModelOption convertInterfaceModelOptionFromString(std::string s);

// Interface type
enum class interfaceType
{
    fluid_fluid,
    solid_solid,
    fluid_solid,
};

interfaceType convertInterfaceTypeFromString(std::string s);

std::string toString(interfaceType type);
#endif /* HAS_INTERFACE */

// Wall-function type
enum class wallFunctionType
{
    standard,
    scalable,
    automatic
};

wallFunctionType convertWallFunctionTypeFromString(std::string s);

// Zone type
enum class domainType
{
    undefined,
    fluid,
    solid
};

domainType convertDomainTypeFromString(std::string s);

// Interpolation scheme
enum class interpolationSchemeType
{
    trilinear,
    linearLinear,
};

interpolationSchemeType convertInterpolationSchemeTypeFromString(std::string s);

enum class timeInterpolationSchemeType
{
    closest,
    piecewiseLinear,
    BSpline,
};

timeInterpolationSchemeType
convertTimeInterpolationSchemeTypeFromString(std::string s);

// Residual type
enum class residualType
{
    RMS
};

residualType convertResidualTypeFromString(std::string s);

// Pressure level information specification
enum class pressureLevelInformationSpecification
{
    automatic,
    cartesianCoordinates
};

pressureLevelInformationSpecification
convertPressureLevelInformationSpecificationFromString(std::string s);

// Mesh motion type
enum class meshMotionType
{
    stationary,
    translating,
    rotating
};

meshMotionType convertMeshMotionTypeFromString(std::string s);

// Mesh deformation option
enum class meshDeformationSpecificationType
{
    none,
    regionsOfMotionSpecified,
    inherent
};

meshDeformationSpecificationType
convertMeshDeformationSpecificationFromString(std::string s);

// Mesh deformation model
enum class meshDeformationModel
{
    displacementDiffusion
};

meshDeformationModel convertMeshDeformationModelTypeFromString(std::string s);

// Mesh stiffness specification
enum class meshStiffnessSpecificationType
{
    value,
    increaseNearSmallVolumes,
    increaseNearBoundaries,
    blendedDistanceAndSmallVolumes
};

meshStiffnessSpecificationType
convertMeshStiffnessSpecificationTypeFromString(std::string s);

// Boundary relative frame type
enum class boundaryRelativeFrameType
{
    absolute,
    relative
};

boundaryRelativeFrameType
convertBoundaryRelativeFrameTypeFromString(std::string s);

// Kinematic formulation type
enum class kinematicFormulationType
{
    totalLagrangian,
    updatedLagrangian
};

kinematicFormulationType
convertKinematicFormulationTypeFromString(std::string s);

// Post process type
enum class postProcessType
{
    reduction,
    force,
    probe
};

postProcessType convertPostProcessTypeFromString(std::string s);

// Monitor type
enum class reductionType
{
    sum,
    average,
    areaAverage
};

reductionType convertStatisticsTypeFromString(std::string s);

} // namespace accel

#endif // TYPES_H
