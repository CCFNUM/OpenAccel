# File       : CMakeLists.txt
# Created    : Thu Oct 16 2025 08:17:36 (+0200)
# Author     : Fabian Wermelinger
# Description: CMake build definitions
# Copyright (c) 2025 CCFNUM, Lucerne University of Applied Sciences and Arts.
# SPDX-License-Identifier: BSD-3-Clause

set(DIR ${CMAKE_CURRENT_SOURCE_DIR})

# project sources
target_sources(main_obj PRIVATE
    ${DIR}/elementField/scalar/elementScalarField.cpp
    ${DIR}/elementField/scalar/heatFlowRate/heatFlowRate.cpp
    ${DIR}/elementField/scalar/massFlowRate/massFlowRate.cpp
    ############################################################################
    ${DIR}/nodeField/scalar/nodeScalarField.cpp
    ${DIR}/nodeField/scalar/property/flow/density.cpp
    ${DIR}/nodeField/scalar/property/flow/dynamicViscosity.cpp
    ${DIR}/nodeField/scalar/property/heatTransfer/compressibility.cpp
    ${DIR}/nodeField/scalar/property/heatTransfer/specificHeatCapacity.cpp
    ${DIR}/nodeField/scalar/property/heatTransfer/thermalConductivity.cpp
    ${DIR}/nodeField/scalar/property/heatTransfer/thermalExpansivity.cpp
    ${DIR}/nodeField/scalar/property/property.cpp
    ${DIR}/nodeField/scalar/property/turbulence/turbulentThermalConductivity.cpp
    ${DIR}/nodeField/scalar/property/turbulence/turbulentViscosity.cpp
    ${DIR}/nodeField/scalar/transport/flow/pressure.cpp
    ${DIR}/nodeField/scalar/transport/heatTransfer/specificEnthalpy.cpp
    ${DIR}/nodeField/scalar/transport/heatTransfer/specificTotalEnthalpy.cpp
    ${DIR}/nodeField/scalar/transport/heatTransfer/temperature.cpp
    ${DIR}/nodeField/scalar/transport/multiphase/volumeFraction.cpp
    ${DIR}/nodeField/scalar/transport/turbulence/transitionOnsetReynoldsNumber.cpp
    ${DIR}/nodeField/scalar/transport/turbulence/turbulentDissipationRate.cpp
    ${DIR}/nodeField/scalar/transport/turbulence/turbulentEddyFrequency.cpp
    ${DIR}/nodeField/scalar/transport/turbulence/turbulentIntermittency.cpp
    ${DIR}/nodeField/scalar/transport/turbulence/turbulentKineticEnergy.cpp
    ${DIR}/nodeField/scalar/property/solidMechanics/youngModulus.cpp
    ${DIR}/nodeField/scalar/property/solidMechanics/poissonRatio.cpp
    ${DIR}/nodeField/tensor/nodeTensorField.cpp
    ${DIR}/nodeField/vector/nodeVectorField.cpp
    ${DIR}/nodeField/vector/transport/flow/velocity.cpp
    ${DIR}/nodeField/vector/addons/displacement.cpp
    ${DIR}/nodeSideField/nodeSideField.cpp
    ############################################################################
    ${DIR}/sideField/sideField.cpp
)

# define header directories (required since #include statements do not quote
# full path and location of headers is inside src)
target_include_directories(main_obj PUBLIC
    ${DIR}
    ############################################################################
    ${DIR}/elementField
    ${DIR}/elementField/scalar
    ${DIR}/elementField/scalar/heatFlowRate
    ${DIR}/elementField/scalar/massFlowRate
    ############################################################################
    ${DIR}/nodeField
    ${DIR}/nodeField/scalar
    ${DIR}/nodeField/scalar/property
    ${DIR}/nodeField/scalar/property/flow
    ${DIR}/nodeField/scalar/property/heatTransfer
    ${DIR}/nodeField/scalar/property/turbulence
    ${DIR}/nodeField/scalar/property/solidMechanics
    ${DIR}/nodeField/scalar/transport/flow
    ${DIR}/nodeField/scalar/transport/heatTransfer
    ${DIR}/nodeField/scalar/transport/multiphase
    ${DIR}/nodeField/scalar/transport/turbulence
    ${DIR}/nodeField/tensor
    ${DIR}/nodeField/vector
    ${DIR}/nodeField/vector/transport/flow
    ${DIR}/nodeField/vector/addons
    ${DIR}/nodeSideField
    ############################################################################
    ${DIR}/sideField
)
