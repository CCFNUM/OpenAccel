# File       : CMakeLists.txt
# Created    : Thu Oct 16 2025 08:17:36 (+0200)
# Author     : Fabian Wermelinger
# Description: CMake build definitions
# Copyright (c) 2025 CCFNUM, Lucerne University of Applied Sciences and Arts.
# SPDX-License-Identifier: BSD-3-Clause

set(DIR ${CMAKE_CURRENT_SOURCE_DIR})

# project sources
target_sources(main_obj PRIVATE
    ${DIR}/model.cpp
    ############################################################################
    ${DIR}/constitutive/thermo/idealGasModel.cpp
    ${DIR}/constitutive/thermo/thermoModel.cpp
    ${DIR}/constitutive/transport/kineticTheoryModel.cpp
    ${DIR}/constitutive/transport/sutherlandsFormulaModel.cpp
    ${DIR}/constitutive/transport/transportModel.cpp
    ############################################################################
    ${DIR}/flow/flowModel.cpp
    ############################################################################
    ${DIR}/heatTransfer/heatTransferModel.cpp
    ############################################################################
    ${DIR}/solidMechanics/solidMechanicsModel.cpp
    ############################################################################
    ${DIR}/multiphase/freeSurfaceFlow/freeSurfaceFlowModel.cpp
    ${DIR}/multiphase/multiphaseModel.cpp
    ############################################################################
    ${DIR}/turbulence/RANS/RANSModel.cpp
    ${DIR}/turbulence/RANS/kEpsilon/kEpsilonModel.cpp
    ${DIR}/turbulence/RANS/kEpsilon/kEpsilonModelWallFunctions.cpp
    ${DIR}/turbulence/RANS/shearStressTransport/shearStressTransportModel.cpp
    ${DIR}/turbulence/RANS/shearStressTransport/shearStressTransportModelWallFunctions.cpp
    ${DIR}/turbulence/RANS/shearStressTransport/transitionShearStressTransport/transitionShearStressTransportModel.cpp
    ${DIR}/turbulence/RANS/shearStressTransport/correlationTransitionShearStressTransport/correlationTransitionShearStressTransportModel.cpp
    ${DIR}/turbulence/turbulenceModel.cpp
    ${DIR}/turbulence/turbulenceModelWallFunctions.cpp
)

# define header directories (required since #include statements do not quote
# full path and location of headers is inside src)
target_include_directories(main_obj PUBLIC
    ${DIR}
    ############################################################################
    ${DIR}/constitutive/thermo
    ${DIR}/constitutive/transport
    ############################################################################
    ${DIR}/flow
    ############################################################################
    ${DIR}/heatTransfer
    ############################################################################
    ${DIR}/solidMechanics
    ############################################################################
    ${DIR}/multiphase
    ${DIR}/multiphase/freeSurfaceFlow
    ############################################################################
    ${DIR}/turbulence
    ${DIR}/turbulence/RANS
    ${DIR}/turbulence/RANS/kEpsilon
    ${DIR}/turbulence/RANS/shearStressTransport
    ${DIR}/turbulence/RANS/shearStressTransport/transitionShearStressTransport
    ${DIR}/turbulence/RANS/shearStressTransport/correlationTransitionShearStressTransport
)
