# File       : CMakeLists.txt
# Created    : Thu Oct 16 2025 08:17:36 (+0200)
# Author     : Fabian Wermelinger
# Description: CMake build definitions
# Copyright (c) 2025 CCFNUM, Lucerne University of Applied Sciences and Arts.
# SPDX-License-Identifier: BSD-3-Clause

set(DIR ${CMAKE_CURRENT_SOURCE_DIR})

# project sources
target_sources(main_obj PRIVATE
    ${DIR}/equation.cpp
    ${DIR}/auxiliary/convergenceAcceleration.cpp
    ############################################################################
    ${DIR}/flow/navierStokes/navierStokesEquation.cpp
    ${DIR}/flow/pressureCorrection/pressureCorrectionEquation.cpp
    ${DIR}/flow/segregatedFlowEquations.cpp
    ############################################################################
    ${DIR}/heatTransfer/thermalTemperature/thermalTemperatureEquation.cpp
    ${DIR}/heatTransfer/thermalEnergy/thermalEnergyEquation.cpp
    ${DIR}/heatTransfer/totalEnergy/totalEnergyEquation.cpp
    ############################################################################
    ${DIR}/multiphase/freeSurfaceFlow/bulkNavierStokes/bulkNavierStokesEquation.cpp
    ${DIR}/multiphase/freeSurfaceFlow/bulkPressureCorrection/bulkPressureCorrectionEquation.cpp
    ${DIR}/multiphase/freeSurfaceFlow/segregatedFreeSurfaceFlowEquations.cpp
    ${DIR}/multiphase/freeSurfaceFlow/volumeFraction/volumeFractionEquation.cpp
    ############################################################################
    ${DIR}/turbulence/RANS/kEpsilon/segregatedKEpsilonEquations.cpp
    ${DIR}/turbulence/RANS/kEpsilon/turbulentDissipationRate/turbulentDissipationRateEquation.cpp
    ${DIR}/turbulence/RANS/kEpsilon/turbulentKineticEnergy/turbulentKineticEnergyKEpsilonEquation.cpp
    ${DIR}/turbulence/RANS/shearStressTransport/segregatedShearStressTransportEquations.cpp
    ${DIR}/turbulence/RANS/shearStressTransport/transitionShearStressTransport/segregatedTransitionShearStressTransportEquations.cpp
    ${DIR}/turbulence/RANS/shearStressTransport/transitionShearStressTransport/transitionOnsetReynoldsNumber/transitionOnsetReynoldsNumberTransitionSSTEquation.cpp
    ${DIR}/turbulence/RANS/shearStressTransport/transitionShearStressTransport/turbulentIntermittency/turbulentIntermittencyTransitionSSTEquation.cpp
    ${DIR}/turbulence/RANS/shearStressTransport/transitionShearStressTransport/turbulentKineticEnergy/turbulentKineticEnergyTransitionSSTEquation.cpp
    ${DIR}/turbulence/RANS/shearStressTransport/correlationTransitionShearStressTransport/segregatedCorrelationTransitionShearStressTransportEquations.cpp
    ${DIR}/turbulence/RANS/shearStressTransport/correlationTransitionShearStressTransport/turbulentIntermittency/turbulentIntermittencyCorrelationTransitionSSTEquation.cpp
    ${DIR}/turbulence/RANS/shearStressTransport/correlationTransitionShearStressTransport/turbulentKineticEnergy/turbulentKineticEnergyCorrelationTransitionSSTEquation.cpp
    ${DIR}/turbulence/RANS/shearStressTransport/turbulentEddyFrequency/turbulentEddyFrequencySSTEquation.cpp
    ${DIR}/turbulence/RANS/shearStressTransport/turbulentKineticEnergy/turbulentKineticEnergySSTEquation.cpp
    ${DIR}/turbulence/turbulentEddyFrequency/turbulentEddyFrequencyEquation.cpp
    ${DIR}/turbulence/turbulentKineticEnergy/turbulentKineticEnergyEquation.cpp
    ############################################################################
    ${DIR}/solidMechanics/solidDisplacement/solidDisplacementEquation.cpp
    ##############################################################
)

# define header directories (required since #include statements do not quote
# full path and location of headers is inside src)
target_include_directories(main_obj PUBLIC
    ${DIR}
    ${DIR}/auxiliary
    ############################################################################
    ${DIR}/flow
    ${DIR}/flow/navierStokes
    ${DIR}/flow/pressureCorrection
    ############################################################################
    ${DIR}/heatTransfer/thermalTemperature
    ${DIR}/heatTransfer/thermalEnergy
    ${DIR}/heatTransfer/totalEnergy
    ############################################################################
    ${DIR}/solidMechanics/solidDisplacement
    ############################################################################
    ${DIR}/multiphase/freeSurfaceFlow
    ${DIR}/multiphase/freeSurfaceFlow/bulkNavierStokes
    ${DIR}/multiphase/freeSurfaceFlow/bulkPressureCorrection
    ${DIR}/multiphase/freeSurfaceFlow/volumeFraction
    ############################################################################
    ${DIR}/turbulence/RANS/kEpsilon
    ${DIR}/turbulence/RANS/kEpsilon/turbulentDissipationRate
    ${DIR}/turbulence/RANS/kEpsilon/turbulentKineticEnergy
    ${DIR}/turbulence/RANS/shearStressTransport
    ${DIR}/turbulence/RANS/shearStressTransport/transitionShearStressTransport
    ${DIR}/turbulence/RANS/shearStressTransport/transitionShearStressTransport/transitionOnsetReynoldsNumber
    ${DIR}/turbulence/RANS/shearStressTransport/transitionShearStressTransport/turbulentIntermittency
    ${DIR}/turbulence/RANS/shearStressTransport/transitionShearStressTransport/turbulentKineticEnergy
    ${DIR}/turbulence/RANS/shearStressTransport/correlationTransitionShearStressTransport
    ${DIR}/turbulence/RANS/shearStressTransport/correlationTransitionShearStressTransport/turbulentIntermittency
    ${DIR}/turbulence/RANS/shearStressTransport/correlationTransitionShearStressTransport/turbulentKineticEnergy
    ${DIR}/turbulence/RANS/shearStressTransport/turbulentEddyFrequency
    ${DIR}/turbulence/RANS/shearStressTransport/turbulentKineticEnergy
    ${DIR}/turbulence/turbulentEddyFrequency
    ${DIR}/turbulence/turbulentKineticEnergy
)
