# File       : CMakeLists.txt
# Created    : Thu Oct 16 2025 08:17:36 (+0200)
# Author     : Fabian Wermelinger
# Description: CMake build definitions
# Copyright (c) 2025 CCFNUM, Lucerne University of Applied Sciences and Arts.
# SPDX-License-Identifier: BSD-3-Clause

set(DIR ${CMAKE_CURRENT_SOURCE_DIR})

# project sources
target_sources(main_obj PRIVATE
    ${DIR}/flow/segregatedFlow/navierStokesAssembler/navierStokesAssembler.cpp
    ${DIR}/flow/segregatedFlow/navierStokesAssembler/navierStokesAssemblerElemBoundaryConditions.cpp
    ${DIR}/flow/segregatedFlow/navierStokesAssembler/navierStokesAssemblerElemTerms.cpp
    ${DIR}/flow/segregatedFlow/navierStokesAssembler/navierStokesAssemblerNodeTerms.cpp
    ${DIR}/flow/segregatedFlow/pressureCorrectionAssembler/pressureCorrectionAssembler.cpp
    ${DIR}/flow/segregatedFlow/pressureCorrectionAssembler/pressureCorrectionAssemblerElemBoundaryConditions.cpp
    ${DIR}/flow/segregatedFlow/pressureCorrectionAssembler/pressureCorrectionAssemblerElemTerms.cpp
    ${DIR}/flow/segregatedFlow/pressureCorrectionAssembler/pressureCorrectionAssemblerNodeTerms.cpp
    ############################################################################
    ${DIR}/heatTransfer/thermalTemperature/thermalTemperatureAssemblerElemBoundaryConditions.cpp
    ${DIR}/heatTransfer/thermalTemperature/thermalTemperatureAssemblerElemTerms.cpp
    ${DIR}/heatTransfer/thermalTemperature/thermalTemperatureAssemblerNodeTerms.cpp
    ${DIR}/heatTransfer/thermalEnergy/thermalEnergyAssemblerElemBoundaryConditions.cpp
    ${DIR}/heatTransfer/thermalEnergy/thermalEnergyAssemblerNodeTerms.cpp
    ${DIR}/heatTransfer/totalEnergy/totalEnergyAssemblerElemBoundaryConditions.cpp
    ${DIR}/heatTransfer/totalEnergy/totalEnergyAssemblerElemTerms.cpp
    ${DIR}/heatTransfer/totalEnergy/totalEnergyAssemblerNodeTerms.cpp
    ############################################################################
    ${DIR}/multiphase/freeSurfaceFlow/segregatedBulkFlow/bulkPressureCorrection/bulkPressureCorrectionAssembler.cpp
    ${DIR}/multiphase/freeSurfaceFlow/segregatedBulkFlow/bulkPressureCorrection/bulkPressureCorrectionAssemblerElemBoundaryConditions.cpp
    ${DIR}/multiphase/freeSurfaceFlow/segregatedBulkFlow/bulkPressureCorrection/bulkPressureCorrectionAssemblerElemTerms.cpp
    ${DIR}/multiphase/freeSurfaceFlow/segregatedBulkFlow/bulkPressureCorrection/bulkPressureCorrectionAssemblerNodeTerms.cpp
    ${DIR}/multiphase/freeSurfaceFlow/volumeFraction/volumeFractionAssembler.cpp
    ${DIR}/multiphase/freeSurfaceFlow/volumeFraction/volumeFractionAssemblerElemTerms.cpp
    ############################################################################
    ${DIR}/turbulence/RANS/kEpsilon/segregatedKEpsilon/turbulentDissipationRate/turbulentDissipationRateAssemblerElemBoundaryConditions.cpp
    ${DIR}/turbulence/RANS/kEpsilon/segregatedKEpsilon/turbulentDissipationRate/turbulentDissipationRateAssemblerNodeTerms.cpp
    ${DIR}/turbulence/RANS/kEpsilon/segregatedKEpsilon/turbulentKineticEnergy/turbulentKineticEnergyKEpsilonAssemblerElemBoundaryConditions.cpp
    ${DIR}/turbulence/RANS/kEpsilon/segregatedKEpsilon/turbulentKineticEnergy/turbulentKineticEnergyKEpsilonAssemblerNodeTerms.cpp
    ${DIR}/turbulence/RANS/shearStressTransport/segregatedShearStressTransport/turbulentEddyFrequency/turbulentEddyFrequencySSTAssemblerElemBoundaryConditions.cpp
    ${DIR}/turbulence/RANS/shearStressTransport/segregatedShearStressTransport/turbulentEddyFrequency/turbulentEddyFrequencySSTAssemblerNodeTerms.cpp
    ${DIR}/turbulence/RANS/shearStressTransport/segregatedShearStressTransport/turbulentKineticEnergy/turbulentKineticEnergySSTAssemblerElemBoundaryConditions.cpp
    ${DIR}/turbulence/RANS/shearStressTransport/segregatedShearStressTransport/turbulentKineticEnergy/turbulentKineticEnergySSTAssemblerNodeTerms.cpp
    ${DIR}/turbulence/RANS/shearStressTransport/segregatedTransitionShearStressTransport/transitionOnsetReynoldsNumberTransitionSST/transitionOnsetReynoldsNumberTransitionSSTAssemblerNodeTerms.cpp
    ${DIR}/turbulence/RANS/shearStressTransport/segregatedTransitionShearStressTransport/turbulentIntermittencyTransitionSST/turbulentIntermittencyTransitionSSTAssemblerNodeTerms.cpp
    ${DIR}/turbulence/RANS/shearStressTransport/segregatedTransitionShearStressTransport/turbulentKineticEnergyTransitionSST/turbulentKineticEnergyTransitionSSTAssemblerElemBoundaryConditions.cpp
    ${DIR}/turbulence/RANS/shearStressTransport/segregatedTransitionShearStressTransport/turbulentKineticEnergyTransitionSST/turbulentKineticEnergyTransitionSSTAssemblerNodeTerms.cpp
    ${DIR}/turbulence/RANS/shearStressTransport/segregatedCorrelationTransitionShearStressTransport/turbulentIntermittencyCorrelationTransitionSST/turbulentIntermittencyCorrelationTransitionSSTAssemblerNodeTerms.cpp
    ${DIR}/turbulence/RANS/shearStressTransport/segregatedCorrelationTransitionShearStressTransport/turbulentKineticEnergyCorrelationTransitionSST/turbulentKineticEnergyCorrelationTransitionSSTAssemblerNodeTerms.cpp
    ${DIR}/turbulence/RANS/shearStressTransport/segregatedCorrelationTransitionShearStressTransport/turbulentKineticEnergyCorrelationTransitionSST/turbulentKineticEnergyCorrelationTransitionSSTAssemblerElemBoundaryConditions.cpp
    ############################################################################
    ${DIR}/solidMechanics/solidDisplacement/solidDisplacementAssembler.cpp
    ${DIR}/solidMechanics/solidDisplacement/solidDisplacementAssemblerElemBoundaryConditions.cpp
    ${DIR}/solidMechanics/solidDisplacement/solidDisplacementAssemblerElemTerms.cpp
    ${DIR}/solidMechanics/solidDisplacement/solidDisplacementAssemblerNodeTerms.cpp
)

# define header directories (required since #include statements do not quote
# full path and location of headers is inside src)
target_include_directories(main_obj PUBLIC
    ${DIR}
    ############################################################################
    ${DIR}/flow/segregatedFlow/navierStokesAssembler
    ${DIR}/flow/segregatedFlow/pressureCorrectionAssembler
    ############################################################################
    ${DIR}/heatTransfer/thermalTemperature
    ${DIR}/heatTransfer/thermalEnergy
    ${DIR}/heatTransfer/totalEnergy
    ############################################################################
    ${DIR}/solidMechanics/solidDisplacement
    ############################################################################
    ${DIR}/multiphase/freeSurfaceFlow/segregatedBulkFlow/bulkPressureCorrection
    ${DIR}/multiphase/freeSurfaceFlow/volumeFraction
    ############################################################################
    ${DIR}/phiAssembler
    ############################################################################
    ${DIR}/turbulence/RANS/kEpsilon/segregatedKEpsilon/turbulentDissipationRate
    ${DIR}/turbulence/RANS/kEpsilon/segregatedKEpsilon/turbulentKineticEnergy
    ${DIR}/turbulence/RANS/shearStressTransport/segregatedShearStressTransport/turbulentEddyFrequency
    ${DIR}/turbulence/RANS/shearStressTransport/segregatedShearStressTransport/turbulentKineticEnergy
    ${DIR}/turbulence/RANS/shearStressTransport/segregatedTransitionShearStressTransport/transitionOnsetReynoldsNumberTransitionSST
    ${DIR}/turbulence/RANS/shearStressTransport/segregatedTransitionShearStressTransport/turbulentIntermittencyTransitionSST
    ${DIR}/turbulence/RANS/shearStressTransport/segregatedTransitionShearStressTransport/turbulentKineticEnergyTransitionSST
    ${DIR}/turbulence/RANS/shearStressTransport/segregatedCorrelationTransitionShearStressTransport/turbulentIntermittencyCorrelationTransitionSST
    ${DIR}/turbulence/RANS/shearStressTransport/segregatedCorrelationTransitionShearStressTransport/turbulentKineticEnergyCorrelationTransitionSST
    ${DIR}/turbulence/RANS/turbulentEddyFrequency
    ${DIR}/turbulence/RANS/turbulentKineticEnergy
)
