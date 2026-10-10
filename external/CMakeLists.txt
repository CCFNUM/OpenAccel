# File       : CMakeLists.txt
# Created    : Thu Oct 16 2025 08:17:36 (+0200)
# Author     : Fabian Wermelinger
# Description: CMake build definitions
# Copyright 2025 CCFNUM HSLU T&A. All Rights Reserved.

target_include_directories(main_obj SYSTEM PUBLIC
    ${CMAKE_CURRENT_SOURCE_DIR}/core/eigen
    ${CMAKE_CURRENT_SOURCE_DIR}/core/exprtk
    ${CMAKE_CURRENT_SOURCE_DIR}/core/nanoflann/include
    ${CMAKE_CURRENT_SOURCE_DIR}/core/nalu/include
)

target_include_directories(main_obj SYSTEM PUBLIC
    ${CMAKE_CURRENT_SOURCE_DIR}/core/gplotpp
)

# SFEM solid mechanics backend (solid_assembler_type: sfem), WITH_SFEM only.
# Uses an installed SFEM when ACCEL_USE_SYSTEM_SFEM=ON (point SFEM_DIR or
# CMAKE_PREFIX_PATH at it), otherwise builds the external/core/sfem submodule.
# MatrixIO is provided by smesh (SFEM_ENABLE_MPI=ON).
if(WITH_SFEM)
    option(ACCEL_USE_SYSTEM_SFEM "Use an installed SFEM instead of external/sfem" OFF)
    if(ACCEL_USE_SYSTEM_SFEM)
        find_package(SFEM CONFIG REQUIRED)
    else()
        set(SFEM_ENABLE_SUBMODULE ON CACHE BOOL "SFEM is built as a subproject")
        set(SFEM_ENABLE_MINIMAL ON CACHE BOOL "Build only the SFEM library (no drivers, no tests)")
        set(SFEM_ENABLE_MPI ON CACHE BOOL "Enable MPI support")
        set(SFEM_ENABLE_AMG OFF CACHE BOOL "Enable AMG solver")
        set(SFEM_ENABLE_RYAML OFF CACHE BOOL "Enable YAML input files with RapidYAML")
        set(SFEM_ENABLE_TRACE OFF CACHE BOOL "Enable trace facilities")
        set(SFEM_ENABLE_MEM_DIAGNOSTICS OFF CACHE BOOL "Enable mem diagnostics")
        set(SFEM_LOG_LEVEL 0 CACHE STRING "Set the default log level for SFEM")
        # LibLinSolve turns BUILD_SHARED_LIBS on globally; SFEM is built static
        set(BUILD_SHARED_LIBS OFF)
        add_subdirectory(core/sfem EXCLUDE_FROM_ALL)
        add_library(SFEM::sfem ALIAS sfem)
    endif()

    target_link_libraries(main_obj PUBLIC SFEM::sfem)
endif()

# Nalu master_element sources
target_sources(main_obj PRIVATE
    ${CMAKE_CURRENT_SOURCE_DIR}/core/nalu/src/master_element/Hex27CVFEM.C
    ${CMAKE_CURRENT_SOURCE_DIR}/core/nalu/src/master_element/Hex8CVFEM.C
    ${CMAKE_CURRENT_SOURCE_DIR}/core/nalu/src/master_element/Hex8FEM.C
    ${CMAKE_CURRENT_SOURCE_DIR}/core/nalu/src/master_element/MasterElement.C
    ${CMAKE_CURRENT_SOURCE_DIR}/core/nalu/src/master_element/MasterElementFactory.C
    ${CMAKE_CURRENT_SOURCE_DIR}/core/nalu/src/master_element/Pyr5CVFEM.C
    ${CMAKE_CURRENT_SOURCE_DIR}/core/nalu/src/master_element/Quad42DCVFEM.C
    ${CMAKE_CURRENT_SOURCE_DIR}/core/nalu/src/master_element/Quad43DCVFEM.C
    ${CMAKE_CURRENT_SOURCE_DIR}/core/nalu/src/master_element/Quad92DCVFEM.C
    ${CMAKE_CURRENT_SOURCE_DIR}/core/nalu/src/master_element/Tet10FEM.C
    ${CMAKE_CURRENT_SOURCE_DIR}/core/nalu/src/master_element/Tet4CVFEM.C
    ${CMAKE_CURRENT_SOURCE_DIR}/core/nalu/src/master_element/Tri32DCVFEM.C
    ${CMAKE_CURRENT_SOURCE_DIR}/core/nalu/src/master_element/Tri6FEM.C
    ${CMAKE_CURRENT_SOURCE_DIR}/core/nalu/src/master_element/Wed6CVFEM.C
    ${CMAKE_CURRENT_SOURCE_DIR}/core/nalu/src/master_element/MasterElementWork.F
    ${CMAKE_CURRENT_SOURCE_DIR}/core/nalu/src/master_element/CalculateBentHexVolume.F
    ${CMAKE_CURRENT_SOURCE_DIR}/core/nalu/src/master_element/CalculateHexVolume.F
    ${CMAKE_CURRENT_SOURCE_DIR}/core/nalu/src/master_element/CalculateOctohedronVolume.F
    ${CMAKE_CURRENT_SOURCE_DIR}/core/nalu/src/master_element/CalculatePolyhedralArea.F
    ${CMAKE_CURRENT_SOURCE_DIR}/core/nalu/src/master_element/CalculatePolyhedralVolume.F
    ${CMAKE_CURRENT_SOURCE_DIR}/core/nalu/src/master_element/CalculateQuadArea.F
)

# Nalu tabular_props sources
target_sources(main_obj PRIVATE
    ${CMAKE_CURRENT_SOURCE_DIR}/core/nalu/src/tabular_props/BSpline.C
    ${CMAKE_CURRENT_SOURCE_DIR}/core/nalu/src/tabular_props/Converter.C
    ${CMAKE_CURRENT_SOURCE_DIR}/core/nalu/src/tabular_props/Functions.C
    ${CMAKE_CURRENT_SOURCE_DIR}/core/nalu/src/tabular_props/H5IO.C
    ${CMAKE_CURRENT_SOURCE_DIR}/core/nalu/src/tabular_props/HDF5FilePtr.C
    ${CMAKE_CURRENT_SOURCE_DIR}/core/nalu/src/tabular_props/HDF5Table.C
    ${CMAKE_CURRENT_SOURCE_DIR}/core/nalu/src/tabular_props/LU.C
)

# Nalu others
target_sources(main_obj PRIVATE
    ${CMAKE_CURRENT_SOURCE_DIR}/core/nalu/src/NaluEnv.C
)
