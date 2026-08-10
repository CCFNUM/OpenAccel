# File       : CMakeLists.txt
# Created    : Thu Oct 16 2025 08:17:36 (+0200)
# Author     : Fabian Wermelinger
# Description: CMake build definitions
# Copyright 2025 CCFNUM HSLU T&A. All Rights Reserved.

target_include_directories(main_obj SYSTEM PUBLIC
    ${CMAKE_CURRENT_SOURCE_DIR}/eigen
    ${CMAKE_CURRENT_SOURCE_DIR}/exprtk
    ${CMAKE_CURRENT_SOURCE_DIR}/nanoflann/include
    ${CMAKE_CURRENT_SOURCE_DIR}/nalu/include
)

target_include_directories(main_obj SYSTEM PUBLIC
    ${CMAKE_CURRENT_SOURCE_DIR}/gplotpp
)

# Nalu master_element sources
target_sources(main_obj PRIVATE
    ${CMAKE_CURRENT_SOURCE_DIR}/nalu/src/master_element/Hex27CVFEM.C
    ${CMAKE_CURRENT_SOURCE_DIR}/nalu/src/master_element/Hex8CVFEM.C
    ${CMAKE_CURRENT_SOURCE_DIR}/nalu/src/master_element/Hex8FEM.C
    ${CMAKE_CURRENT_SOURCE_DIR}/nalu/src/master_element/MasterElement.C
    ${CMAKE_CURRENT_SOURCE_DIR}/nalu/src/master_element/MasterElementFactory.C
    ${CMAKE_CURRENT_SOURCE_DIR}/nalu/src/master_element/Pyr5CVFEM.C
    ${CMAKE_CURRENT_SOURCE_DIR}/nalu/src/master_element/Quad42DCVFEM.C
    ${CMAKE_CURRENT_SOURCE_DIR}/nalu/src/master_element/Quad43DCVFEM.C
    ${CMAKE_CURRENT_SOURCE_DIR}/nalu/src/master_element/Quad92DCVFEM.C
    ${CMAKE_CURRENT_SOURCE_DIR}/nalu/src/master_element/Tet10FEM.C
    ${CMAKE_CURRENT_SOURCE_DIR}/nalu/src/master_element/Tet4CVFEM.C
    ${CMAKE_CURRENT_SOURCE_DIR}/nalu/src/master_element/Tri32DCVFEM.C
    ${CMAKE_CURRENT_SOURCE_DIR}/nalu/src/master_element/Tri6FEM.C
    ${CMAKE_CURRENT_SOURCE_DIR}/nalu/src/master_element/Wed6CVFEM.C
    ${CMAKE_CURRENT_SOURCE_DIR}/nalu/src/master_element/MasterElementWork.F
    ${CMAKE_CURRENT_SOURCE_DIR}/nalu/src/master_element/CalculateBentHexVolume.F
    ${CMAKE_CURRENT_SOURCE_DIR}/nalu/src/master_element/CalculateHexVolume.F
    ${CMAKE_CURRENT_SOURCE_DIR}/nalu/src/master_element/CalculateOctohedronVolume.F
    ${CMAKE_CURRENT_SOURCE_DIR}/nalu/src/master_element/CalculatePolyhedralArea.F
    ${CMAKE_CURRENT_SOURCE_DIR}/nalu/src/master_element/CalculatePolyhedralVolume.F
    ${CMAKE_CURRENT_SOURCE_DIR}/nalu/src/master_element/CalculateQuadArea.F
)

# Nalu tabular_props sources
target_sources(main_obj PRIVATE
    ${CMAKE_CURRENT_SOURCE_DIR}/nalu/src/tabular_props/BSpline.C
    ${CMAKE_CURRENT_SOURCE_DIR}/nalu/src/tabular_props/Converter.C
    ${CMAKE_CURRENT_SOURCE_DIR}/nalu/src/tabular_props/Functions.C
    ${CMAKE_CURRENT_SOURCE_DIR}/nalu/src/tabular_props/H5IO.C
    ${CMAKE_CURRENT_SOURCE_DIR}/nalu/src/tabular_props/HDF5FilePtr.C
    ${CMAKE_CURRENT_SOURCE_DIR}/nalu/src/tabular_props/HDF5Table.C
    ${CMAKE_CURRENT_SOURCE_DIR}/nalu/src/tabular_props/LU.C
)

# Nalu others
target_sources(main_obj PRIVATE
    ${CMAKE_CURRENT_SOURCE_DIR}/nalu/src/NaluEnv.C
)

# SSDF (signed distance function library used for wall distance)
if(NOT EXISTS ${CMAKE_CURRENT_SOURCE_DIR}/ssdf/CMakeLists.txt)
    message(FATAL_ERROR "SSDF submodule is missing, run: git submodule update --init --recursive")
endif()

# SSDF defaults for OpenAccel (no FORCE, still overridable with -D on the command line)
set(SSDF_ENABLE_MPI ON CACHE BOOL "Enable MPI support in SSDF")
set(SSDF_ENABLE_OPENMP OFF CACHE BOOL "Enable OpenMP support in SSDF")
set(SSDF_ENABLE_CUBIQL OFF CACHE BOOL "Enable cuBQL (CUDA BVH) support in SSDF")

# SSDF would otherwise turn this ON for the whole project
set(BUILD_SHARED_LIBS OFF CACHE BOOL "Build libraries as shared libraries")

# EXCLUDE_FROM_ALL drops the SSDF drivers and tests from the build, libssdf
# itself is still built as a link dependency of main_obj. Header search paths
# (serial/ and the generated ssdf_config.hpp) come with the target.
add_subdirectory(ssdf EXCLUDE_FROM_ALL)
target_link_libraries(main_obj PUBLIC ssdf)
