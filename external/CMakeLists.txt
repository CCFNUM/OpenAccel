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

if(NOT WITH_CVFEM_SOLID_MECHANICS)
    set(SMESH_ENABLE_DEMO OFF CACHE BOOL "" FORCE)
    set(SMESH_ENABLE_DEV_MODE OFF CACHE BOOL "" FORCE)
    set(SMESH_ENABLE_MPI OFF CACHE BOOL "" FORCE)
    set(SMESH_ENABLE_TESTING OFF CACHE BOOL "" FORCE)
    set(SMESH_ENABLE_BENCHMARK OFF CACHE BOOL "" FORCE)
    set(SMESH_ENABLE_RYAML OFF CACHE BOOL "" FORCE)
    add_subdirectory(sfem/external/smesh EXCLUDE_FROM_ALL)

    set(SFEM_LOG_LEVEL 0)
    set(SFEM_TET4_CUDA 0)
    set(SFEM_TET10_CUDA 0)
    set(SFEM_TET10_WENO 1)
    set(SFEM_CUDA_MEMORY_MODEL 0)
    configure_file(
        ${CMAKE_CURRENT_SOURCE_DIR}/sfem/configuration/sfem_config.h.in
        ${CMAKE_CURRENT_BINARY_DIR}/sfem_config.h
    )

    add_library(openaccel_sfem_elasticity STATIC
        ${CMAKE_CURRENT_SOURCE_DIR}/sfem/algebra/sortreduce.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/sfem/base/sfem_base.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/sfem/operators/linear_elasticity.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/sfem/operators/tri3/tri3_linear_elasticity.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/sfem/operators/tet4/tet4_linear_elasticity.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/sfem/operators/tet10/tet10_linear_elasticity.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/sfem/operators/macro_tet4/macro_tet4_linear_elasticity.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/sfem/operators/hex8/hex8_linear_elasticity.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/sfem/operators/sshex8/sshex8_linear_elasticity.cpp
    )
    target_include_directories(openaccel_sfem_elasticity PUBLIC
        ${CMAKE_CURRENT_BINARY_DIR}
        ${CMAKE_CURRENT_SOURCE_DIR}/sfem/base
        ${CMAKE_CURRENT_SOURCE_DIR}/sfem/algebra
        ${CMAKE_CURRENT_SOURCE_DIR}/sfem/operators
        ${CMAKE_CURRENT_SOURCE_DIR}/sfem/operators/beam2
        ${CMAKE_CURRENT_SOURCE_DIR}/sfem/operators/hex8
        ${CMAKE_CURRENT_SOURCE_DIR}/sfem/operators/macro_tet4
        ${CMAKE_CURRENT_SOURCE_DIR}/sfem/operators/sshex8
        ${CMAKE_CURRENT_SOURCE_DIR}/sfem/operators/tet10
        ${CMAKE_CURRENT_SOURCE_DIR}/sfem/operators/tet4
        ${CMAKE_CURRENT_SOURCE_DIR}/sfem/operators/tri3
        $<TARGET_PROPERTY:smesh,INTERFACE_INCLUDE_DIRECTORIES>
    )
    target_link_libraries(openaccel_sfem_elasticity PUBLIC MPI::MPI_C)
    target_link_libraries(main_obj PUBLIC openaccel_sfem_elasticity)
endif()

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
