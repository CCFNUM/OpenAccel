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

# SFEM external kernels (linked unconditionally: SFEM and CVFEM solid
# mechanics are both always compiled in, with runtime dispatch between them
# -- see solid_assembler_type expert parameter).
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

# Unlike neohookean_ogden's element_api.hpp (header-only templates,
# instantiated at the call site -- no .cpp needed), the generated
# modified_mooney_rivlin *_isoparametric_mesh_soa functions are plain
# extern "C" functions declared in sfem_GeneratedModifiedMooneyRivlin_
# c_abi.hpp but defined in separate .cpp files (dispatch + per-topology
# kernels) under frontend/ops/generated/modified_mooney_rivlin/. Those
# must be compiled into this library or the symbols stay undefined at
# link time.
#
# Listed explicitly rather than globbed: a plain glob over that
# directory also pulls in sfem_GeneratedModifiedMooneyRivlin.cpp and
# sfem_GeneratedModifiedMooneyRivlin_registration.cpp, which are the
# sfem::Op-framework wrapper/registration for this technology and
# transitively #include sfem_Op.hpp (via sfem_GeneratedModifiedMooney
# Rivlin.hpp) -- a framework OpenAccel does not build against. Only the
# dispatch entry points OpenAccel actually calls
# (modified_mooney_rivlin_gradient_2d_isoparametric_mesh_soa,
# modified_mooney_rivlin_hessian_bsr_2d_isoparametric_mesh_soa, both
# defined in the dispatch .cpp below) plus every per-topology kernel
# their switch statements dispatch to are listed.
#
# This means both the 2D (QUAD4/TRI3, TRI6, PROTEUS_QUAD4) *and* 3D
# (HEX8/HEX27/TET4/TET10/PROTEUS_HEX8/27/64) topology kernels are
# included, even though OpenAccel's SPATIAL_DIM=2 build only ever
# calls the 2D dispatcher: the dispatch .cpp is one translation unit
# that defines both the 2D and 3D dispatchers together, and this
# build has no -ffunction-sections/--gc-sections, so once the linker
# pulls that object file in for the 2D symbols, the 3D dispatcher's
# own unresolved references (to the 3D per-topology symbols) must
# also resolve from somewhere -- hence they're listed here too.
# QUAD4/HEX8/HEX27 forward to their PROTEUS_* counterparts; TRI3,
# TRI6, TET4, TET10, and the PROTEUS_* kernels are self-contained.
# None of the twelve files below, nor anything they #include,
# references sfem_Op.hpp.
set(SFEM_MODIFIED_MOONEY_RIVLIN_SOURCES
    ${CMAKE_CURRENT_SOURCE_DIR}/sfem/frontend/ops/generated/modified_mooney_rivlin/op/sfem_GeneratedModifiedMooneyRivlin_isoparametric_dispatch.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/sfem/frontend/ops/generated/modified_mooney_rivlin/d2/proteus_quad4/modified_mooney_rivlin_proteus_quad4_operator.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/sfem/frontend/ops/generated/modified_mooney_rivlin/d2/quad4/modified_mooney_rivlin_quad4_operator.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/sfem/frontend/ops/generated/modified_mooney_rivlin/d2/tri3/modified_mooney_rivlin_tri3_operator.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/sfem/frontend/ops/generated/modified_mooney_rivlin/d2/tri6/modified_mooney_rivlin_tri6_operator.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/sfem/frontend/ops/generated/modified_mooney_rivlin/d3/hex27/modified_mooney_rivlin_hex27_operator.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/sfem/frontend/ops/generated/modified_mooney_rivlin/d3/hex8/modified_mooney_rivlin_hex8_operator.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/sfem/frontend/ops/generated/modified_mooney_rivlin/d3/proteus_hex27/modified_mooney_rivlin_proteus_hex27_operator.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/sfem/frontend/ops/generated/modified_mooney_rivlin/d3/proteus_hex64/modified_mooney_rivlin_proteus_hex64_operator.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/sfem/frontend/ops/generated/modified_mooney_rivlin/d3/proteus_hex8/modified_mooney_rivlin_proteus_hex8_operator.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/sfem/frontend/ops/generated/modified_mooney_rivlin/d3/tet10/modified_mooney_rivlin_tet10_operator.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/sfem/frontend/ops/generated/modified_mooney_rivlin/d3/tet4/modified_mooney_rivlin_tet4_operator.cpp
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
    ${SFEM_MODIFIED_MOONEY_RIVLIN_SOURCES}
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
    ${CMAKE_CURRENT_SOURCE_DIR}/sfem/frontend/ops/generated
    ${CMAKE_CURRENT_SOURCE_DIR}/sfem/frontend/ops/generated/neohookean_ogden/op
    ${CMAKE_CURRENT_SOURCE_DIR}/sfem/frontend/ops/generated/modified_mooney_rivlin/op
    $<TARGET_PROPERTY:smesh,INTERFACE_INCLUDE_DIRECTORIES>
)
target_link_libraries(openaccel_sfem_elasticity PUBLIC MPI::MPI_C)
target_link_libraries(main_obj PUBLIC openaccel_sfem_elasticity)

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
