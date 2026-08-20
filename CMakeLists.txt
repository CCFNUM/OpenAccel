cmake_minimum_required(VERSION 3.10)
cmake_policy(SET CMP0079 NEW)

project(OpenAccel)

# Set C++ standard (optional but recommended)
set(CMAKE_CXX_STANDARD 20)

# options
option(WITH_AMGSOLVER "Build with AMG solver" OFF)
option(WITH_CANONICAL_SUFFIX "Enable canonical executable suffix" OFF)
option(WITH_CVFEM_SOLID_MECHANICS "Use CVFEM for solid mechanics" ON)
option(WITH_DEBUG_LIB "Build debugging utilities library" OFF)
option(WITH_FASTMATH "Use fast math optimization" OFF)
option(WITH_THERMAL_TEMPERATURE "Enable thermal temperature energy equation" OFF)
option(WITH_TRILINOS_SOLVER "Enable Trilinos solver support" ON)
option(WITH_WARNINGS "Enable compiler warnings" OFF)
option(WITH_VECTORIZATION "Enable x86_64 AVX2/FMA vectorization flags (-mavx2 -mfma)" OFF)

# Get Git hash
include(cmake/CheckGit.cmake)
CheckGitSetup()

# Find MPI package
option(MPI_CXX_SKIP_MPICXX "Skip MPI C++ bindings" OFF)
if(DEFINED PETSC_DIR)
    get_filename_component(PETSC_PARENT_DIR "${PETSC_DIR}" DIRECTORY)
    file(GLOB SPACK_OPENMPI_PREFIXES "${PETSC_PARENT_DIR}/openmpi-*")
    list(LENGTH SPACK_OPENMPI_PREFIXES SPACK_OPENMPI_PREFIX_COUNT)
    if(SPACK_OPENMPI_PREFIX_COUNT EQUAL 1)
        list(GET SPACK_OPENMPI_PREFIXES 0 SPACK_OPENMPI_PREFIX)
        if(EXISTS "${SPACK_OPENMPI_PREFIX}/bin/mpicc" AND EXISTS "${SPACK_OPENMPI_PREFIX}/bin/mpicxx")
            message("-- Using OpenMPI from PETSc sibling prefix: ${SPACK_OPENMPI_PREFIX}")
            set(MPI_C_COMPILER "${SPACK_OPENMPI_PREFIX}/bin/mpicc" CACHE FILEPATH "MPI C compiler" FORCE)
            set(MPI_CXX_COMPILER "${SPACK_OPENMPI_PREFIX}/bin/mpicxx" CACHE FILEPATH "MPI CXX compiler" FORCE)
            include_directories(BEFORE "${SPACK_OPENMPI_PREFIX}/include")
            link_directories("${SPACK_OPENMPI_PREFIX}/lib" "${SPACK_OPENMPI_PREFIX}/lib64")
            list(APPEND CMAKE_PREFIX_PATH "${SPACK_OPENMPI_PREFIX}/lib/cmake" "${SPACK_OPENMPI_PREFIX}/lib64/cmake")
        endif()
    elseif(SPACK_OPENMPI_PREFIX_COUNT GREATER 1)
        message(WARNING "Multiple OpenMPI prefixes found next to PETSC_DIR; using default MPI detection")
    endif()
endif()
find_package(MPI REQUIRED COMPONENTS C CXX)

# Set a default value of 3 for SPATIAL_DIM if not provided
if(NOT DEFINED SPATIAL_DIM)
    set(SPATIAL_DIM 3)
endif()

# Enforce SPATIAL_DIM ∈ {2, 3}
if(SPATIAL_DIM LESS 2 OR SPATIAL_DIM GREATER 3)
    message(FATAL_ERROR "SPATIAL_DIM must be either 2 or 3, but got ${SPATIAL_DIM}")
endif()

# Set the executable name
set(exec_suffix "-${SPATIAL_DIM}D")
if(WITH_CANONICAL_SUFFIX)
    set(exec_suffix "-stk")
endif()
set(exec "openaccel${exec_suffix}.exe")
if(CMAKE_BUILD_TYPE STREQUAL DEBUG OR CMAKE_BUILD_TYPE STREQUAL Debug)
   set(exec "openaccel${exec_suffix}-debug.exe")
endif()

# NOTE [faw 2026-05-04]: silence trilinos stk_topology related deprecation
# warnings
add_compile_options(
    -Wno-deprecated-declarations
)

if(WITH_WARNINGS)
add_compile_options(
    -Wall
    -Wextra
    -Wpedantic
    # -Werror
)
endif()

# Spatial dimension
add_compile_definitions(SPATIAL_DIM=${SPATIAL_DIM})

# various system and dependencies for code in external/
set(THREADS_PREFER_PTHREAD_FLAG ON)
find_package(PkgConfig)
find_package(Threads REQUIRED)
find_package(HDF5 REQUIRED)

# For Trilinos (required for STK, etc.)
set(CMAKE_PREFIX_PATH ${Trilinos_DIR} ${CMAKE_PREFIX_PATH})
find_package(Trilinos NO_MODULE REQUIRED)
message("-- Found Trilinos = ${Trilinos_DIR}")
message("   Trilinos version: ${Trilinos_VERSION}")

set(TRILINOS_INCLUDE_HINTS)
set(TRILINOS_LIBRARY_HINTS)
if(DEFINED ENV{TRILINOS_ROOT})
    list(APPEND TRILINOS_INCLUDE_HINTS "$ENV{TRILINOS_ROOT}/include")
    list(APPEND TRILINOS_LIBRARY_HINTS
         "$ENV{TRILINOS_ROOT}/lib" "$ENV{TRILINOS_ROOT}/lib64")
else()
    # A package-manager installation normally places TrilinosConfig.cmake in
    # <prefix>/lib*/cmake/Trilinos. Derive the prefix from the package that
    # find_package selected instead of embedding a developer-machine path.
    get_filename_component(TRILINOS_CMAKE_DIR "${Trilinos_DIR}" DIRECTORY)
    get_filename_component(TRILINOS_LIB_DIR "${TRILINOS_CMAKE_DIR}" DIRECTORY)
    get_filename_component(TRILINOS_PREFIX "${TRILINOS_LIB_DIR}" DIRECTORY)
    if(EXISTS "${TRILINOS_PREFIX}/include")
        list(APPEND TRILINOS_INCLUDE_HINTS "${TRILINOS_PREFIX}/include")
        list(APPEND TRILINOS_LIBRARY_HINTS
             "${TRILINOS_PREFIX}/lib" "${TRILINOS_PREFIX}/lib64")
    elseif(EXISTS "/usr/include/trilinos")
        list(APPEND TRILINOS_INCLUDE_HINTS "/usr/include/trilinos")
    endif()
endif()
list(APPEND TRILINOS_INCLUDE_HINTS ${Trilinos_INCLUDE_DIRS})
list(APPEND TRILINOS_LIBRARY_HINTS ${Trilinos_LIBRARY_DIRS})
list(REMOVE_DUPLICATES TRILINOS_INCLUDE_HINTS)
list(REMOVE_DUPLICATES TRILINOS_LIBRARY_HINTS)
include_directories(BEFORE SYSTEM ${TRILINOS_INCLUDE_HINTS})
link_directories(${TRILINOS_LIBRARY_HINTS})

# Check if Trilinos linear solver packages are available (Belos, Ifpack2, Tpetra)
set(TRILINOS_LINEAR_SOLVER_AVAILABLE FALSE)
if(WITH_TRILINOS_SOLVER)
    set(TRILINOS_LINEAR_SOLVER_AVAILABLE TRUE)
    foreach(pkg Belos Ifpack2 Tpetra)
        list(FIND Trilinos_PACKAGE_LIST ${pkg} pkg_index)
        if(pkg_index EQUAL -1)
            set(TRILINOS_LINEAR_SOLVER_AVAILABLE FALSE)
            message("   Trilinos package ${pkg}: NOT FOUND")
        else()
            message("   Trilinos package ${pkg}: found")
        endif()
    endforeach()
endif()

if(TRILINOS_LINEAR_SOLVER_AVAILABLE)
    add_compile_definitions(HAS_TRILINOS)
    message("-- Trilinos linear solver support: ENABLED")
else()
    message("-- Trilinos linear solver support: DISABLED (missing required packages)")
endif()

# For YAML
set(CMAKE_PREFIX_PATH ${YAML_DIR}/lib/cmake ${YAML_DIR}/lib64/cmake
    $ENV{YAML_CPP_DIR}/lib/cmake $ENV{YAML_CPP_DIR}/lib64/cmake ${CMAKE_PREFIX_PATH})
find_package(YAML-CPP)
if(YAML-CPP_FOUND)
  # YAML system install is used
  include_directories(${YAML_CPP_INCLUDE_DIR})
  link_directories(${YAML_CPP_LIBRARY_DIR})
else()
  # YAML 0.5.3 is used
  find_library(YAML_CPP_LIBRARIES NAMES yaml-cpp PATHS ${YAML_DIR}/lib)
  find_path(YAML_CPP_INCLUDE_DIR yaml.h PATHS ${YAML_DIR}/include/yaml-cpp)
  if((DEFINED YAML_CPP_LIBRARIES) AND (DEFINED YAML_CPP_INCLUDE_DIR))
    include_directories(SYSTEM ${YAML_CPP_INCLUDE_DIR}/..)
    set(YAML-CPP_FOUND TRUE)
  endif()
endif()
if(YAML-CPP_FOUND)
  message("-- Found YAML-CPP = ${YAML_DIR}")
else()
  message(FATAL_ERROR "YAML-CPP NOT FOUND")
endif()

# PETSc
if(DEFINED PETSC_DIR)
    if(NOT IS_DIRECTORY "${PETSC_DIR}")
        message(FATAL_ERROR "PETSC_DIR is defined but does not point to a directory: ${PETSC_DIR}")
    endif()

    if(NOT EXISTS "${PETSC_DIR}/include/petscksp.h" AND
       NOT EXISTS "${PETSC_DIR}/${PETSC_ARCH}/include/petscksp.h" AND
       NOT EXISTS "${PETSC_DIR}/include/petsc/finclude/petscksp.h")
        message(FATAL_ERROR "Cannot find petscksp.h in PETSC_DIR. Please verify PETSC_DIR and PETSC_ARCH.\n  PETSC_DIR=${PETSC_DIR}\n  PETSC_ARCH=${PETSC_ARCH}")
    endif()

    message("-- PETSC_DIR (${PETSC_DIR}) is provided. PETSc functionalities are enabled.")
    add_compile_definitions(HAS_PETSC)
    include_directories(BEFORE ${PETSC_DIR}/include)

    if(DEFINED PETSC_ARCH)
        include_directories(BEFORE ${PETSC_DIR}/${PETSC_ARCH}/include)
    endif()
else()
    if (PKG_CONFIG_FOUND)
        # https://gitlab.com/petsc/petsc/-/merge_requests/5992
        pkg_check_modules(PETSC QUIET petsc)

        if(PETSC_FOUND)
            add_compile_definitions(HAS_PETSC)
            include_directories(BEFORE ${PETSC_INCLUDE_DIRS})
            link_directories(${PETSC_LIBRARY_DIRS})
        else()
            message("-- PETSc functionalities are DISABLED")
        endif()
    endif()
endif()

# HYPRE
find_package(HYPRE QUIET)
if(HYPRE_FOUND)
    message("-- HYPRE found and enabled.")
    add_compile_definitions(HAS_HYPRE)
    include_directories(${HYPRE_INCLUDE_DIRS})
    link_directories(${HYPRE_LIBRARY_DIRS})
else()
    message("-- HYPRE functionalities are DISABLED")
endif()

# Solid mechanics assembler
if(WITH_CVFEM_SOLID_MECHANICS)
    add_compile_definitions(USE_CVFEM_SOLID_MECHANICS)
    message("-- Solid mechanics assembler: CVFEM")
else()
    message("-- Solid mechanics assembler: FEM")
endif()

# Which thermal energy equation to use?
if(WITH_THERMAL_TEMPERATURE)
    add_compile_definitions(HAS_THERMAL_TEMPERATURE)
endif()

# linear solver library
add_compile_definitions(NO_LINSOLVE_CONF)
add_subdirectory(src/solver)
if(WITH_AMGSOLVER)
    get_filename_component(libAMG_PATH src/solver/libamg REALPATH)
    if(NOT EXISTS ${libAMG_PATH})
        message(FATAL_ERROR "WITH_AMGSOLVER set but library source does not exist")
    endif()
    message("-- Compiling sources with AMGSolver")
endif()

# Enable Fortran language support (if not already done)
enable_language(Fortran)

# Add main object library on which targets below depend on
add_library(main_obj OBJECT)

# Add target sources and header directories
add_subdirectory(src)
add_subdirectory(external)

# third-party includes
target_include_directories(main_obj PUBLIC
    ${libLinSolve_INC}
    ${Trilinos_INCLUDE_DIRS}
    ${Trilinos_TPL_INCLUDE_DIRS}
)

# linear solver library
target_link_libraries(main_obj PUBLIC ${libLinSolve})

# Link your executable with Trilinos libraries
target_link_libraries(main_obj PUBLIC ${Trilinos_LIBRARIES})

# Link your executable with MPI libraries
target_link_libraries(main_obj PUBLIC MPI::MPI_C MPI::MPI_CXX)

# Link to Git version object and various dependencies for external/
target_link_libraries(main_obj PRIVATE git_version)
target_link_libraries(main_obj PRIVATE Threads::Threads)
target_link_libraries(main_obj PRIVATE HDF5::HDF5)

# Link your executable with YAML libraries
if(EXISTS ${YAML_DIR}/lib/libyaml-cpp.so)
    # Link to the library directly if it exists
    target_link_libraries(main_obj PRIVATE ${YAML_DIR}/lib/libyaml-cpp.so)
else()
    # Fallback to using YAML_CPP_LIBRARIES if the direct path does not work
    target_link_libraries(main_obj PRIVATE ${YAML_CPP_LIBRARIES})
endif()

# Link to PETSc libraries
if(PETSC_FOUND)
    if(${CMAKE_SYSTEM_NAME} MATCHES Darwin)
        target_link_libraries(main_obj PUBLIC ${PETSC_INCLUDE_DIRS}/../lib/libpetsc.3.20.1.dylib)
    else()
        target_link_libraries(main_obj PUBLIC ${PETSC_LIBRARIES})
    endif()
else()
    if(DEFINED PETSC_DIR)
        if(DEFINED PETSC_ARCH)
            target_link_libraries(main_obj PUBLIC ${PETSC_DIR}/${PETSC_ARCH}/lib/libpetsc.so)
        else()
            target_link_libraries(main_obj PUBLIC petsc)
            target_link_directories(main_obj PUBLIC ${PETSC_DIR}/lib ${PETSC_DIR}/lib64)
        endif()
    endif()
endif()

# Link to HYPRE libraries
if(HYPRE_FOUND)
    target_link_libraries(main_obj PUBLIC HYPRE::HYPRE)
endif()

if(WITH_VECTORIZATION)
    include(CheckCXXCompilerFlag)
    check_cxx_compiler_flag("-mavx2" COMPILER_SUPPORTS_AVX2)
    check_cxx_compiler_flag("-mfma" COMPILER_SUPPORTS_FMA)
    if(COMPILER_SUPPORTS_AVX2 AND COMPILER_SUPPORTS_FMA)
        string(APPEND CMAKE_CXX_FLAGS_RELEASE " -mavx2 -mfma")
        string(APPEND CMAKE_CXX_FLAGS_RELWITHDEBINFO " -mavx2 -mfma")
        message("-- Using AVX2/FMA vectorization flags")
    else()
        message("-- AVX2/FMA not supported by compiler, vectorization disabled")
    endif()
endif()

if(WITH_FASTMATH)
    string(APPEND CMAKE_CXX_FLAGS_RELEASE " -Ofast -mavx2 -mfma")
    string(APPEND CMAKE_CXX_FLAGS_RELWITHDEBINFO " -Ofast -mavx2 -mfma")
    message("-- Using fast math optimizations (NOT compliant with IEEE)")
endif()

add_executable(${exec} main.cpp)
target_link_libraries(${exec} PRIVATE main_obj)

# debugging utilities
if(WITH_DEBUG_LIB)
    add_library(debug OBJECT
        debug/debug.cpp
    )
    target_compile_options(debug PRIVATE -O0 -g -UNDEBUG) # always build this object library in debug mode
    target_link_libraries(debug PRIVATE main_obj)
endif()
