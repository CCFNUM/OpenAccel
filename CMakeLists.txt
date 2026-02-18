cmake_minimum_required(VERSION 3.10)

project(Accel VERSION 1.0 LANGUAGES C CXX Fortran)

# Set C++ standard (optional but recommended)
set(CMAKE_CXX_STANDARD 20)

# Get Git hash
include(cmake/CheckGit.cmake)
CheckGitSetup()

# Find MPI package
find_package(MPI REQUIRED)

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
if(DEFINED CANONICAL_EXEC_SUFFIX AND CANONICAL_EXEC_SUFFIX)
    set(exec_suffix "-stk")
endif()
set(exec "accel${exec_suffix}.exe")
if(CMAKE_BUILD_TYPE STREQUAL DEBUG OR CMAKE_BUILD_TYPE STREQUAL Debug)
   set(exec "accel${exec_suffix}-debug.exe")
endif()

if(DEFINED ACCEL_WARNING AND ACCEL_WARNING)
add_compile_options(
    -Wall
    -Wextra
    -Wpedantic
    # -Werror
)
endif()

# Spatial dimension
add_compile_definitions(SPATIAL_DIM=${SPATIAL_DIM})

# Detect gnuplot for residual plotting support.
# If GNUPLOT_DIR is provided, search only within that path; otherwise search
# standard system locations via find_program.
set(GNUPLOT_DIR "" CACHE PATH "Path to gnuplot installation (prefix or bin directory)")
if(GNUPLOT_DIR)
    if(NOT IS_DIRECTORY "${GNUPLOT_DIR}")
        message(FATAL_ERROR
            "GNUPLOT_DIR is set to '${GNUPLOT_DIR}', but this directory does not exist.")
    endif()

    set(GNUPLOT_EXECUTABLE "")
    set(_gnuplot_candidates
        "${GNUPLOT_DIR}/gnuplot"
        "${GNUPLOT_DIR}/pgnuplot"
        "${GNUPLOT_DIR}/gnuplot.exe"
        "${GNUPLOT_DIR}/pgnuplot.exe"
        "${GNUPLOT_DIR}/bin/gnuplot"
        "${GNUPLOT_DIR}/bin/pgnuplot"
        "${GNUPLOT_DIR}/bin/gnuplot.exe"
        "${GNUPLOT_DIR}/bin/pgnuplot.exe"
    )
    foreach(_gnuplot_candidate IN LISTS _gnuplot_candidates)
        if(EXISTS "${_gnuplot_candidate}")
            set(GNUPLOT_EXECUTABLE "${_gnuplot_candidate}")
            break()
        endif()
    endforeach()

    if(NOT GNUPLOT_EXECUTABLE)
        message(FATAL_ERROR
            "GNUPLOT_DIR is set to '${GNUPLOT_DIR}', but no gnuplot executable "
            "was found. Expected one of: gnuplot, pgnuplot, gnuplot.exe, "
            "pgnuplot.exe (directly in GNUPLOT_DIR or in GNUPLOT_DIR/bin).")
    endif()
else()
    # No GNUPLOT_DIR provided: search standard system paths.
    # Unset the cache entry so every cmake run performs a fresh search and
    # does not reuse a stale NOTFOUND result from a previous invocation.
    unset(GNUPLOT_EXECUTABLE CACHE)
    find_program(GNUPLOT_EXECUTABLE
        NAMES gnuplot gnuplot.exe pgnuplot pgnuplot.exe
    )
endif()

if(GNUPLOT_EXECUTABLE)
    set(HAS_GNUPLOT TRUE)
    add_compile_definitions(HAS_GNUPLOT)
    message("-- Gnuplot residual plotting: ENABLED (${GNUPLOT_EXECUTABLE})")
else()
    set(HAS_GNUPLOT FALSE)
    message("-- Gnuplot residual plotting: DISABLED (gnuplot not found)")
endif()

# PkgConfig
find_package(PkgConfig)

# For Trilinos (required for STK, etc.)
set(CMAKE_PREFIX_PATH ${Trilinos_DIR} ${CMAKE_PREFIX_PATH})
find_package(Trilinos NO_MODULE REQUIRED)
message("-- Found Trilinos = ${Trilinos_DIR}")
message("   Trilinos version: ${Trilinos_VERSION}")

# Check if Trilinos linear solver packages are available (Belos, Ifpack2, Tpetra)
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
  # YAML master branch is used
  include_directories(SYSTEM ${YAML_CPP_INCLUDE_DIR})
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
    message("PETSC_DIR (${PETSC_DIR}) is provided. PETSc functionalities are enabled.")
    add_compile_definitions(HAS_PETSC)
    include_directories(${PETSC_DIR}/include)

    if(DEFINED PETSC_ARCH)
        include_directories(${PETSC_DIR}/${PETSC_ARCH}/include)
    endif()
else()
    if (PKG_CONFIG_FOUND)
        # https://gitlab.com/petsc/petsc/-/merge_requests/5992
        pkg_check_modules(PETSC PETSc)

        if(PETSC_FOUND)
            message("PETSc functionalities are enabled.")
            add_compile_definitions(HAS_PETSC)
            include_directories(${PETSC_INCLUDE_DIRS})
        else()
            message("PETSc functionalities are disabled.")
        endif()
    endif()
endif()

# HYPRE
if(DEFINED HYPRE_DIR)
    message("HYPRE_DIR (${HYPRE_DIR}) is provided. HYPRE functionalities are enabled.")
    find_library(HYPRE_LIBRARIES NAMES HYPRE PATHS ${HYPRE_DIR}/lib
    ${HYPRE_DIR}/lib64)
    if(DEFINED HYPRE_LIBRARIES)
        include_directories(${HYPRE_DIR}/include)
        set(HYPRE_FOUND TRUE)
    endif()
endif()
if(HYPRE_FOUND)
    add_compile_definitions(HAS_HYPRE)
    message("-- Found HYPRE = ${HYPRE_DIR}")
endif()

# Solid mechanics assembler
option(USE_CVFEM_SOLID_MECHANICS "Use CVFEM for solid mechanics" ON)
if(USE_CVFEM_SOLID_MECHANICS)
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

# Enable Fortran language support (if not already done)
enable_language(Fortran)

# Add main object library on which targets below depend on
add_library(main_obj OBJECT)

# Add target sources and header directories
add_subdirectory(src)
add_subdirectory(external)

# third-party includes
target_include_directories(main_obj PUBLIC
    ${MPI_C_INCLUDE_DIRS}
    ${libLinSolve_INC}
    ${Trilinos_INCLUDE_DIRS}
    ${Trilinos_TPL_INCLUDE_DIRS}
)

# Link your executable with liblinsolve and Trilinos libraries
target_link_libraries(main_obj PUBLIC ${libLinSolve})
target_link_libraries(main_obj PUBLIC ${Trilinos_LIBRARIES})

# Link your executable with MPI libraries
target_link_libraries(main_obj PUBLIC ${MPI_C_LIBRARIES} ${MPI_CXX_LIBRARIES})

# Link to Git hash object
target_link_libraries(main_obj PRIVATE git_revision)

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
if(DEFINED HYPRE_DIR)
    if(DEFINED HYPRE_ARCH)
        target_link_libraries(main_obj PUBLIC ${HYPRE_DIR}/${HYPRE_ARCH}/lib/libHYPRE.so)
    else()
        target_link_libraries(main_obj PUBLIC HYPRE)
        target_link_directories(main_obj PUBLIC ${HYPRE_DIR}/lib ${HYPRE_DIR}/lib64)
    endif()
else()
    if(HYPRE_FOUND)
        if(${CMAKE_SYSTEM_NAME} MATCHES Darwin)
            target_link_libraries(main_obj PUBLIC ${HYPRE_INCLUDE_DIRS}/../lib/libHYPRE.3.20.1.dylib)
        else()
            target_link_libraries(main_obj PUBLIC ${HYPRE_LIBRARIES})
        endif()
    endif()
endif()

if(DEFINED WITH_FASTMATH AND WITH_FASTMATH)
    string(APPEND CMAKE_CXX_FLAGS_RELEASE " -Ofast -mavx2 -mfma")
    string(APPEND CMAKE_CXX_FLAGS_RELWITHDEBINFO " -Ofast -mavx2 -mfma")
    message("-- Using fast math optimizations (NOT compliant with IEEE)")
endif()

add_executable(${exec} main.cpp)
target_link_libraries(${exec} PRIVATE main_obj)
