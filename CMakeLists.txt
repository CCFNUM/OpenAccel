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

# HYPRE (Strict Version Check)
set(HYPRE_FOUND FALSE)

if(DEFINED HYPRE_DIR)
    set(HYPRE_SEARCH_PATHS ${HYPRE_DIR}/include)
else()
    set(HYPRE_SEARCH_PATHS /usr/include/hypre /usr/local/include/hypre /usr/include)
endif()

# 1. Look for the header and library
find_path(HYPRE_INCLUDE_DIR NAMES HYPRE.h hypre.h PATHS ${HYPRE_SEARCH_PATHS} NO_DEFAULT_PATH)

if(DEFINED HYPRE_DIR)
    find_library(HYPRE_LIB NAMES HYPRE hypre
                 PATHS ${HYPRE_DIR}/lib ${HYPRE_DIR}/lib64
                 NO_DEFAULT_PATH)
else()
    find_library(HYPRE_LIB NAMES HYPRE hypre
                 PATHS /usr/lib/x86_64-linux-gnu /usr/local/lib /usr/lib /usr/lib64)
endif()

if(HYPRE_INCLUDE_DIR AND HYPRE_LIB)
    # 2. Extract version from HYPRE_config.h for compatibility check
    find_path(HYPRE_CONFIG_PATH NAMES HYPRE_config.h PATHS ${HYPRE_INCLUDE_DIR})
    
    if(HYPRE_CONFIG_PATH)
        file(STRINGS "${HYPRE_CONFIG_PATH}/HYPRE_config.h" HYPRE_VER_LINE REGEX "HYPRE_RELEASE_NUMBER")
        if(HYPRE_VER_LINE MATCHES "([0-9]+)")
            set(HYPRE_VERSION_VAL ${CMAKE_MATCH_1})

            # 3. Enforce Version 2.29.0 (22900) or higher
            if(${HYPRE_VERSION_VAL} LESS 22900)
                message(WARNING "HYPRE version ${HYPRE_VERSION_VAL} is too old. Required: 2.29.0+. Disabling HYPRE.")
                set(HYPRE_FOUND FALSE)
            else()
                set(HYPRE_FOUND TRUE)
                set(HYPRE_INCLUDE_DIRS ${HYPRE_INCLUDE_DIR})
                set(HYPRE_LIBRARIES ${HYPRE_LIB})
            endif()
        endif()
    else()
        message(WARNING "Found HYPRE headers but could not find HYPRE_config.h to verify version.")
    endif()
endif()

# 4. Final Activation Logic
if(HYPRE_FOUND)
    message("-- HYPRE version ${HYPRE_VERSION_VAL} found and enabled.")
    add_compile_definitions(HAS_HYPRE)
    include_directories(${HYPRE_INCLUDE_DIRS})
else()
    message("-- HYPRE functionalities are DISABLED (Version mismatch or not found).")
endif()

# Interface support (enabled by default)
option(HAS_INTERFACE "Enable interface support" ON)
if(HAS_INTERFACE)
    add_compile_definitions(HAS_INTERFACE)
    message("-- Interface support: ENABLED")
else()
    message("-- Interface support: DISABLED")
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
if(HYPRE_FOUND)
    # If we found it via the custom HYPRE_DIR path
    if(DEFINED HYPRE_DIR)
        # Find the actual library file in the custom directory
        find_library(HYPRE_CUSTOM_LIB NAMES HYPRE hypre 
                     PATHS ${HYPRE_DIR}/lib ${HYPRE_DIR}/lib64 
                     NO_DEFAULT_PATH)
        target_link_libraries(main_obj PUBLIC ${HYPRE_CUSTOM_LIB})
        message("-- Linking against custom HYPRE: ${HYPRE_CUSTOM_LIB}")
    else()
        # Fallback for system-wide/pkg-config installs
        target_link_libraries(main_obj PUBLIC ${HYPRE_LIBRARIES})
    endif()
else()
    # ... your existing 'disabled' logic ...
endif()

if(DEFINED WITH_FASTMATH AND WITH_FASTMATH)
    string(APPEND CMAKE_CXX_FLAGS_RELEASE " -Ofast -mavx2 -mfma")
    string(APPEND CMAKE_CXX_FLAGS_RELWITHDEBINFO " -Ofast -mavx2 -mfma")
    message("-- Using fast math optimizations (NOT compliant with IEEE)")
endif()

add_executable(${exec} main.cpp)
target_link_libraries(${exec} PRIVATE main_obj)
