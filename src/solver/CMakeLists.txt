# Linear solver wrappers (PETSc/HYPRE interfaces and sparse matrix utils)
set(DIR ${CMAKE_CURRENT_SOURCE_DIR})

# project sources
target_sources(main_obj PRIVATE
    ${DIR}/matrix/CRSNodeGraph.cpp
    ${DIR}/utils/Profiler.cpp
)

# Optional Trilinos linear solver support
# Set Trilinos_DIR in do-configure.sh or via -DTrilinos_DIR=<path>
# Requires Belos, Ifpack2, and Tpetra packages to be enabled in Trilinos
find_package(Trilinos QUIET)
if(Trilinos_FOUND)
    message(STATUS "Found Trilinos: ${Trilinos_DIR}")
    message(STATUS "Trilinos version: ${Trilinos_VERSION}")
    # Check if linear solver packages are available
    set(TRILINOS_LINEAR_SOLVER_AVAILABLE TRUE)
    foreach(pkg Belos Ifpack2 Tpetra)
        list(FIND Trilinos_PACKAGE_LIST ${pkg} pkg_index)
        if(pkg_index EQUAL -1)
            set(TRILINOS_LINEAR_SOLVER_AVAILABLE FALSE)
            message(STATUS "  Trilinos package ${pkg}: NOT FOUND")
        else()
            message(STATUS "  Trilinos package ${pkg}: found")
        endif()
    endforeach()
    if(TRILINOS_LINEAR_SOLVER_AVAILABLE)
        add_compile_definitions(HAS_TRILINOS)
        message(STATUS "Trilinos linear solver support: ENABLED")
    else()
        message(STATUS "Trilinos linear solver support: DISABLED (missing required packages)")
    endif()
else()
    message(STATUS "Trilinos not found - Trilinos linear solver support disabled")
endif()

# header search paths
target_include_directories(main_obj PUBLIC
    ${DIR}
    ${DIR}/HYPRE
    ${DIR}/PETSc
    ${DIR}/Trilinos
    ${DIR}/matrix
    ${DIR}/matrix/layout
    ${DIR}/matrix/operators/blockMatrix
    ${DIR}/matrix/operators/distributed
    ${DIR}/utils
)

if(TRILINOS_LINEAR_SOLVER_AVAILABLE)
    set(HAS_TRILINOS TRUE PARENT_SCOPE)
    set(Trilinos_INCLUDE_DIRS ${Trilinos_INCLUDE_DIRS} PARENT_SCOPE)
    set(Trilinos_LIBRARIES ${Trilinos_LIBRARIES} PARENT_SCOPE)
endif()
