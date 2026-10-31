include(FindPackageHandleStandardArgs)

find_path(HYPRE_INCLUDE_DIR
    NAMES HYPRE.h
    HINTS ${HYPRE_ROOT} $ENV{HYPRE_ROOT}
    PATH_SUFFIXES include
)

find_library(HYPRE_LIBRARY
    NAMES HYPRE
    HINTS ${HYPRE_ROOT} $ENV{HYPRE_ROOT}
    PATH_SUFFIXES lib lib64
)

find_package_handle_standard_args(HYPRE
    REQUIRED_VARS HYPRE_INCLUDE_DIR HYPRE_LIBRARY
)

if(HYPRE_FOUND)
    set(HYPRE_INCLUDE_DIRS ${HYPRE_INCLUDE_DIR})
    set(HYPRE_LIBRARIES ${HYPRE_LIBRARY})
    get_filename_component(HYPRE_LIBRARY_DIRS ${HYPRE_LIBRARY} DIRECTORY)

    if(NOT TARGET HYPRE::HYPRE)
        add_library(HYPRE::HYPRE UNKNOWN IMPORTED)
        set_target_properties(HYPRE::HYPRE PROPERTIES
            IMPORTED_LOCATION "${HYPRE_LIBRARY}"
            INTERFACE_INCLUDE_DIRECTORIES "${HYPRE_INCLUDE_DIR}"
        )
    endif()
endif()

mark_as_advanced(HYPRE_INCLUDE_DIR HYPRE_LIBRARY)
