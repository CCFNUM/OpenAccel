# CheckGit.cmake — generates git_revision.h with the current commit hash

function(CheckGitSetup)
    # Find git
    find_package(Git QUIET)

    # Generate git_revision.h at configure time
    set(GIT_REVISION "unknown")
    if(GIT_FOUND)
        execute_process(
            COMMAND ${GIT_EXECUTABLE} rev-parse --short HEAD
            WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
            OUTPUT_VARIABLE GIT_REVISION
            OUTPUT_STRIP_TRAILING_WHITESPACE
            ERROR_QUIET
        )
    endif()

    # Write the header
    set(GIT_REVISION_H "${CMAKE_BINARY_DIR}/git_revision.h")
    file(WRITE ${GIT_REVISION_H}
        "#ifndef GIT_REVISION_H\n"
        "#define GIT_REVISION_H\n"
        "namespace accel {\n"
        "  constexpr const char* git_revision = \"${GIT_REVISION}\";\n"
        "}\n"
        "#endif\n"
    )

    # Create a library target so we can link against it
    add_library(git_revision INTERFACE)
    target_include_directories(git_revision INTERFACE ${CMAKE_BINARY_DIR})
endfunction()
