set(CURRENT_LIST_DIR ${CMAKE_CURRENT_LIST_DIR})
if (NOT DEFINED pre_configure_dir)
    set(pre_configure_dir ${CMAKE_CURRENT_LIST_DIR})
endif()

if (NOT DEFINED post_configure_dir)
    set(post_configure_dir ${CMAKE_BINARY_DIR}/generated)
endif()

set(pre_configure_file ${pre_configure_dir}/version.cpp.in)
set(post_configure_file ${post_configure_dir}/version.cpp)

function(CheckGitWrite git_hash)
    file(WRITE ${CMAKE_BINARY_DIR}/git-state.txt ${git_hash})
endfunction()

function(CheckGitRead git_hash)
    if (EXISTS ${CMAKE_BINARY_DIR}/git-state.txt)
        file(STRINGS ${CMAKE_BINARY_DIR}/git-state.txt CONTENT)
        LIST(GET CONTENT 0 var)

        set(${git_hash} ${var} PARENT_SCOPE)
    endif()
endfunction()

function(CheckGitVersion)
    if (GIT_REPO_BUILD)
        find_package(Git REQUIRED QUIET)
        execute_process(
            COMMAND ${GIT_EXECUTABLE} describe --tags --always --long --dirty --broken
            WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
            OUTPUT_VARIABLE GIT_DESCRIBE
            OUTPUT_STRIP_TRAILING_WHITESPACE
            )
        execute_process(
            COMMAND ${GIT_EXECUTABLE} rev-parse HEAD
            WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
            OUTPUT_VARIABLE GIT_COMMIT_HASH
            OUTPUT_STRIP_TRAILING_WHITESPACE
        )
    else()
        set(GIT_DESCRIBE "out-of-repo-build")
        set(GIT_COMMIT_HASH "out-of-repo-build")
    endif()

    # the following region is special, do not remove the marker comments below
    # %GIT_TAG%
    # %EXTERNAL_REPLACE_START%
    if (NOT GIT_REPO_BUILD)
        message(FATAL_ERROR "Cannot determine version: this is not a proper release tarball, neither a Git repository/worktree")
    endif()

    execute_process(
        COMMAND ${GIT_EXECUTABLE} describe --always
        WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
        OUTPUT_VARIABLE GIT_TAG
        OUTPUT_STRIP_TRAILING_WHITESPACE
        RESULT_VARIABLE GIT_RESULT
    )

    if (NOT GIT_RESULT EQUAL 0)
        message(FATAL_ERROR "Failed to get version from git tag")
    endif()
    # %EXTERNAL_REPLACE_END%

    string(REGEX REPLACE "^v-?" "" CLEAN_TAG "${GIT_TAG}")
    if(CLEAN_TAG MATCHES "^([0-9]+)\\.([0-9]+)\\.([0-9]+)(-.*)?$")
        set(PROJECT_VERSION_MAJOR ${CMAKE_MATCH_1})
        set(PROJECT_VERSION_MINOR ${CMAKE_MATCH_2})
        set(PROJECT_VERSION_PATCH ${CMAKE_MATCH_3})

        set(GIT_HASH "${GIT_COMMIT_HASH}")
        set(GIT_DESCRIBE "${GIT_DESCRIBE}")
        set(PROJECT_VERSION "${CMAKE_MATCH_1}.${CMAKE_MATCH_2}.${CMAKE_MATCH_3}") 
    else()
        message(WARNING "Tag '${CLEAN_TAG}' does not match semver format")
        return()
    endif()

    # caching
    CheckGitRead(GIT_HASH_CACHE)
    if (NOT EXISTS ${post_configure_dir})
        file(MAKE_DIRECTORY ${post_configure_dir})
    endif()

    if (NOT EXISTS ${post_configure_dir}/version.h)
        file(COPY ${pre_configure_dir}/version.h DESTINATION ${post_configure_dir})
    endif()

    if (NOT DEFINED GIT_HASH_CACHE)
        set(GIT_HASH_CACHE "INVALID")
    endif()

    # Only update the version.cpp if the hash has changed. This will
    # prevent us from rebuilding the project more than we need to.
    if (NOT ${GIT_DESCRIBE} STREQUAL ${GIT_HASH_CACHE} OR NOT EXISTS ${post_configure_file})
        # Set the GIT_HASH_CACHE variable the next build won't have
        # to regenerate the source file.
        CheckGitWrite(${GIT_DESCRIBE})
        configure_file(${pre_configure_file} ${post_configure_file} @ONLY)
    endif()
endfunction()

function(CheckGitSetup)
    if (NOT DEFINED GIT_REPO_BUILD)
        # check if this is a git repository or a git worktree
        if (IS_DIRECTORY ${PROJECT_SOURCE_DIR}/.git OR EXISTS ${PROJECT_SOURCE_DIR}/.git)
            set(GIT_REPO_BUILD TRUE)
        else()
            set(GIT_REPO_BUILD FALSE)
        endif()
    endif()

    add_custom_target(AlwaysCheckGit COMMAND ${CMAKE_COMMAND}
        -DRUN_CHECK_GIT_REVISION=1
        -DGIT_REPO_BUILD=${GIT_REPO_BUILD}
        -Dpre_configure_dir=${pre_configure_dir}
        -Dpost_configure_file=${post_configure_dir}
        -DGIT_HASH_CACHE=${GIT_HASH_CACHE}
        -P ${CURRENT_LIST_DIR}/CheckGit.cmake
        BYPRODUCTS ${post_configure_file}
    )

    add_library(git_version ${CMAKE_BINARY_DIR}/generated/version.cpp)
    target_include_directories(git_version PUBLIC ${CMAKE_BINARY_DIR}/generated)
    add_dependencies(git_version AlwaysCheckGit)

    CheckGitVersion()
endfunction()

# This is used to run this function from an external cmake process.
if (RUN_CHECK_GIT_REVISION)
    CheckGitVersion()
endif()
