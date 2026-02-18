# File       : CMakeLists.txt
# Created    : Mon Feb 09 2026 08:42:10 (+0100)
# Author     : Mhamad Mahdi Alloush
# Description:
# Copyright (c) 2026 CCFNUM, Lucerne University of Applied Sciences and Arts.
# SPDX-License-Identifier: BSD-3-Clause

set(DIR ${CMAKE_CURRENT_SOURCE_DIR})

# project sources
target_sources(main_obj PRIVATE
    ${DIR}/overrides.cpp
)

# define header directories (required since #include statements do not quote
# full path and location of headers is inside src)
target_include_directories(main_obj PUBLIC
    ${DIR}
)
