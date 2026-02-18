# File       : CMakeLists.txt
# Created    : Thu Oct 16 2025 08:17:36 (+0200)
# Author     : Fabian Wermelinger
# Description: CMake build definitions
# Copyright (c) 2025 CCFNUM, Lucerne University of Applied Sciences and Arts.
# SPDX-License-Identifier: BSD-3-Clause

set(DIR ${CMAKE_CURRENT_SOURCE_DIR})

# project sources
target_sources(main_obj PRIVATE
    ${DIR}/graph/nodeGraph.cpp
    ############################################################################
    ${DIR}/interface/dgInfo.cpp
    ${DIR}/interface/interface.cpp
    ${DIR}/interface/interfaceIO.cpp
    ${DIR}/interface/interfaceSideInfo.cpp
    ${DIR}/interface/transfer/dataTransfer.cpp
    ############################################################################
    ${DIR}/mesh.cpp
    ${DIR}/meshGeometry.cpp
    ${DIR}/meshIO.cpp
    ${DIR}/meshNodeGraph.cpp
    ############################################################################
    ${DIR}/validation/elementValidator.cpp
    ############################################################################
    ${DIR}/zone/boundary/boundary.cpp
    ${DIR}/zone/motion/zoneDeformation.cpp
    ${DIR}/zone/motion/zoneTransformation.cpp
    ${DIR}/zone/zone.cpp
    ${DIR}/zone/zoneIO.cpp
)

# define header directories (required since #include statements do not quote
# full path and location of headers is inside src)
target_include_directories(main_obj PUBLIC
    ${DIR}
    ############################################################################
    ${DIR}/graph
    ############################################################################
    ${DIR}/interface
    ${DIR}/interface/transfer
    ############################################################################
    ${DIR}/validation
    ############################################################################
    ${DIR}/zone
    ${DIR}/zone/boundary
    ${DIR}/zone/motion
)
