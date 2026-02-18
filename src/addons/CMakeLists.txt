# File       : CMakeLists.txt
# Created    : Thu Oct 16 2025 08:17:36 (+0200)
# Author     : Fabian Wermelinger
# Description: CMake build definitions
# Copyright (c) 2025 CCFNUM, Lucerne University of Applied Sciences and Arts.
# SPDX-License-Identifier: BSD-3-Clause

set(DIR ${CMAKE_CURRENT_SOURCE_DIR})

# project sources
target_sources(main_obj PRIVATE
    ${DIR}/meshMotion/deformation/deformation.cpp
    ${DIR}/meshMotion/deformation/displacementDiffusion/displacementDiffusionAssembler.cpp
    ${DIR}/meshMotion/deformation/displacementDiffusion/displacementDiffusionAssemblerElemBoundaryConditions.cpp
    ${DIR}/meshMotion/deformation/displacementDiffusion/displacementDiffusionEquation.cpp
    ${DIR}/meshMotion/deformation/displacementDiffusion/displacementDiffusionModel.cpp
    ${DIR}/meshMotion/meshMotion.cpp
    ${DIR}/meshMotion/transformation/transformation.cpp
    ${DIR}/rigidBody/rigidBody.cpp
    ${DIR}/wallDistance/wallDistance.cpp
    ${DIR}/wallDistance/wallScale/wallScale.cpp
    ${DIR}/wallDistance/wallScale/wallScaleDiffusion/wallScaleDiffusionModel.cpp
    ${DIR}/wallDistance/wallScale/wallScaleDiffusion/wallScaleDiffusionEquation.cpp
    ${DIR}/wallDistance/wallScale/wallScaleDiffusion/wallScaleDiffusionAssemblerElemBoundaryConditions.cpp
    ${DIR}/wallDistance/wallScale/wallScaleDiffusion/wallScaleDiffusionAssemblerElemTerms.cpp
    ${DIR}/wallDistance/wallScale/wallScaleDiffusion/wallScaleDiffusionAssemblerNodeTerms.cpp
)

# define header directories (required since #include statements do not quote
# full path and location of headers is inside src)
target_include_directories(main_obj PUBLIC
    ${DIR}/meshMotion
    ${DIR}/meshMotion/deformation
    ${DIR}/meshMotion/deformation/displacementDiffusion
    ${DIR}/meshMotion/transformation
    ${DIR}/rigidBody
    ${DIR}/wallDistance
    ${DIR}/wallDistance/wallScale
    ${DIR}/wallDistance/wallScale/wallScaleDiffusion
)
