// File       : dgInfo.h
// Created    : Fri Aug 25 2023 12:55:24 (+0100)
// Author     : Mhamad Mahdi Alloush
// Description: A storage class for an ip info on an interface side
// Copyright (c) 2023 CCFNUM, Lucerne University of Applied Sciences and Arts.
// SPDX-License-Identifier: BSD-3-Clause

#ifndef DGINFO_H
#define DGINFO_H

#ifdef HAS_INTERFACE

// code
#include "types.h"

namespace accel
{

class MasterElement;

class dgInfo
{
public:
    // constructor and destructor
    dgInfo(label parallelRank,
           uint64_t globalFaceId,
           uint64_t localGaussPointId,
           label currentGaussPointId,
           stk::mesh::Entity currentFace,
           stk::mesh::Entity currentElement,
           const label currentFaceOrdinal,
           MasterElement* meFCCurrent,
           MasterElement* meSCSurrent,
           stk::topology currentElementTopo,
           scalar searchTolerance);

    ~dgInfo();

    void dumpInfo();

    label parallelRank_;

    uint64_t globalFaceId_;

    uint64_t localGaussPointId_;

    label currentGaussPointId_;

    stk::mesh::Entity currentFace_;

    stk::mesh::Entity currentElement_;

    const label currentFaceOrdinal_;

    MasterElement* meFCCurrent_;

    MasterElement* meSCSCurrent_;

    stk::topology currentElementTopo_;

    const scalar bestXRef_;

    scalar bestX_;

    scalar nearestDistance_;

    const scalar nearestDistanceSafety_;

    bool gaussPointExposed_;

    label opposingFaceIsGhosted_;

    // search provides opposing face
    stk::mesh::Entity opposingFace_;

    // face:element relations provide connected element to opposing face
    stk::mesh::Entity opposingElement_;

    // opposing element topo
    stk::topology opposingElementTopo_;

    // for the opposing face, what is its ordinal?
    label opposingFaceOrdinal_;

    // master element for opposing face
    MasterElement* meFCOpposing_;

    // master element for opposing face connected element
    MasterElement* meSCSOpposing_;

    // coordinates of gauss points on current face
    std::vector<scalar> currentGaussPointCoords_;

    // iso-parametric coordinates for gauss point on current face (-1:1)
    std::vector<scalar> currentIsoParCoords_;

    // iso-parametric coordinates for gauss point on opposing face (-1:1)
    // the coordinates will not change upon any rigid motion of the face
    std::vector<scalar> opposingIsoParCoords_;

    // opposing gauss point id (set only for conformal interfaces)
    label opposingGaussPointId_ = -1;

    // possible reuse
    std::vector<uint64_t> allOpposingFaceIds_;
};

} // namespace accel

#endif /* HAS_INTERFACE */
#endif
