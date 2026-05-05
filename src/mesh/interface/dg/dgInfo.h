// File       : dgInfo.h
// Created    : Fri Aug 25 2023 12:55:24 (+0100)
// Author     : Mhamad Mahdi Alloush
// Description: A storage class for an ip info on a DG interface side. Inherits
//              from ipInfo so that DG and GGI per-IP records live in a unified
//              ipInfoVec_ on the base interfaceSideInfo.
// Copyright 2023 CCFNUM HSLU T&A. All Rights Reserved.

#ifndef DGINFO_H
#define DGINFO_H

#ifdef HAS_INTERFACE

// code
#include "ipInfo.h"

namespace accel
{

class dgInfo : public ipInfo
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

    // DG-specific fields (everything else lives in ipInfo)

    uint64_t localGaussPointId_;

    const scalar bestXRef_;

    scalar bestX_;

    scalar nearestDistance_;

    const scalar nearestDistanceSafety_;

    // True when this gauss point has no opposing face. Mirrors
    // ipInfo::isExposed_ but kept as a separate field so existing DG consumer
    // call sites (which read `dgInfo->gaussPointExposed_`) continue to work
    // until they migrate to the unified `isExposed_`.
    bool gaussPointExposed_;

    // possible reuse
    std::vector<uint64_t> allOpposingFaceIds_;
};

} // namespace accel

#endif /* HAS_INTERFACE */
#endif
