// File       : ggiInfo.h
// Created    : Thu Feb 26 2026
// Author     : Mhamad Mahdi Alloush
// Description: A storage class for a single CS connection on a GGI interface
//              side: one current-side IP connected to one opposing face with
//              its rasterisation-computed area fraction. Inherits from ipInfo
//              so that DG and GGI per-IP records live in a unified ipInfoVec_
//              on the base interfaceSideInfo. One current IP can produce
//              MULTIPLE ggiInfo objects (one per overlapping opposing face).
// Copyright 2026 CCFNUM HSLU T&A. All Rights Reserved.

#ifndef GGIINFO_H
#define GGIINFO_H

#ifdef HAS_INTERFACE

// code
#include "ipInfo.h"

namespace accel
{

class ggiInfo : public ipInfo
{
public:
    ggiInfo(label parallelRank,
            uint64_t globalFaceId,
            label currentGaussPointId,
            stk::mesh::Entity currentFace,
            stk::mesh::Entity currentElement,
            label currentFaceOrdinal,
            MasterElement* meFCCurrent,
            MasterElement* meSCSCurrent,
            stk::topology currentElementTopo,
            stk::mesh::Entity opposingFace,
            stk::mesh::Entity opposingElement,
            label opposingFaceOrdinal,
            label opposingGaussPointId,
            stk::topology opposingElementTopo,
            MasterElement* meFCOpposing,
            MasterElement* meSCSOpposing,
            scalar areaFraction,
            label opposingFaceIsGhosted);

    ~ggiInfo() = default;
};

} // namespace accel

#endif /* HAS_INTERFACE */
#endif
