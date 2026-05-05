// File       : ggiInfo.cpp
// Created    : Thu Feb 26 2026
// Author     : Mhamad Mahdi Alloush
// Description:
// Copyright 2026 CCFNUM HSLU T&A. All Rights Reserved.

#ifdef HAS_INTERFACE

#include "ggiInfo.h"

namespace accel
{

ggiInfo::ggiInfo(label parallelRank,
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
                 label opposingFaceIsGhosted)
{
    // populate inherited (ipInfo) fields
    parallelRank_ = parallelRank;
    globalFaceId_ = globalFaceId;
    currentGaussPointId_ = currentGaussPointId;
    currentFace_ = currentFace;
    currentElement_ = currentElement;
    currentFaceOrdinal_ = currentFaceOrdinal;
    meFCCurrent_ = meFCCurrent;
    meSCSCurrent_ = meSCSCurrent;
    currentElementTopo_ = currentElementTopo;
    opposingFace_ = opposingFace;
    opposingElement_ = opposingElement;
    opposingFaceOrdinal_ = opposingFaceOrdinal;
    opposingGaussPointId_ = opposingGaussPointId;
    opposingElementTopo_ = opposingElementTopo;
    meFCOpposing_ = meFCOpposing;
    meSCSOpposing_ = meSCSOpposing;
    areaFraction_ = areaFraction;
    opposingFaceIsGhosted_ = opposingFaceIsGhosted;

    currentGaussPointCoords_.resize(SPATIAL_DIM, 0.0);
    currentIsoParCoords_.resize(SPATIAL_DIM, 0.0);
    opposingIsoParCoords_.resize(SPATIAL_DIM, 0.0);
}

} // namespace accel

#endif /* HAS_INTERFACE */
