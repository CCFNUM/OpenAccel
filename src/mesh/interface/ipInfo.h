// File       : ipInfo.h
// Created    : 2026
// Author     : Mhamad Mahdi Alloush
// Description: Common base class for per-integration-point info shared by DG
//              (Discontinuous Galerkin) and GGI (General Grid Interface)
//              non-conformal interface methods.
// Copyright 2026 CCFNUM HSLU T&A. All Rights Reserved.

#ifndef IPINFO_H
#define IPINFO_H

// code
#include "types.h"

namespace accel
{

class ipInfo
{
public:
    virtual ~ipInfo() = default;

    // Current-side identification
    label parallelRank_ = 0;
    uint64_t globalFaceId_ = 0;
    label currentGaussPointId_ = 0;
    stk::mesh::Entity currentFace_;
    stk::mesh::Entity currentElement_;
    label currentFaceOrdinal_ = 0;
    MasterElement* meFCCurrent_ = nullptr;
    MasterElement* meSCSCurrent_ = nullptr;
    stk::topology currentElementTopo_;

    // IP coordinates and parametric coordinates (sized SPATIAL_DIM by derived
    // ctors)
    std::vector<scalar> currentGaussPointCoords_;
    std::vector<scalar> currentIsoParCoords_;
    std::vector<scalar> opposingIsoParCoords_;

    // Opposing-side identification
    stk::mesh::Entity opposingFace_;
    stk::mesh::Entity opposingElement_;
    label opposingFaceOrdinal_ = 0;
    label opposingGaussPointId_ = -1;
    stk::topology opposingElementTopo_;
    MasterElement* meFCOpposing_ = nullptr;
    MasterElement* meSCSOpposing_ = nullptr;
    label opposingFaceIsGhosted_ = 0;

    // Surface area fraction. GGI fills with the rasterised value (0..1); DG
    // keeps the default of 1.0 so unified assembler loops collapse to the
    // same math when SCS surfaces are fully covered.
    scalar areaFraction_ = 1.0;

    // True when this IP has no opposing connection (non-overlap region).
    // GGI sets this directly; for DG the corresponding flag is mirrored on
    // dgInfo::gaussPointExposed_ until consumer loops are unified.
    bool isExposed_ = false;
};

} // namespace accel

#endif
