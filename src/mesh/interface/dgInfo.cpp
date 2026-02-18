// File       : dgInfo.cpp
// Created    : Fri Aug 25 2023 12:55:24 (+0100)
// Author     : Mhamad Mahdi Alloush
// Description:
// Copyright (c) 2023 CCFNUM, Lucerne University of Applied Sciences and Arts.
// SPDX-License-Identifier: BSD-3-Clause

#ifdef HAS_INTERFACE

// #code
#include "dgInfo.h"
#include "master_element/MasterElement.h"

namespace accel
{

dgInfo::dgInfo(label parallelRank,
               uint64_t globalFaceId,
               uint64_t localGaussPointId,
               label currentGaussPointId,
               stk::mesh::Entity currentFace,
               stk::mesh::Entity currentElement,
               const label currentFaceOrdinal,
               MasterElement* meFCCurrent,
               MasterElement* meSCSCurrent,
               stk::topology currentElementTopo,
               const scalar searchTolerance)
    : parallelRank_(parallelRank), globalFaceId_(globalFaceId),
      localGaussPointId_(localGaussPointId),
      currentGaussPointId_(currentGaussPointId), currentFace_(currentFace),
      currentElement_(currentElement), currentFaceOrdinal_(currentFaceOrdinal),
      meFCCurrent_(meFCCurrent), meSCSCurrent_(meSCSCurrent),
      currentElementTopo_(currentElementTopo), bestXRef_(1.0e16),
      bestX_(bestXRef_), nearestDistance_(searchTolerance),
      nearestDistanceSafety_(2.0), gaussPointExposed_(false),
      opposingFaceIsGhosted_(0)
{
    // resize internal vectors
    currentGaussPointCoords_.resize(SPATIAL_DIM);

    // isoPar coords will map to full volume element
    currentIsoParCoords_.resize(SPATIAL_DIM);
    opposingIsoParCoords_.resize(SPATIAL_DIM);
}

dgInfo::~dgInfo()
{
    // nothing to delete
}

void dgInfo::dumpInfo()
{
    std::cout << "------------------------------------------------- "
              << std::endl;
    std::cout << "DGInfo::dump_info() for localGaussPointId_ "
              << localGaussPointId_ << " On Rank " << parallelRank_
              << std::endl;
    std::cout << "parallelRank_ " << parallelRank_ << std::endl;
    std::cout << "globalFaceId_ " << globalFaceId_ << std::endl;
    std::cout << "currentGaussPointId_ " << currentGaussPointId_ << std::endl;
    std::cout << "currentFace_ " << currentFace_ << std::endl;
    std::cout << "currentElement_ " << currentElement_ << std::endl;
    std::cout << "currentElementTopo_ " << currentElementTopo_ << std::endl;
    std::cout << "nDim " << SPATIAL_DIM << std::endl;
    std::cout << "bestXRef_ " << bestXRef_ << std::endl;
    std::cout << "bestX_ " << bestX_ << std::endl;
    std::cout << "nearestDistance_ " << nearestDistance_ << std::endl;
    std::cout << "opposingFaceIsGhosted_ " << opposingFaceIsGhosted_
              << std::endl;
    std::cout << "opposingFace_ " << opposingFace_ << std::endl;
    std::cout << "opposingElement_ " << opposingElement_ << std::endl;
    std::cout << "opposingElementTopo_ " << opposingElementTopo_ << std::endl;
    std::cout << "opposingFaceOrdinal_ " << opposingFaceOrdinal_ << std::endl;
    std::cout << "meFCOpposing_ " << meFCOpposing_ << std::endl;
    std::cout << "meSCSOpposing_ " << meSCSOpposing_ << std::endl;
    std::cout << "currentGaussPointCoords_ " << std::endl;
    for (size_t k = 0; k < currentGaussPointCoords_.size(); ++k)
        std::cout << currentGaussPointCoords_[k] << std::endl;
    std::cout << "currentIsoParCoords_ " << std::endl;
    for (size_t k = 0; k < currentIsoParCoords_.size(); ++k)
        std::cout << currentIsoParCoords_[k] << std::endl;
    std::cout << "opposingIsoParCoords_ " << std::endl;
    for (size_t k = 0; k < opposingIsoParCoords_.size(); ++k)
        std::cout << opposingIsoParCoords_[k] << std::endl;
    std::cout << "allOpposingFaceIds_ " << std::endl;
    for (size_t k = 0; k < allOpposingFaceIds_.size(); ++k)
        std::cout << allOpposingFaceIds_[k] << std::endl;
    std::cout << "gaussPointExposed_ " << gaussPointExposed_ << std::endl;
    std::cout << "------------------------------------------------- "
              << std::endl;
    std::cout << std::endl;
}

} // namespace accel

#endif /* HAS_INTERFACE */
