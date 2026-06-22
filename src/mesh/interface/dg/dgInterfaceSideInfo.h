// File       : dgInterfaceSideInfo.h
// Created    : Fri Aug 25 2023 12:55:24 (+0100)
// Author     : Mhamad Mahdi Alloush
// Description: DG (Discontinuous Galerkin) implementation of interfaceSideInfo
// Copyright 2023 CCFNUM HSLU T&A. All Rights Reserved.

#ifndef DGINTERFACESIDEINFO_H
#define DGINTERFACESIDEINFO_H

// code
#include "interfaceSideInfo.h"

namespace accel
{

class dgInfo;

class dgInterfaceSideInfo : public interfaceSideInfo
{
public:
    // constructor and destructor
    dgInterfaceSideInfo(interface* interfPtr,
                        bool isMasterSide,
                        const stk::mesh::PartVector currentPartVec,
                        const stk::mesh::PartVector opposingPartVec,
                        interfaceModelOption option,
                        const scalar expandBoxPercentage,
                        const std::string& searchMethodName,
                        const bool clipIsoParametricCoords,
                        const scalar searchTolerance,
                        const bool dynamicSearchTolAlg,
                        const bool useShifted,
                        const std::string name);

    ~dgInterfaceSideInfo();

    // delete dg info
    void deleteDgInfo();

    // virtual overrides
    void setup() override;
    void initialize() override;
    void update() override;
    void completeSearch() override;
    void determineElemsToGhost() override;
    void provideDiagnosis() override;
    size_t errorCheck() override;

    // DG-specific methods

    void reset();

    void constructDgInfo();

    void constructBoundingPoints();

    void constructBoundingBoxes();

    // For conformal interfaces, determine the opposing gauss point id
    // for each dgInfo by matching opposingIsoParCoords_ to the
    // integration point parametric locations of the opposing face
    void determineOpposingGaussPointIds() override;

    void dumpSearchResults();

    void dumpFaceToFaceResults();

private:
    /* expand search box */
    scalar expandBoxPercentage_;

    /* clip isoparametric coordinates if they are out of bounds */
    const bool clipIsoParametricCoords_;

    /* allow for some finite search tolerance for bounding box */
    const scalar searchTolerance_ = 1e-3;

    /* allow for dynamic search tolerance algorithm where search tolerance is
     * used as point radius from isInElem */
    const bool dynamicSearchTolAlg_ = false;

    bool useShifted_ = true;

    // Methods

    void deleteRangePointsFound_(
        std::vector<boundingSphere>& boundingSphereVec,
        const std::vector<std::pair<theKey, theKey>>& searchKeyPair) const;

    void repeatSearchIfNeeded_(
        const std::vector<boundingSphere>& boundingSphereVec,
        std::vector<std::pair<theKey, theKey>>& searchKeyPair) const;

    void doSearch_();
};

} // namespace accel

#endif
