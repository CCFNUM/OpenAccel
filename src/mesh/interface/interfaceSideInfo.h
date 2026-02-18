// File       : interfaceSideInfo.h
// Created    : Fri Aug 25 2023 12:55:24 (+0100)
// Author     : Mhamad Mahdi Alloush
// Description: Operations at an interface side (single side of a pair)
// Copyright (c) 2023 CCFNUM, Lucerne University of Applied Sciences and Arts.
// SPDX-License-Identifier: BSD-3-Clause

#ifndef INTERFACESIDEINFO_H
#define INTERFACESIDEINFO_H

#ifdef HAS_INTERFACE

// code
#include "dataHandler.h"
#include "vectorUtils.h"

namespace accel
{

class mesh;
class dgInfo;
class interface;
class zone;

typedef stk::search::IdentProc<uint64_t, int> theKey;
typedef stk::search::Point<scalar> Point;
typedef stk::search::Sphere<scalar> Sphere;
typedef stk::search::Box<scalar> Box;
typedef std::pair<Sphere, theKey> boundingSphere;
typedef std::pair<Box, theKey> boundingElementBox;
typedef stk::search::IdentProc<stk::mesh::EntityKey> SearchId;
typedef std::vector<std::pair<Sphere, SearchId>> SphereIdVector;

class interfaceSideInfo
{
public:
    // constructor and destructor
    interfaceSideInfo(interface* interfPtr,
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
                      const std::string debugName);

    ~interfaceSideInfo();

    // delete dg info
    void deleteDgInfo();

    void setup();

    // perform initialization such as dgInfoVec creation and search point/boxes
    void initialize();

    // perform another search due to mesh motion
    void update();

    // reset all search containers
    void reset();

    // perform the 'new'
    void constructDgInfo();

    void constructBoundingPoints();

    void constructBoundingBoxes();

    void determineElemsToGhost();

    void completeSearch();

    // For conformal interfaces, determine the opposing gauss point id
    // for each dgInfo by matching opposingIsoParCoords_ to the
    // integration point parametric locations of the opposing face
    void determineOpposingGaussPointIds();

    void provideDiagnosis();

    void dumpSearchResults();

    void dumpFaceToFaceResults();

    size_t errorCheck();

    void transformCoordinateList(std::vector<scalar>& coordsList,
                                 label npts) const;

    // rotate a vector list where the storage is, i.e. xyzxyzxyz
    template <size_t N>
    void rotateVectorList(std::vector<scalar>& vectorList, label nvecs) const
    {
    }

    // rotate a vector list where the storage is, i.e. xxxyyyzzz
    template <size_t N>
    void rotateVectorListCompact(std::vector<scalar>& vectorList,
                                 label nvecs) const
    {
    }

    // rotate a vector
    template <size_t N>
    void rotateVector(std::vector<scalar>& vector) const
    {
    }

    // reverse rotate a vector list where the storage is, i.e. xyzxyzxyz
    template <size_t N>
    void reverseRotateVectorList(std::vector<scalar>& vectorList,
                                 label nvecs) const
    {
    }

    // reverse rotate a vector list where the storage is, i.e. xxxyyyzzz
    template <size_t N>
    void reverseRotateVectorListCompact(std::vector<scalar>& vectorList,
                                        label nvecs) const
    {
    }

    // reverse rotate a vector
    template <size_t N>
    void reverseRotateVector(std::vector<scalar>& vector) const
    {
    }

    /* vector of DgInfo */
    std::vector<std::vector<dgInfo*>> dgInfoVec_;

    // monarch subject parts; subject part can be subsetted while monarch is
    // not..
    const stk::mesh::PartVector currentPartVec_;

    const stk::mesh::PartVector opposingPartVec_;

    // When current and opposing parts do not fully overlap
    bool hasNonoverlap_ = false;

    // Transformation data for periodicity

    utils::matrix rotationMatrix_;

    utils::vector translationVector_;

    // Access

    std::string name() const
    {
        return name_;
    };

    bool isMasterSide() const
    {
        return isMasterSide_;
    };

    const interface* interfPtr() const;

    const zone* zonePtr() const;

    dataHandler& dataHandlerRef();

    const dataHandler& dataHandlerRef() const;

    // Setters and Getters

    void setDomainType(domainType type)
    {
        parentDomainType_ = type;
    }

    domainType parentDomainType() const
    {
        return parentDomainType_;
    };

private:
    interface* interfPtr_;

    const std::string name_;

    bool isMasterSide_;

    domainType parentDomainType_ = domainType::fluid;

    /* expand search box */
    scalar expandBoxPercentage_;

    const stk::search::SearchMethod searchMethod_;

    /* clip isoparametric coordinates if they are out of bounds */
    const bool clipIsoParametricCoords_;

    /* allow for some finite search tolereance for bounding box */
    const scalar searchTolerance_ = 1e-3;

    /* allow for dynamic search tolerance algorithm where search tolerance is
     * used as point radius from isInElem */
    const bool dynamicSearchTolAlg_ = false;

    bool useShifted_ = true;

    // other possible data handling
    std::unique_ptr<dataHandler> dataHandler_ = nullptr;

    /* bounding box data types for stk_search */
    std::vector<boundingSphere> boundingSphereVec_;

    std::vector<boundingElementBox> boundingFaceElementBoxVec_;

    /* save off product of search */
    std::vector<std::pair<theKey, theKey>> searchKeyPair_;

    interfaceModelOption interfaceModelOption_;

    // Methods

    void deleteRangePointsFound_(
        std::vector<boundingSphere>& boundingSphereVec,
        const std::vector<std::pair<theKey, theKey>>& searchKeyPair) const;

    void repeatSearchIfNeeded_(
        const std::vector<boundingSphere>& boundingSphereVec,
        std::vector<std::pair<theKey, theKey>>& searchKeyPair) const;

    void doSearch_();
};

// Specializations

// rotate a vector list where the storage is, i.e. xyzxyzxyz
template <>
void interfaceSideInfo::rotateVectorList<1>(std::vector<scalar>& vectorList,
                                            label nvecs) const;

// rotate a vector list where the storage is, i.e. xxxyyyzzz
template <>
void interfaceSideInfo::rotateVectorListCompact<1>(
    std::vector<scalar>& vectorList,
    label nvecs) const;

// rotate a vector list where the storage is, i.e. xxxyyyzzz
template <>
void interfaceSideInfo::rotateVector<1>(std::vector<scalar>& vector) const;

// rotate a vector list where the storage is, i.e. xyzxyzxyz
template <>
void interfaceSideInfo::rotateVectorList<SPATIAL_DIM>(
    std::vector<scalar>& vectorList,
    label nvecs) const;

// rotate a vector list where the storage is, i.e. xxxyyyzzz
template <>
void interfaceSideInfo::rotateVectorListCompact<SPATIAL_DIM>(
    std::vector<scalar>& vectorList,
    label nvecs) const;

// reverse rotate a vector list where the storage is, i.e. xyzxyzxyz
template <>
void interfaceSideInfo::reverseRotateVectorList<1>(
    std::vector<scalar>& vectorList,
    label nvecs) const;

// reverse rotate a vector list where the storage is, i.e. xxxyyyzzz
template <>
void interfaceSideInfo::reverseRotateVectorListCompact<1>(
    std::vector<scalar>& vectorList,
    label nvecs) const;

template <>
void interfaceSideInfo::reverseRotateVector<1>(
    std::vector<scalar>& vector) const;

// reverse rotate a vector list where the storage is, i.e. xyzxyzxyz
template <>
void interfaceSideInfo::reverseRotateVectorList<SPATIAL_DIM>(
    std::vector<scalar>& vectorList,
    label nvecs) const;

// reverse rotate a vector list where the storage is, i.e. xxxyyyzzz
template <>
void interfaceSideInfo::reverseRotateVectorListCompact<SPATIAL_DIM>(
    std::vector<scalar>& vectorList,
    label nvecs) const;

template <>
void interfaceSideInfo::reverseRotateVector<SPATIAL_DIM>(
    std::vector<scalar>& vector) const;

} // namespace accel

#endif /* HAS_INTERFACE */
#endif
