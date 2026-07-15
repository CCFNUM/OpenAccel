// File       : interfaceSideInfo.h
// Created    : Fri Aug 25 2023 12:55:24 (+0100)
// Author     : Mhamad Mahdi Alloush
// Description: Abstract base class for a single interface side
// Copyright 2023 CCFNUM HSLU T&A. All Rights Reserved.

#ifndef INTERFACESIDEINFO_H
#define INTERFACESIDEINFO_H

// code
#include "dataHandler.h"
#include "vectorUtils.h"

namespace accel
{

class interface;
class zone;
class ipInfo;

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
    interfaceSideInfo(interface* interfPtr,
                      bool isMasterSide,
                      stk::mesh::PartVector currentPartVec,
                      stk::mesh::PartVector opposingPartVec,
                      interfaceModelOption option,
                      std::string name);

    virtual ~interfaceSideInfo();

    // Core lifecycle methods (differ per non-conformal method)
    virtual void setup() = 0;
    virtual void initialize() = 0;
    virtual void update() = 0;
    virtual void completeSearch() = 0;
    virtual void determineElemsToGhost() = 0;
    virtual void provideDiagnosis() = 0;
    virtual size_t errorCheck() = 0;

    // Common accessors — concrete because identical for all derived methods
    bool isMasterSide() const
    {
        return isMasterSide_;
    }

    bool hasNonoverlap() const
    {
        return hasNonoverlap_;
    }

    std::string name() const
    {
        return name_;
    }

    const interface* interfPtr() const
    {
        return interfPtr_;
    }

    const zone* zonePtr() const;

    dataHandler& dataHandlerRef()
    {
        return *dataHandler_;
    }

    const dataHandler& dataHandlerRef() const
    {
        return *dataHandler_;
    }

    void setDomainType(domainType type)
    {
        parentDomainType_ = type;
    }

    domainType parentDomainType() const
    {
        return parentDomainType_;
    }

    const stk::mesh::PartVector& currentPartVec() const
    {
        return currentPartVec_;
    }

    const stk::mesh::PartVector& opposingPartVec() const
    {
        return opposingPartVec_;
    }

    // Unified per-IP info storage (concrete ipInfo subtypes).  Outer
    // index = current face index; inner = IP / connection index.  Owned by
    // this base class — populated by derived classes during initialize() /
    // update() and freed by the base destructor (via clearIpInfo_()).
    const std::vector<std::vector<ipInfo*>>& ipInfoVec() const
    {
        return ipInfoVec_;
    }

    // For conformal interfaces only: match opposing gauss point ids.
    // Only DG needs this and overrides it; the default is a no-op.
    virtual void determineOpposingGaussPointIds()
    {
    }

    // Part vectors — populated by the base constructor and accessed
    // directly by the assembly layer (historically as public data members).
    stk::mesh::PartVector currentPartVec_;
    stk::mesh::PartVector opposingPartVec_;

    // Coordinate transformation for periodic interfaces.  Defaults do
    // nothing; periodic DG specialises via dgInterfaceSideInfo (callers
    // needing rotation/translation downcast or know the concrete type).
    // The transformation data below is common to all methods.
    utils::matrix rotationMatrix_;
    utils::vector translationVector_;

    void transformCoordinateList(std::vector<scalar>& coordsList,
                                 label npts) const;

    template <size_t N>
    void rotateVectorList(std::vector<scalar>& vectorList, label nvecs) const
    {
    }

    template <size_t N>
    void rotateVectorListCompact(std::vector<scalar>& vectorList,
                                 label nvecs) const
    {
    }

    template <size_t N>
    void rotateVector(std::vector<scalar>& vector) const
    {
    }

    template <size_t N>
    void reverseRotateVectorList(std::vector<scalar>& vectorList,
                                 label nvecs) const
    {
    }

    template <size_t N>
    void reverseRotateVectorListCompact(std::vector<scalar>& vectorList,
                                        label nvecs) const
    {
    }

    template <size_t N>
    void reverseRotateVector(std::vector<scalar>& vector) const
    {
    }

protected:
    // Frees and clears ipInfoVec_.  Called by ~interfaceSideInfo() and may
    // also be invoked by derived classes from their reset()/update() paths.
    void clearIpInfo_();

    // Polymorphic per-IP info storage.  Populated by derived classes with
    // their concrete ipInfo records.  Ownership lives here.
    std::vector<std::vector<ipInfo*>> ipInfoVec_;

    // Members shared by all derived methods; accessed directly.
    interface* interfPtr_;
    std::string name_;
    bool isMasterSide_;
    bool hasNonoverlap_ = false;
    domainType parentDomainType_ = domainType::fluid;
    interfaceModelOption interfaceModelOption_;

    std::unique_ptr<dataHandler> dataHandler_;
    std::vector<boundingSphere> boundingSphereVec_;
    std::vector<boundingElementBox> boundingFaceElementBoxVec_;
    std::vector<std::pair<theKey, theKey>> searchKeyPair_;
    stk::search::SearchMethod searchMethod_ = stk::search::KDTREE;
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

// rotate a vector
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

#endif
