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

    // Core lifecycle methods (differ between DG and GGI)
    virtual void setup() = 0;
    virtual void initialize() = 0;
    virtual void update() = 0;
    virtual void completeSearch() = 0;
    virtual void determineElemsToGhost() = 0;
    virtual void provideDiagnosis() = 0;
    virtual size_t errorCheck() = 0;

    // Common accessors — concrete (non-virtual) implementations because both
    // DG and GGI use identical bodies that read base-class state.
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

    // Unified per-IP info storage (DG: dgInfo*, GGI: ggiInfo*).  Outer
    // index = current face index; inner = IP / connection index.  Owned by
    // this base class — populated by derived classes during initialize() /
    // update() and freed by the base destructor (via clearIpInfo_()).
    const std::vector<std::vector<ipInfo*>>& ipInfoVec() const
    {
        return ipInfoVec_;
    }

    // For conformal interfaces only: match opposing gauss point ids.
    // GGI does not need this; DG overrides with its full implementation.
    virtual void determineOpposingGaussPointIds()
    {
    }

    // Part vectors — populated by the base constructor and accessed
    // directly by the assembly layer (historically as public data members).
    stk::mesh::PartVector currentPartVec_;
    stk::mesh::PartVector opposingPartVec_;

    // Coordinate transformation for periodic interfaces.
    // Default implementations do nothing (GGI with generalConnection),
    // specialisations are provided for periodic DG cases via the
    // concrete dgInterfaceSideInfo class.  External callers that use
    // rotation/translation must downcast to dgInterfaceSideInfo when
    // GGI is not relevant, or rely on the concrete type known at the
    // call site.  The transformation data below is common to both
    // methods so it lives here.
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
    // either dgInfo* (from dgInterfaceSideInfo) or ggiInfo* (from
    // ggiInterfaceSideInfo).  Ownership lives here.
    std::vector<std::vector<ipInfo*>> ipInfoVec_;

    // Members shared between DG and GGI (formerly duplicated in each
    // child).  Derived classes access these directly.
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
