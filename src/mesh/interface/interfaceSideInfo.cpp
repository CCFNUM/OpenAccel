// File       : interfaceSideInfo.cpp
// Created    : Fri Aug 25 2023
// Author     : Mhamad Mahdi Alloush
// Description: Coordinate transformation and rotation specializations for
//              interfaceSideInfo (shared by DG and GGI).
// Copyright 2023 CCFNUM HSLU T&A. All Rights Reserved.

// code
#include "interfaceSideInfo.h"
#include "interface.h"
#include "ipInfo.h"
#include "mesh.h"

namespace accel
{

interfaceSideInfo::interfaceSideInfo(interface* interfPtr,
                                     bool isMasterSide,
                                     stk::mesh::PartVector currentPartVec,
                                     stk::mesh::PartVector opposingPartVec,
                                     interfaceModelOption option,
                                     std::string name)
    : currentPartVec_(std::move(currentPartVec)),
      opposingPartVec_(std::move(opposingPartVec)), interfPtr_(interfPtr),
      name_(std::move(name)), isMasterSide_(isMasterSide),
      interfaceModelOption_(option), dataHandler_(new dataHandler)
{
}

interfaceSideInfo::~interfaceSideInfo()
{
    clearIpInfo_();
}

const zone* interfaceSideInfo::zonePtr() const
{
    const label zoneIndex = isMasterSide_ ? interfPtr_->masterZoneIndex()
                                          : interfPtr_->slaveZoneIndex();
    return interfPtr_->meshPtr()->zonePtr(zoneIndex);
}

void interfaceSideInfo::clearIpInfo_()
{
    for (auto& faceVec : ipInfoVec_)
    {
        for (ipInfo* p : faceVec)
        {
            delete p;
        }
        faceVec.clear();
    }
    ipInfoVec_.clear();
}

void interfaceSideInfo::transformCoordinateList(std::vector<scalar>& coordsList,
                                                label npts) const
{
    assert(coordsList.size() / SPATIAL_DIM == npts);

    for (label ni = 0; ni < npts; ++ni)
    {
        const utils::vectorViewC coords(&coordsList[ni * SPATIAL_DIM]);
        utils::vector new_coords =
            utils::transformVector(rotationMatrix_, coords);
        new_coords += translationVector_;

        for (label j = 0; j < SPATIAL_DIM; ++j)
        {
            coordsList[ni * SPATIAL_DIM + j] = new_coords(j);
        }
    }
}

// ---------------------------------------------------------------------------
// Forward rotation specializations
// ---------------------------------------------------------------------------

template <>
void interfaceSideInfo::rotateVectorList<1>(std::vector<scalar>& vectorList,
                                            label nvecs) const
{
    // do nothing
}

template <>
void interfaceSideInfo::rotateVectorListCompact<1>(
    std::vector<scalar>& vectorList,
    label nvecs) const
{
    // do nothing
}

template <>
void interfaceSideInfo::rotateVector<1>(std::vector<scalar>& vector) const
{
    // do nothing
}

template <>
void interfaceSideInfo::rotateVectorList<SPATIAL_DIM>(
    std::vector<scalar>& vectorList,
    label nvecs) const
{
    assert(vectorList.size() / SPATIAL_DIM == nvecs);

    for (label ni = 0; ni < nvecs; ++ni)
    {
        utils::vectorView vector(&vectorList[ni * SPATIAL_DIM]);
        utils::transformVectorInPlace(rotationMatrix_, vector);
    }
}

template <>
void interfaceSideInfo::rotateVectorListCompact<SPATIAL_DIM>(
    std::vector<scalar>& vectorList,
    label nvecs) const
{
    assert(vectorList.size() / SPATIAL_DIM == nvecs);

    for (label ni = 0; ni < nvecs; ++ni)
    {
        // clone
        utils::vector newVec;
        for (label j = 0; j < SPATIAL_DIM; ++j)
        {
            label offset = j * nvecs + ni;
            newVec(j) = vectorList[offset];
        }

        // rotate
        utils::transformVectorInPlace(rotationMatrix_, newVec);

        // copy back
        for (label j = 0; j < SPATIAL_DIM; ++j)
        {
            label offset = j * nvecs + ni;
            vectorList[offset] = newVec(j);
        }
    }
}

template <>
void interfaceSideInfo::rotateVector<SPATIAL_DIM>(
    std::vector<scalar>& vector) const
{
    assert(vector.size() == SPATIAL_DIM);
    utils::vectorView vec(vector.data());
    utils::transformVectorInPlace(rotationMatrix_, vec);
}

// ---------------------------------------------------------------------------
// Reverse rotation specializations
// ---------------------------------------------------------------------------

template <>
void interfaceSideInfo::reverseRotateVectorList<1>(
    std::vector<scalar>& vectorList,
    label nvecs) const
{
    // do nothing
}

template <>
void interfaceSideInfo::reverseRotateVectorListCompact<1>(
    std::vector<scalar>& vectorList,
    label nvecs) const
{
    // do nothing
}

template <>
void interfaceSideInfo::reverseRotateVector<1>(
    std::vector<scalar>& vector) const
{
    // do nothing
}

template <>
void interfaceSideInfo::reverseRotateVectorList<SPATIAL_DIM>(
    std::vector<scalar>& vectorList,
    label nvecs) const
{
    assert(vectorList.size() / SPATIAL_DIM == nvecs);

    for (label ni = 0; ni < nvecs; ++ni)
    {
        utils::vectorView vector(&vectorList[ni * SPATIAL_DIM]);
        utils::transformVectorInPlace(rotationMatrix_.transpose(), vector);
    }
}

template <>
void interfaceSideInfo::reverseRotateVectorListCompact<SPATIAL_DIM>(
    std::vector<scalar>& vectorList,
    label nvecs) const
{
    assert(vectorList.size() / SPATIAL_DIM == nvecs);

    for (label ni = 0; ni < nvecs; ++ni)
    {
        // clone
        utils::vector newVec;
        for (label j = 0; j < SPATIAL_DIM; ++j)
        {
            label offset = j * nvecs + ni;
            newVec(j) = vectorList[offset];
        }

        // rotate
        utils::transformVectorInPlace(rotationMatrix_.transpose(), newVec);

        // copy back
        for (label j = 0; j < SPATIAL_DIM; ++j)
        {
            label offset = j * nvecs + ni;
            vectorList[offset] = newVec(j);
        }
    }
}

template <>
void interfaceSideInfo::reverseRotateVector<SPATIAL_DIM>(
    std::vector<scalar>& vector) const
{
    assert(vector.size() == SPATIAL_DIM);
    utils::vectorView vec(vector.data());
    utils::transformVectorInPlace(rotationMatrix_.transpose(), vec);
}

} // namespace accel
