// File       : dataTransfer.cpp
// Created    : Fri Nov 21 2025 14:01:11 (+0100)
// Author     : Mhamad Mahdi Alloush
// Description:
// Copyright (c) 2025 CCFNUM, Lucerne University of Applied Sciences and Arts.
// SPDX-License-Identifier: BSD-3-Clause

// code
#ifdef HAS_INTERFACE
#include "dataTransfer.h"

namespace accel
{

// Conformal data transfer

void conformalDataTransfer::setup()
{
}

void conformalDataTransfer::initialize()
{
}

void conformalDataTransfer::update()
{
    // get the vector of matching node pairs
    const std::vector<std::pair<stk::mesh::Entity, stk::mesh::Entity>>&
        matchingNodes = interfacePtr_->matchingNodePairVector();

    // get volume field if required
    const STKScalarField* dualNodalVolumeSTKFieldPtr =
        (type_ == dataTransferType::volumeAverage ||
         type_ == dataTransferType::customized)
            ? meta_.template get_field<scalar>(stk::topology::NODE_RANK,
                                               mesh::dual_nodal_volume_ID)
            : nullptr;

    // get rotation tensor (identity in case of translation periodicity or
    // general connection)
    // the following rotation matrix applies for a quantity on the slave
    // side to be rotated to the master side
    const utils::matrix& rotMat =
        noRotation_ ? utils::matrix::Identity()
                    : interfacePtr_->masterInfoPtr()->rotationMatrix_;

    for (auto fieldPairName : fieldPairNames_)
    {
        stk::mesh::Field<scalar>* fieldPtr1 = meta_.get_field<scalar>(
            stk::topology::NODE_RANK, fieldPairName.first);
        stk::mesh::Field<scalar>* fieldPtr2 = meta_.get_field<scalar>(
            stk::topology::NODE_RANK, fieldPairName.second);

        const label N = fieldPtr1->max_size();

        if (N == 1)
        {
            label iPair = 0;
            for (const auto& nodePair : matchingNodes)
            {
                // retreive values of master and slave nodes
                scalar* value1 =
                    stk::mesh::field_data(*fieldPtr1, nodePair.first);
                scalar* value2 =
                    stk::mesh::field_data(*fieldPtr2, nodePair.second);

                switch (type_)
                {
                    case dataTransferType::copy:
                        {
                            if (reverse_)
                            {
                                // copy value2 to value1
                                value1[0] = value2[0];
                            }
                            else
                            {
                                // copy value1 to value2
                                value2[0] = value1[0];
                            }
                        }
                        break;

                    case dataTransferType::add:
                        {
                            if (reverse_)
                            {
                                // add value2 to value1
                                value1[0] += value2[0];
                            }
                            else
                            {
                                // add value1 to value2
                                value2[0] += value1[0];
                            }
                        }
                        break;

                    case dataTransferType::subtract:
                        {
                            if (reverse_)
                            {
                                // subtract value2 from value1
                                value1[0] -= value2[0];
                            }
                            else
                            {
                                // subtract value1 from value2
                                value2[0] -= value1[0];
                            }
                        }
                        break;

                    case dataTransferType::gather:
                        {
                            // add value2 to value1
                            value1[0] += value2[0];

                            // copy value1 to value2
                            value2[0] = value1[0];
                        }
                        break;

                    case dataTransferType::average:
                        {
                            // store the average in value1
                            value1[0] = (value1[0] + value2[0]) * 0.5;

                            // copy value1 to value2
                            value2[0] = value1[0];
                        }
                        break;

                    case dataTransferType::volumeAverage:
                        {
                            // retreive values of master and slave nodes
                            scalar vol1 = *stk::mesh::field_data(
                                *dualNodalVolumeSTKFieldPtr, nodePair.first);
                            scalar vol2 = *stk::mesh::field_data(
                                *dualNodalVolumeSTKFieldPtr, nodePair.second);

                            // get volume ratio
                            scalar r = vol1 / (vol1 + vol2);

                            // store the average in value1
                            value1[0] = r * value1[0] + (1.0 - r) * value2[0];

                            // copy value1 to value2
                            value2[0] = value1[0];
                        }
                        break;

                    case dataTransferType::move:
                        {
                            if (reverse_)
                            {
                                // move value2 to value1
                                value1[0] += value2[0];

                                // 0 value2
                                value2[0] = 0;
                            }
                            else
                            {
                                // move value1 to value2
                                value2[0] += value1[0];

                                // 0 value1
                                value1[0] = 0;
                            }
                        }
                        break;

                    case dataTransferType::customized:
                        {
                            // retreive values of master and slave nodes
                            scalar vol1 = *stk::mesh::field_data(
                                *dualNodalVolumeSTKFieldPtr, nodePair.first);
                            scalar vol2 = *stk::mesh::field_data(
                                *dualNodalVolumeSTKFieldPtr, nodePair.second);

                            // get volume ratio
                            scalar r = vol1 / (vol1 + vol2);

                            value1[0] /= r;
                            value2[0] /= (1.0 - r);
                        }
                        break;
                }

                // increment
                iPair++;
            }
        }
        else if (N == SPATIAL_DIM)
        {
            label iPair = 0;
            for (const auto& nodePair : matchingNodes)
            {
                // retreive values of master and slave nodes
                scalar* value1 =
                    stk::mesh::field_data(*fieldPtr1, nodePair.first);
                scalar* value2 =
                    stk::mesh::field_data(*fieldPtr2, nodePair.second);

                switch (type_)
                {
                    case dataTransferType::copy:
                        {
                            if (reverse_)
                            {
                                // copy value2 to value1
                                for (label i = 0; i < SPATIAL_DIM; ++i)
                                {
                                    value1[i] = 0;
                                    for (label j = 0; j < SPATIAL_DIM; ++j)
                                    {
                                        value1[i] += rotMat(i, j) * value2[j];
                                    }
                                }
                            }
                            else
                            {
                                // copy value1 to value2
                                for (label i = 0; i < SPATIAL_DIM; ++i)
                                {
                                    value2[i] = 0;
                                    for (label j = 0; j < SPATIAL_DIM; ++j)
                                    {
                                        value2[i] += rotMat(j, i) * value1[j];
                                    }
                                }
                            }
                        }
                        break;

                    case dataTransferType::add:
                        {
                            if (reverse_)
                            {
                                // add value2 to value1
                                for (label i = 0; i < SPATIAL_DIM; ++i)
                                {
                                    for (label j = 0; j < SPATIAL_DIM; ++j)
                                    {
                                        value1[i] += rotMat(i, j) * value2[j];
                                    }
                                }
                            }
                            else
                            {
                                // add value1 to value2
                                for (label i = 0; i < SPATIAL_DIM; ++i)
                                {
                                    for (label j = 0; j < SPATIAL_DIM; ++j)
                                    {
                                        value2[i] += rotMat(j, i) * value1[j];
                                    }
                                }
                            }
                        }
                        break;

                    case dataTransferType::subtract:
                        {
                            if (reverse_)
                            {
                                // subtract value2 from value1
                                for (label i = 0; i < SPATIAL_DIM; ++i)
                                {
                                    for (label j = 0; j < SPATIAL_DIM; ++j)
                                    {
                                        value1[i] -= rotMat(i, j) * value2[j];
                                    }
                                }
                            }
                            else
                            {
                                // subtract value1 from value2
                                for (label i = 0; i < SPATIAL_DIM; ++i)
                                {
                                    for (label j = 0; j < SPATIAL_DIM; ++j)
                                    {
                                        value2[i] -= rotMat(j, i) * value1[j];
                                    }
                                }
                            }
                        }
                        break;

                    case dataTransferType::gather:
                        {
                            // add value2 to value1
                            for (label i = 0; i < SPATIAL_DIM; ++i)
                            {
                                for (label j = 0; j < SPATIAL_DIM; ++j)
                                {
                                    value1[i] += rotMat(i, j) * value2[j];
                                }
                            }

                            // copy value1 to value2
                            for (label i = 0; i < SPATIAL_DIM; ++i)
                            {
                                value2[i] = 0;
                                for (label j = 0; j < SPATIAL_DIM; ++j)
                                {
                                    value2[i] += rotMat(j, i) * value1[j];
                                }
                            }
                        }
                        break;

                    case dataTransferType::average:
                        {
                            // add value2 to value1
                            for (label i = 0; i < SPATIAL_DIM; ++i)
                            {
                                for (label j = 0; j < SPATIAL_DIM; ++j)
                                {
                                    value1[i] += rotMat(i, j) * value2[j];
                                }
                            }

                            // half value1
                            for (label i = 0; i < SPATIAL_DIM; ++i)
                            {
                                value1[i] *= 0.5;
                            }

                            // copy value1 to value2
                            for (label i = 0; i < SPATIAL_DIM; ++i)
                            {
                                value2[i] = 0;
                                for (label j = 0; j < SPATIAL_DIM; ++j)
                                {
                                    value2[i] += rotMat(j, i) * value1[j];
                                }
                            }
                        }
                        break;

                    case dataTransferType::volumeAverage:
                        {
                            // retreive values of master and slave nodes
                            scalar vol1 = *stk::mesh::field_data(
                                *dualNodalVolumeSTKFieldPtr, nodePair.first);
                            scalar vol2 = *stk::mesh::field_data(
                                *dualNodalVolumeSTKFieldPtr, nodePair.second);

                            // get volume ratio
                            scalar r = vol1 / (vol1 + vol2);

                            // multiply value1 by volume ratio
                            for (label i = 0; i < SPATIAL_DIM; ++i)
                            {
                                value1[i] *= r;
                            }

                            // add value2 to value1: v1 = (1 - r) * R v2
                            for (label i = 0; i < SPATIAL_DIM; ++i)
                            {
                                for (label j = 0; j < SPATIAL_DIM; ++j)
                                {
                                    value1[i] +=
                                        (1.0 - r) * rotMat(i, j) * value2[j];
                                }
                            }

                            // copy value1 to value2
                            for (label i = 0; i < SPATIAL_DIM; ++i)
                            {
                                value2[i] = 0;
                                for (label j = 0; j < SPATIAL_DIM; ++j)
                                {
                                    value2[i] += rotMat(j, i) * value1[j];
                                }
                            }
                        }
                        break;

                    case dataTransferType::move:
                        {
                            if (reverse_)
                            {
                                // move value2 to value1
                                for (label i = 0; i < SPATIAL_DIM; ++i)
                                {
                                    for (label j = 0; j < SPATIAL_DIM; ++j)
                                    {
                                        value1[i] += rotMat(i, j) * value2[j];
                                    }
                                }

                                // 0 value2
                                for (label i = 0; i < SPATIAL_DIM; ++i)
                                {
                                    value2[i] = 0;
                                }
                            }
                            else
                            {
                                // move value1 to value2
                                for (label i = 0; i < SPATIAL_DIM; ++i)
                                {
                                    for (label j = 0; j < SPATIAL_DIM; ++j)
                                    {
                                        value2[i] += rotMat(j, i) * value1[j];
                                    }
                                }

                                // 0 value1
                                for (label i = 0; i < SPATIAL_DIM; ++i)
                                {
                                    value1[i] = 0;
                                }
                            }
                        }
                        break;

                    case dataTransferType::customized:
                        {
                            // retreive values of master and slave nodes
                            scalar vol1 = *stk::mesh::field_data(
                                *dualNodalVolumeSTKFieldPtr, nodePair.first);
                            scalar vol2 = *stk::mesh::field_data(
                                *dualNodalVolumeSTKFieldPtr, nodePair.second);

                            // get volume ratio
                            scalar r = vol1 / (vol1 + vol2);

                            for (label i = 0; i < SPATIAL_DIM; ++i)
                            {
                                value1[i] /= r;
                                value2[i] /= (1.0 - r);
                            }
                        }
                        break;
                }

                // increment
                iPair++;
            }
        }
        else if (N == SPATIAL_DIM * SPATIAL_DIM)
        {
            label iPair = 0;
            for (const auto& nodePair : matchingNodes)
            {
                // retreive values of master and slave nodes
                scalar* value1 =
                    stk::mesh::field_data(*fieldPtr1, nodePair.first);
                scalar* value2 =
                    stk::mesh::field_data(*fieldPtr2, nodePair.second);

                switch (type_)
                {
                    case dataTransferType::copy:
                        {
                            if (reverse_)
                            {
                                // copy value2 to value1: value1 = R value2
                                // R^T
                                for (label i = 0; i < SPATIAL_DIM; ++i)
                                {
                                    for (label j = 0; j < SPATIAL_DIM; ++j)
                                    {
                                        scalar sum = 0.0;
                                        for (label k = 0; k < SPATIAL_DIM; ++k)
                                        {
                                            for (label l = 0; l < SPATIAL_DIM;
                                                 ++l)
                                            {
                                                // sum += R(i,k) * T(k,l) *
                                                // R(j,l)
                                                sum += rotMat(i, k) *
                                                       value2[k * SPATIAL_DIM +
                                                              l] *
                                                       rotMat(j, l);
                                            }
                                        }
                                        value1[i * SPATIAL_DIM + j] = sum;
                                    }
                                }
                            }
                            else
                            {
                                // copy value1 to value2: value2 = R^T
                                // value1 R
                                for (label i = 0; i < SPATIAL_DIM; ++i)
                                {
                                    for (label j = 0; j < SPATIAL_DIM; ++j)
                                    {
                                        scalar sum = 0.0;
                                        for (label k = 0; k < SPATIAL_DIM; ++k)
                                        {
                                            for (label l = 0; l < SPATIAL_DIM;
                                                 ++l)
                                            {
                                                // sum += R(k,i) * T(k,l) *
                                                // R(l,j)
                                                sum += rotMat(k, i) *
                                                       value1[k * SPATIAL_DIM +
                                                              l] *
                                                       rotMat(l, j);
                                            }
                                        }
                                        value2[i * SPATIAL_DIM + j] = sum;
                                    }
                                }
                            }
                        }
                        break;

                    case dataTransferType::add:
                        {
                            if (reverse_)
                            {
                                // add value2 to value1: value1 += R value2
                                // R^T
                                for (label i = 0; i < SPATIAL_DIM; ++i)
                                {
                                    for (label j = 0; j < SPATIAL_DIM; ++j)
                                    {
                                        for (label k = 0; k < SPATIAL_DIM; ++k)
                                        {
                                            for (label l = 0; l < SPATIAL_DIM;
                                                 ++l)
                                            {
                                                // sum += R(i,k) * T(k,l) *
                                                // R(j,l)
                                                value1[i * SPATIAL_DIM + j] +=
                                                    rotMat(i, k) *
                                                    value2[k * SPATIAL_DIM +
                                                           l] *
                                                    rotMat(j, l);
                                            }
                                        }
                                    }
                                }
                            }
                            else
                            {
                                // add value1 to value2: value2 += R^T
                                // value1 R
                                for (label i = 0; i < SPATIAL_DIM; ++i)
                                {
                                    for (label j = 0; j < SPATIAL_DIM; ++j)
                                    {
                                        for (label k = 0; k < SPATIAL_DIM; ++k)
                                        {
                                            for (label l = 0; l < SPATIAL_DIM;
                                                 ++l)
                                            {
                                                // sum += R(k,i) * T(k,l) *
                                                // R(l,j)
                                                value2[i * SPATIAL_DIM + j] +=
                                                    rotMat(k, i) *
                                                    value1[k * SPATIAL_DIM +
                                                           l] *
                                                    rotMat(l, j);
                                            }
                                        }
                                    }
                                }
                            }
                        }
                        break;

                    case dataTransferType::subtract:
                        {
                            if (reverse_)
                            {
                                // subtract value2 from value1: value1 -= R
                                // value2 R^T
                                for (label i = 0; i < SPATIAL_DIM; ++i)
                                {
                                    for (label j = 0; j < SPATIAL_DIM; ++j)
                                    {
                                        for (label k = 0; k < SPATIAL_DIM; ++k)
                                        {
                                            for (label l = 0; l < SPATIAL_DIM;
                                                 ++l)
                                            {
                                                // sum += R(i,k) * T(k,l) *
                                                // R(j,l)
                                                value1[i * SPATIAL_DIM + j] -=
                                                    rotMat(i, k) *
                                                    value2[k * SPATIAL_DIM +
                                                           l] *
                                                    rotMat(j, l);
                                            }
                                        }
                                    }
                                }
                            }
                            else
                            {
                                // subtract value1 from value2: value2 -=
                                // R^T value1 R
                                for (label i = 0; i < SPATIAL_DIM; ++i)
                                {
                                    for (label j = 0; j < SPATIAL_DIM; ++j)
                                    {
                                        for (label k = 0; k < SPATIAL_DIM; ++k)
                                        {
                                            for (label l = 0; l < SPATIAL_DIM;
                                                 ++l)
                                            {
                                                // sum += R(k,i) * T(k,l) *
                                                // R(l,j)
                                                value2[i * SPATIAL_DIM + j] -=
                                                    rotMat(k, i) *
                                                    value1[k * SPATIAL_DIM +
                                                           l] *
                                                    rotMat(l, j);
                                            }
                                        }
                                    }
                                }
                            }
                        }
                        break;

                    case dataTransferType::gather:
                        {
                            // add value2 to value1: value1 += R value2 R^T
                            for (label i = 0; i < SPATIAL_DIM; ++i)
                            {
                                for (label j = 0; j < SPATIAL_DIM; ++j)
                                {
                                    for (label k = 0; k < SPATIAL_DIM; ++k)
                                    {
                                        for (label l = 0; l < SPATIAL_DIM; ++l)
                                        {
                                            // sum += R(i,k) * T(k,l) * R(j,l)
                                            value1[i * SPATIAL_DIM + j] +=
                                                rotMat(i, k) *
                                                value2[k * SPATIAL_DIM + l] *
                                                rotMat(j, l);
                                        }
                                    }
                                }
                            }

                            // copy value1 to value2: value2 = R^T value1 R
                            for (label i = 0; i < SPATIAL_DIM; ++i)
                            {
                                for (label j = 0; j < SPATIAL_DIM; ++j)
                                {
                                    scalar sum = 0.0;
                                    for (label k = 0; k < SPATIAL_DIM; ++k)
                                    {
                                        for (label l = 0; l < SPATIAL_DIM; ++l)
                                        {
                                            // sum += R(k,i) * T(k,l) * R(l,j)
                                            sum += rotMat(k, i) *
                                                   value1[k * SPATIAL_DIM + l] *
                                                   rotMat(l, j);
                                        }
                                    }
                                    value2[i * SPATIAL_DIM + j] = sum;
                                }
                            }
                        }
                        break;

                    case dataTransferType::average:
                        {
                            // add value2 to value1: value1 += R value2 R^T
                            for (label i = 0; i < SPATIAL_DIM; ++i)
                            {
                                for (label j = 0; j < SPATIAL_DIM; ++j)
                                {
                                    for (label k = 0; k < SPATIAL_DIM; ++k)
                                    {
                                        for (label l = 0; l < SPATIAL_DIM; ++l)
                                        {
                                            // sum += R(i,k) * T(k,l) * R(j,l)
                                            value1[i * SPATIAL_DIM + j] +=
                                                rotMat(i, k) *
                                                value2[k * SPATIAL_DIM + l] *
                                                rotMat(j, l);
                                        }
                                    }
                                }
                            }

                            // half value1
                            for (label i = 0; i < SPATIAL_DIM; ++i)
                            {
                                for (label j = 0; j < SPATIAL_DIM; ++j)
                                {
                                    value1[i * SPATIAL_DIM + j] *= 0.5;
                                }
                            }

                            // copy value1 to value2: value2 = R^T value1 R
                            for (label i = 0; i < SPATIAL_DIM; ++i)
                            {
                                for (label j = 0; j < SPATIAL_DIM; ++j)
                                {
                                    scalar sum = 0.0;
                                    for (label k = 0; k < SPATIAL_DIM; ++k)
                                    {
                                        for (label l = 0; l < SPATIAL_DIM; ++l)
                                        {
                                            // sum += R(k,i) * T(k,l) * R(l,j)
                                            sum += rotMat(k, i) *
                                                   value1[k * SPATIAL_DIM + l] *
                                                   rotMat(l, j);
                                        }
                                    }
                                    value2[i * SPATIAL_DIM + j] = sum;
                                }
                            }
                        }
                        break;

                    case dataTransferType::volumeAverage:
                        {
                            // retreive values of master and slave nodes
                            scalar vol1 = *stk::mesh::field_data(
                                *dualNodalVolumeSTKFieldPtr, nodePair.first);
                            scalar vol2 = *stk::mesh::field_data(
                                *dualNodalVolumeSTKFieldPtr, nodePair.second);

                            // get volume ratio
                            scalar r = vol1 / (vol1 + vol2);

                            // multiply value1 by volume ratio
                            for (label i = 0; i < SPATIAL_DIM; ++i)
                            {
                                for (label j = 0; j < SPATIAL_DIM; ++j)
                                {
                                    value1[i * SPATIAL_DIM + j] *= r;
                                }
                            }

                            // add value2 to value1: value1 += (1-r) R value2
                            // R^T
                            for (label i = 0; i < SPATIAL_DIM; ++i)
                            {
                                for (label j = 0; j < SPATIAL_DIM; ++j)
                                {
                                    for (label k = 0; k < SPATIAL_DIM; ++k)
                                    {
                                        for (label l = 0; l < SPATIAL_DIM; ++l)
                                        {
                                            // sum += (1-r) * R(i,k) * T(k,l) *
                                            // R(j,l)
                                            value1[i * SPATIAL_DIM + j] +=
                                                (1.0 - r) * rotMat(i, k) *
                                                value2[k * SPATIAL_DIM + l] *
                                                rotMat(j, l);
                                        }
                                    }
                                }
                            }

                            // copy value1 to value2: value2 = R^T value1 R
                            for (label i = 0; i < SPATIAL_DIM; ++i)
                            {
                                for (label j = 0; j < SPATIAL_DIM; ++j)
                                {
                                    scalar sum = 0.0;
                                    for (label k = 0; k < SPATIAL_DIM; ++k)
                                    {
                                        for (label l = 0; l < SPATIAL_DIM; ++l)
                                        {
                                            // sum += R(k,i) * T(k,l) * R(l,j)
                                            sum += rotMat(k, i) *
                                                   value1[k * SPATIAL_DIM + l] *
                                                   rotMat(l, j);
                                        }
                                    }
                                    value2[i * SPATIAL_DIM + j] = sum;
                                }
                            }
                        }
                        break;

                    case dataTransferType::move:
                        {
                            if (reverse_)
                            {
                                // move value2 to value1: value1 += R value2
                                // R^T
                                for (label i = 0; i < SPATIAL_DIM; ++i)
                                {
                                    for (label j = 0; j < SPATIAL_DIM; ++j)
                                    {
                                        for (label k = 0; k < SPATIAL_DIM; ++k)
                                        {
                                            for (label l = 0; l < SPATIAL_DIM;
                                                 ++l)
                                            {
                                                // sum += R(i,k) * T(k,l) *
                                                // R(j,l)
                                                value1[i * SPATIAL_DIM + j] +=
                                                    rotMat(i, k) *
                                                    value2[k * SPATIAL_DIM +
                                                           l] *
                                                    rotMat(j, l);
                                            }
                                        }
                                    }
                                }

                                // 0 value2
                                for (label i = 0; i < SPATIAL_DIM; ++i)
                                {
                                    for (label j = 0; j < SPATIAL_DIM; ++j)
                                    {
                                        value2[i * SPATIAL_DIM + j] = 0;
                                    }
                                }
                            }
                            else
                            {
                                // move value1 to value2: value2 += R^T
                                // value1 R
                                for (label i = 0; i < SPATIAL_DIM; ++i)
                                {
                                    for (label j = 0; j < SPATIAL_DIM; ++j)
                                    {
                                        for (label k = 0; k < SPATIAL_DIM; ++k)
                                        {
                                            for (label l = 0; l < SPATIAL_DIM;
                                                 ++l)
                                            {
                                                // sum += R(k,i) * T(k,l) *
                                                // R(l,j)
                                                value2[i * SPATIAL_DIM + j] +=
                                                    rotMat(k, i) *
                                                    value1[k * SPATIAL_DIM +
                                                           l] *
                                                    rotMat(l, j);
                                            }
                                        }
                                    }
                                }

                                // 0 value1
                                for (label i = 0; i < SPATIAL_DIM; ++i)
                                {
                                    for (label j = 0; j < SPATIAL_DIM; ++j)
                                    {
                                        value1[i * SPATIAL_DIM + j] = 0;
                                    }
                                }
                            }
                        }
                        break;

                    case dataTransferType::customized:
                        {
                            // retreive values of master and slave nodes
                            scalar vol1 = *stk::mesh::field_data(
                                *dualNodalVolumeSTKFieldPtr, nodePair.first);
                            scalar vol2 = *stk::mesh::field_data(
                                *dualNodalVolumeSTKFieldPtr, nodePair.second);

                            // get volume ratio
                            scalar r = vol1 / (vol1 + vol2);

                            for (label i = 0; i < SPATIAL_DIM; ++i)
                            {
                                for (label j = 0; j < SPATIAL_DIM; ++j)
                                {
                                    value1[i * SPATIAL_DIM + j] /= r;
                                    value2[i * SPATIAL_DIM + j] /= (1.0 - r);
                                }
                            }
                        }
                        break;
                }

                // increment
                iPair++;
            }
        }
    }
}

// Non-conformal data transfer

void nonconformalDataTransfer::setup()
{
    std::shared_ptr<fromMesh> fm(
        new fromMesh(meta_,
                     bulk_,
                     meta_.coordinate_field_name(),
                     fieldPairNames_,
                     reverse_ ? interfacePtr_->masterInfoPtr()->opposingPartVec_
                              : interfacePtr_->masterInfoPtr()->currentPartVec_,
                     bulk_.parallel()));

    std::shared_ptr<toMesh> tm(
        new toMesh(meta_,
                   bulk_,
                   meta_.coordinate_field_name(),
                   fieldPairNames_,
                   reverse_ ? interfacePtr_->masterInfoPtr()->currentPartVec_
                            : interfacePtr_->masterInfoPtr()->opposingPartVec_,
                   bulk_.parallel(),
                   searchTolerance_,
                   clipMap_));

    typedef stk::transfer::GeometricTransfer<
        linearInterpolation<fromMesh, toMesh>>
        STKTransfer;

    // extract search type; at present, only one supported
    const stk::search::SearchMethod searchMethod = stk::search::KDTREE;
    if (searchMethodName_ != "stk_kdtree")
    {
        std::cout << "error" << std::endl;
        exit(1);
    }

    STKTransfer_.reset(
        new STKTransfer(fm, tm, name_, searchExpansionFactor_, searchMethod));
}

void nonconformalDataTransfer::initialize()
{
    // do the search
    STKTransfer_->coarse_search();

    // ghost_from_elements
    {
        bulk_.modification_begin();

        typedef stk::transfer::GeometricTransfer<
            linearInterpolation<fromMesh, toMesh>>
            STKTransfer;

        const std::shared_ptr<STKTransfer> transferPtr =
            std::dynamic_pointer_cast<STKTransfer>(STKTransfer_);
        typename STKTransfer::MeshA::EntityProcVec entity_keys;
        transferPtr->determine_entities_to_copy(entity_keys);

        const std::shared_ptr<typename STKTransfer::MeshA> mesha =
            transferPtr->meshA();
        mesha->update_ghosting(entity_keys);

        bulk_.modification_end();
    }

    // local search
    STKTransfer_->local_search();
}

void nonconformalDataTransfer::update()
{
    STKTransfer_->apply();
}

} // namespace accel

#endif /* HAS_INTERFACE */
