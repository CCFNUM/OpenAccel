// File : macros.h
// Created : Tue Apr 30 2024 12:55:24 (+0100)
// Author : Mhamad Mahdi Alloush
// Description: Common field operations, IO helpers, and STK mesh utility
// functions
// Copyright (c) 2024 CCFNUM, Lucerne University of Applied Sciences
// and Arts. SPDX-License-Identifier: BSD-3-Clause

#ifndef MACROS_H
#define MACROS_H

// std
#include <iostream>
#include <span>
#include <string>
#include <vector>

// stk
#include <stk_mesh/base/BulkData.hpp>
#include <stk_mesh/base/MetaData.hpp>
#include <stk_util/environment/memory_util.hpp>
#include <stk_util/parallel/CommSparse.hpp>
#include <stk_util/parallel/ParallelReduce.hpp>
#include <stk_util/util/SortAndUnique.hpp>

namespace accel
{

// print out std::vector
template <class T>
std::ostream& operator<<(std::ostream& os, const std::vector<T>& data)
{
    os << data.size() << std::endl << "(\n";
    for (int i = 0; i < data.size(); i++)
    {
        os << data[i] << std::endl;
    }
    os << ")\n";

    return os;
}

// print out std::span
template <class T>
std::ostream& operator<<(std::ostream& os, const std::span<T>& data)
{
    os << data.size() << std::endl << "(\n";
    for (int i = 0; i < data.size(); i++)
    {
        os << data[i] << std::endl;
    }
    os << ")\n";

    return os;
}

// IO

void tolower(std::string& s);

void errorMsg(std::string msg);

void warningMsg(std::string msg);

void infoMsg(std::string msg);

namespace ops
{

// Field

template <class T>
void zero(stk::mesh::Field<T>* fldPtr, stk::mesh::PartVector parts = {})
{
    auto& metaData = fldPtr->mesh_meta_data();
    auto& bulkData = fldPtr->mesh_meta_data().mesh_bulk_data();

    stk::mesh::Selector selAllEntities =
        parts.empty()
            ? metaData.universal_part() & stk::mesh::selectField(*fldPtr)
            : metaData.universal_part() & stk::mesh::selectField(*fldPtr) &
                  stk::mesh::selectUnion(parts);

    stk::mesh::BucketVector const& buckets =
        bulkData.get_buckets(fldPtr->entity_rank(), selAllEntities);

    auto fldSize = fldPtr->max_size();

    for (stk::mesh::BucketVector::const_iterator ib = buckets.begin();
         ib != buckets.end();
         ++ib)
    {
        stk::mesh::Bucket& bucket = **ib;
        const stk::mesh::Bucket::size_type nEntitiesPerBucket = bucket.size();

        for (stk::mesh::Bucket::size_type iEntity = 0;
             iEntity < nEntitiesPerBucket;
             ++iEntity)
        {
            stk::mesh::Entity entity = bucket[iEntity];

            T* value = stk::mesh::field_data(*fldPtr, entity);

            for (auto i = 0; i < fldSize; i++)
            {
                value[i] = static_cast<T>(0);
            }
        }
    }
}

template <class T>
void setValue(stk::mesh::Field<T>* fldPtr,
              T* val,
              stk::mesh::PartVector parts = {})
{
    auto& metaData = fldPtr->mesh_meta_data();
    auto& bulkData = fldPtr->mesh_meta_data().mesh_bulk_data();

    stk::mesh::Selector selAllEntities =
        parts.empty()
            ? metaData.universal_part() & stk::mesh::selectField(*fldPtr)
            : metaData.universal_part() & stk::mesh::selectField(*fldPtr) &
                  stk::mesh::selectUnion(parts);

    stk::mesh::BucketVector const& buckets =
        bulkData.get_buckets(fldPtr->entity_rank(), selAllEntities);

    auto fldSize = fldPtr->max_size();

    for (stk::mesh::BucketVector::const_iterator ib = buckets.begin();
         ib != buckets.end();
         ++ib)
    {
        stk::mesh::Bucket& bucket = **ib;
        const stk::mesh::Bucket::size_type nEntitiesPerBucket = bucket.size();

        for (stk::mesh::Bucket::size_type iEntity = 0;
             iEntity < nEntitiesPerBucket;
             ++iEntity)
        {
            stk::mesh::Entity entity = bucket[iEntity];

            T* value = stk::mesh::field_data(*fldPtr, entity);

            for (auto i = 0; i < fldSize; i++)
            {
                value[i] = val[i];
            }
        }
    }
}

template <class T>
void copy(const stk::mesh::Field<T>* srcFldPtr,
          stk::mesh::Field<T>* dstFldPtr,
          stk::mesh::PartVector parts = {})
{
    auto& metaData = srcFldPtr->mesh_meta_data();
    auto& bulkData = srcFldPtr->mesh_meta_data().mesh_bulk_data();

    stk::mesh::Selector selAllEntities =
        parts.empty()
            ? metaData.universal_part() & stk::mesh::selectField(*srcFldPtr)
            : metaData.universal_part() & stk::mesh::selectField(*srcFldPtr) &
                  stk::mesh::selectUnion(parts);

    stk::mesh::BucketVector const& buckets =
        bulkData.get_buckets(srcFldPtr->entity_rank(), selAllEntities);

    auto srcFldSize = srcFldPtr->max_size();
    auto dstFldSize = dstFldPtr->max_size();

    assert(srcFldSize == dstFldSize);

    for (stk::mesh::BucketVector::const_iterator ib = buckets.begin();
         ib != buckets.end();
         ++ib)
    {
        stk::mesh::Bucket& bucket = **ib;
        const stk::mesh::Bucket::size_type nEntitiesPerBucket = bucket.size();

        for (stk::mesh::Bucket::size_type iEntity = 0;
             iEntity < nEntitiesPerBucket;
             ++iEntity)
        {
            stk::mesh::Entity entity = bucket[iEntity];

            const T* srcValue = stk::mesh::field_data(*srcFldPtr, entity);
            T* dstValue = stk::mesh::field_data(*dstFldPtr, entity);

            assert(srcValue);
            assert(dstValue);

            for (auto i = 0; i < srcFldSize; i++)
            {
                dstValue[i] = srcValue[i];
            }
        }
    }
}

template <class T>
void print(const stk::mesh::Field<T>* fldPtr, stk::mesh::PartVector parts = {})
{
    auto& metaData = fldPtr->mesh_meta_data();
    auto& bulkData = fldPtr->mesh_meta_data().mesh_bulk_data();

    stk::mesh::Selector selAllEntities =
        parts.empty()
            ? metaData.universal_part() & stk::mesh::selectField(*fldPtr)
            : metaData.universal_part() & stk::mesh::selectField(*fldPtr) &
                  stk::mesh::selectUnion(parts);

    stk::mesh::BucketVector const& buckets =
        bulkData.get_buckets(fldPtr->entity_rank(), selAllEntities);

    auto fldSize = fldPtr->max_size();

    if (bulkData.parallel_rank() == 0)
    {
        std::cout << fldPtr->name() << std::endl;
    }

    for (stk::mesh::BucketVector::const_iterator ib = buckets.begin();
         ib != buckets.end();
         ++ib)
    {
        stk::mesh::Bucket& bucket = **ib;
        const stk::mesh::Bucket::size_type nEntitiesPerBucket = bucket.size();

        for (stk::mesh::Bucket::size_type iEntity = 0;
             iEntity < nEntitiesPerBucket;
             ++iEntity)
        {
            stk::mesh::Entity entity = bucket[iEntity];

            const T* value = stk::mesh::field_data(*fldPtr, entity);

            for (auto i = 0; i < fldSize; i++)
            {
                std::cout << value[i] << " ";
            }
            std::cout << std::endl;
        }
    }
}

// ghosting

void populateGhostCommProcs(const stk::mesh::BulkData& bulk_data,
                            stk::mesh::Ghosting& ghosting,
                            std::vector<int>& ghostCommProcs);

// Memory Diagnostics

void printMemoryDiag(const stk::mesh::BulkData& bulk);

void printHwmMemoryDiag(const stk::mesh::BulkData& bulk);

}; // namespace ops

} // namespace accel

#endif // MACROS_H
