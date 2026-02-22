// File       : dataTransfer.h
// Created    : Fri Nov 21 2025 14:01:11 (+0100)
// Author     : Mhamad Mahdi Alloush
// Description:
// Copyright (c) 2025 CCFNUM, Lucerne University of Applied Sciences and Arts.
// SPDX-License-Identifier: BSD-3-Clause

#ifndef DATATRANSFER_H
#define DATATRANSFER_H

#ifdef HAS_INTERFACE

// code
#include "interface.h"

namespace accel
{

enum class dataTransferType
{
    copy,
    add,
    subtract,
    gather,
    average,
    volumeAverage,
    move,
    customized
};

class dataTransfer
{
protected:
    interface* interfacePtr_ = nullptr;

    std::string name_;

    stk::mesh::BulkData& bulk_;

    stk::mesh::MetaData& meta_;

    std::vector<std::pair<std::string, std::string>> fieldPairNames_;

    dataTransferType type_ = dataTransferType::copy;

    bool reverse_ = false;

    std::map<std::string, std::pair<scalar, scalar>> clipMap_;

public:
    dataTransfer(interface* interfacePtr,
                 std::string name,
                 std::vector<std::pair<std::string, std::string>> fields,
                 dataTransferType type,
                 bool reverse = false,
                 std::vector<std::pair<scalar, scalar>> minMaxClipVector = {})
        : interfacePtr_(interfacePtr), name_(name),
          bulk_(interfacePtr->meshPtr()->bulkDataRef()),
          meta_(interfacePtr->meshPtr()->metaDataRef()),
          fieldPairNames_(fields), type_(type), reverse_(reverse)
    {
        // default clipping vector
        if (minMaxClipVector.empty())
        {
            minMaxClipVector.resize(fields.size());
            for (label iField = 0; iField < fields.size(); iField++)
            {
                minMaxClipVector[iField] = std::pair<scalar, scalar>(-BIG, BIG);
            }
        }
#ifndef NDEBUG
        // diagnostics
        for (auto fieldPairName : fieldPairNames_)
        {
            stk::mesh::Field<scalar>* fieldPtr1 = meta_.get_field<scalar>(
                stk::topology::NODE_RANK, fieldPairName.first);
            stk::mesh::Field<scalar>* fieldPtr2 = meta_.get_field<scalar>(
                stk::topology::NODE_RANK, fieldPairName.second);

            if (fieldPtr1 == nullptr)
            {
                errorMsg("Field 1 passed to interface data transfer not found");
            }

            if (fieldPtr2 == nullptr)
            {
                errorMsg("Field 2 passed to interface data transfer not found");
            }

            if (fieldPtr1->entity_rank() != stk::topology::NODE_RANK)
            {
                errorMsg("Field 1 passed to interface data transfer is not "
                         "defined over nodes");
            }

            if (fieldPtr2->entity_rank() != stk::topology::NODE_RANK)
            {
                errorMsg("Field 2 passed to interface data transfer is not "
                         "defined over nodes");
            }

            if (fieldPtr1->max_size() != fieldPtr2->max_size())
            {
                errorMsg("transferring data between fields of different size "
                         "is not allowed");
            }
        }

        if (fieldPairNames_.size() != minMaxClipVector.size())
        {
            errorMsg("Inconsistency in fields and bounds");
        }
#endif /* NDEBUG */

        // populate clip map
        for (label i = 0; i < fieldPairNames_.size(); i++)
        {
            clipMap_.insert({fieldPairNames_[i].first, minMaxClipVector[i]});
        }
    }

    void setDataTransferType(dataTransferType type)
    {
        type_ = type;
    }

    void setReversedDataTransfer(bool reverse)
    {
        reverse_ = reverse;
    }

    virtual void setup() = 0;

    virtual void initialize() = 0;

    virtual void update() = 0;
};

class conformalDataTransfer : public dataTransfer
{
protected:
    // force no rotation
    bool noRotation_ = false;

public:
    conformalDataTransfer(
        interface* interfacePtr,
        std::string name,
        std::vector<std::pair<std::string, std::string>> fields,
        dataTransferType type,
        bool reverse = false,
        std::vector<std::pair<scalar, scalar>> minMaxClipVector = {},
        bool noRotation = false)
        : dataTransfer(interfacePtr,
                       name,
                       fields,
                       type,
                       reverse,
                       minMaxClipVector),
          noRotation_(noRotation)
    {
    }

    void setup() override;

    void initialize() override;

    void update() override;
};

class nonconformalDataTransfer : public dataTransfer
{
private:
    std::shared_ptr<stk::transfer::TransferBase> STKTransfer_;

    scalar searchTolerance_ = 1.0e-4;

    scalar searchExpansionFactor_ = 1.5;

    std::string searchMethodName_ = "stk_kdtree";

    // private types

    class fromMesh
    {
    public:
        typedef stk::mesh::Entity Entity;
        typedef std::vector<Entity> EntityVec;
        typedef stk::mesh::EntityKey EntityKey;
        typedef std::set<EntityKey> EntityKeySet;
        typedef stk::search::IdentProc<EntityKey, unsigned> EntityProc;
        typedef std::vector<EntityProc> EntityProcVec;

        typedef stk::search::Point<scalar> Point;
        typedef stk::search::Box<scalar> Box;
        typedef std::pair<Box, EntityProc> BoundingBox;

        enum
        {
            Dimension = 3
        };

        typedef std::vector<std::pair<std::string, std::string>> PairNames;

        std::vector<const stk::mesh::FieldBase*>
        getFields(const stk::mesh::MetaData& fromMetaData,
                  const PairNames& VarPairName)
        {
            // will want to check that all is well with field registration
            bool allFieldsAreFine = true;
            std::vector<const stk::mesh::FieldBase*> fromFieldVec;
            // provide field names
            for (PairNames::const_iterator i = VarPairName.begin();
                 i != VarPairName.end();
                 ++i)
            {
                const std::string& fieldName = i->first;
                const stk::mesh::FieldBase* fromfield =
                    stk::mesh::get_field_by_name(fieldName, fromMetaData);
                if (NULL == fromfield)
                {
                    allFieldsAreFine = false;
                }
                else
                {
                    // always push back; check for errors below
                    fromFieldVec.push_back(fromfield);

                    // check that the field is defined on **all** parts
                    stk::mesh::Selector fieldSelector =
                        stk::mesh::selectField(*fromfield);
                    for (size_t k = 0; k < fromPartVec_.size(); ++k)
                    {
                        stk::mesh::BucketVector const& partBuckets =
                            fromBulkData_.get_buckets(
                                stk::topology::NODE_RANK,
                                stk::mesh::Selector(*fromPartVec_[k]));

                        bool fieldIsFine = true;
                        for (stk::mesh::BucketVector::const_iterator ib =
                                 partBuckets.begin();
                             ib != partBuckets.end();
                             ++ib)
                        {
                            stk::mesh::Bucket& b = **ib;
                            fieldIsFine &= fieldSelector(b);
                        }

                        // local check to make sure that the field is somewhere
                        // (delay the throw)
                        if (!fieldIsFine)
                        {
                            allFieldsAreFine = false;
                        }
                    }
                }
            }

            // final error check; only return when all is well
            if (allFieldsAreFine)
            {
                return fromFieldVec;
            }
            else
            {
                throw std::runtime_error("Xfer::FromMesh:Error field "
                                         "registration on desired parts of "
                                         "the mesh is not complete");
            }
        }

        fromMesh(const stk::mesh::MetaData& fromMetaData,
                 stk::mesh::BulkData& fromBulkData,
                 const std::string& coordinates_name,
                 const PairNames& VarPairName,
                 const stk::mesh::PartVector& fromPartVec,
                 const stk::ParallelMachine comm)
            : fromMetaData_(fromMetaData), fromBulkData_(fromBulkData),
              fromcoordinates_(
                  fromMetaData.get_field<scalar>(stk::topology::NODE_RANK,
                                                 coordinates_name)),
              fromPartVec_(fromPartVec),
              fromFieldVec_(getFields(fromMetaData, VarPairName)), comm_(comm),
              mesh_modified_(false), ghosting_(0), ghosting_map_()
        {
            // nothing to do
        }

        ~fromMesh()
        {
        }

        struct BoundingBoxCompare
        {
            bool operator()(const BoundingBox& a, const BoundingBox& b) const
            {
                return a.second.id() < b.second.id();
            }
        };

        // Needed for STK Transfer
        stk::ParallelMachine comm() const
        {
            return comm_;
        }

        void bounding_boxes(std::vector<BoundingBox>& v) const
        {
            Point min_corner, max_corner;

            stk::mesh::Selector s_locally_owned_union =
                fromMetaData_.locally_owned_part() &
                stk::mesh::selectUnion(fromPartVec_);

            // determine entity rank for the part served up; should be
            // homogeneous
            stk::mesh::EntityRank partEntityRank =
                fromPartVec_[0]->primary_entity_rank();

            stk::mesh::BucketVector const& entity_buckets =
                fromBulkData_.get_buckets(partEntityRank,
                                          s_locally_owned_union);
            for (stk::mesh::BucketVector::const_iterator ib =
                     entity_buckets.begin();
                 ib != entity_buckets.end();
                 ++ib)
            {
                stk::mesh::Bucket& b = **ib;

                const stk::mesh::Bucket::size_type length = b.size();

                for (stk::mesh::Bucket::size_type k = 0; k < length; ++k)
                {
                    // get entity
                    stk::mesh::Entity theEntity = b[k];

                    // initialize max and min
                    for (unsigned j = 0; j < SPATIAL_DIM; ++j)
                    {
                        min_corner[j] = +1.0e16;
                        max_corner[j] = -1.0e16;
                    }

                    stk::mesh::Entity const* entity_node_rels =
                        fromBulkData_.begin_nodes(theEntity);
                    label num_entity_nodes = fromBulkData_.num_nodes(theEntity);
                    for (label ni = 0; ni < num_entity_nodes; ++ni)
                    {
                        stk::mesh::Entity node = entity_node_rels[ni];

                        scalar* coords =
                            stk::mesh::field_data(*fromcoordinates_, node);
                        for (unsigned j = 0; j < SPATIAL_DIM; ++j)
                        {
                            min_corner[j] = std::min(min_corner[j], coords[j]);
                            max_corner[j] = std::max(max_corner[j], coords[j]);
                        }
                    }

                    // setup ident
                    fromMesh::EntityProc theIdent(
                        fromBulkData_.entity_key(theEntity),
                        fromBulkData_.parallel_rank());

                    v.push_back(
                        BoundingBox(Box(min_corner, max_corner), theIdent));
                }
            }
            std::sort(v.begin(), v.end(), BoundingBoxCompare());
        }

        void update_ghosting(const EntityProcVec& entity_keys)
        {
            ghosting_map_.resize(entity_keys.size());
            for (size_t i = 0; i < entity_keys.size(); ++i)
            {
                // convert from EntityProc based on EntityKey to EntityProc
                // based on raw Entity.
                const EntityProc& key_proc = entity_keys[i];
                const stk::mesh::EntityKey key = key_proc.id();
                const unsigned proc = key_proc.proc();
                const stk::mesh::Entity e = entity(key);
                const stk::mesh::EntityProc ep(e, proc);
                ghosting_map_[i] = ep;
            }

            unsigned s = !ghosting_map_.empty();
            stk::all_reduce(comm_, stk::ReduceSum<1>(&s));

            if (s)
            {
                std::sort(ghosting_map_.begin(), ghosting_map_.end());
                stk::mesh::EntityProcVec::iterator del =
                    std::unique(ghosting_map_.begin(), ghosting_map_.end());
                ghosting_map_.resize(std::distance(ghosting_map_.begin(), del));

                std::string theGhostName = "accel_transfer_ghosting";
                for (unsigned i = 0; i != fromFieldVec_.size(); ++i)
                    theGhostName += "_" + fromFieldVec_[i]->name();
                ghosting_ = &fromBulkData_.create_ghosting(theGhostName);
                fromBulkData_.change_ghosting(*ghosting_, ghosting_map_);
                mesh_modified_ = true;
            }
        }

        void update_values()
        {
            if (ghosting_)
            {
                std::vector<const stk::mesh::FieldBase*> fields(
                    fromFieldVec_.begin(), fromFieldVec_.end());
                if (mesh_modified_)
                {
                    // Copy coordinates to the newly ghosted nodes
                    mesh_modified_ = false;
                    fields.push_back(fromcoordinates_);
                }
                stk::mesh::communicate_field_data(*ghosting_, fields);
                stk::mesh::copy_owned_to_shared(fromBulkData_, fields);
            }
        }

        Entity entity(const EntityKey k) const
        {
            return fromBulkData_.get_entity(k);
        }

        const stk::mesh::MetaData& fromMetaData_;
        stk::mesh::BulkData& fromBulkData_;
        const stk::mesh::Field<scalar>* fromcoordinates_;
        const stk::mesh::PartVector fromPartVec_;
        const std::vector<const stk::mesh::FieldBase*> fromFieldVec_;
        const stk::ParallelMachine comm_;

        bool mesh_modified_;
        stk::mesh::Ghosting* ghosting_;
        stk::mesh::EntityProcVec ghosting_map_;
    };

    class toMesh
    {
    public:
        typedef stk::mesh::Entity Entity;
        typedef std::vector<Entity> EntityVec;
        typedef stk::mesh::EntityKey EntityKey;
        typedef std::set<EntityKey> EntityKeySet;
        typedef stk::search::IdentProc<EntityKey, unsigned> EntityProc;
        typedef std::vector<EntityProc> EntityProcVec;

        typedef stk::search::Point<scalar> Point;
        typedef stk::search::Sphere<scalar> Sphere;
        typedef std::pair<Sphere, EntityProc> BoundingBox;

        enum
        {
            Dimension = 3
        };

        typedef std::vector<std::pair<std::string, std::string>> PairNames;

        std::vector<const stk::mesh::FieldBase*>
        getFields(const stk::mesh::MetaData& toMetaData,
                  const PairNames& VarPairName)
        {
            // will want to check that all is well with field registration
            bool allFieldsAreFine = true;
            std::vector<const stk::mesh::FieldBase*> toFieldVec;
            // provide field names
            for (PairNames::const_iterator i = VarPairName.begin();
                 i != VarPairName.end();
                 ++i)
            {
                const std::string& fieldName = i->second;
                const stk::mesh::FieldBase* tofield =
                    stk::mesh::get_field_by_name(fieldName, toMetaData);
                if (NULL == tofield)
                {
                    allFieldsAreFine = false;
                }
                else
                {
                    // always push back; check for errors below
                    toFieldVec.push_back(tofield);

                    // check that the field is defined on **all** parts
                    stk::mesh::Selector fieldSelector =
                        stk::mesh::selectField(*tofield);
                    for (size_t k = 0; k < toPartVec_.size(); ++k)
                    {
                        const stk::mesh::BucketVector& partBuckets =
                            toBulkData_.get_buckets(
                                stk::topology::NODE_RANK,
                                stk::mesh::Selector(*toPartVec_[k]));

                        bool fieldIsFine = true;
                        for (stk::mesh::BucketVector::const_iterator ib =
                                 partBuckets.begin();
                             ib != partBuckets.end();
                             ++ib)
                        {
                            stk::mesh::Bucket& b = **ib;
                            fieldIsFine &= fieldSelector(b);
                        }

                        // local check to make sure that the field is somewhere
                        // (delay the throw)
                        if (!fieldIsFine)
                        {
                            allFieldsAreFine = false;
                        }
                    }
                }
            }

            // final error check; only return when all is well
            if (allFieldsAreFine)
            {
                return toFieldVec;
            }
            else
            {
                throw std::runtime_error("Xfer::ToMesh:Error field "
                                         "registration on desired parts of the "
                                         "mesh is not complete");
            }
        }

        toMesh(stk::mesh::MetaData& toMetaData,
               stk::mesh::BulkData& toBulkData,
               const std::string& coordinates_name,
               const PairNames& VarPairName,
               const stk::mesh::PartVector& toPartVec,
               const stk::ParallelMachine comm,
               const scalar radius,
               const std::map<std::string, std::pair<scalar, scalar>> clipMap)
            : toMetaData_(toMetaData), toBulkData_(toBulkData),
              tocoordinates_(
                  toMetaData.get_field<scalar>(stk::topology::NODE_RANK,
                                               coordinates_name)),
              toPartVec_(toPartVec),
              toFieldVec_(getFields(toMetaData, VarPairName)), comm_(comm),
              radius_(radius), clipMap_(clipMap)
        {
            // nothing to do
        }

        ~toMesh(){};

        struct BoundingBoxCompare
        {
            bool operator()(const BoundingBox& a, const BoundingBox& b) const
            {
                return a.second.id() < b.second.id();
            }
        };

        // Needed for STK Transfer
        stk::ParallelMachine comm() const
        {
            return comm_;
        }

        void bounding_boxes(std::vector<BoundingBox>& v) const
        {
            const unsigned spatial_dimension = toMetaData_.spatial_dimension();

            Point center;

            stk::mesh::Selector s_locally_owned_union =
                toMetaData_.locally_owned_part() &
                stk::mesh::selectUnion(toPartVec_);

            stk::mesh::BucketVector const& node_buckets =
                toBulkData_.get_buckets(stk::topology::NODE_RANK,
                                        s_locally_owned_union);
            for (stk::mesh::BucketVector::const_iterator ib =
                     node_buckets.begin();
                 ib != node_buckets.end();
                 ++ib)
            {
                stk::mesh::Bucket& b = **ib;

                const stk::mesh::Bucket::size_type length = b.size();

                for (stk::mesh::Bucket::size_type k = 0; k < length; ++k)
                {
                    // get node
                    stk::mesh::Entity node = b[k];

                    scalar* coord =
                        stk::mesh::field_data(*tocoordinates_, node);
                    for (unsigned i = 0; i < spatial_dimension; ++i)
                    {
                        center[i] = coord[i];
                    }

                    // setup ident
                    toMesh::EntityProc theIdent(toBulkData_.entity_key(node),
                                                toBulkData_.parallel_rank());

                    toMesh::BoundingBox theBox(Sphere(center, radius_),
                                               theIdent);
                    v.push_back(theBox);
                }
            }
            std::sort(v.begin(), v.end(), BoundingBoxCompare());
        }

        void update_values()
        {
            std::vector<const stk::mesh::FieldBase*> fields(toFieldVec_.begin(),
                                                            toFieldVec_.end());
            stk::mesh::copy_owned_to_shared(toBulkData_, fields);
        }

        stk::mesh::MetaData& toMetaData_;
        stk::mesh::BulkData& toBulkData_;
        const stk::mesh::Field<scalar>* tocoordinates_;
        const stk::mesh::PartVector toPartVec_;
        const std::vector<const stk::mesh::FieldBase*> toFieldVec_;
        const stk::ParallelMachine comm_;
        const scalar radius_;
        std::map<std::string, std::pair<scalar, scalar>> clipMap_;
        typedef std::map<stk::mesh::EntityKey, std::vector<scalar>>
            TransferInfo;
        TransferInfo TransferInfo_;
    };

    template <class FROM, class TO>
    class linearInterpolation
    {
    public:
        typedef FROM MeshA;
        typedef TO MeshB;
        typedef typename MeshA::EntityKey EntityKeyA;
        typedef typename MeshB::EntityKey EntityKeyB;
        typedef typename MeshA::EntityProc EntityProcA;
        typedef typename MeshB::EntityProc EntityProcB;

        typedef std::pair<EntityProcB, EntityProcA> EntityProcRelation;
        typedef std::vector<EntityProcRelation> EntityProcRelationVec;

        typedef std::multimap<EntityKeyB, EntityKeyA> EntityKeyMap;

        enum
        {
            Dimension = MeshA::Dimension
        };

        static void filter_to_nearest(EntityKeyMap& rangeToDomain,
                                      const MeshA& fromElem,
                                      MeshB& toPoints);

        static void apply(MeshB& toPoints,
                          const MeshA& fromElem,
                          const EntityKeyMap& rangeToDomain);
    };

public:
    nonconformalDataTransfer(
        interface* interfacePtr,
        std::string name,
        std::vector<std::pair<std::string, std::string>> fields,
        dataTransferType type,
        bool reverse,
        std::vector<std::pair<scalar, scalar>> minMaxClipVector,
        scalar searchTolerance,
        scalar searchExpansionFactor)
        : dataTransfer(interfacePtr,
                       name,
                       fields,
                       type,
                       reverse,
                       minMaxClipVector),
          searchTolerance_(searchTolerance),
          searchExpansionFactor_(searchExpansionFactor)
    {
        // only applies for a copy interface data transfer type
        assert(type_ == dataTransferType::copy);
    }

    void setup() override;

    void initialize() override;

    void update() override;
};

} // namespace accel

#include "dataTransfer.hpp"

#endif /* HAS_INTERFACE */
#endif
