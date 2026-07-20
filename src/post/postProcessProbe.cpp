// File       : postProcessProbe.cpp
// Created    : Tue Aug 05 2025 19:49:24 (+0100)
// Author     : Mhamad Mahdi Alloush
// Description: Point field probe; interpolates a nodal field to a fixed point
// Copyright 2025 CCFNUM HSLU T&A. All Rights Reserved.

#include "messager.h"
#include "postProcess.h"

namespace accel
{

probeObject::probeObject(postProcess* postProcessPtr,
                         std::string name,
                         postProcessType type,
                         std::vector<std::string> location,
                         label frequency,
                         bool writeToFile,
                         std::array<scalar, SPATIAL_DIM> probeLocation,
                         std::string field)
    : postProcessObject(postProcessPtr,
                        name,
                        type,
                        location,
                        frequency,
                        writeToFile),
      probeLocation_(probeLocation), field_(field)
{
    const stk::mesh::MetaData& metaData =
        postProcessPtr_->meshRef().metaDataRef();

    // validate the field name only; elements are not populated yet, so locating
    // the encapsulating element is deferred to the first update()
    const STKScalarField* STKFieldPtr =
        metaData.get_field<scalar>(stk::topology::NODE_RANK, field_);

    if (!STKFieldPtr)
    {
        std::string msg = field_ + " not valid for post-process. "
                                   "Available fields:\n";
        const auto& fields = metaData.get_fields();
        for (auto f : fields)
        {
            if (f->entity_rank() == stk::topology::NODE_RANK)
            {
                msg += f->name() + "\n";
            }
        }

        errorMsg(msg);
    }

    nComp_ = STKFieldPtr->max_size();

    if (writeToFile_ && messager::master())
    {
        if (fs::is_regular_file(postProcessPtr_->directory() / name_))
        {
            fs::remove_all(postProcessPtr_->directory() / name_);
        }
        std::string fileName((postProcessPtr_->directory() / name_).string());
        std::ofstream file(fileName);

        file << postProcessPtr_->instanceHeader();

        // scalar -> field name; spatial vector -> field_x/_y/_z; else numbered
        if (nComp_ == 1)
        {
            file << "\t" << field_;
        }
        else if (nComp_ == SPATIAL_DIM)
        {
            const char* axis[3] = {"x", "y", "z"};
            for (label i = 0; i < nComp_; i++)
            {
                file << "\t" << field_ << "_" << axis[i];
            }
        }
        else
        {
            for (label i = 0; i < nComp_; i++)
            {
                file << "\t" << field_ << "_" << std::to_string(i);
            }
        }

        file << "\n";
        file.close();
    }
}

void probeObject::update()
{
    const stk::mesh::MetaData& metaData =
        postProcessPtr_->meshRef().metaDataRef();
    const stk::mesh::BulkData& bulkData =
        postProcessPtr_->meshRef().bulkDataRef();

    // lazy one-time location; every rank searches its own elements and the
    // ranks then agree on a single owner by global argmin of the distance
    if (!located_)
    {
        const scalar isoTol = 1.0e-3;

        const STKScalarField* coordField = metaData.get_field<scalar>(
            stk::topology::NODE_RANK,
            postProcessPtr_->meshRef().coordinates_ID);

        const stk::mesh::Selector ownedSel = metaData.locally_owned_part();
        const stk::mesh::BucketVector& elemBuckets =
            bulkData.get_buckets(stk::topology::ELEMENT_RANK, ownedSel);

        scalar bestDist = std::numeric_limits<scalar>::max();
        std::array<scalar, SPATIAL_DIM> bestIso = {0.0};
        stk::mesh::EntityId bestId = stk::mesh::InvalidEntityId;

        for (const stk::mesh::Bucket* bptr : elemBuckets)
        {
            const stk::mesh::Bucket& bkt = *bptr;

            if (bkt.topology().rank() != stk::topology::ELEMENT_RANK)
            {
                continue;
            }

            MasterElement* meSCS =
                MasterElementRepo::get_surface_master_element(bkt.topology());

            if (!meSCS)
            {
                continue;
            }

            const label nodesPerElement = meSCS->nodesPerElement_;

            std::vector<scalar> elemCoords(SPATIAL_DIM * nodesPerElement);
            std::array<scalar, 3> iso = {0.0, 0.0, 0.0};

            for (size_t ie = 0; ie < bkt.size(); ++ie)
            {
                const stk::mesh::Entity elem = bkt[ie];

                stk::mesh::Entity const* elem_nodes =
                    bulkData.begin_nodes(elem);
                const label num_nodes = bulkData.num_nodes(elem);

                for (label ni = 0; ni < num_nodes; ++ni)
                {
                    const scalar* xc =
                        stk::mesh::field_data(*coordField, elem_nodes[ni]);
                    for (label j = 0; j < SPATIAL_DIM; ++j)
                    {
                        elemCoords[j * nodesPerElement + ni] = xc[j];
                    }
                }

                const scalar dist = meSCS->isInElement(
                    &elemCoords[0], &probeLocation_[0], &iso[0]);

                if (dist < bestDist)
                {
                    bestDist = dist;
                    for (label j = 0; j < SPATIAL_DIM; ++j)
                    {
                        bestIso[j] = iso[j];
                    }
                    bestId = bulkData.identifier(elem);
                }
            }
        }

        // a rank is a candidate only if its best element actually contains the
        // point; the parametric distance is unity on the element boundary
        const bool localContains = (bestId != stk::mesh::InvalidEntityId) &&
                                   (bestDist <= 1.0 + isoTol);

        // global argmin over the containing candidates: (distance, rank)
        stk::ParallelMachine comm = bulkData.parallel();
        const int myRank = stk::parallel_machine_rank(comm);

        struct DistRank
        {
            double dist;
            int rank;
        } locMin, glbMin;

        locMin.dist = localContains ? static_cast<double>(bestDist)
                                    : std::numeric_limits<double>::max();
        locMin.rank = myRank;

        MPI_Allreduce(&locMin, &glbMin, 1, MPI_DOUBLE_INT, MPI_MINLOC, comm);

        // no rank contains the point; fatal, decided identically on every rank
        if (glbMin.dist == std::numeric_limits<double>::max())
        {
            std::string ptStr;
            for (label j = 0; j < SPATIAL_DIM; ++j)
            {
                ptStr += (j ? ", " : "") + std::to_string(probeLocation_[j]);
            }

            errorMsg("Probe '" + name_ + "': point [" + ptStr +
                     "] was not found inside any element (outside the mesh "
                     "domain).");
        }

        owner_ = (myRank == glbMin.rank);

        if (owner_)
        {
            encapsulatingElementIdent_ = bestId;
            isoParCoords_ = bestIso;
        }

        located_ = true;
    }

    // the owner interpolates at the cached iso-coords; its contribution is
    // reduced to master, which writes one row, and non-owners send zeros
    std::vector<scalar> l_result(nComp_, 0.0);
    std::vector<scalar> g_result(nComp_, 0.0);

    if (owner_)
    {
        const STKScalarField* fieldPtr =
            metaData.get_field<scalar>(stk::topology::NODE_RANK, field_);

        const stk::mesh::Entity elem = bulkData.get_entity(
            stk::topology::ELEMENT_RANK, encapsulatingElementIdent_);

        const stk::mesh::Bucket& bkt = bulkData.bucket(elem);
        MasterElement* meSCS =
            MasterElementRepo::get_surface_master_element(bkt.topology());

        const label nodesPerElement = meSCS->nodesPerElement_;

        stk::mesh::Entity const* elem_nodes = bulkData.begin_nodes(elem);
        const label num_nodes = bulkData.num_nodes(elem);

        std::vector<scalar> coeff(nComp_ * nodesPerElement, 0.0);
        for (label ni = 0; ni < num_nodes; ++ni)
        {
            const scalar* fv = stk::mesh::field_data(*fieldPtr, elem_nodes[ni]);
            for (label j = 0; j < nComp_; ++j)
            {
                coeff[j * nodesPerElement + ni] = fv[j];
            }
        }

        const int nc = static_cast<int>(nComp_);
        meSCS->interpolatePoint(nc, &isoParCoords_[0], &coeff[0], &l_result[0]);
    }

    stk::all_reduce_sum(
        bulkData.parallel(), &l_result[0], &g_result[0], nComp_);

    if (messager::master() && writeToFile_)
    {
        std::string fileName((postProcessPtr_->directory() / name_).string());
        std::ofstream file(fileName, std::ios::app);

        file << postProcessPtr_->instance();
        for (label j = 0; j < nComp_; j++)
        {
            file << "\t" << g_result[j];
        }
        file << "\n";
    }
}

} // namespace accel
