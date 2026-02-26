// File       : equation.h
// Created    : Fri Jan 26 2024 09:07:47 (+0100)
// Author     : Fabian Wermelinger
// Description: Abstract base class for a physics equation
// Copyright (c) 2024 CCFNUM, Lucerne University of Applied Sciences and Arts.
// SPDX-License-Identifier: BSD-3-Clause

#ifndef EQUATION_H
#define EQUATION_H

// code
#include "convergenceAcceleration.h"
#include "domain.h"
#include "mesh.h"
#include "types.h"
#include "zone.h"

#include <memory>

namespace accel
{

#define FOREACH_DOMAIN(FCN, ...)                                               \
    do                                                                         \
    {                                                                          \
        for (const auto& domain : this->domainVector_)                         \
        {                                                                      \
            FCN(domain, ##__VA_ARGS__);                                        \
        }                                                                      \
    } while (0)

#define FOREACH_DOMAIN_PTR(FCN, ...)                                           \
    do                                                                         \
    {                                                                          \
        for (const auto& domain : this->domainVector_)                         \
        {                                                                      \
            FCN(domain.get(), ##__VA_ARGS__);                                  \
        }                                                                      \
    } while (0)

#define FOREACH_DOMAIN_IF(FCN, CONDITION, ...)                                 \
    do                                                                         \
    {                                                                          \
        for (const auto& domain : this->domainVector_)                         \
        {                                                                      \
            if (CONDITION)                                                     \
            {                                                                  \
                FCN(domain, ##__VA_ARGS__);                                    \
            }                                                                  \
        }                                                                      \
    } while (0)

#define FOREACH_DOMAIN_PTR_IF(FCN, CONDITION, ...)                             \
    do                                                                         \
    {                                                                          \
        for (const auto& domain : this->domainVector_)                         \
        {                                                                      \
            if (CONDITION)                                                     \
            {                                                                  \
                FCN(domain.get(), ##__VA_ARGS__);                              \
            }                                                                  \
        }                                                                      \
    } while (0)

#define FOREACH_DOMAIN_RAW(BODY)                                               \
    do                                                                         \
    {                                                                          \
        for (const auto& domain : this->domainVector_)                         \
        {                                                                      \
            BODY                                                               \
        }                                                                      \
    } while (0)

class equation
{
public:
    equation(const std::string name, bool sub = false, label subIters = 1)
        : name_(name), sub_(sub), subIters_(subIters), isInitialized_(false)
    {
    }

    equation() = delete;

    virtual ~equation()
    {
    }

    // public API
    const std::string& name() const
    {
        return name_;
    }

    bool isCreated() const
    {
        return isCreated_;
    }

    bool isInitialized() const
    {
        return isInitialized_;
    }

    bool sub() const
    {
        return sub_;
    }

    label subIters() const
    {
        return subIters_;
    }

    virtual void addDomain(std::shared_ptr<domain> domain)
    {
        checkDomain(domain);

        domainVector_.push_back(domain);
    }

    virtual void checkDomain(const std::shared_ptr<domain> domain)
    {
    }

    // pure methods
    virtual bool isConverged() const = 0;

    virtual void setup() = 0;

    // initialization of equation (called once before entering main loop)
    virtual void initialize() = 0;

    virtual void postInitialize() = 0;

    // timestep-specific
    virtual void preTimeStep() = 0;

    // most equations do not require this explicitly
    virtual void postTimeStep()
    {
    }

    // iter-specific
    virtual void preSolve() = 0;

    virtual void solve() = 0;

    // most equations do not require this explicitly
    virtual void postSolve()
    {
    }

    virtual void printScales()
    {
    }

    // operation necessary for matrix manipulation purposes
    virtual stk::mesh::PartVector collectInactiveInteriorParts();

    virtual stk::mesh::PartVector collectInteriorParts();

    virtual stk::mesh::PartVector collectBoundaryParts();

    virtual stk::mesh::PartVector collectStationaryParts();

    // unique equation identifier
    virtual equationID getID()
    {
        return equationID::noID;
    }

protected:
    const std::string name_;

    bool isCreated_;

    bool isInitialized_;

    // is the equation a sub equation (part of a system)?
    bool sub_ = false;

    // sub-iterations within the coefficient loop
    label subIters_ = 1;

    std::vector<std::shared_ptr<domain>> domainVector_;

private:
    std::unique_ptr<convergenceAcceleration> accelerationPtr_;
    std::vector<Vector> accelerationScratch_;
    const Vector* lastAccelCorrectionPtr_ = nullptr;
    label lastAccelIter_ = -1;
    label lastAccelTimeStep_ = -1;
    scalar lastAccelRelaxValue_ = 1.0;
    bool lastAccelUsesScratch_ = false;
    bool accelerationInitialized_ = false;

    void initializeAcceleration_();

    const Vector& applyAcceleration_(const Vector& correction,
                                     const scalar relaxValue,
                                     scalar& outRelaxValue);

protected:
    template <int BLOCKSIZE,
              int FIELD_DIM = 1,
              int STRIDE = 0,
              int CLIP = 0,
              int OFFSET = 0>
    void correctField_(const domain* domain,
                       const Vector& correction,
                       const stk::mesh::EntityRank entityRank,
                       STKScalarField& stk_dst,
                       const scalar relaxValue = 1.0,
                       const scalar lowerBoundValue = 0,
                       const scalar upperBoundValue = BIG,
                       const scalar clipFactor = 1.0,
                       const scalar offset = 0.0)
    {
        static_assert(BLOCKSIZE > 0,
                      "equation::correctField_: BLOCKSIZE must be > 0");
        static_assert(FIELD_DIM > 0,
                      "equation::correctField_: FIELD_DIM must be > 0");
        static_assert(
            FIELD_DIM + STRIDE <= BLOCKSIZE,
            "equation::correctField_: ill-formed FIELD_DIM/STRIDE combination");

        using Bucket = stk::mesh::Bucket;
        using BucketVec = stk::mesh::BucketVector;

        const auto& mesh = domain->meshRef();
        const stk::mesh::MetaData& metaData = mesh.metaDataRef();
        const stk::mesh::BulkData& bulkData = mesh.bulkDataRef();

        const stk::mesh::PartVector& domain_parts =
            domain->zonePtr()->interiorParts();
        const stk::mesh::Selector selection =
            metaData.locally_owned_part() &
            stk::mesh::selectUnion(domain_parts);

        const BucketVec& buckets = bulkData.get_buckets(entityRank, selection);

        scalar effectiveRelaxValue = relaxValue;
        const Vector& effectiveCorrection =
            applyAcceleration_(correction, relaxValue, effectiveRelaxValue);

        for (size_t ib = 0; ib < buckets.size(); ib++)
        {
            const Bucket& bucket = *buckets[ib];
            const Bucket::size_type n_entities = bucket.size();
            scalar* fieldVal = stk::mesh::field_data(stk_dst, bucket);
            for (Bucket::size_type i = 0; i < n_entities; ++i)
            {
                const stk::mesh::Entity entity = bucket[i];
                const auto id = bulkData.local_id(entity);
                for (int k = 0; k < FIELD_DIM; k++)
                {
                    scalar newVal =
                        fieldVal[i * FIELD_DIM + k] +
                        effectiveRelaxValue *
                            (effectiveCorrection[id * BLOCKSIZE + STRIDE + k] +
                             OFFSET * offset);

                    if (CLIP)
                    {
                        if (newVal > upperBoundValue)
                        {
                            newVal = fieldVal[i * FIELD_DIM + k] +
                                     clipFactor * (upperBoundValue -
                                                   fieldVal[i * FIELD_DIM + k]);
                        }
                        else if (newVal < lowerBoundValue)
                        {
                            newVal = fieldVal[i * FIELD_DIM + k] +
                                     clipFactor * (lowerBoundValue -
                                                   fieldVal[i * FIELD_DIM + k]);
                        }
                    }

                    fieldVal[i * FIELD_DIM + k] = newVal;
                }
            }
        }

        // any dependent nodes to be updated?
        applyDependencyUpdates_(domain, entityRank, stk_dst);
    }

    virtual void applyDependencyUpdates_(const domain* domain,
                                         const stk::mesh::EntityRank entityRank,
                                         STKScalarField& stk_dst)
    {
    }
};

} /* namespace accel */

#endif // EQUATION_H
