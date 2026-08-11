// File       : pressureCorrectionEquation.h
// Created    : Thu Mar 14 2024 12:50:04 (+0100)
// Author     : Fabian Wermelinger
// Description: Pressure correction (continuity) equation
// Copyright 2024 CCFNUM HSLU T&A. All Rights Reserved.

#ifndef PRESSURECORRECTIONEQUATION_H
#define PRESSURECORRECTIONEQUATION_H

// code
#include "equation.h"
#include "flowModel.h"
#include "linearSystem.h"
#include "pressureCorrectionAssembler.h"

namespace accel
{

class pressureCorrectionEquation : public equation, public linearSystem<1>
{
private:
    flowModel* model_;

    using Assembler = pressureCorrectionAssembler;

public:
    static constexpr equationID ID = equationID::pressureCorrection;

    pressureCorrectionEquation(realm* realm, flowModel* model);

    void checkDomain(const std::shared_ptr<domain> domain) override;

    bool isConverged() const override;

    void setup() override;

    void initialize() override;

    void postInitialize() override;

    void preSolve() override;

    void solve() override;

    void applyOverrides() override;

    void postSolve() override;

    void preTimeStep() override;

    void printScales() override;

    equationID getID() override
    {
        return ID;
    }

protected:
    void setResidualScales_() override;

    // Pressure-correction-specific override of the generic field corrector.
    template <int BLOCKSIZE,
              int FIELD_DIM = 1,
              int STRIDE = 0,
              int CLIP = 0,
              int OFFSET = 0>
    void correctField_(const domain* domain,
                       ::linearSolver::coefficients<BLOCKSIZE>* coeffs,
                       const stk::mesh::EntityRank entityRank,
                       STKScalarField& stk_dst,
                       const scalar relaxValue = 1.0,
                       const scalar lowerBoundValue = 0,
                       const scalar upperBoundValue = BIG,
                       const scalar clipFactor = 1.0,
                       const scalar offset = 0.0)
    {
        // 1) base correction: under-relaxes and updates the pressure field
        equation::correctField_<BLOCKSIZE, FIELD_DIM, STRIDE, CLIP, OFFSET>(
            domain,
            coeffs,
            entityRank,
            stk_dst,
            relaxValue,
            lowerBoundValue,
            upperBoundValue,
            clipFactor,
            offset);

        // 2) store the full (un-relaxed) pressure correction p'
        pressureCorrection& pCorr = model_->pCorrRef();
        STKScalarField& pCorrSTKFieldRef = pCorr.stkFieldRef();
        const Vector& correction = coeffs->getXVector();

        const mesh& meshObj = domain->meshRef();
        const stk::mesh::MetaData& metaData = meshObj.metaDataRef();
        const stk::mesh::BulkData& bulkData = meshObj.bulkDataRef();

        const stk::mesh::Selector selection =
            metaData.locally_owned_part() &
            stk::mesh::selectUnion(domain->zonePtr()->interiorParts());
        const stk::mesh::BucketVector& buckets =
            bulkData.get_buckets(entityRank, selection);
        for (size_t ib = 0; ib < buckets.size(); ++ib)
        {
            const stk::mesh::Bucket& bucket = *buckets[ib];
            const stk::mesh::Bucket::size_type nEntities = bucket.size();
            scalar* pCorrVal = stk::mesh::field_data(pCorrSTKFieldRef, bucket);
            for (stk::mesh::Bucket::size_type i = 0; i < nEntities; ++i)
            {
                const int64_t row = coeffs->getGraph()->localToRow(
                    bulkData.local_id(bucket[i]));
                if (row < 0) // node not part of this (subset) system
                {
                    pCorrVal[i] = 0.0;
                    continue;
                }
                pCorrVal[i] = correction[row * BLOCKSIZE + STRIDE];
            }
        }

        // 3) refresh ghosts, then build grad(p') (gradient relaxation == 1)
        pCorr.synchronizeGhostedEntities(domain->index());
        pCorr.updateGradientField(domain->index());
    }

    std::unique_ptr<Assembler> assembler_;
};

} /* namespace accel */

#endif // PRESSURECORRECTIONEQUATION_H
