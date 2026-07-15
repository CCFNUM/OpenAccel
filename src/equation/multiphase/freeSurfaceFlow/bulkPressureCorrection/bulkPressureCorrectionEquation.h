// File       : bulkPressureCorrectionEquation.h
// Created    : Sun Feb 02 2025 20:44:56 (+0100)
// Author     : Mhamad Mahdi Alloush
// Description: Pressure correction equation for bulk free-surface flow
// Copyright 2025 CCFNUM HSLU T&A. All Rights Reserved.

#ifndef BULKPRESSURECORRECTIONEQUATION_H
#define BULKPRESSURECORRECTIONEQUATION_H

// code
#include "bulkPressureCorrectionAssembler.h"
#include "equation.h"
#include "freeSurfaceFlowModel.h"
#include "linearSystem.h"

namespace accel
{

class bulkPressureCorrectionEquation : public equation, public linearSystem<1>
{
private:
    freeSurfaceFlowModel* model_;

    using Assembler = bulkPressureCorrectionAssembler;

public:
    bulkPressureCorrectionEquation(realm* realm, freeSurfaceFlowModel* model);

    void checkDomain(const std::shared_ptr<domain> domain) override;

    bool isConverged() const override;

    void setup() override;

    void initialize() override;

    void postInitialize() override;

    void preSolve() override;

    void solve() override;

    void postSolve() override;

    void preTimeStep() override;

    void printScales() override;

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
        const Vector& correction = coeffs->getXVector();
        pressureCorrection& pCorr = model_->pCorrRef();
        STKScalarField& pCorrSTKFieldRef = pCorr.stkFieldRef();

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

private:
    std::unique_ptr<Assembler> assembler_;
};

} /* namespace accel */

#endif // BULKPRESSURECORRECTIONEQUATION_H
