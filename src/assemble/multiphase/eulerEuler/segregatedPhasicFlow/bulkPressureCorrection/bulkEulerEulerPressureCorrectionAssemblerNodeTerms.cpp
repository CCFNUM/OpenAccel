// File       : bulkEulerEulerPressureCorrectionAssemblerNodeTerms.cpp
// Description: Nodal terms for common Euler-Euler pressure correction

#include "bulkEulerEulerPressureCorrectionAssembler.h"

namespace accel
{

void bulkEulerEulerPressureCorrectionAssembler::assembleNodeTermsFused_(
    const domain* domain,
    Context* ctx)
{
    if (!model_->controlsRef().isTransient())
    {
        return;
    }

    const auto& mesh = model_->meshRef();
    const auto& metaData = mesh.metaDataRef();
    const auto& bulkData = mesh.bulkDataRef();
    const auto& volumeField = *metaData.get_field<scalar>(
        stk::topology::NODE_RANK, mesh::dual_nodal_volume_ID);
    const scalar dt = model_->controlsRef().getTimestep();
    const auto scheme = model_->controlsRef()
                            .solverRef()
                            .solverControl_.basicSettings_.transientScheme_;
    std::array<scalar, 3> coefficients = {1.0, -1.0, 0.0};
    label numberOfStates = 2;
    if (scheme == transientSchemeType::secondOrderBackwardEuler)
    {
        coefficients =
            BDF2::coeff(dt, model_->controlsRef().getTimestep(-1));
        numberOfStates = 3;
    }
    else
    {
        const auto firstOrder = BDF1::coeff();
        coefficients[0] = firstOrder[0];
        coefficients[1] = firstOrder[1];
    }

    Matrix& matrix = ctx->getAMatrix();
    Vector& rhsVector = ctx->getBVector();
    const auto selection =
        metaData.locally_owned_part() &
        stk::mesh::selectUnion(domain->zonePtr()->interiorParts());
    const auto& buckets =
        bulkData.get_buckets(stk::topology::NODE_RANK, selection);
    std::vector<scalar> lhs(1);
    std::vector<scalar> rhs(1);
    std::vector<label> scratchIds(1);
    std::vector<scalar> scratchValues(1);
    std::vector<stk::mesh::Entity> connectedNodes(1);

    for (const stk::mesh::Bucket* bucketPtr : buckets)
    {
        const auto& bucket = *bucketPtr;
        const scalar* volume = stk::mesh::field_data(volumeField, bucket);
        for (stk::mesh::Bucket::size_type nodeIndex = 0;
             nodeIndex < bucket.size();
             ++nodeIndex)
        {
            lhs[0] = 0.0;
            rhs[0] = 0.0;
            connectedNodes[0] = bucket[nodeIndex];
            for (const label phaseIndex : model_->activePhaseIndices(domain))
            {
                const label localPhase =
                    domain->globalToLocalMaterialIndex(phaseIndex);
                const auto& alpha = model_->alphaRef(phaseIndex);
                const auto& rho = model_->rhoRef(phaseIndex);
                const scalar alphaCurrent = *stk::mesh::field_data(
                    alpha.stkFieldRef(), bucket, nodeIndex);
                const scalar rhoCurrent = *stk::mesh::field_data(
                    rho.stkFieldRef(), bucket, nodeIndex);
                scalar transientResidual =
                    coefficients[0] * alphaCurrent * rhoCurrent;
                const scalar alphaOld = *stk::mesh::field_data(
                    alpha.prevTimeRef().stkFieldRef(), bucket, nodeIndex);
                const scalar rhoOld = *stk::mesh::field_data(
                    rho.prevTimeRef().stkFieldRef(), bucket, nodeIndex);
                transientResidual +=
                    coefficients[1] * alphaOld * rhoOld;
                if (numberOfStates == 3)
                {
                    const scalar alphaOldOld = *stk::mesh::field_data(
                        alpha.prevTimeRef()
                            .prevTimeRef()
                            .stkFieldRef(),
                        bucket,
                        nodeIndex);
                    const scalar rhoOldOld = *stk::mesh::field_data(
                        rho.prevTimeRef()
                            .prevTimeRef()
                            .stkFieldRef(),
                        bucket,
                        nodeIndex);
                    transientResidual +=
                        coefficients[2] * alphaOldOld * rhoOldOld;
                }
                rhs[0] -= transientResidual / (rhoCurrent + SMALL) *
                          volume[nodeIndex] / dt;

                if (domain->isMaterialCompressible(localPhase))
                {
                    const scalar psi = *stk::mesh::field_data(
                        model_->psiRef(phaseIndex).stkFieldRef(),
                        bucket,
                        nodeIndex);
                    lhs[0] += coefficients[0] * alphaCurrent * psi /
                              (rhoCurrent + SMALL) * volume[nodeIndex] / dt;
                }
            }

            Base::applyCoeff_(matrix,
                              rhsVector,
                              connectedNodes,
                              scratchIds,
                              scratchValues,
                              rhs,
                              lhs);
        }
    }
}

} // namespace accel
