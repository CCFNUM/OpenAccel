// File       : phasicNavierStokesAssemblerNodeTerms.cpp
// Description: Nodal sources for one Euler-Euler phasic momentum equation

#include "phasicNavierStokesAssembler.h"

#include "zoneTransformation.h"

namespace accel
{

void phasicNavierStokesAssembler::assembleNodeTermsFused_(
    const domain* domain,
    Context* ctx)
{
    Base::assembleNodeTermsFused_(domain, ctx);
    assembleSources_(domain, ctx);
}

void phasicNavierStokesAssembler::assembleSources_(const domain* domain,
                                                    Context* ctx)
{
    Matrix& matrix = ctx->getAMatrix();
    Vector& rhsVector = ctx->getBVector();
    const auto& mesh = model_->meshRef();
    const auto& metaData = mesh.metaDataRef();
    const auto& bulkData = mesh.bulkDataRef();

    const auto& alphaField = model_->alphaRef(phaseIndex_).stkFieldRef();
    const auto& rhoField = model_->rhoRef(phaseIndex_).stkFieldRef();
    const auto& velocityField = model_->URef(phaseIndex_).stkFieldRef();
    const auto& pressureGradientField =
        model_->pRef().gradRef().stkFieldRef();
    const auto& dragSourceField =
        model_->interphaseMomentumSourceRef(phaseIndex_).stkFieldRef();
    const auto& dragDiagonalField =
        model_->dragDiagonalRef(phaseIndex_).stkFieldRef();
    const auto& volumeField = *metaData.get_field<scalar>(
        stk::topology::NODE_RANK, mesh::dual_nodal_volume_ID);

    const auto coriolis =
        domain->zonePtr()->frameRotating()
            ? domain->zonePtr()->transformationRef().rotation().coriolisMatrix_
            : utils::matrix::Zero();
    const scalar* coriolisData = coriolis.data();
    const bool buoyant =
        domain->buoyancy_.option_ != buoyancyOption::nonBuoyant;
    // OpenAccel solves for the reduced (excess) pressure: the body force is
    // (rho - rho_ref) g (see flowModel::computeBodyForces), the pressure is
    // initialized hydrostatically from rho_ref, and the pressure level is
    // shifted by rho_ref g.(x - x_ref). eulerEulerModel::updatePhaseBodyForces
    // builds alpha_k (rho_k - rho_ref) g per phase and redistributes it with
    // the shared flowModel kernel, so that the body force sits on the same
    // discrete footing as the compact pressure gradient. Read that field here
    // rather than applying gravity directly.
    const auto& phaseBodyForceField =
        model_->phaseBodyForceRef(phaseIndex_).stkFieldRef();
    const auto& generalSource = domain->generalMomentumSourceRef().value_;

    const auto selection =
        metaData.locally_owned_part() &
        stk::mesh::selectUnion(domain->zonePtr()->interiorParts());
    const auto& buckets =
        bulkData.get_buckets(stk::topology::NODE_RANK, selection);

    std::vector<scalar> lhs(SPATIAL_DIM * SPATIAL_DIM, 0.0);
    std::vector<scalar> rhs(SPATIAL_DIM, 0.0);
    std::vector<label> scratchIds(SPATIAL_DIM);
    std::vector<scalar> scratchValues(SPATIAL_DIM);
    std::vector<stk::mesh::Entity> connectedNodes(1);

    for (const stk::mesh::Bucket* bucketPtr : buckets)
    {
        const auto& bucket = *bucketPtr;
        const scalar* alpha = stk::mesh::field_data(alphaField, bucket);
        const scalar* rho = stk::mesh::field_data(rhoField, bucket);
        const scalar* velocity =
            stk::mesh::field_data(velocityField, bucket);
        const scalar* pressureGradient =
            stk::mesh::field_data(pressureGradientField, bucket);
        const scalar* dragSource =
            stk::mesh::field_data(dragSourceField, bucket);
        const scalar* dragDiagonal =
            stk::mesh::field_data(dragDiagonalField, bucket);
        const scalar* volume = stk::mesh::field_data(volumeField, bucket);
        const scalar* phaseBodyForce =
            stk::mesh::field_data(phaseBodyForceField, bucket);

        for (stk::mesh::Bucket::size_type nodeIndex = 0;
             nodeIndex < bucket.size();
             ++nodeIndex)
        {
            std::fill(lhs.begin(), lhs.end(), 0.0);
            std::fill(rhs.begin(), rhs.end(), 0.0);
            connectedNodes[0] = bucket[nodeIndex];
            const scalar alphaRho = alpha[nodeIndex] * rho[nodeIndex];
            const scalar nodalVolume = volume[nodeIndex];

            for (label row = 0; row < SPATIAL_DIM; ++row)
            {
                rhs[row] -=
                    alpha[nodeIndex] *
                    pressureGradient[SPATIAL_DIM * nodeIndex + row] *
                    nodalVolume;
                rhs[row] +=
                    dragSource[SPATIAL_DIM * nodeIndex + row] * nodalVolume;
                rhs[row] +=
                    alpha[nodeIndex] * generalSource[row] * nodalVolume;
                if (buoyant)
                {
                    rhs[row] +=
                        phaseBodyForce[SPATIAL_DIM * nodeIndex + row] *
                        nodalVolume;
                }

                lhs[row * SPATIAL_DIM + row] +=
                    dragDiagonal[nodeIndex] * nodalVolume;
                for (label column = 0; column < SPATIAL_DIM; ++column)
                {
                    const scalar coefficient =
                        alphaRho *
                        coriolisData[row * SPATIAL_DIM + column] *
                        nodalVolume;
                    lhs[row * SPATIAL_DIM + column] += coefficient;
                    rhs[row] -=
                        coefficient *
                        velocity[SPATIAL_DIM * nodeIndex + column];
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
