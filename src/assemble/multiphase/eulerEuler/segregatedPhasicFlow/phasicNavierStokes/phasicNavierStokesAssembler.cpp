// File       : phasicNavierStokesAssembler.cpp
// Description: Assembler for one Euler-Euler phasic momentum equation

#include "phasicNavierStokesAssembler.h"

namespace accel
{


phasicNavierStokesAssembler::phasicNavierStokesAssembler(
    eulerEulerModel* model,
    label phaseIndex)
    : Base(model), model_(model), phaseIndex_(phaseIndex)
{
}

void phasicNavierStokesAssembler::setup(
    const std::vector<std::shared_ptr<domain>>& domains)
{
    typename Base::GammaFunction viscosity =
        [this](const domain* domain, STKScalarField& field)
    { model_->fillPhaseEffectiveViscosity(domain, phaseIndex_, field); };
    Base::setup(&model_->URef(phaseIndex_),
                advectionDiffusion,
                domains,
                viscosity);
}

void phasicNavierStokesAssembler::computeDUCoefficients_(
    const domain* domain,
    Context* ctx)
{
    auto& matrix = ctx->getAMatrix();
    const auto& mesh = model_->meshRef();
    const auto& metaData = mesh.metaDataRef();
    const auto& bulkData = mesh.bulkDataRef();
    auto& duField = model_->duRef(phaseIndex_).stkFieldRef();
    auto& duTildeField = model_->duTildeRef(phaseIndex_).stkFieldRef();
    const auto& alphaField = model_->alphaRef(phaseIndex_).stkFieldRef();
    const auto& rhoField = model_->rhoRef(phaseIndex_).stkFieldRef();
    const auto& volumeField = *metaData.get_field<scalar>(
        stk::topology::NODE_RANK, mesh::dual_nodal_volume_ID);

    const bool transient = model_->controlsRef().isTransient();
    const scalar dt = transient ? model_->controlsRef().getTimestep()
                                : model_->controlsRef().getPhysicalTimescale();
    scalar gamma = 1.0;
    if (transient)
    {
        const auto scheme = model_->controlsRef()
                                .solverRef()
                                .solverControl_.basicSettings_
                                .transientScheme_;
        if (scheme == transientSchemeType::firstOrderBackwardEuler)
        {
            gamma = BDF1::coeff()[0];
        }
        else if (scheme == transientSchemeType::secondOrderBackwardEuler)
        {
            gamma = BDF2::coeff(
                        dt, model_->controlsRef().getTimestep(-1))[0];
        }
    }
    const bool fractionalStep =
        model_->controlsRef()
            .solverRef()
            .solverControl_.expertParameters_.fractionalStepMethod_;
    const bool consistent =
        model_->controlsRef()
            .solverRef()
            .solverControl_.expertParameters_.consistent_;

    const auto selection =
        metaData.locally_owned_part() &
        stk::mesh::selectUnion(domain->zonePtr()->interiorParts());
    const auto& buckets =
        bulkData.get_buckets(stk::topology::NODE_RANK, selection);
    for (const stk::mesh::Bucket* bucketPtr : buckets)
    {
        const auto& bucket = *bucketPtr;
        scalar* du = stk::mesh::field_data(duField, bucket);
        scalar* duTilde = stk::mesh::field_data(duTildeField, bucket);
        const scalar* alpha = stk::mesh::field_data(alphaField, bucket);
        const scalar* rho = stk::mesh::field_data(rhoField, bucket);
        const scalar* volume = stk::mesh::field_data(volumeField, bucket);
        for (stk::mesh::Bucket::size_type nodeIndex = 0;
             nodeIndex < bucket.size();
             ++nodeIndex)
        {
            const int64_t row = matrix.getGraph()->localToRow(
                bulkData.local_id(bucket[nodeIndex]));
            if (row < 0)
            {
                continue;
            }
            const auto rowValues = matrix.rowVals(row);
            const auto rowColumns = matrix.rowCols(row);
            const auto& diagonalOffsets = matrix.diagOffsetRef();
            for (label component = 0; component < SPATIAL_DIM; ++component)
            {
                const scalar diagonal = matrix.dofDiag(row, component);
                const scalar pressureVolume =
                    alpha[nodeIndex] * volume[nodeIndex];
                const label offset = SPATIAL_DIM * nodeIndex + component;
                du[offset] = fractionalStep
                                 ? dt / (gamma * rho[nodeIndex] + SMALL)
                                 : pressureVolume / (diagonal + SMALL);

                scalar offDiagonal = 0.0;
                if (consistent)
                {
                    for (label column = 0;
                         column < static_cast<label>(rowColumns.size());
                         ++column)
                    {
                        if (column == diagonalOffsets[row])
                        {
                            continue;
                        }
                        offDiagonal +=
                            rowValues[SPATIAL_DIM * SPATIAL_DIM * column +
                                      SPATIAL_DIM * component + component];
                    }
                }
                duTilde[offset] =
                    pressureVolume / (diagonal + offDiagonal + SMALL);
            }
        }
    }

    model_->duRef(phaseIndex_).synchronizeGhostedEntities(domain->index());
    model_->duTildeRef(phaseIndex_)
        .synchronizeGhostedEntities(domain->index());
}

void phasicNavierStokesAssembler::postAssemble(const domain* domain,
                                                Context* ctx)
{
    Base::postAssemble(domain, ctx);
    computeDUCoefficients_(domain, ctx);
    assembleBoundaryRelaxation_(domain, ctx, 0.75);
    applySymmetryConditions_(domain, ctx);
}

} // namespace accel
