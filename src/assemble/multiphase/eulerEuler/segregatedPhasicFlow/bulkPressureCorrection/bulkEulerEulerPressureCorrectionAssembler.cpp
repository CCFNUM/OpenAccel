// File       : bulkEulerEulerPressureCorrectionAssembler.cpp
// Description: Common pressure-correction assembler setup and constraints

#include "bulkEulerEulerPressureCorrectionAssembler.h"

namespace accel
{

bulkEulerEulerPressureCorrectionAssembler::
    bulkEulerEulerPressureCorrectionAssembler(eulerEulerModel* model)
    : Base(model), model_(model)
{
}

void bulkEulerEulerPressureCorrectionAssembler::postAssemble(
    const domain* domain,
    Context* ctx)
{
    this->applyConstraints(domain, ctx);
}

void bulkEulerEulerPressureCorrectionAssembler::adjustMatrixForPressureReference(
    const domain* domain,
    Context* ctx)
{
    if (domain->associatedPartitionRankForPressureLevelNode() !=
        messager::myProcNo())
    {
        return;
    }

    Matrix& matrix = ctx->getAMatrix();
    Vector& rhs = ctx->getBVector();
    const auto referenceNode = model_->meshRef().bulkDataRef().get_entity(
        stk::topology::NODE_RANK, domain->pressureLevelNodeId());
    const int64_t row = matrix.getGraph()->localToRow(
        model_->meshRef().bulkDataRef().local_id(referenceNode));
    if (row < 0)
    {
        return;
    }

    const scalar diagonal = matrix.dofDiag(row);
    auto values = matrix.rowVals(row);
    std::fill(values.begin(), values.end(), 0.0);
    matrix.dofDiag(row) = diagonal;
    rhs[row] = 0.0;

}

} // namespace accel
