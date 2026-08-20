// File       : bulkEulerEulerPressureCorrectionAssemblerElemTerms.cpp
// Description: Interior terms for common Euler-Euler pressure correction

#include "bulkEulerEulerPressureCorrectionAssembler.h"

namespace accel
{

void bulkEulerEulerPressureCorrectionAssembler::assembleElemTermsInterior_(
    const domain* domain,
    Context* ctx)
{
    const auto phaseIndices = model_->activePhaseIndices(domain);
    const auto& mesh = model_->meshRef();
    const auto& metaData = mesh.metaDataRef();
    const auto& bulkData = mesh.bulkDataRef();
    const auto& coordinatesField = *metaData.get_field<scalar>(
        stk::topology::NODE_RANK, mesh::coordinates_ID);
    const bool consistent = model_->controlsRef()
                                .solverRef()
                                .solverControl_.expertParameters_.consistent_;
    Matrix& matrix = ctx->getAMatrix();
    Vector& rhsVector = ctx->getBVector();

    const auto selection =
        metaData.universal_part() &
        stk::mesh::selectUnion(domain->zonePtr()->interiorParts());
    const auto& buckets =
        bulkData.get_buckets(stk::topology::ELEMENT_RANK, selection);
    std::vector<scalar> lhs;
    std::vector<scalar> rhs;
    std::vector<label> scratchIds;
    std::vector<scalar> scratchValues;
    std::vector<stk::mesh::Entity> connectedNodes;
    std::vector<scalar> coordinates;
    std::vector<scalar> shapeFunction;
    std::vector<scalar> areaVector;
    std::vector<scalar> dndx;
    std::vector<scalar> deriv;
    std::vector<scalar> detJ;

    for (const stk::mesh::Bucket* bucketPtr : buckets)
    {
        const auto& bucket = *bucketPtr;
        MasterElement* meSCS =
            MasterElementRepo::get_surface_master_element(bucket.topology());
        const label nodesPerElement = meSCS->nodesPerElement_;
        const label numScsIp = meSCS->numIntPoints_;
        const label* adjacentNodes = meSCS->adjacentNodes();
        lhs.resize(nodesPerElement * nodesPerElement);
        rhs.resize(nodesPerElement);
        scratchIds.resize(nodesPerElement);
        scratchValues.resize(nodesPerElement);
        connectedNodes.resize(nodesPerElement);
        coordinates.resize(nodesPerElement * SPATIAL_DIM);
        shapeFunction.resize(numScsIp * nodesPerElement);
        areaVector.resize(numScsIp * SPATIAL_DIM);
        dndx.resize(numScsIp * nodesPerElement * SPATIAL_DIM);
        deriv.resize(numScsIp * nodesPerElement * SPATIAL_DIM);
        detJ.resize(numScsIp);
        meSCS->shape_fcn(shapeFunction.data());

        for (stk::mesh::Bucket::size_type elementIndex = 0;
             elementIndex < bucket.size();
             ++elementIndex)
        {
            std::fill(lhs.begin(), lhs.end(), 0.0);
            std::fill(rhs.begin(), rhs.end(), 0.0);
            const auto* nodes = bucket.begin_nodes(elementIndex);
            for (label nodeIndex = 0; nodeIndex < nodesPerElement;
                 ++nodeIndex)
            {
                connectedNodes[nodeIndex] = nodes[nodeIndex];
                const scalar* nodeCoordinates =
                    stk::mesh::field_data(coordinatesField, nodes[nodeIndex]);
                for (label component = 0; component < SPATIAL_DIM;
                     ++component)
                {
                    coordinates[SPATIAL_DIM * nodeIndex + component] =
                        nodeCoordinates[component];
                }
            }
            scalar error = 0.0;
            meSCS->determinant(
                1, coordinates.data(), areaVector.data(), &error);
            meSCS->grad_op(1,
                           coordinates.data(),
                           dndx.data(),
                           deriv.data(),
                           detJ.data(),
                           &error);

            for (label ip = 0; ip < numScsIp; ++ip)
            {
                const label left = adjacentNodes[2 * ip];
                const label right = adjacentNodes[2 * ip + 1];
                scalar volumeFlux = 0.0;
                for (const label phaseIndex : phaseIndices)
                {
                    const scalar* phaseFlux = stk::mesh::field_data(
                        model_->mDotRef(phaseIndex).stkFieldRef(),
                        bucket,
                        elementIndex);
                    scalar rhoIp = 0.0;
                    scalar alphaIp = 0.0;
                    std::array<scalar, SPATIAL_DIM> duIp{};
                    const auto& rhoField =
                        model_->rhoRef(phaseIndex).stkFieldRef();
                    const auto& alphaField =
                        model_->alphaRef(phaseIndex).stkFieldRef();
                    const auto& duField =
                        consistent
                            ? model_->duTildeRef(phaseIndex).stkFieldRef()
                            : model_->duRef(phaseIndex).stkFieldRef();
                    for (label nodeIndex = 0;
                         nodeIndex < nodesPerElement;
                         ++nodeIndex)
                    {
                        const scalar shape =
                            shapeFunction[ip * nodesPerElement + nodeIndex];
                        rhoIp += shape * *stk::mesh::field_data(
                                             rhoField, nodes[nodeIndex]);
                        alphaIp += shape * *stk::mesh::field_data(
                                               alphaField, nodes[nodeIndex]);
                        const scalar* nodeDu =
                            stk::mesh::field_data(duField, nodes[nodeIndex]);
                        for (label component = 0;
                             component < SPATIAL_DIM;
                             ++component)
                        {
                            duIp[component] += shape * nodeDu[component];
                        }
                    }
                    volumeFlux += phaseFlux[ip] / (rhoIp + SMALL);
                    for (label column = 0;
                         column < nodesPerElement;
                         ++column)
                    {
                        scalar coefficient = 0.0;
                        for (label component = 0;
                             component < SPATIAL_DIM;
                             ++component)
                        {
                            coefficient -=
                                alphaIp * duIp[component] *
                                dndx[(ip * nodesPerElement + column) *
                                          SPATIAL_DIM +
                                      component] *
                                areaVector[ip * SPATIAL_DIM + component];
                        }
                        lhs[left * nodesPerElement + column] += coefficient;
                        lhs[right * nodesPerElement + column] -= coefficient;
                    }
                }
                rhs[left] -= volumeFlux;
                rhs[right] += volumeFlux;
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
