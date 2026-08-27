// File       : bulkEulerEulerPressureCorrectionAssemblerElemBoundaryConditions.cpp
// Description: Boundary flux terms for common Euler-Euler pressure correction

#include "bulkEulerEulerPressureCorrectionAssembler.h"

namespace accel
{

void bulkEulerEulerPressureCorrectionAssembler::assembleSideFlux_(
    const domain* domain,
    const stk::mesh::PartVector& parts,
    Context* ctx)
{
    if (parts.empty())
    {
        return;
    }
    const auto phaseIndices = model_->activePhaseIndices(domain);
    const auto& mesh = model_->meshRef();
    const auto& metaData = mesh.metaDataRef();
    const auto& bulkData = mesh.bulkDataRef();
    Matrix& matrix = ctx->getAMatrix();
    Vector& rhsVector = ctx->getBVector();
    const auto selection =
        metaData.universal_part() & stk::mesh::selectUnion(parts);
    const auto& buckets =
        bulkData.get_buckets(metaData.side_rank(), selection);
    std::vector<scalar> lhs(1, 0.0);
    std::vector<scalar> rhs(1, 0.0);
    std::vector<label> scratchIds(1);
    std::vector<scalar> scratchValues(1);
    std::vector<stk::mesh::Entity> connectedNodes(1);
    std::vector<scalar> shapeFunction;
    std::vector<scalar> shiftedShapeFunction;

    for (const stk::mesh::Bucket* bucketPtr : buckets)
    {
        const auto& bucket = *bucketPtr;
        MasterElement* meFC =
            MasterElementRepo::get_surface_master_element(bucket.topology());
        const label nodesPerSide = bucket.topology().num_nodes();
        const label numScsIp = meFC->numIntPoints_;
        const label* ipNodeMap = meFC->ipNodeMap();
        shapeFunction.resize(numScsIp * nodesPerSide);
        shiftedShapeFunction.resize(numScsIp * nodesPerSide);
        meFC->shape_fcn(shapeFunction.data());
        meFC->shifted_shape_fcn(shiftedShapeFunction.data());
        for (stk::mesh::Bucket::size_type sideIndex = 0;
             sideIndex < bucket.size();
             ++sideIndex)
        {
            const auto side = bucket[sideIndex];
            const auto* nodes = bucket.begin_nodes(sideIndex);
            for (label ip = 0; ip < numScsIp; ++ip)
            {
                scalar volumeFlux = 0.0;
                bool hasFlux = false;
                for (const label phaseIndex : phaseIndices)
                {
                    if (!model_->mDotRef(phaseIndex).sideFieldPtr())
                    {
                        continue;
                    }
                    const scalar* flux = stk::mesh::field_data(
                        model_->mDotRef(phaseIndex)
                            .sideFieldRef()
                            .stkFieldRef(),
                        side);
                    if (!flux)
                    {
                        continue;
                    }
                    hasFlux = true;
                    scalar rho = 0.0;
                    const auto& rhoField =
                        model_->rhoRef(phaseIndex).stkFieldRef();
                    const auto& phaseShape =
                        model_->URef(phaseIndex).isShifted()
                            ? shiftedShapeFunction
                            : shapeFunction;
                    for (label node = 0; node < nodesPerSide; ++node)
                    {
                        rho += phaseShape[ip * nodesPerSide + node] *
                               *stk::mesh::field_data(rhoField, nodes[node]);
                    }
                    volumeFlux += flux[ip] / (rho + SMALL);
                }
                if (!hasFlux)
                {
                    continue;
                }
                lhs[0] = 0.0;
                rhs[0] = -volumeFlux;
                connectedNodes[0] = nodes[ipNodeMap[ip]];
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
}

void bulkEulerEulerPressureCorrectionAssembler::assembleElemTermsBoundary_(
    const domain* domain,
    Context* ctx)
{
    for (label boundaryIndex = 0;
         boundaryIndex < domain->zonePtr()->nBoundaries();
         ++boundaryIndex)
    {
        const auto& boundary = domain->zonePtr()->boundaryRef(boundaryIndex);
        assembleSideFlux_(domain, boundary.parts(), ctx);

        // A boundary whose pressure is specified has p' = 0 on the face. Its
        // residual still depends on the pressure correction in the adjacent
        // element; without this sensitivity the Euler-Euler pressure system
        // receives a flux residual but has no pressure response there. This
        // covers both a static-pressure outlet and an opening, whose total
        // pressure is likewise prescribed (the opening side pressure itself is
        // built in eulerEulerModel::updatePhasicPressureSideFieldOpening_).
        if (boundary.type() == boundaryPhysicalType::outlet ||
            boundary.type() == boundaryPhysicalType::opening)
        {
            assembleSpecifiedPressureCorrection_(
                domain, &boundary, ctx);
        }
    }
}

void bulkEulerEulerPressureCorrectionAssembler::
    assembleSpecifiedPressureCorrection_(
        const domain* domain,
        const boundary* boundary,
        Context* ctx)
{
    const auto& parts = boundary->parts();
    if (parts.empty())
    {
        return;
    }

    const auto phaseIndices = model_->activePhaseIndices(domain);
    const auto& mesh = model_->meshRef();
    const auto& metaData = mesh.metaDataRef();
    const auto& bulkData = mesh.bulkDataRef();
    const auto selection =
        metaData.universal_part() & stk::mesh::selectUnion(parts);
    const auto& buckets =
        bulkData.get_buckets(metaData.side_rank(), selection);
    const auto& areaField = *metaData.get_field<scalar>(
        metaData.side_rank(), mesh::exposed_area_vector_ID);
    const bool consistent = model_->controlsRef()
                                .solverRef()
                                .solverControl_.expertParameters_.consistent_;

    std::vector<stk::topology> parentTopologies;
    std::vector<scalar> coordinates;
    std::vector<scalar> shapeFunction;
    std::vector<scalar> shiftedShapeFunction;
    std::vector<scalar> dndx;
    std::vector<scalar> detJ;
    std::vector<scalar> lhs;
    std::vector<scalar> rhs;
    std::vector<label> scratchIds;
    std::vector<scalar> scratchValues;
    std::vector<stk::mesh::Entity> connectedNodes;

    for (const stk::mesh::Bucket* bucketPtr : buckets)
    {
        const auto& bucket = *bucketPtr;
        MasterElement* meFC =
            MasterElementRepo::get_surface_master_element(bucket.topology());
        const label nodesPerSide = bucket.topology().num_nodes();
        const label numFaceIps = meFC->numIntPoints_;
        shapeFunction.resize(numFaceIps * nodesPerSide);
        shiftedShapeFunction.resize(numFaceIps * nodesPerSide);
        meFC->shape_fcn(shapeFunction.data());
        meFC->shifted_shape_fcn(shiftedShapeFunction.data());

        bucket.parent_topology(stk::topology::ELEMENT_RANK,
                               parentTopologies);
        STK_ThrowAssert(parentTopologies.size() == 1);
        MasterElement* meSCS = MasterElementRepo::get_surface_master_element(
            parentTopologies[0]);
        const label nodesPerElement = meSCS->nodesPerElement_;

        coordinates.resize(nodesPerElement * SPATIAL_DIM);
        dndx.resize(numFaceIps * nodesPerElement * SPATIAL_DIM);
        detJ.resize(numFaceIps);
        lhs.resize(nodesPerElement * nodesPerElement);
        rhs.resize(nodesPerElement);
        scratchIds.resize(nodesPerElement);
        scratchValues.resize(nodesPerElement);
        connectedNodes.resize(nodesPerElement);

        for (stk::mesh::Bucket::size_type sideIndex = 0;
             sideIndex < bucket.size();
             ++sideIndex)
        {
            const auto side = bucket[sideIndex];
            const auto* elementRelations = bulkData.begin_elements(side);
            STK_ThrowAssert(bulkData.num_elements(side) == 1);
            const auto element = elementRelations[0];
            const auto* elementOrdinals =
                bulkData.begin_element_ordinals(side);
            const label faceOrdinal = elementOrdinals[0];
            const label* faceNodeOrdinals =
                meSCS->side_node_ordinals(faceOrdinal);
            const label* scsIpNodeMap = meSCS->ipNodeMap(faceOrdinal);
            const auto* elementNodes = bulkData.begin_nodes(element);

            for (label node = 0; node < nodesPerElement; ++node)
            {
                connectedNodes[node] = elementNodes[node];
                const scalar* nodeCoordinates =
                    stk::mesh::field_data(
                        *metaData.get_field<scalar>(
                            stk::topology::NODE_RANK, mesh::coordinates_ID),
                        elementNodes[node]);
                for (label component = 0; component < SPATIAL_DIM;
                     ++component)
                {
                    coordinates[node * SPATIAL_DIM + component] =
                        nodeCoordinates[component];
                }
            }

            scalar error = 0.0;
            meSCS->face_grad_op(1,
                                faceOrdinal,
                                coordinates.data(),
                                dndx.data(),
                                detJ.data(),
                                &error);

            const auto* sideNodes = bulkData.begin_nodes(side);
            for (label ip = 0; ip < numFaceIps; ++ip)
            {
                scalar alphaDu[SPATIAL_DIM] = {};
                for (const label phaseIndex : phaseIndices)
                {
                    const auto* reversalField =
                        model_->URef(phaseIndex).reversalFlagPtr();
                    const label* reversed =
                        reversalField
                            ? stk::mesh::field_data(
                                  reversalField->stkFieldRef(), side)
                            : nullptr;
                    if (boundary->type() == boundaryPhysicalType::outlet &&
                        reversed && reversed[ip] == 1)
                    {
                        continue;
                    }
                    if (model_->URef(phaseIndex)
                            .boundaryConditionRef(domain->index(),
                                                  boundary->index())
                            .type() == boundaryConditionType::specifiedValue)
                    {
                        continue;
                    }
                    const auto& phaseAlphaField =
                        model_->alphaRef(phaseIndex).stkFieldRef();
                    const auto& phaseDuField =
                        (consistent
                             ? model_->duTildeRef(phaseIndex).stkFieldRef()
                             : model_->duRef(phaseIndex).stkFieldRef());
                    const auto* sideAlphaField =
                        model_->alphaRef(phaseIndex).sideFieldPtr()
                            ? model_->alphaRef(phaseIndex)
                                  .sideFieldRef()
                                  .stkFieldPtr()
                            : nullptr;
                    const auto* intrinsicFluxField =
                        model_->intrinsicMDotRef(phaseIndex).sideFieldPtr()
                            ? &model_->intrinsicMDotRef(phaseIndex)
                                   .sideFieldRef()
                                   .stkFieldRef()
                            : nullptr;
                    const scalar* specifiedAlpha =
                        sideAlphaField
                            ? stk::mesh::field_data(*sideAlphaField, side)
                            : nullptr;
                    const scalar* intrinsicFlux =
                        intrinsicFluxField
                            ? stk::mesh::field_data(*intrinsicFluxField, side)
                            : nullptr;
                    const auto& phaseShape =
                        model_->URef(phaseIndex).isShifted()
                            ? shiftedShapeFunction
                            : shapeFunction;

                    scalar phaseAlphaIp = 0.0;
                    scalar phaseDuIp[SPATIAL_DIM] = {};
                    for (label sideNode = 0; sideNode < nodesPerSide;
                         ++sideNode)
                    {
                        const scalar shape =
                            phaseShape[ip * nodesPerSide + sideNode];
                        phaseAlphaIp +=
                            shape * *stk::mesh::field_data(
                                        phaseAlphaField,
                                        sideNodes[sideNode]);
                        const scalar* phaseDu = stk::mesh::field_data(
                            phaseDuField, sideNodes[sideNode]);
                        for (label component = 0;
                             component < SPATIAL_DIM;
                             ++component)
                        {
                            phaseDuIp[component] +=
                                shape * phaseDu[component];
                        }
                    }
                    if (specifiedAlpha && intrinsicFlux &&
                        intrinsicFlux[ip] < 0.0)
                    {
                        phaseAlphaIp = specifiedAlpha[ip];
                    }
                    for (label component = 0; component < SPATIAL_DIM;
                         ++component)
                    {
                        alphaDu[component] +=
                            phaseAlphaIp * phaseDuIp[component];
                    }
                }

                std::fill(lhs.begin(), lhs.end(), 0.0);
                std::fill(rhs.begin(), rhs.end(), 0.0);
                const scalar* area = stk::mesh::field_data(areaField, side);
                const label rowNode = scsIpNodeMap[ip];
                for (label column = 0; column < nodesPerElement; ++column)
                {
                    bool boundaryNode = false;
                    for (label sideNode = 0; sideNode < nodesPerSide;
                         ++sideNode)
                    {
                        boundaryNode = boundaryNode ||
                                       faceNodeOrdinals[sideNode] == column;
                    }
                    if (boundaryNode)
                    {
                        continue;
                    }

                    scalar coefficient = 0.0;
                    for (label component = 0; component < SPATIAL_DIM;
                         ++component)
                    {
                        coefficient -=
                            alphaDu[component] *
                            dndx[(ip * nodesPerElement + column) *
                                     SPATIAL_DIM +
                                 component] *
                            area[ip * SPATIAL_DIM + component];
                    }
                    lhs[rowNode * nodesPerElement + column] += coefficient;
                }

                Base::applyCoeff_(ctx->getAMatrix(),
                                  ctx->getBVector(),
                                  connectedNodes,
                                  scratchIds,
                                  scratchValues,
                                  rhs,
                                  lhs);
            }
        }
    }
}

} // namespace accel
