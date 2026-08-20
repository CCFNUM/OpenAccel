// File       : phasicNavierStokesAssemblerElemBoundaryConditions.cpp
// Description: Boundary conditions for one Euler-Euler phasic momentum
//              equation. The phase momentum equation is a vector transport
//              equation, so the generic phiAssembler boundary switch (which
//              only knows the scalar Dirichlet/Neumann types) cannot dispatch
//              the velocity boundary types `noSlip` and `slip`. This file maps
//              the momentum boundary types onto the generic kernels and adds
//              the vector symmetry treatment.
//
//              Diffusivity at the boundary is the phase effective viscosity
//              alpha_k * mu_k supplied through the GammaFunction in setup(),
//              and the transient/convective coefficients come from the
//              rhoRef()/mDotRef() overrides, so the generic kernels already
//              produce alpha-weighted wall shear and boundary fluxes.

#include "phasicNavierStokesAssembler.h"

#include <cmath>
#include <unordered_set>

namespace accel
{

// These mirror the momentum kernels in
// navierStokesAssemblerElemBoundaryConditions.cpp, where they are file-local.
// They are duplicated (rather than shared) to keep this file's assembly
// self-contained; the phase weighting enters through tmDot and muEffBip.
#define PHASIC_IP_EXPLICIT_ADVECTIVE_FLUX__(flux)                              \
    do                                                                         \
    {                                                                          \
        for (label i = 0; i < SPATIAL_DIM; ++i)                                \
        {                                                                      \
            const label indexR = nearestNode * SPATIAL_DIM + i;                \
            p_rhs[indexR] -= (flux);                                           \
        }                                                                      \
    } while (0)

#define PHASIC_IP_IMPLICIT_UPWIND_ADVECTIVE_FLUX__(flux)                       \
    do                                                                         \
    {                                                                          \
        for (label i = 0; i < SPATIAL_DIM; ++i)                                \
        {                                                                      \
            const label indexR = nearestNode * SPATIAL_DIM + i;                \
            const label rowR = indexR * nodesPerElement * SPATIAL_DIM;         \
            p_rhs[indexR] -= (flux);                                           \
            p_lhs[rowR + nearestNode * SPATIAL_DIM + i] += tmDot;              \
        }                                                                      \
    } while (0)

#define PHASIC_IP_ZERO_NORMAL_STRESS__()                                       \
    do                                                                         \
    {                                                                          \
        for (label ic = 0; ic < nodesPerElement; ++ic)                         \
        {                                                                      \
            const label offSetDnDx =                                           \
                SPATIAL_DIM * nodesPerElement * ip + ic * SPATIAL_DIM;         \
            for (label j = 0; j < SPATIAL_DIM; ++j)                            \
            {                                                                  \
                const scalar axj = areaVec[faceOffSet + j];                    \
                const scalar dndxj = p_dndx[offSetDnDx + j];                   \
                const scalar uxj = p_U[ic * SPATIAL_DIM + j];                  \
                const scalar divUstress =                                      \
                    2.0 / 3.0 * muEffBip * dndxj * uxj * axj * comp;           \
                for (label i = 0; i < SPATIAL_DIM; ++i)                        \
                {                                                              \
                    label indexR = nearestNode * SPATIAL_DIM + i;              \
                    label rowR = indexR * nodesPerElement * SPATIAL_DIM;       \
                    const scalar dndxi = p_dndx[offSetDnDx + i];               \
                    const scalar uxi = p_U[ic * SPATIAL_DIM + i];              \
                    const scalar nxi = p_nx[i];                                \
                    const scalar om_nxinxi = 1.0 - nxi * nxi;                  \
                                                                               \
                    scalar lhsfac = -muEffBip * dndxj * axj * om_nxinxi;       \
                    p_lhs[rowR + ic * SPATIAL_DIM + i] += lhsfac;              \
                    p_rhs[indexR] -= lhsfac * uxi + divUstress * om_nxinxi;    \
                                                                               \
                    lhsfac = -muEffBip * dndxi * axj * om_nxinxi;              \
                    p_lhs[rowR + ic * SPATIAL_DIM + j] += lhsfac;              \
                    p_rhs[indexR] -= lhsfac * uxj;                             \
                                                                               \
                    for (label l = 0; l < SPATIAL_DIM; ++l)                    \
                    {                                                          \
                        if (i != l)                                            \
                        {                                                      \
                            const scalar nxinxl = nxi * p_nx[l];               \
                            const scalar uxl = p_U[ic * SPATIAL_DIM + l];      \
                            const scalar dndxl = p_dndx[offSetDnDx + l];       \
                                                                               \
                            lhsfac = muEffBip * dndxj * axj * nxinxl;          \
                            p_lhs[rowR + ic * SPATIAL_DIM + l] += lhsfac;      \
                            p_rhs[indexR] -=                                   \
                                lhsfac * uxl + divUstress * nxinxl;            \
                                                                               \
                            lhsfac = muEffBip * dndxl * axj * nxinxl;          \
                            p_lhs[rowR + ic * SPATIAL_DIM + j] += lhsfac;      \
                            p_rhs[indexR] -= lhsfac * uxj;                     \
                        }                                                      \
                    }                                                          \
                }                                                              \
            }                                                                  \
        }                                                                      \
    } while (0)

#define PHASIC_IP_ZERO_NORMAL_STRESS_NO_GRADU_TRANSPOSE__()                    \
    do                                                                         \
    {                                                                          \
        for (label ic = 0; ic < nodesPerElement; ++ic)                         \
        {                                                                      \
            const label offSetDnDx =                                           \
                SPATIAL_DIM * nodesPerElement * ip + ic * SPATIAL_DIM;         \
            for (label j = 0; j < SPATIAL_DIM; ++j)                            \
            {                                                                  \
                const scalar axj = areaVec[faceOffSet + j];                    \
                const scalar dndxj = p_dndx[offSetDnDx + j];                   \
                const scalar uxj = p_U[ic * SPATIAL_DIM + j];                  \
                const scalar divUstress =                                      \
                    2.0 / 3.0 * muEffBip * dndxj * uxj * axj * comp;           \
                for (label i = 0; i < SPATIAL_DIM; ++i)                        \
                {                                                              \
                    label indexR = nearestNode * SPATIAL_DIM + i;              \
                    label rowR = indexR * nodesPerElement * SPATIAL_DIM;       \
                    const scalar uxi = p_U[ic * SPATIAL_DIM + i];              \
                    const scalar nxi = p_nx[i];                                \
                    const scalar om_nxinxi = 1.0 - nxi * nxi;                  \
                                                                               \
                    scalar lhsfac = -muEffBip * dndxj * axj * om_nxinxi;       \
                    p_lhs[rowR + ic * SPATIAL_DIM + i] += lhsfac;              \
                    p_rhs[indexR] -= lhsfac * uxi + divUstress * om_nxinxi;    \
                                                                               \
                    for (label l = 0; l < SPATIAL_DIM; ++l)                    \
                    {                                                          \
                        if (i != l)                                            \
                        {                                                      \
                            const scalar nxinxl = nxi * p_nx[l];               \
                            const scalar uxl = p_U[ic * SPATIAL_DIM + l];      \
                                                                               \
                            lhsfac = muEffBip * dndxj * axj * nxinxl;          \
                            p_lhs[rowR + ic * SPATIAL_DIM + l] += lhsfac;      \
                            p_rhs[indexR] -=                                   \
                                lhsfac * uxl + divUstress * nxinxl;            \
                        }                                                      \
                    }                                                          \
                }                                                              \
            }                                                                  \
        }                                                                      \
    } while (0)


void phasicNavierStokesAssembler::assembleElemTermsBoundary_(
    const domain* domain,
    Context* ctx)
{
    const zone* zonePtr = domain->zonePtr();

    for (label iBoundary = 0; iBoundary < zonePtr->nBoundaries(); iBoundary++)
    {
        const boundary* boundary = zonePtr->boundaryPtr(iBoundary);

        const boundaryPhysicalType type = boundary->type();

        const boundaryConditionType UBCType =
            model_->URef(phaseIndex_)
                .boundaryConditionRef(domain->index(), iBoundary)
                .type();

        switch (type)
        {
            case boundaryPhysicalType::symmetry:
                {
                    // Zero shear and zero normal flux. The normal momentum
                    // component is decoupled from the system in
                    // applySymmetryConditions_ once the row is assembled.
                    Base::assembleElemTermsBoundarySymmetry_(
                        domain, boundary, ctx);
                }
                break;

            case boundaryPhysicalType::wall:
                {
                    switch (UBCType)
                    {
                        case boundaryConditionType::noSlip:
                        case boundaryConditionType::specifiedValue:
                            {
                                // Dirichlet U_k = U_wall (zero for no-slip).
                                // The fixed-value kernel evaluates the viscous
                                // flux with the phase diffusivity, which is the
                                // alpha_k-weighted wall shear on this phase.
                                Base::assembleElemTermsBoundaryWallFixedValue_(
                                    domain, boundary, ctx);
                            }
                            break;

                        case boundaryConditionType::slip:
                            {
                                // A free-slip wall carries no shear: identical
                                // treatment to a symmetry plane.
                                Base::assembleElemTermsBoundarySymmetry_(
                                    domain, boundary, ctx);
                            }
                            break;

                        case boundaryConditionType::zeroGradient:
                            {
                                Base::assembleElemTermsBoundaryWallZeroGradient_(
                                    domain, boundary, ctx);
                            }
                            break;

                        default:
                            errorMsg(
                                "phasic momentum: invalid velocity boundary "
                                "condition at wall `" +
                                boundary->name() + "` for phase index " +
                                std::to_string(phaseIndex_));
                    }
                }
                break;

            case boundaryPhysicalType::inlet:
                {
                    switch (UBCType)
                    {
                        case boundaryConditionType::specifiedValue:
                            {
                                Base::assembleElemTermsBoundaryInletFixedValue_(
                                    domain, boundary, ctx);
                            }
                            break;

                        case boundaryConditionType::zeroGradient:
                            {
                                Base::
                                    assembleElemTermsBoundaryOutletZeroGradient_(
                                        domain, boundary, ctx);
                            }
                            break;

                        default:
                            errorMsg(
                                "phasic momentum: invalid velocity boundary "
                                "condition at inlet `" +
                                boundary->name() + "` for phase index " +
                                std::to_string(phaseIndex_));
                    }
                }
                break;

            case boundaryPhysicalType::outlet:
                {
                    switch (UBCType)
                    {
                        case boundaryConditionType::zeroGradient:
                            {
                                Base::
                                    assembleElemTermsBoundaryOutletZeroGradient_(
                                        domain, boundary, ctx);
                            }
                            break;

                        case boundaryConditionType::specifiedValue:
                            {
                                Base::assembleElemTermsBoundaryInletFixedValue_(
                                    domain, boundary, ctx);
                            }
                            break;

                        default:
                            errorMsg(
                                "phasic momentum: invalid velocity boundary "
                                "condition at outlet `" +
                                boundary->name() + "` for phase index " +
                                std::to_string(phaseIndex_));
                    }
                }
                break;

            case boundaryPhysicalType::opening:
                {
                    // specialized below; see the header for why the generic
                    // scalar opening kernel cannot be reused here
                    assembleElemTermsBoundaryOpening_(domain, boundary, ctx);
                }
                break;

            default:
                break;
        }
    }
}

void phasicNavierStokesAssembler::applySymmetryConditions_(const domain* domain,
                                                           Context* ctx)
{
    Matrix& A = ctx->getAMatrix();
    Vector& b = ctx->getBVector();

    const auto& mesh = model_->meshRef();
    const stk::mesh::MetaData& metaData = mesh.metaDataRef();
    const stk::mesh::BulkData& bulkData = mesh.bulkDataRef();
    const zone* zonePtr = domain->zonePtr();

    // collect symmetry boundary parts (slip walls are handled as symmetry too)
    stk::mesh::PartVector partVec;
    for (label iB = 0; iB < zonePtr->nBoundaries(); ++iB)
    {
        const auto& boundaryRef = zonePtr->boundaryRef(iB);
        const bool isSymmetry =
            boundaryRef.type() == boundaryPhysicalType::symmetry;
        const bool isSlipWall =
            boundaryRef.type() == boundaryPhysicalType::wall &&
            model_->URef(phaseIndex_)
                    .boundaryConditionRef(domain->index(), iB)
                    .type() == boundaryConditionType::slip;
        if (isSymmetry || isSlipWall)
        {
            for (auto* pp : boundaryRef.parts())
            {
                partVec.push_back(pp);
            }
        }
    }
    if (partVec.empty())
    {
        return;
    }

    const auto& symmArea = *metaData.template get_field<scalar>(
        stk::topology::NODE_RANK, mesh::assembled_symm_area_ID);

    // interface slave nodes shared with a symmetry boundary must not have their
    // normal residual cancelled twice
    std::unordered_set<std::size_t> slaveIfcNodes;
    for (const interface* interf : domain->interfacesRef())
    {
        if (!interf->isConformalTreatment())
        {
            continue;
        }
        for (const auto& np : interf->matchingNodePairVector())
        {
            slaveIfcNodes.insert(np.second.local_offset());
        }
    }

    scalar n[SPATIAL_DIM];

    const auto& diagOffsets = A.diagOffsetRef();

    stk::mesh::Selector sel =
        metaData.locally_owned_part() & stk::mesh::selectUnion(partVec);
    for (const stk::mesh::Bucket* bkt :
         bulkData.get_buckets(stk::topology::NODE_RANK, sel))
    {
        for (size_t iN = 0; iN < bkt->size(); ++iN)
        {
            const stk::mesh::Entity node = (*bkt)[iN];
            if (slaveIfcNodes.count(node.local_offset()))
            {
                continue;
            }

            const int64_t row =
                A.getGraph()->localToRow(bulkData.local_id(node));
            if (row < 0) // node not part of this (subset) system
            {
                continue;
            }

            // unit symmetry normal
            const scalar* aarea = stk::mesh::field_data(symmArea, node);
            scalar asq = 0.0;
            for (label i = 0; i < SPATIAL_DIM; ++i)
            {
                n[i] = aarea[i];
                asq += n[i] * n[i];
            }
            if (asq < SMALL)
            {
                continue;
            }
            const scalar amag = std::sqrt(asq);
            for (label i = 0; i < SPATIAL_DIM; ++i)
            {
                n[i] /= amag;
            }

            auto vals = A.rowVals(row);
            const label nBlocks =
                static_cast<label>(vals.size()) / (SPATIAL_DIM * SPATIAL_DIM);
            const label diagBk = diagOffsets[row];

            // normal stiffness scale = n^T D n (keeps du physical after decoup)
            scalar* Dblk = &vals[SPATIAL_DIM * SPATIAL_DIM * diagBk];
            scalar scale = 0.0;
            for (label i = 0; i < SPATIAL_DIM; ++i)
            {
                for (label j = 0; j < SPATIAL_DIM; ++j)
                {
                    scale += n[i] * Dblk[i * SPATIAL_DIM + j] * n[j];
                }
            }
            if (std::abs(scale) < SMALL)
            {
                scale = amag;
            }

            // For every block B in the row: B <- T B (zero its normal row), and
            // for the diagonal block restore the normal coefficient (+ scale
            // N). Result: the normal equation is scale * (n.U) = 0 with no
            // tangential or neighbour coupling, i.e. U.n = 0 fully decoupled.
            for (label bk = 0; bk < nBlocks; ++bk)
            {
                scalar* B = &vals[SPATIAL_DIM * SPATIAL_DIM * bk];
                scalar colProj[SPATIAL_DIM];
                for (label j = 0; j < SPATIAL_DIM; ++j)
                {
                    scalar s = 0.0;
                    for (label k = 0; k < SPATIAL_DIM; ++k)
                    {
                        s += n[k] * B[k * SPATIAL_DIM + j];
                    }
                    colProj[j] = s;
                }
                for (label i = 0; i < SPATIAL_DIM; ++i)
                {
                    for (label j = 0; j < SPATIAL_DIM; ++j)
                    {
                        B[i * SPATIAL_DIM + j] -= n[i] * colProj[j];
                    }
                }
                if (bk == diagBk)
                {
                    for (label i = 0; i < SPATIAL_DIM; ++i)
                    {
                        for (label j = 0; j < SPATIAL_DIM; ++j)
                        {
                            B[i * SPATIAL_DIM + j] += scale * n[i] * n[j];
                        }
                    }
                }
            }

            // RHS: remove the normal component (normal equation rhs = 0)
            scalar* rhs = &b[SPATIAL_DIM * row];
            scalar bproj = 0.0;
            for (label i = 0; i < SPATIAL_DIM; ++i)
            {
                bproj += n[i] * rhs[i];
            }
            for (label i = 0; i < SPATIAL_DIM; ++i)
            {
                rhs[i] -= n[i] * bproj;
            }
        }
    }
}


// Momentum opening kernel for one phase. This mirrors
// navierStokesAssembler::assembleElemTermsBoundaryOpening_: outflow is an
// implicit upwind advective flux with zero normal stress; inflow entrains
// along the opening flow direction with the normal speed recovered from the
// interior. The phase mass flux mDot_k = alpha_k*rho_k*U_k.A and the phase
// effective viscosity alpha_k*mu_k (the assembler Gamma) already carry the
// alpha_k weighting, so no extra scaling is applied here.
void phasicNavierStokesAssembler::assembleElemTermsBoundaryOpening_(
    const domain* domain,
    const boundary* boundary,
    Context* ctx)
{
    auto& mesh = model_->meshRef();

    Matrix& A = ctx->getAMatrix();
    Vector& b = ctx->getBVector();

    const label localPhaseIndex =
        domain->globalToLocalMaterialIndex(phaseIndex_);
    const scalar comp =
        domain->isMaterialCompressible(localPhaseIndex) ? 1.0 : 0.0;

    stk::mesh::BulkData& bulkData = mesh.bulkDataRef();
    stk::mesh::MetaData& metaData = mesh.metaDataRef();

    std::vector<scalar> lhs;
    std::vector<scalar> rhs;
    std::vector<label> scratchIds;
    std::vector<scalar> scratchVals;
    std::vector<stk::mesh::Entity> connectedNodes;

    std::vector<scalar> uBip(SPATIAL_DIM);
    std::vector<scalar> nx(SPATIAL_DIM);
    scalar* p_uBip = &uBip[0];
    scalar* p_nx = &nx[0];

    // nodal fields to gather
    std::vector<scalar> ws_U;
    std::vector<scalar> ws_coords;
    std::vector<scalar> ws_muEff;

    // master element
    std::vector<scalar> ws_velocity_face_shape_function;
    std::vector<scalar> ws_dndx;
    std::vector<scalar> ws_det_j;

    auto& phaseVelocity = model_->URef(phaseIndex_);
    const auto& USTKFieldRef = phaseVelocity.stkFieldRef();

    // alpha_k*mu_k, filled by the GammaFunction installed in setup()
    STK_ThrowRequireMsg(GammaSTKFieldPtr_ != nullptr,
                        "phasic momentum opening `"
                            << boundary->name()
                            << "`: phase effective viscosity is not set up");
    const auto& muEffSTKFieldRef = *GammaSTKFieldPtr_;

    // mDot_k side field; registered for inlet/outlet/opening patches in
    // fieldBroker::setupMassFlowRate
    STK_ThrowRequireMsg(mDotSideSTKFieldPtr_ != nullptr,
                        "phasic momentum opening `"
                            << boundary->name()
                            << "`: phase mass flux side field is missing");
    const auto& mDotSideSTKFieldRef = *mDotSideSTKFieldPtr_;

    // opening flow direction; registered in fieldBroker::setupVelocity for
    // every phase at an opening. Fall back to the inward face normal.
    const bool hasSideFlowDirection =
        phaseVelocity.hasSideFlowDirectionFields();
    const auto* sideFlowDirectionSTKFieldPtr =
        hasSideFlowDirection
            ? &phaseVelocity.sideFlowDirectionFieldRef().stkFieldRef()
            : nullptr;

    // Get geometric fields
    const auto& exposedAreaVecSTKFieldRef = *metaData.get_field<scalar>(
        metaData.side_rank(), this->getExposedAreaVectorID_(domain));
    const auto& coordsSTKFieldRef = *metaData.get_field<scalar>(
        stk::topology::NODE_RANK, this->getCoordinatesID_(domain));

    std::vector<stk::topology> parentTopo;

    stk::mesh::Selector selAllSides =
        metaData.universal_part() & stk::mesh::selectUnion(boundary->parts());

    const bool isUShifted = phaseVelocity.isShifted();
    const bool isUGradientShifted = phaseVelocity.isGradientShifted();

    stk::mesh::BucketVector const& sideBuckets =
        bulkData.get_buckets(metaData.side_rank(), selAllSides);

    for (stk::mesh::BucketVector::const_iterator ib = sideBuckets.begin();
         ib != sideBuckets.end();
         ++ib)
    {
        stk::mesh::Bucket& sideBucket = **ib;

        sideBucket.parent_topology(stk::topology::ELEMENT_RANK, parentTopo);
        STK_ThrowAssert(parentTopo.size() == 1);
        stk::topology theElemTopo = parentTopo[0];

        MasterElement* meSCS =
            MasterElementRepo::get_surface_master_element(theElemTopo);
        const label nodesPerElement = meSCS->nodesPerElement_;

        MasterElement* meFC = MasterElementRepo::get_surface_master_element(
            sideBucket.topology());
        const label nodesPerSide = meFC->nodesPerElement_;
        const label numScsBip = meFC->numIntPoints_;

        const label lhsSize =
            nodesPerElement * SPATIAL_DIM * nodesPerElement * SPATIAL_DIM;
        const label rhsSize = nodesPerElement * SPATIAL_DIM;
        lhs.resize(lhsSize);
        rhs.resize(rhsSize);
        scratchIds.resize(rhsSize);
        scratchVals.resize(rhsSize);
        connectedNodes.resize(nodesPerElement);

        ws_U.resize(nodesPerElement * SPATIAL_DIM);
        ws_coords.resize(nodesPerElement * SPATIAL_DIM);
        ws_muEff.resize(nodesPerSide);
        ws_velocity_face_shape_function.resize(numScsBip * nodesPerSide);
        ws_dndx.resize(SPATIAL_DIM * numScsBip * nodesPerElement);
        ws_det_j.resize(numScsBip);

        scalar* p_lhs = &lhs[0];
        scalar* p_rhs = &rhs[0];
        scalar* p_U = &ws_U[0];
        scalar* p_coords = &ws_coords[0];
        scalar* p_muEff = &ws_muEff[0];
        scalar* p_velocity_face_shape_function =
            &ws_velocity_face_shape_function[0];
        scalar* p_dndx = &ws_dndx[0];

        if (isUShifted)
        {
            meFC->shifted_shape_fcn(&p_velocity_face_shape_function[0]);
        }
        else
        {
            meFC->shape_fcn(&p_velocity_face_shape_function[0]);
        }

        const stk::mesh::Bucket::size_type nSidesPerBucket = sideBucket.size();

        for (stk::mesh::Bucket::size_type iSide = 0; iSide < nSidesPerBucket;
             ++iSide)
        {
            for (label p = 0; p < lhsSize; ++p)
            {
                p_lhs[p] = 0.0;
            }
            for (label p = 0; p < rhsSize; ++p)
            {
                p_rhs[p] = 0.0;
            }

            stk::mesh::Entity side = sideBucket[iSide];

            const scalar* mDot =
                stk::mesh::field_data(mDotSideSTKFieldRef, side);
            const scalar* dir = sideFlowDirectionSTKFieldPtr
                                    ? stk::mesh::field_data(
                                          *sideFlowDirectionSTKFieldPtr, side)
                                    : nullptr;
            const scalar* areaVec =
                stk::mesh::field_data(exposedAreaVecSTKFieldRef, side);

            stk::mesh::Entity const* faceElemRels =
                bulkData.begin_elements(side);
            STK_ThrowAssert(bulkData.num_elements(side) == 1);

            stk::mesh::Entity element = faceElemRels[0];
            const label faceOrdinal = bulkData.begin_element_ordinals(side)[0];

            const label* ipNodeMap = meSCS->ipNodeMap(faceOrdinal);
            const label* faceNodeOrdinals =
                meSCS->side_node_ordinals(faceOrdinal);

            //==========================================
            // gather nodal data off of element
            //==========================================
            stk::mesh::Entity const* elemNodeRels =
                bulkData.begin_nodes(element);
            label numNodes = bulkData.num_nodes(element);

            STK_ThrowAssert(numNodes == nodesPerElement);
            for (label ni = 0; ni < numNodes; ++ni)
            {
                stk::mesh::Entity node = elemNodeRels[ni];

                connectedNodes[ni] = node;

                const scalar* U = stk::mesh::field_data(USTKFieldRef, node);
                const scalar* coords =
                    stk::mesh::field_data(coordsSTKFieldRef, node);
                const label offSet = ni * SPATIAL_DIM;
                for (label j = 0; j < SPATIAL_DIM; ++j)
                {
                    p_U[offSet + j] = U[j];
                    p_coords[offSet + j] = coords[j];
                }
            }

            //======================================
            // gather nodal data off of face
            //======================================
            stk::mesh::Entity const* sideNodeRels = bulkData.begin_nodes(side);
            label numSideNodes = bulkData.num_nodes(side);

            STK_ThrowAssert(numSideNodes == nodesPerSide);
            for (label ni = 0; ni < numSideNodes; ++ni)
            {
                stk::mesh::Entity node = sideNodeRels[ni];
                p_muEff[ni] = *stk::mesh::field_data(muEffSTKFieldRef, node);
            }

            // compute dndx
            scalar scs_error = 0.0;
            if (isUGradientShifted)
            {
                meSCS->shifted_face_grad_op(1,
                                            faceOrdinal,
                                            &p_coords[0],
                                            &p_dndx[0],
                                            &ws_det_j[0],
                                            &scs_error);
            }
            else
            {
                meSCS->face_grad_op(1,
                                    faceOrdinal,
                                    &p_coords[0],
                                    &p_dndx[0],
                                    &ws_det_j[0],
                                    &scs_error);
            }

            // loop over boundary ips
            for (label ip = 0; ip < numScsBip; ++ip)
            {
                const label nearestNode = ipNodeMap[ip];

                const label faceOffSet = ip * SPATIAL_DIM;
                const label offSetSF = ip * nodesPerSide;

                scalar asq = 0.0;
                for (label j = 0; j < SPATIAL_DIM; ++j)
                {
                    const scalar axj = areaVec[faceOffSet + j];
                    asq += axj * axj;
                }
                const scalar amag = std::sqrt(asq);

                for (label j = 0; j < SPATIAL_DIM; ++j)
                {
                    p_uBip[j] = 0.0;
                }

                scalar muEffBip = 0.0;
                for (label ic = 0; ic < nodesPerSide; ++ic)
                {
                    const label inn = faceNodeOrdinals[ic];

                    const scalar r_vel =
                        p_velocity_face_shape_function[offSetSF + ic];

                    muEffBip += r_vel * p_muEff[ic];

                    const label icNdim = inn * SPATIAL_DIM;
                    for (label j = 0; j < SPATIAL_DIM; ++j)
                    {
                        p_uBip[j] += r_vel * p_U[icNdim + j];
                    }
                }

                for (label i = 0; i < SPATIAL_DIM; ++i)
                {
                    p_nx[i] = areaVec[faceOffSet + i] / amag;
                }

                const scalar tmDot = mDot[ip];

                if (tmDot > 0.0) // outflow
                {
                    PHASIC_IP_IMPLICIT_UPWIND_ADVECTIVE_FLUX__(
                        tmDot * p_U[SPATIAL_DIM * nearestNode + i]);
                    PHASIC_IP_ZERO_NORMAL_STRESS__();
                }
                else // inflow: entrain along the opening direction
                {
                    scalar num = 0.0;
                    scalar den = 0.0;
                    for (label i = 0; i < SPATIAL_DIM; ++i)
                    {
                        // Opening directions point into the domain. Without a
                        // phase-specific direction field, use the inward face
                        // normal.
                        const scalar d =
                            dir ? dir[SPATIAL_DIM * ip + i] : -p_nx[i];
                        num += p_uBip[i] * p_nx[i];
                        den += d * p_nx[i];
                    }
                    const scalar Iu_ipI = num / (den + SMALL);
                    PHASIC_IP_EXPLICIT_ADVECTIVE_FLUX__(
                        tmDot * Iu_ipI *
                        (dir ? dir[SPATIAL_DIM * ip + i] : -p_nx[i]));
                    PHASIC_IP_ZERO_NORMAL_STRESS_NO_GRADU_TRANSPOSE__();
                }
            }

            this->applyCoeff_(
                A, b, connectedNodes, scratchIds, scratchVals, rhs, lhs);
        }
    }
}

} // namespace accel
