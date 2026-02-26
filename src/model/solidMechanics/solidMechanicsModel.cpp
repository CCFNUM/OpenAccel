// File       : solidMechanicsModel.cpp
// Created    : Thu Dec 04 2025 08:42:10 (+0100)
// Author     : Mhamad Mahdi Alloush
// Description:
// Copyright (c) 2025 CCFNUM, Lucerne University of Applied Sciences and Arts.
// SPDX-License-Identifier: BSD-3-Clause

#include "solidMechanicsModel.h"
#include "domain.h"
#include "initialConditions.h"
#include "realm.h"
#include "simulation.h"

namespace accel
{

solidMechanicsModel::solidMechanicsModel(realm* realm) : model(realm)
{
    // create field instances
    DRef();
    rhoRef();
    ERef();
    nuRef();

    // create post-processing field instances
    stressRef();
    strainRef();
}

void solidMechanicsModel::setupDisplacement(
    const std::shared_ptr<domain> domain)
{
    if (DRef().isZoneUnset(domain->index()))
    {
        DRef().setZone(domain->index());

        // set initial conditions from input
        initialCondition::setupFieldInitializationOverDomainFromInput(
            DRef(), realm::D_ID, domain);

#ifdef HAS_INTERFACE
        // In the case of fluid-structure interaction, that is, the presence of
        // fluid domains adjacent to the current solid domain (provided that
        // this solid domain is naturally deforming), therefore, the fluid-solid
        // interface must store traction information coming from fluid
        for (interface* interf : domain->interfacesRef())
        {
            if (interf->isInternal())
            {
                // do nothing
            }
            else
            {
                if (interf->isFluidSolidType())
                {
                    // this is a fluid-structure interface, and displacement
                    // flux field must be defined on this side
                    DRef().registerSideFluxFieldForInterfaceSide(
                        interf->index(), interf->isMasterZone(domain->index()));

                    if (!domain->zonePtr()->meshDeforming())
                    {
                        // if no mesh deformation is set (simple physics), then
                        // the side flux field on the fluid side must be created
                        // here, because the displacement diffusion equation on
                        // the fluid zone is not created
                        DRef().registerSideFluxFieldForInterfaceSide(
                            interf->index(),
                            !interf->isMasterZone(domain->index()));
                    }
                }
                else
                {
                    // do nothing
                }
            }
        }
#endif /* HAS_INTERFACE */

        // boundary conditions for this domain
        setupBoundaryConditions_(
            domain,
            // anonymous function to set boundary conditions for this model
            [this](const ::accel::domain* domain,
                   const label iBoundary,
                   const boundaryPhysicalType bc_type,
                   const YAML::Node boundaryDetailsNode)
        {
            auto& bc = DRef().boundaryConditionRef(domain->index(), iBoundary);

#ifndef NDEBUG
            if (messager::master())
            {
                // clang-format off
		std::cout << "Setting boundary conditions:\n";
		std::cout << "\tdomain name: " << domain->name() << "\n";
		std::cout << "\tdomain index: " << domain->index() << "\n";
		std::cout << "\tpatch index: " << iBoundary << "\n";
		std::cout << "\tBoundary type: " << ::accel::toString(bc_type) << "\n";
		std::cout << "\tYAML values:\n" << boundaryDetailsNode << "\n\n";
                // clang-format on
            }
#endif /* NDEBUG */

            switch (bc_type)
            {
                case boundaryPhysicalType::wall:
                    {
                        if (boundaryDetailsNode["solid_mechanics"])
                        {
                            const auto& node =
                                boundaryDetailsNode["solid_mechanics"];
                            std::string option =
                                node["option"].template as<std::string>();

                            if (option == "fixed")
                            {
                                bc.setType(
                                    boundaryConditionType::specifiedValue);
                                bc.setConstantValue<SPATIAL_DIM>(
                                    "value",
                                    std::vector<scalar>(SPATIAL_DIM, 0));
                                DRef().registerSideFields(domain->index(),
                                                          iBoundary);
                            }
                            else if (option == "prescribed")
                            {
                                // Check if user specified which directions to
                                // prescribe
                                if (node["fixed_directions"])
                                {
                                    // Partial prescription: fix only specified
                                    // directions, zero-gradient on others
                                    // (roller-like BC)
                                    auto fixedDirs =
                                        node["fixed_directions"]
                                            .as<std::vector<std::string>>();

                                    // Parse displacement values
                                    std::vector<scalar> dispValues(SPATIAL_DIM,
                                                                   0.0);
                                    if (node["displacement"])
                                    {
                                        dispValues =
                                            node["displacement"]
                                                .as<std::vector<scalar>>();
                                    }

                                    // Create component flags
                                    std::vector<scalar> fvFlag(SPATIAL_DIM,
                                                               0.0);
                                    std::vector<std::string> dirNames = {
                                        "x", "y", "z"};

                                    // Mark which directions are fixed
                                    for (const auto& dir : fixedDirs)
                                    {
                                        for (label i = 0; i < SPATIAL_DIM; ++i)
                                        {
                                            if (dir == dirNames[i])
                                            {
                                                fvFlag[i] = 1.0; // Fixed
                                                break;
                                            }
                                        }
                                    }

                                    // Check if all directions are fixed
                                    bool allFixed = true;
                                    for (label i = 0; i < SPATIAL_DIM; ++i)
                                    {
                                        if (fvFlag[i] < 0.5)
                                        {
                                            allFixed = false;
                                            break;
                                        }
                                    }

                                    if (allFixed)
                                    {
                                        // All directions fixed - pure Dirichlet
                                        bc.setType(boundaryConditionType::
                                                       specifiedValue);
                                        bc.setConstantValue<SPATIAL_DIM>(
                                            "value", dispValues);
                                        DRef().registerSideFields(
                                            domain->index(), iBoundary);
                                    }
                                    else
                                    {
                                        // Mixed BC: some fixed, others
                                        // zero-gradient
                                        bc.setType(
                                            boundaryConditionType::mixed);
                                        bc.setConstantValue<SPATIAL_DIM>(
                                            "value", dispValues);
                                        bc.setConstantValue<SPATIAL_DIM>(
                                            "flux",
                                            std::vector<scalar>(SPATIAL_DIM,
                                                                0.0));
                                        bc.setConstantValue<SPATIAL_DIM>(
                                            "fixed_value_flag", fvFlag);
                                        DRef().registerSideFields(
                                            domain->index(), iBoundary);
                                        DRef().registerSideFluxField(
                                            domain->index(), iBoundary);
                                    }
                                }
                                else
                                {
                                    // Original behavior: fix all directions
                                    bc.setType(
                                        boundaryConditionType::specifiedValue);
                                    bc.query<SPATIAL_DIM>(
                                        node, "value", "displacement");
                                    DRef().registerSideFields(domain->index(),
                                                              iBoundary);
                                }
                            }
                            else if (option == "roller")
                            {
                                // Roller BC: constrain displacement in
                                // specified direction(s), free (zero-gradient)
                                // in other directions Common in structural
                                // mechanics for supports

                                if (!node["direction"])
                                {
                                    errorMsg("direction not specified for "
                                             "roller BC");
                                }

                                auto constrainedDir =
                                    node["direction"]
                                        .as<std::vector<std::string>>();

                                // Default: zero displacement in constrained
                                // directions
                                std::vector<scalar> dispValues(SPATIAL_DIM,
                                                               0.0);
                                if (node["displacement"])
                                {
                                    dispValues = node["displacement"]
                                                     .as<std::vector<scalar>>();
                                }

                                // Create component flags
                                std::vector<scalar> fvFlag(SPATIAL_DIM, 0.0);
                                std::vector<std::string> dirNames = {
                                    "x", "y", "z"};

                                // Mark constrained directions
                                for (const auto& dir : constrainedDir)
                                {
                                    for (label i = 0; i < SPATIAL_DIM; ++i)
                                    {
                                        if (dir == dirNames[i])
                                        {
                                            fvFlag[i] = 1.0; // Constrained
                                            break;
                                        }
                                    }
                                }

                                // Mixed BC: constrained in specified
                                // directions, free in others
                                bc.setType(boundaryConditionType::mixed);
                                bc.setConstantValue<SPATIAL_DIM>("value",
                                                                 dispValues);
                                bc.setConstantValue<SPATIAL_DIM>(
                                    "flux",
                                    std::vector<scalar>(SPATIAL_DIM, 0.0));
                                bc.setConstantValue<SPATIAL_DIM>(
                                    "fixed_value_flag", fvFlag);
                                DRef().registerSideFields(domain->index(),
                                                          iBoundary);
                                DRef().registerSideFluxField(domain->index(),
                                                             iBoundary);
                            }
                            else if (option == "traction")
                            {
                                bc.setType(
                                    boundaryConditionType::specifiedFlux);
                                bc.query<1>(node, "pressure", "pressure");
                                bc.query<SPATIAL_DIM>(node, "shear", "shear");
                                DRef().registerSideFluxField(domain->index(),
                                                             iBoundary);
                            }
                            else if (option == "mixed")
                            {
#if SPATIAL_DIM == 3
                                // Read component specifications
                                std::string xspec, yspec, zspec;
                                if (node["x_specification"])
                                {
                                    xspec = node["x_specification"]["option"]
                                                .template as<std::string>();
                                }
                                else
                                {
                                    errorMsg("x_specification not provided in "
                                             "solid_mechanics mixed BC");
                                }

                                if (node["y_specification"])
                                {
                                    yspec = node["y_specification"]["option"]
                                                .template as<std::string>();
                                }
                                else
                                {
                                    errorMsg("y_specification not provided in "
                                             "solid_mechanics mixed BC");
                                }

                                if (node["z_specification"])
                                {
                                    zspec = node["z_specification"]["option"]
                                                .template as<std::string>();
                                }
                                else
                                {
                                    errorMsg("z_specification not provided in "
                                             "solid_mechanics mixed BC");
                                }

                                // Component multiplier for mixed BC
                                // fvFlag[i] = 1 means fixed (Dirichlet)
                                // fvFlag[i] = 0 means flux (Neumann)
                                std::vector<scalar> fvFlag(SPATIAL_DIM, 0.0);

                                // Process each component
                                scalar xval, yval, zval;
                                scalar xflux, yflux, zflux;

                                // X-component
                                if (xspec == "fixed")
                                {
                                    xval = node["x_specification"]["value"]
                                               .template as<scalar>();
                                    fvFlag[0] = 1.0;
                                    xflux = 0.0; // dummy
                                }
                                else if (xspec == "fixed_flux")
                                {
                                    xflux = node["x_specification"]["value"]
                                                .template as<scalar>();
                                    xval = 0.0; // dummy
                                }
                                else if (xspec == "zero_flux")
                                {
                                    xflux = 0.0;
                                    xval = 0.0; // dummy
                                }
                                else
                                {
                                    errorMsg("invalid specification for x in "
                                             "solid_mechanics mixed BC");
                                }

                                // Y-component
                                if (yspec == "fixed")
                                {
                                    yval = node["y_specification"]["value"]
                                               .template as<scalar>();
                                    fvFlag[1] = 1.0;
                                    yflux = 0.0; // dummy
                                }
                                else if (yspec == "fixed_flux")
                                {
                                    yflux = node["y_specification"]["value"]
                                                .template as<scalar>();
                                    yval = 0.0; // dummy
                                }
                                else if (yspec == "zero_flux")
                                {
                                    yflux = 0.0;
                                    yval = 0.0; // dummy
                                }
                                else
                                {
                                    errorMsg("invalid specification for y in "
                                             "solid_mechanics mixed BC");
                                }

                                // Z-component
                                if (zspec == "fixed")
                                {
                                    zval = node["z_specification"]["value"]
                                               .template as<scalar>();
                                    fvFlag[2] = 1.0;
                                    zflux = 0.0; // dummy
                                }
                                else if (zspec == "fixed_flux")
                                {
                                    zflux = node["z_specification"]["value"]
                                                .template as<scalar>();
                                    zval = 0.0; // dummy
                                }
                                else if (zspec == "zero_flux")
                                {
                                    zflux = 0.0;
                                    zval = 0.0; // dummy
                                }
                                else
                                {
                                    errorMsg("invalid specification for z in "
                                             "solid_mechanics mixed BC");
                                }

                                // Determine BC type based on specifications
                                if (xspec == "fixed" && yspec == "fixed" &&
                                    zspec == "fixed")
                                {
                                    // All fixed - pure Dirichlet
                                    bc.setType(
                                        boundaryConditionType::specifiedValue);
                                    bc.setConstantValue<SPATIAL_DIM>(
                                        "value", {xval, yval, zval});
                                    DRef().registerSideFields(domain->index(),
                                                              iBoundary);
                                }
                                else if (xspec == "fixed_flux" &&
                                         yspec == "fixed_flux" &&
                                         zspec == "fixed_flux")
                                {
                                    // All flux - pure Neumann
                                    bc.setType(
                                        boundaryConditionType::specifiedFlux);
                                    bc.setConstantValue<SPATIAL_DIM>(
                                        "flux", {xflux, yflux, zflux});
                                    DRef().registerSideFluxField(
                                        domain->index(), iBoundary);
                                }
                                else if (xspec == "zero_flux" &&
                                         yspec == "zero_flux" &&
                                         zspec == "zero_flux")
                                {
                                    // All zero flux
                                    bc.setType(
                                        boundaryConditionType::zeroGradient);
                                }
                                else
                                {
                                    // Mixed BC - combination of Dirichlet and
                                    // Neumann
                                    bc.setType(boundaryConditionType::mixed);
                                    bc.setConstantValue<SPATIAL_DIM>(
                                        "value", {xval, yval, zval});
                                    bc.setConstantValue<SPATIAL_DIM>(
                                        "flux", {xflux, yflux, zflux});
                                    bc.setConstantValue<SPATIAL_DIM>(
                                        "fixed_value_flag", fvFlag);
                                    DRef().registerSideFields(domain->index(),
                                                              iBoundary);
                                    DRef().registerSideFluxField(
                                        domain->index(), iBoundary);
                                }
#elif SPATIAL_DIM == 2
                                // 2D version
                                std::string xspec, yspec;
                                if (node["x_specification"])
                                {
                                    xspec = node["x_specification"]["option"]
                                                .template as<std::string>();
                                }
                                else
                                {
                                    errorMsg("x_specification not provided in "
                                             "solid_mechanics mixed BC");
                                }

                                if (node["y_specification"])
                                {
                                    yspec = node["y_specification"]["option"]
                                                .template as<std::string>();
                                }
                                else
                                {
                                    errorMsg("y_specification not provided in "
                                             "solid_mechanics mixed BC");
                                }

                                std::vector<scalar> fvFlag(SPATIAL_DIM, 0.0);
                                scalar xval, yval;
                                scalar xflux, yflux;

                                // X-component
                                if (xspec == "fixed")
                                {
                                    xval = node["x_specification"]["value"]
                                               .template as<scalar>();
                                    fvFlag[0] = 1.0;
                                    xflux = 0.0;
                                }
                                else if (xspec == "fixed_flux")
                                {
                                    xflux = node["x_specification"]["value"]
                                                .template as<scalar>();
                                    xval = 0.0;
                                }
                                else if (xspec == "zero_flux")
                                {
                                    xflux = 0.0;
                                    xval = 0.0;
                                }
                                else
                                {
                                    errorMsg("invalid specification for x in "
                                             "solid_mechanics mixed BC");
                                }

                                // Y-component
                                if (yspec == "fixed")
                                {
                                    yval = node["y_specification"]["value"]
                                               .template as<scalar>();
                                    fvFlag[1] = 1.0;
                                    yflux = 0.0;
                                }
                                else if (yspec == "fixed_flux")
                                {
                                    yflux = node["y_specification"]["value"]
                                                .template as<scalar>();
                                    yval = 0.0;
                                }
                                else if (yspec == "zero_flux")
                                {
                                    yflux = 0.0;
                                    yval = 0.0;
                                }
                                else
                                {
                                    errorMsg("invalid specification for y in "
                                             "solid_mechanics mixed BC");
                                }

                                // Determine BC type
                                if (xspec == "fixed" && yspec == "fixed")
                                {
                                    bc.setType(
                                        boundaryConditionType::specifiedValue);
                                    bc.setConstantValue<SPATIAL_DIM>(
                                        "value", {xval, yval});
                                    DRef().registerSideFields(domain->index(),
                                                              iBoundary);
                                }
                                else if (xspec == "fixed_flux" &&
                                         yspec == "fixed_flux")
                                {
                                    bc.setType(
                                        boundaryConditionType::specifiedFlux);
                                    bc.setConstantValue<SPATIAL_DIM>(
                                        "flux", {xflux, yflux});
                                    DRef().registerSideFluxField(
                                        domain->index(), iBoundary);
                                }
                                else if (xspec == "zero_flux" &&
                                         yspec == "zero_flux")
                                {
                                    bc.setType(
                                        boundaryConditionType::zeroGradient);
                                }
                                else
                                {
                                    bc.setType(boundaryConditionType::mixed);
                                    bc.setConstantValue<SPATIAL_DIM>(
                                        "value", {xval, yval});
                                    bc.setConstantValue<SPATIAL_DIM>(
                                        "flux", {xflux, yflux});
                                    bc.setConstantValue<SPATIAL_DIM>(
                                        "fixed_value_flag", fvFlag);
                                    DRef().registerSideFields(domain->index(),
                                                              iBoundary);
                                    DRef().registerSideFluxField(
                                        domain->index(), iBoundary);
                                }
#endif
                            }
                            else
                            {
                                errorMsg(std::string("option for ") +
                                         "solid_mechanics" +
                                         std::string(" not provided at wall"));
                            }
                        }
                        else
                        {
                            bc.setType(boundaryConditionType::zeroGradient);
                        }
                    }
                    break;

                case boundaryPhysicalType::symmetry:
                    bc.setType(boundaryConditionType::symmetry);
                    break;

                default:
                    break;
            }
        });
    }
}

void solidMechanicsModel::updateDisplacement(
    const std::shared_ptr<domain> domain)
{
    // raw update
    fieldBroker::updateDisplacement(domain);

    // model-based side updates
    updateDisplacementSideFields_(domain);
}

void solidMechanicsModel::updateDisplacementSideFields_(
    const std::shared_ptr<domain> domain)
{
#ifdef HAS_INTERFACE
    // Interface
    for (const interface* interf : domain->interfacesRef())
    {
        updateDisplacementInterfaceSideFieldTraction_(
            domain, interf->interfaceSideInfoPtr(domain->index()));
    }
#endif /* HAS_INTERFACE */

    // Boundary
    for (label iBoundary = 0;
         iBoundary < this->meshRef().zonePtr(domain->index())->nBoundaries();
         iBoundary++)
    {
        const auto* boundary =
            this->meshRef().zonePtr(domain->index())->boundaryPtr(iBoundary);

        boundaryPhysicalType physicalType = boundary->type();
        const auto& bcType =
            DRef().boundaryConditionRef(domain->index(), iBoundary).type();

        switch (physicalType)
        {
            case boundaryPhysicalType::wall:
                {
                    switch (bcType)
                    {
                        case boundaryConditionType::specifiedFlux:
                            {
                                updateDisplacementBoundarySideFieldTraction_(
                                    domain, boundary);
                            }
                            break;

                        default:
                            break;
                    }
                }
                break;

            default:
                break;
        }
    }
}

#ifdef HAS_INTERFACE
void solidMechanicsModel::updateDisplacementInterfaceSideFieldTraction_(
    const std::shared_ptr<domain> domain,
    const interfaceSideInfo* interfaceSideInfoPtr)
{
    // update side flux field on the fluid side from pressure and wall
    // shear stress information, and interpolate that to the solid side
    // as traction
    if (!interfaceSideInfoPtr->interfPtr()->isFluidSolidType())
        return;

    // Consider only if a side field is defined on the boundary.
    if ((!this->DRef().sideFluxFieldRef().definedOn(
            interfaceSideInfoPtr->currentPartVec_)))
        errorMsg("Side flux field not defined on interface side");
    if ((!this->DRef().sideFluxFieldRef().definedOn(
            interfaceSideInfoPtr->opposingPartVec_)))
        errorMsg("Side flux field not defined on interface side");

    // Common mesh data references
    stk::mesh::BulkData& bulkData = this->meshRef().bulkDataRef();
    stk::mesh::MetaData& metaData = this->meshRef().metaDataRef();

    // nodal fields to gather
    std::vector<scalar> ws_p;

    // master element
    std::vector<scalar> ws_face_shape_function;

    // Get field references and mesh data (only once)
    const auto& pSTKFieldRef = this->pRef().stkFieldRef();
    const auto& wallShearStressSTKFieldRef =
        this->wallShearStressRef().stkFieldRef();
    auto& sideFluxSTKFieldRef = this->DRef().sideFluxFieldRef().stkFieldRef();

    // Area vector field needed only for pressure (normal traction). We always
    // use current area for fluid side assuming it is deforming
    const STKScalarField* exposedAreaVecSTKFieldPtr =
        metaData.get_field<scalar>(metaData.side_rank(),
                                   mesh::exposed_area_vector_ID);

    // update traction on fluid side: need to negate the fluxes
    stk::mesh::Selector selOwnedSides =
        metaData.locally_owned_part() &
        stk::mesh::selectUnion(interfaceSideInfoPtr->opposingPartVec_);

    // shifted ip's for field?
    const bool isPShifted = this->pRef().isShifted();

    stk::mesh::BucketVector const& sideBuckets =
        bulkData.get_buckets(metaData.side_rank(), selOwnedSides);

    for (const stk::mesh::Bucket* bucket : sideBuckets)
    {
        MasterElement* meFC =
            MasterElementRepo::get_surface_master_element(bucket->topology());
        const label nodesPerSide = meFC->nodesPerElement_;
        const label numScsBip = meFC->numIntPoints_;

        // algorithm related; element
        ws_p.resize(nodesPerSide);
        ws_face_shape_function.resize(numScsBip * nodesPerSide);

        // pointers
        scalar* p_p = &ws_p[0];
        scalar* p_face_shape_function = &ws_face_shape_function[0];

        // shape functions
        if (isPShifted)
        {
            meFC->shifted_shape_fcn(&p_face_shape_function[0]);
        }
        else
        {
            meFC->shape_fcn(&p_face_shape_function[0]);
        }

        for (const stk::mesh::Entity side : *bucket)
        {
            // face node relations
            stk::mesh::Entity const* sideNodeRels = bulkData.begin_nodes(side);

            //======================================
            // gather nodal data off of side
            //======================================
            for (label ni = 0; ni < nodesPerSide; ++ni)
            {
                stk::mesh::Entity node = sideNodeRels[ni];

                // gather scalars
                p_p[ni] = *stk::mesh::field_data(pSTKFieldRef, node);
            }

            // Get traction values at integration points
            scalar* tractionValues =
                stk::mesh::field_data(sideFluxSTKFieldRef, side);

            // Get area vector if needed for pressure
            const scalar* areaVec =
                stk::mesh::field_data(*exposedAreaVecSTKFieldPtr, side);

            const scalar* wallShearStressBip =
                stk::mesh::field_data(wallShearStressSTKFieldRef, side);

            // Loop over integration points
            for (label ip = 0; ip < numScsBip; ++ip)
            {
                // offsets
                const label offSetSF_face = ip * nodesPerSide;

                // Apply pressure contribution (normal force on solid)
                // Force on solid = +p * n, where n is the fluid's outward
                // normal (pointing into solid). Positive pressure causes
                // compression.
                // Compute magnitude of area vector to normalize
                scalar aMag = 0.0;
                for (label j = 0; j < SPATIAL_DIM; ++j)
                {
                    const scalar axj = areaVec[SPATIAL_DIM * ip + j];
                    aMag += axj * axj;
                }

                // Inverse of area magnitude for normalization
                const scalar invAMag = 1.0 / std::sqrt(aMag);

                // interpolate to bip
                scalar pBip = 0.0;
                for (label ic = 0; ic < nodesPerSide; ++ic)
                {
                    const scalar r = p_face_shape_function[offSetSF_face + ic];
                    pBip += r * p_p[ic];
                }

                // Set traction as pressure * unit normal
                for (label j = 0; j < SPATIAL_DIM; ++j)
                {
                    tractionValues[SPATIAL_DIM * ip + j] =
                        pBip * areaVec[SPATIAL_DIM * ip + j] * invAMag;
                }

                // Accumulate shear contribution (tangential traction)
                // Note: If only shear is present, this adds to existing
                // values
                for (label j = 0; j < SPATIAL_DIM; ++j)
                {
                    tractionValues[SPATIAL_DIM * ip + j] +=
                        wallShearStressBip[SPATIAL_DIM * ip + j];
                }
            }
        }
    }

    // interpolate data from fluid side to solid side
    this->DRef().sideFluxFieldRef().transfer(
        interfaceSideInfoPtr->interfPtr()->index(),
        interfaceSideInfoPtr->isMasterSide(),
        DRef().isShifted());

#ifndef NDEBUG
    // ========================================================================
    // Debug: Compute force statistics on FLUID side
    // ========================================================================
    std::array<scalar, SPATIAL_DIM> fluidForceMin;
    std::array<scalar, SPATIAL_DIM> fluidForceMax;
    std::array<scalar, SPATIAL_DIM> fluidForceSum;
    fluidForceMin.fill(std::numeric_limits<scalar>::max());
    fluidForceMax.fill(std::numeric_limits<scalar>::lowest());
    fluidForceSum.fill(0.0);
    label fluidTotalIps = 0;

    for (const stk::mesh::Bucket* bucket : sideBuckets)
    {
        MasterElement* meFC =
            MasterElementRepo::get_surface_master_element(bucket->topology());
        const label numScsBip = meFC->numIntPoints_;

        for (const stk::mesh::Entity side : *bucket)
        {
            const scalar* tractionValues =
                stk::mesh::field_data(sideFluxSTKFieldRef, side);
            const scalar* areaVec =
                stk::mesh::field_data(*exposedAreaVecSTKFieldPtr, side);

            for (label ip = 0; ip < numScsBip; ++ip)
            {
                // Compute area magnitude for this ip
                scalar aMag = 0.0;
                for (label j = 0; j < SPATIAL_DIM; ++j)
                {
                    const scalar axj = areaVec[SPATIAL_DIM * ip + j];
                    aMag += axj * axj;
                }
                aMag = std::sqrt(aMag);

                for (label j = 0; j < SPATIAL_DIM; ++j)
                {
                    // Force = traction * area
                    const scalar force =
                        tractionValues[SPATIAL_DIM * ip + j] * aMag;
                    fluidForceMin[j] = std::min(fluidForceMin[j], force);
                    fluidForceMax[j] = std::max(fluidForceMax[j], force);
                    fluidForceSum[j] += force;
                }
                fluidTotalIps++;
            }
        }
    }

    // MPI reduction for fluid side statistics
    std::array<scalar, SPATIAL_DIM> globalFluidForceMin;
    std::array<scalar, SPATIAL_DIM> globalFluidForceMax;
    std::array<scalar, SPATIAL_DIM> globalFluidForceSum;
    label globalFluidTotalIps = 0;

    stk::all_reduce_min(bulkData.parallel(),
                        fluidForceMin.data(),
                        globalFluidForceMin.data(),
                        SPATIAL_DIM);
    stk::all_reduce_max(bulkData.parallel(),
                        fluidForceMax.data(),
                        globalFluidForceMax.data(),
                        SPATIAL_DIM);
    stk::all_reduce_sum(bulkData.parallel(),
                        fluidForceSum.data(),
                        globalFluidForceSum.data(),
                        SPATIAL_DIM);
    stk::all_reduce_sum(
        bulkData.parallel(), &fluidTotalIps, &globalFluidTotalIps, 1);

    // ========================================================================
    // Debug: Compute force statistics on SOLID side
    // ========================================================================
    std::array<scalar, SPATIAL_DIM> solidForceMin;
    std::array<scalar, SPATIAL_DIM> solidForceMax;
    std::array<scalar, SPATIAL_DIM> solidForceSum;
    solidForceMin.fill(std::numeric_limits<scalar>::max());
    solidForceMax.fill(std::numeric_limits<scalar>::lowest());
    solidForceSum.fill(0.0);
    label solidTotalIps = 0;

    // Select solid side (current side)
    stk::mesh::Selector selOwnedSidesSolid =
        metaData.locally_owned_part() &
        stk::mesh::selectUnion(interfaceSideInfoPtr->currentPartVec_);

    // Area vector on solid side
    const STKScalarField* solidExposedAreaVecPtr = metaData.get_field<scalar>(
        metaData.side_rank(), mesh::exposed_area_vector_ID);

    stk::mesh::BucketVector const& solidSideBuckets =
        bulkData.get_buckets(metaData.side_rank(), selOwnedSidesSolid);

    for (const stk::mesh::Bucket* bucket : solidSideBuckets)
    {
        MasterElement* meFC =
            MasterElementRepo::get_surface_master_element(bucket->topology());
        const label numScsBip = meFC->numIntPoints_;

        for (const stk::mesh::Entity side : *bucket)
        {
            const scalar* tractionValues =
                stk::mesh::field_data(sideFluxSTKFieldRef, side);
            const scalar* areaVec =
                stk::mesh::field_data(*solidExposedAreaVecPtr, side);

            for (label ip = 0; ip < numScsBip; ++ip)
            {
                // Compute area magnitude for this ip
                scalar aMag = 0.0;
                for (label j = 0; j < SPATIAL_DIM; ++j)
                {
                    const scalar axj = areaVec[SPATIAL_DIM * ip + j];
                    aMag += axj * axj;
                }
                aMag = std::sqrt(aMag);

                for (label j = 0; j < SPATIAL_DIM; ++j)
                {
                    // Force = traction * area
                    const scalar force =
                        tractionValues[SPATIAL_DIM * ip + j] * aMag;
                    solidForceMin[j] = std::min(solidForceMin[j], force);
                    solidForceMax[j] = std::max(solidForceMax[j], force);
                    solidForceSum[j] += force;
                }
                solidTotalIps++;
            }
        }
    }

    // MPI reduction for solid side statistics
    std::array<scalar, SPATIAL_DIM> globalSolidForceMin;
    std::array<scalar, SPATIAL_DIM> globalSolidForceMax;
    std::array<scalar, SPATIAL_DIM> globalSolidForceSum;
    label globalSolidTotalIps = 0;

    stk::all_reduce_min(bulkData.parallel(),
                        solidForceMin.data(),
                        globalSolidForceMin.data(),
                        SPATIAL_DIM);
    stk::all_reduce_max(bulkData.parallel(),
                        solidForceMax.data(),
                        globalSolidForceMax.data(),
                        SPATIAL_DIM);
    stk::all_reduce_sum(bulkData.parallel(),
                        solidForceSum.data(),
                        globalSolidForceSum.data(),
                        SPATIAL_DIM);
    stk::all_reduce_sum(
        bulkData.parallel(), &solidTotalIps, &globalSolidTotalIps, 1);

    // Print debug information
    if (messager::master())
    {
        std::vector<std::string> dirNames = {"x", "y", "z"};

        // Save original flags to restore them later (good practice)
        std::ios_base::fmtflags oldFlags = std::cout.flags();

        std::cout << "\n===== FSI Force Transfer Debug (Interface: "
                  << interfaceSideInfoPtr->interfPtr()->name()
                  << ") =====" << std::endl;

        // Set scientific notation and high precision
        std::cout << std::scientific << std::setprecision(8);

        std::cout << "FLUID side (source):\n";
        std::cout << "  Total IPs: " << globalFluidTotalIps << "\n";
        for (label j = 0; j < SPATIAL_DIM; ++j)
        {
            scalar mean = (globalFluidTotalIps > 0)
                              ? globalFluidForceSum[j] / globalFluidTotalIps
                              : 0.0;
            std::cout << "  F_" << dirNames[j] << ": min=" << std::setw(15)
                      << globalFluidForceMin[j] << ", max=" << std::setw(15)
                      << globalFluidForceMax[j] << ", mean=" << std::setw(15)
                      << mean << ", total=" << std::setw(15)
                      << globalFluidForceSum[j] << "\n";
        }

        std::cout << "SOLID side (target, after transfer):\n";
        std::cout << "  Total IPs: " << globalSolidTotalIps << "\n";
        for (label j = 0; j < SPATIAL_DIM; ++j)
        {
            scalar mean = (globalSolidTotalIps > 0)
                              ? globalSolidForceSum[j] / globalSolidTotalIps
                              : 0.0;
            std::cout << "  F_" << dirNames[j] << ": min=" << std::setw(15)
                      << globalSolidForceMin[j] << ", max=" << std::setw(15)
                      << globalSolidForceMax[j] << ", mean=" << std::setw(15)
                      << mean << ", total=" << std::setw(15)
                      << globalSolidForceSum[j] << "\n";
        }

        // Conservation check
        std::cout << "Force conservation (solid - fluid):\n";
        for (label j = 0; j < SPATIAL_DIM; ++j)
        {
            scalar diff = globalSolidForceSum[j] - globalFluidForceSum[j];
            scalar relErr =
                (std::abs(globalFluidForceSum[j]) > 1e-12)
                    ? std::abs(diff / globalFluidForceSum[j]) * 100.0
                    : 0.0;

            // Use fixed for the percentage but keep scientific for the delta
            std::cout << "  Delta_F_" << dirNames[j] << ": " << std::scientific
                      << diff << " (rel. error: " << std::fixed
                      << std::setprecision(6) << relErr << "%)\n";
            std::cout << std::scientific
                      << std::setprecision(8); // Reset to scientific
        }
        std::cout << "=================================================\n"
                  << std::endl;

        // Restore original formatting
        std::cout.flags(oldFlags);
    }
#endif /* NDEBUG */
}
#endif /* HAS_INTERFACE */

void solidMechanicsModel::updateDisplacementBoundarySideFieldTraction_(
    const std::shared_ptr<domain> domain,
    const boundary* boundary)
{
    // Consider only if a side field is defined on the boundary.
    if ((!this->DRef().sideFluxFieldRef().definedOn(boundary->parts())))
        errorMsg("Side flux field not defined on boundary");

    // Common mesh data references
    stk::mesh::BulkData& bulkData = this->meshRef().bulkDataRef();
    stk::mesh::MetaData& metaData = this->meshRef().metaDataRef();

    auto& bc =
        this->DRef().boundaryConditionRef(domain->index(), boundary->index());

    // Get boundary condition data for both contributions
    auto& pressureData = bc.template data<1>("pressure");
    auto& shearData = bc.template data<SPATIAL_DIM>("shear");

    // Determine which contributions are active
    const bool hasPressure = (pressureData.type() == inputDataType::constant ||
                              pressureData.type() == inputDataType::timeTable);
    const bool hasShear = (shearData.type() == inputDataType::constant ||
                           shearData.type() == inputDataType::timeTable);

    // Check for unsupported input types
    if (pressureData.type() == inputDataType::expression ||
        pressureData.type() == inputDataType::profileData)
    {
        errorMsg("expression and profileData not provided yet for pressure");
    }
    if (shearData.type() == inputDataType::expression ||
        shearData.type() == inputDataType::profileData)
    {
        errorMsg("expression and profileData not provided yet for shear");
    }

    // Early return if no work needed
    if (!hasPressure && !hasShear)
        return;

    // Interpolate/retrieve values for active contributions
    scalar pressureValue = 0.0;
    std::array<scalar, SPATIAL_DIM> shearValue = {0.0};

    if (hasPressure)
    {
        if (pressureData.type() == inputDataType::constant)
        {
            pressureValue = *pressureData.value();
        }
        else // timeTable
        {
            auto inputValue =
                pressureData.interpolate(this->meshRef().controlsRef().time);
            pressureValue = inputValue[0];
        }
    }

    if (hasShear)
    {
        if (shearData.type() == inputDataType::constant)
        {
            std::copy(shearData.value(),
                      shearData.value() + SPATIAL_DIM,
                      shearValue.begin());
        }
        else // timeTable
        {
            shearValue =
                shearData.interpolate(this->meshRef().controlsRef().time);
        }
    }

    // Get field references and mesh data (only once)
    auto& sideFluxSTKFieldRef = this->DRef().sideFluxFieldRef().stkFieldRef();

    stk::mesh::Selector selAllSides =
        metaData.universal_part() & stk::mesh::selectUnion(boundary->parts());

    stk::mesh::BucketVector const& sideBuckets =
        bulkData.get_buckets(metaData.side_rank(), selAllSides);

    // Area vector field needed only for pressure (normal traction)
    const STKScalarField* exposedAreaVecSTKFieldPtr =
        metaData.get_field<scalar>(metaData.side_rank(),
                                   mesh::exposed_area_vector_ID);
    const STKScalarField* originalExposedAreaVecSTKFieldPtr =
        this->meshRef().controlsRef().isTransient()
            ? (metaData.get_field<scalar>(
                  metaData.side_rank(),
                  domain->solidMechanics_.formulation_ ==
                          kinematicFormulationType::updatedLagrangian
                      ? this->getExposedAreaVectorID_(domain)
                      : mesh::original_exposed_area_vector_ID))
            : exposedAreaVecSTKFieldPtr;

    for (const stk::mesh::Bucket* bucket : sideBuckets)
    {
        const MasterElement* meFC =
            MasterElementRepo::get_surface_master_element(bucket->topology());
        const label numScsBip = meFC->numIntPoints_;

        for (const stk::mesh::Entity side : *bucket)
        {
            // Get traction values at integration points
            scalar* tractionValues =
                stk::mesh::field_data(sideFluxSTKFieldRef, side);

            // Get area vector if needed for pressure
            const scalar* areaVec =
                stk::mesh::field_data(*exposedAreaVecSTKFieldPtr, side);
            const scalar* orgAreaVec =
                stk::mesh::field_data(*originalExposedAreaVecSTKFieldPtr, side);

            // Loop over integration points
            for (label ip = 0; ip < numScsBip; ++ip)
            {
                // Apply pressure contribution (normal traction)
                // Traction = -p * n, where n is the unit normal vector
                scalar aMag0 = 0.0;
                for (label j = 0; j < SPATIAL_DIM; ++j)
                {
                    const scalar axj = orgAreaVec[SPATIAL_DIM * ip + j];
                    aMag0 += axj * axj;
                }

                // Inverse of area magnitude for normalization
                const scalar invAMag0 = 1.0 / std::sqrt(aMag0);

                // Set traction as pressure * unit normal
                for (label j = 0; j < SPATIAL_DIM; ++j)
                {
                    tractionValues[SPATIAL_DIM * ip + j] =
                        -pressureValue * areaVec[SPATIAL_DIM * ip + j] *
                        invAMag0;
                }

                // Accumulate shear contribution (tangential traction)
                // Note: If only shear is present, this adds to existing values
                for (label j = 0; j < SPATIAL_DIM; ++j)
                {
                    tractionValues[SPATIAL_DIM * ip + j] += shearValue[j];
                }
            }
        }
    }
}

void solidMechanicsModel::updateStressAndStrain_(
    const std::shared_ptr<domain> domain)
{
    auto& mesh = this->meshRef();
    stk::mesh::BulkData& bulkData = mesh.bulkDataRef();
    stk::mesh::MetaData& metaData = mesh.metaDataRef();

    // Get field pointers
    STKScalarField* sigmaSTKFieldPtr = stressRef().stkFieldPtr();
    STKScalarField* epsilonSTKFieldPtr = strainRef().stkFieldPtr();
    const STKScalarField* DSTKFieldPtr = DRef().stkFieldPtr();
    const STKScalarField* ESTKFieldPtr = ERef().stkFieldPtr();
    const STKScalarField* nuSTKFieldPtr = nuRef().stkFieldPtr();

    // Get geometric fields
    const auto& coordinatesRef = *metaData.get_field<scalar>(
        stk::topology::NODE_RANK, this->getCoordinatesID_(domain));

    // Check if plane stress or plane strain
    const bool planeStress = domain->solidMechanics_.planeStress_;

    // Get interior parts
    const stk::mesh::PartVector& partVec = domain->zonePtr()->interiorParts();
    stk::mesh::Selector selAllElements =
        metaData.universal_part() & stk::mesh::selectUnion(partVec);

    stk::mesh::BucketVector const& elemBuckets =
        bulkData.get_buckets(stk::topology::ELEMENT_RANK, selAllElements);

    // Scratch arrays
    std::vector<scalar> ws_coordinates;
    std::vector<scalar> ws_disp;
    std::vector<scalar> ws_E;
    std::vector<scalar> ws_nu;
    std::vector<scalar> ws_dndx;
    std::vector<scalar> ws_deriv;
    std::vector<scalar> ws_det_j;
    std::vector<scalar> ws_grad_u;
    std::vector<scalar> ws_strain;
    std::vector<scalar> ws_stress;

    for (stk::mesh::BucketVector::const_iterator ib = elemBuckets.begin();
         ib != elemBuckets.end();
         ++ib)
    {
        stk::mesh::Bucket& elemBucket = **ib;
        const stk::mesh::Bucket::size_type length = elemBucket.size();

        MasterElement* meSCV =
            MasterElementRepo::get_volume_master_element(elemBucket.topology());
        const label nodesPerElement = meSCV->nodesPerElement_;
        const label numScvIp = meSCV->numIntPoints_;

        // Resize scratch arrays
        ws_coordinates.resize(nodesPerElement * SPATIAL_DIM);
        ws_disp.resize(nodesPerElement * SPATIAL_DIM);
        ws_E.resize(nodesPerElement);
        ws_nu.resize(nodesPerElement);
        ws_dndx.resize(SPATIAL_DIM * numScvIp * nodesPerElement);
        ws_deriv.resize(SPATIAL_DIM * numScvIp * nodesPerElement);
        ws_det_j.resize(numScvIp);
        ws_grad_u.resize(SPATIAL_DIM * SPATIAL_DIM);
        ws_strain.resize(SPATIAL_DIM * SPATIAL_DIM);
        ws_stress.resize(SPATIAL_DIM * SPATIAL_DIM);

        // Create pointers
        scalar* p_coordinates = &ws_coordinates[0];
        scalar* p_disp = &ws_disp[0];
        scalar* p_E = &ws_E[0];
        scalar* p_nu = &ws_nu[0];
        scalar* p_dndx = &ws_dndx[0];
        scalar* p_grad_u = &ws_grad_u[0];
        scalar* p_strain = &ws_strain[0];
        scalar* p_stress = &ws_stress[0];

        for (stk::mesh::Bucket::size_type k = 0; k < length; ++k)
        {
            stk::mesh::Entity const* elemNodes =
                bulkData.begin_nodes(elemBucket[k]);

            // Gather nodal data
            for (label ni = 0; ni < nodesPerElement; ++ni)
            {
                stk::mesh::Entity node = elemNodes[ni];

                const scalar* coords =
                    stk::mesh::field_data(coordinatesRef, node);
                const scalar* disp = stk::mesh::field_data(*DSTKFieldPtr, node);
                const scalar* E = stk::mesh::field_data(*ESTKFieldPtr, node);
                const scalar* nu = stk::mesh::field_data(*nuSTKFieldPtr, node);

                const label offSet = ni * SPATIAL_DIM;
                for (label j = 0; j < SPATIAL_DIM; ++j)
                {
                    p_coordinates[offSet + j] = coords[j];
                    p_disp[offSet + j] = disp[j];
                }
                p_E[ni] = E[0];
                p_nu[ni] = nu[0];
            }

            // Compute gradients
            scalar scs_error = 0.0;
            meSCV->grad_op(1,
                           &ws_coordinates[0],
                           &ws_dndx[0],
                           &ws_deriv[0],
                           &ws_det_j[0],
                           &scs_error);

            // Loop over nodes to compute stress and strain at each node
            for (label ni = 0; ni < nodesPerElement; ++ni)
            {
                stk::mesh::Entity node = elemNodes[ni];

                // Get material properties at this node
                const scalar E = p_E[ni];
                const scalar nu = p_nu[ni];

                // Compute Lame parameters
                const scalar mu = E / (2.0 * (1.0 + nu));
                scalar lambda;
                if (planeStress)
                {
                    lambda = nu * E / ((1.0 + nu) * (1.0 - nu));
                }
                else
                {
                    lambda = nu * E / ((1.0 + nu) * (1.0 - 2.0 * nu));
                }

                // Compute strain at this node (average over integration points)
                // Zero strain
                for (label i = 0; i < SPATIAL_DIM * SPATIAL_DIM; ++i)
                {
                    p_strain[i] = 0.0;
                }

                for (label ip = 0; ip < numScvIp; ++ip)
                {
                    // Compute displacement gradient at integration point
                    // Zero grad_u
                    for (label i = 0; i < SPATIAL_DIM * SPATIAL_DIM; ++i)
                    {
                        p_grad_u[i] = 0.0;
                    }

                    for (label ic = 0; ic < nodesPerElement; ++ic)
                    {
                        const label offSetDnDx =
                            SPATIAL_DIM * nodesPerElement * ip +
                            ic * SPATIAL_DIM;

                        for (label i = 0; i < SPATIAL_DIM; ++i)
                        {
                            const scalar ui = p_disp[ic * SPATIAL_DIM + i];
                            for (label j = 0; j < SPATIAL_DIM; ++j)
                            {
                                p_grad_u[i * SPATIAL_DIM + j] +=
                                    ui * p_dndx[offSetDnDx + j];
                            }
                        }
                    }

                    // Compute strain: ε_ij = 0.5 * (∂u_i/∂x_j + ∂u_j/∂x_i)
                    for (label i = 0; i < SPATIAL_DIM; ++i)
                    {
                        for (label j = 0; j < SPATIAL_DIM; ++j)
                        {
                            p_strain[i * SPATIAL_DIM + j] +=
                                0.5 *
                                (p_grad_u[i * SPATIAL_DIM + j] +
                                 p_grad_u[j * SPATIAL_DIM + i]) /
                                numScvIp;
                        }
                    }
                }

                // Compute stress: σ_ij = λ * ε_kk * δ_ij + 2μ * ε_ij
                // Zero stress
                for (label i = 0; i < SPATIAL_DIM * SPATIAL_DIM; ++i)
                {
                    p_stress[i] = 0.0;
                }

                scalar trace_strain = 0.0;
                for (label i = 0; i < SPATIAL_DIM; ++i)
                {
                    trace_strain += p_strain[i * SPATIAL_DIM + i];
                }

                for (label i = 0; i < SPATIAL_DIM; ++i)
                {
                    for (label j = 0; j < SPATIAL_DIM; ++j)
                    {
                        scalar delta_ij = (i == j) ? 1.0 : 0.0;
                        p_stress[i * SPATIAL_DIM + j] =
                            lambda * trace_strain * delta_ij +
                            2.0 * mu * p_strain[i * SPATIAL_DIM + j];
                    }
                }

                // Store stress and strain
                scalar* node_stress =
                    stk::mesh::field_data(*sigmaSTKFieldPtr, node);
                scalar* node_strain =
                    stk::mesh::field_data(*epsilonSTKFieldPtr, node);

                for (label i = 0; i < SPATIAL_DIM * SPATIAL_DIM; ++i)
                {
                    node_stress[i] = p_stress[i];
                    node_strain[i] = p_strain[i];
                }
            }
        }
    }
}

} /* namespace accel */
