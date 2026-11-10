// File       : solidMechanicsModel.cpp
// Created    : Thu Dec 04 2025 08:42:10 (+0100)
// Author     : Mhamad Mahdi Alloush
// Description:
// Copyright 2025 CCFNUM HSLU T&A. All Rights Reserved.

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

    // thermal expansion coefficient only for domains with thermal stress
    for (const auto& domain : realm->simulationRef().domainVector())
    {
        if (domain->type() == domainType::solid &&
            domain->solidMechanics_.thermalStress_)
        {
            alphaTRef();
            break;
        }
    }

    // create post-processing field instances
    stressRef();
    stressVonMisesRef();
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

        // raw solid correction X^k (un-relaxed FSI residual)
        if (solidCorrectionRef().isZoneUnset(domain->index()))
        {
            solidCorrectionRef().setZone(domain->index());
            initialCondition::setupFieldInitializationOverDomainFromValues(
                solidCorrectionRef(),
                domain->index(),
                std::vector<scalar>(SPATIAL_DIM, 0.0));
        }

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

                                // pressure_type: dead (default) or follower
                                bool followerPressure = false;
                                if (node["pressure_type"])
                                {
                                    std::string pt =
                                        node["pressure_type"]
                                            .template as<std::string>();
                                    if (pt == "follower")
                                        followerPressure = true;
                                    else if (pt == "dead")
                                        followerPressure = false;
                                    else
                                        errorMsg("pressure_type must be "
                                                 "'follower' or 'dead'");
                                }
                                bc.setConstantValue<1>(
                                    "follower_pressure",
                                    {followerPressure ? 1.0 : 0.0});

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
    // Interface
    for (const interface* interf : domain->interfacesRef())
    {
        updateDisplacementInterfaceSideFieldTraction_(
            domain, interf->interfaceSideInfoPtr(domain->index()));
    }

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

    // current (exposed_area_vector) and reference (original_...) area vectors
    const STKScalarField* exposedAreaVecSTKFieldPtr =
        metaData.get_field<scalar>(metaData.side_rank(),
                                   mesh::exposed_area_vector_ID);
    const STKScalarField* origExposedAreaVecSTKFieldPtr =
        metaData.get_field<scalar>(metaData.side_rank(),
                                   mesh::original_exposed_area_vector_ID);

    // update traction on fluid side: need to negate the fluxes
    stk::mesh::Selector selSides =
        metaData.universal_part() &
        stk::mesh::selectUnion(interfaceSideInfoPtr->opposingPartVec_);

    // follower load p*a/|a0| (default) or dead load p*a0/|a0|
    const bool useDeformedArea =
        interfaceSideInfoPtr->interfPtr()->useDeformedArea();

    std::vector<scalar> deformedAreaVec;
    if (useDeformedArea)
    {
        this->meshRef().computeDeformedExposedAreaVector(
            interfaceSideInfoPtr->opposingPartVec_,
            *this->DRef().stkFieldPtr(),
            deformedAreaVec);
    }
    label ipOffset = 0;

    // shifted ip's for field?
    const bool isPShifted = this->pRef().isShifted();

    stk::mesh::BucketVector const& sideBuckets =
        bulkData.get_buckets(metaData.side_rank(), selSides);

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

            // this side's entries of deformedAreaVec
            const scalar* deformedAreaVecForSide =
                useDeformedArea ? &deformedAreaVec[ipOffset] : nullptr;

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

            // current area vector (see above)
            const scalar* areaVec =
                stk::mesh::field_data(*exposedAreaVecSTKFieldPtr, side);

            // reference area vector; the current one for a fixed fluid zone
            const scalar* origAreaVec =
                origExposedAreaVecSTKFieldPtr
                    ? stk::mesh::field_data(*origExposedAreaVecSTKFieldPtr,
                                            side)
                    : nullptr;
            if (!origAreaVec)
            {
                origAreaVec = areaVec;
            }

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
                // reference area dA0: the solid integrates the traction over it
                scalar aMag0 = 0.0;
                for (label j = 0; j < SPATIAL_DIM; ++j)
                {
                    const scalar axj = origAreaVec[SPATIAL_DIM * ip + j];
                    aMag0 += axj * axj;
                }

                // zero-area ip: no force (aMag0 is squared, hence VSMALL)
                if (aMag0 < VSMALL)
                {
                    for (label j = 0; j < SPATIAL_DIM; ++j)
                    {
                        tractionValues[SPATIAL_DIM * ip + j] = 0.0;
                    }
                    continue;
                }

                const scalar invAMag0 = 1.0 / std::sqrt(aMag0);

                // interpolate pressure to bip
                scalar pBip = 0.0;
                for (label ic = 0; ic < nodesPerSide; ++ic)
                {
                    const scalar r = p_face_shape_function[offSetSF_face + ic];
                    pBip += r * p_p[ic];
                }

                for (label j = 0; j < SPATIAL_DIM; ++j)
                {
                    if (useDeformedArea)
                    {
                        // follower load, Nanson pullback: p * a_current / dA0
                        tractionValues[SPATIAL_DIM * ip + j] =
                            pBip *
                            deformedAreaVecForSide[SPATIAL_DIM * ip + j] *
                            invAMag0;
                    }
                    else
                    {
                        // dead load: p * a_0 / dA0 = p * n_0
                        tractionValues[SPATIAL_DIM * ip + j] =
                            pBip * origAreaVec[SPATIAL_DIM * ip + j] * invAMag0;
                    }
                }

                scalar aMag_current = 0.0;
                if (useDeformedArea)
                {
                    for (label j = 0; j < SPATIAL_DIM; ++j)
                    {
                        const scalar axj =
                            deformedAreaVecForSide[SPATIAL_DIM * ip + j];
                        aMag_current += axj * axj;
                    }
                    aMag_current = std::sqrt(aMag_current);
                }
                else
                {
                    aMag_current = 1.0 / invAMag0; // dA0: dead load
                }

                // shear stored as tau*|a|/dA0, as the pressure term
                for (label j = 0; j < SPATIAL_DIM; ++j)
                {
                    tractionValues[SPATIAL_DIM * ip + j] +=
                        wallShearStressBip[SPATIAL_DIM * ip + j] *
                        aMag_current * invAMag0;
                }
            }

            // Advance offset to next side's deformed area entries
            ipOffset += numScsBip * SPATIAL_DIM;
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

            // traction stored per reference area: force = traction * |a0|
            const scalar* origAreaVec =
                origExposedAreaVecSTKFieldPtr
                    ? stk::mesh::field_data(*origExposedAreaVecSTKFieldPtr,
                                            side)
                    : nullptr;
            if (!origAreaVec)
            {
                origAreaVec = areaVec;
            }

            for (label ip = 0; ip < numScsBip; ++ip)
            {
                // Reference area magnitude for this ip
                scalar aMag = 0.0;
                for (label j = 0; j < SPATIAL_DIM; ++j)
                {
                    const scalar axj = origAreaVec[SPATIAL_DIM * ip + j];
                    aMag += axj * axj;
                }
                aMag = std::sqrt(aMag);

                for (label j = 0; j < SPATIAL_DIM; ++j)
                {
                    // Force = traction * reference area
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

            // traction stored per reference area: force = traction * |a0|
            const scalar* origAreaVec =
                origExposedAreaVecSTKFieldPtr
                    ? stk::mesh::field_data(*origExposedAreaVecSTKFieldPtr,
                                            side)
                    : nullptr;
            if (!origAreaVec)
            {
                origAreaVec = areaVec;
            }

            for (label ip = 0; ip < numScsBip; ++ip)
            {
                // Reference area magnitude for this ip
                scalar aMag = 0.0;
                for (label j = 0; j < SPATIAL_DIM; ++j)
                {
                    const scalar axj = origAreaVec[SPATIAL_DIM * ip + j];
                    aMag += axj * axj;
                }
                aMag = std::sqrt(aMag);

                for (label j = 0; j < SPATIAL_DIM; ++j)
                {
                    // Force = traction * reference area
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

    // guard: callers without follower_pressure get dead load
    auto& followerPressureData = bc.template data<1>("follower_pressure");
    const bool followerPressure =
        (followerPressureData.type() == inputDataType::constant) &&
        (*followerPressureData.value() > 0.5);

    // Determine which contributions are active
    const bool pressureExpression =
        pressureData.type() == inputDataType::expression;
    const bool shearExpression = shearData.type() == inputDataType::expression;
    const bool hasPressure =
        (pressureData.type() == inputDataType::constant ||
         pressureData.type() == inputDataType::timeTable || pressureExpression);
    const bool hasShear =
        (shearData.type() == inputDataType::constant ||
         shearData.type() == inputDataType::timeTable || shearExpression);

    // Check for unsupported input types
    if (pressureData.type() == inputDataType::profileData)
    {
        errorMsg("profileData not provided yet for pressure");
    }
    if (shearData.type() == inputDataType::profileData)
    {
        errorMsg("profileData not provided yet for shear");
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

    // expressions in t, x, y, z, evaluated at the side integration points
    typedef exprtk::symbol_table<scalar> symbol_table_t;
    typedef exprtk::expression<scalar> expression_t;
    typedef exprtk::parser<scalar> parser_t;

    scalar t = this->meshRef().controlsRef().time;
    scalar x = 0.0, y = 0.0, z = 0.0;
    symbol_table_t symbolTable;
    symbolTable.add_constants();
    symbolTable.add_variable("t", t);
    symbolTable.add_variable("x", x);
    symbolTable.add_variable("y", y);
    symbolTable.add_variable("z", z);

    parser_t parser;
    expression_t pressureExpr;
    std::vector<expression_t> shearExpr(shearExpression ? SPATIAL_DIM : 0);
    if (pressureExpression)
    {
        pressureExpr.register_symbol_table(symbolTable);
        if (!parser.compile(pressureData.expression()[0], pressureExpr))
        {
            errorMsg("Error in the expression provided for pressure: " +
                     pressureData.expression()[0]);
        }
    }
    for (label i = 0; i < static_cast<label>(shearExpr.size()); ++i)
    {
        shearExpr[i].register_symbol_table(symbolTable);
        if (!parser.compile(shearData.expression()[i], shearExpr[i]))
        {
            errorMsg("Error in the expression provided for shear: " +
                     shearData.expression()[i]);
        }
    }
    const bool anyExpression = pressureExpression || shearExpression;
    const STKScalarField* coordinatesSTKFieldPtr = metaData.get_field<scalar>(
        stk::topology::NODE_RANK, this->getCoordinatesID_(domain));
    std::vector<scalar> faceShapeFunction;

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

    // follower pressure: area vectors on the current configuration
    std::vector<scalar> deformedAreaVec;
    if (followerPressure && hasPressure)
    {
        this->meshRef().computeDeformedExposedAreaVector(
            boundary->parts(), *this->DRef().stkFieldPtr(), deformedAreaVec);
    }
    label ipOffset = 0;

    for (const stk::mesh::Bucket* bucket : sideBuckets)
    {
        MasterElement* meFC =
            MasterElementRepo::get_surface_master_element(bucket->topology());
        const label numScsBip = meFC->numIntPoints_;
        const label nodesPerSide = meFC->nodesPerElement_;

        if (anyExpression)
        {
            faceShapeFunction.resize(numScsBip * nodesPerSide);
            if (this->DRef().isShifted())
            {
                meFC->shifted_shape_fcn(faceShapeFunction.data());
            }
            else
            {
                meFC->shape_fcn(faceShapeFunction.data());
            }
        }

        for (const stk::mesh::Entity side : *bucket)
        {
            // Get traction values at integration points
            scalar* tractionValues =
                stk::mesh::field_data(sideFluxSTKFieldRef, side);

            // Get area vector if needed for pressure
            const scalar* orgAreaVec =
                stk::mesh::field_data(*originalExposedAreaVecSTKFieldPtr, side);
            const scalar* deformedAreaVecForSide =
                (followerPressure && hasPressure) ? &deformedAreaVec[ipOffset]
                                                  : nullptr;

            const stk::mesh::Entity* sideNodes = bulkData.begin_nodes(side);

            // Loop over integration points
            for (label ip = 0; ip < numScsBip; ++ip)
            {
                // expression values at the integration point coordinates
                if (anyExpression)
                {
                    scalar ipCoordinates[3] = {0.0, 0.0, 0.0};
                    for (label n = 0; n < nodesPerSide; ++n)
                    {
                        const scalar* coords = stk::mesh::field_data(
                            *coordinatesSTKFieldPtr, sideNodes[n]);
                        for (label j = 0; j < SPATIAL_DIM; ++j)
                        {
                            ipCoordinates[j] +=
                                faceShapeFunction[ip * nodesPerSide + n] *
                                coords[j];
                        }
                    }
                    x = ipCoordinates[0];
                    y = ipCoordinates[1];
                    z = ipCoordinates[2];

                    if (pressureExpression)
                    {
                        pressureValue = pressureExpr.value();
                    }
                    for (label j = 0; j < static_cast<label>(shearExpr.size());
                         ++j)
                    {
                        shearValue[j] = shearExpr[j].value();
                    }
                }

                // reference area magnitude dA0
                scalar aMag0 = 0.0;
                for (label j = 0; j < SPATIAL_DIM; ++j)
                {
                    const scalar axj = orgAreaVec[SPATIAL_DIM * ip + j];
                    aMag0 += axj * axj;
                }

                // Inverse of area magnitude for normalization
                const scalar invAMag0 = 1.0 / std::sqrt(aMag0);

                if (followerPressure)
                {
                    // follower: Nanson pullback t0 = -p * a_current / dA0
                    for (label j = 0; j < SPATIAL_DIM; ++j)
                    {
                        tractionValues[SPATIAL_DIM * ip + j] =
                            -pressureValue *
                            deformedAreaVecForSide[SPATIAL_DIM * ip + j] *
                            invAMag0;
                    }
                }
                else
                {
                    // dead load (default): traction = -p * n0
                    for (label j = 0; j < SPATIAL_DIM; ++j)
                    {
                        tractionValues[SPATIAL_DIM * ip + j] =
                            -pressureValue * orgAreaVec[SPATIAL_DIM * ip + j] *
                            invAMag0;
                    }
                }

                // Accumulate shear contribution (tangential traction)
                // Note: If only shear is present, this adds to existing values
                for (label j = 0; j < SPATIAL_DIM; ++j)
                {
                    tractionValues[SPATIAL_DIM * ip + j] += shearValue[j];
                }
            }

            ipOffset += numScsBip * SPATIAL_DIM;
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
    STKScalarField* sigmaVonMisesSTKFieldPtr =
        stressVonMisesRef().stkFieldPtr();
    const STKScalarField* DSTKFieldPtr = DRef().stkFieldPtr();
    const STKScalarField* ESTKFieldPtr = ERef().stkFieldPtr();
    const STKScalarField* nuSTKFieldPtr = nuRef().stkFieldPtr();

    // Get geometric fields
    const auto& coordinatesRef = *metaData.get_field<scalar>(
        stk::topology::NODE_RANK, this->getCoordinatesID_(domain));

    // Check if plane stress or plane strain
    const bool planeStress = domain->solidMechanics_.planeStress_;

    // thermal stress -beta*(T - T_ref)*I (linear elastic only)
    const bool thermalStress = domain->solidMechanics_.thermalStress_;
    const scalar TRefValue = domain->solidMechanics_.referenceTemperature_;
    const STKScalarField* TSTKFieldPtr =
        thermalStress ? TRef().stkFieldPtr() : nullptr;
    const STKScalarField* alphaTSTKFieldPtr =
        thermalStress ? alphaTRef().stkFieldPtr() : nullptr;

    // FEM: post-processed stress/strain must match the solve's model
    bool useNeoHookeanStress = false;
    bool useMooneyRivlinStress = false;
    bool useSaintVenantKirchhoffStress = false;

    // c1/c2/kappa are per-domain constants, fetched once
    scalar mooneyC1 = 0.0;
    scalar mooneyC2 = 0.0;
    scalar mooneyKappa = 0.0;

    if (!this->controlsRef().isCvfemSolidMechanics())
    {
        useNeoHookeanStress = (domain->solidMechanics_.option_ ==
                               solidMechanicsOption::neoHookean);

        // same for the Flory-split Mooney-Rivlin model
        useMooneyRivlinStress = (domain->solidMechanics_.option_ ==
                                 solidMechanicsOption::modifiedMooneyRivlin);

        useSaintVenantKirchhoffStress =
            (domain->solidMechanics_.option_ ==
             solidMechanicsOption::saintVenantKirchhoff);

        if (useMooneyRivlinStress)
        {
            const auto& mechProps = domain->materialRef().mechanicalProperties_;
            mooneyC1 = mechProps.c1_;
            mooneyC2 = mechProps.c2_;
            mooneyKappa = mechProps.kappa_;
        }
    }

    // Get interior parts
    const stk::mesh::PartVector& partVec = domain->zonePtr()->interiorParts();
    stk::mesh::Selector selAllElements =
        metaData.universal_part() & stk::mesh::selectUnion(partVec);

    stk::mesh::BucketVector const& elemBuckets =
        bulkData.get_buckets(stk::topology::ELEMENT_RANK, selAllElements);

    // nodal recovery of the element values: per node [stress, strain,
    // sigma_zz]; sigma_zz is carried in the von Mises field until the end
    constexpr label nTensor = SPATIAL_DIM * SPATIAL_DIM;
    constexpr label nQ = 2 * nTensor + 1;
    const auto nodeQ = [&](stk::mesh::Entity node, label q) -> scalar&
    {
        if (q < nTensor)
        {
            return stk::mesh::field_data(*sigmaSTKFieldPtr, node)[q];
        }
        if (q < 2 * nTensor)
        {
            return stk::mesh::field_data(*epsilonSTKFieldPtr,
                                         node)[q - nTensor];
        }
        return *stk::mesh::field_data(*sigmaVonMisesSTKFieldPtr, node);
    };

    // boundary nodes: nodes of sides with one element in this domain
    std::vector<char> isBoundaryNode(bulkData.get_size_of_entity_index_space(),
                                     0);
    for (const stk::mesh::Bucket* sideBucket :
         bulkData.get_buckets(metaData.side_rank(), metaData.universal_part()))
    {
        for (stk::mesh::Entity side : *sideBucket)
        {
            label nDomainElements = 0;
            const stk::mesh::Entity* sideElements =
                bulkData.begin_elements(side);
            for (unsigned ie = 0; ie < bulkData.num_elements(side); ++ie)
            {
                if (selAllElements(bulkData.bucket(sideElements[ie])))
                {
                    nDomainElements++;
                }
            }

            if (nDomainElements == 1)
            {
                const stk::mesh::Entity* sideNodes = bulkData.begin_nodes(side);
                for (unsigned in = 0; in < bulkData.num_nodes(side); ++in)
                {
                    isBoundaryNode[sideNodes[in].local_offset()] = 1;
                }
            }
        }
    }

    // zero the nodal accumulators of this domain
    std::vector<scalar> nodeWeight(bulkData.get_size_of_entity_index_space(),
                                   0.0);
    stk::mesh::BucketVector const& nodeBuckets =
        bulkData.get_buckets(stk::topology::NODE_RANK, selAllElements);
    for (const stk::mesh::Bucket* nodeBucket : nodeBuckets)
    {
        for (stk::mesh::Entity node : *nodeBucket)
        {
            for (label q = 0; q < nQ; ++q)
            {
                nodeQ(node, q) = 0.0;
            }
        }
    }

    // elements with boundary nodes, kept for the boundary extrapolation
    struct boundaryElement
    {
        std::vector<stk::mesh::Entity> nodes;
        std::vector<scalar> coordinates; // node coordinates
        std::vector<scalar> centroid;
        std::vector<scalar> gradOp; // dN/dx at the centre (SCV mean)
        std::vector<scalar> Q;      // element values seen by each node
        std::vector<scalar> weight; // SCV volume of each node
    };

    std::vector<boundaryElement> boundaryElements;

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
    std::vector<scalar> ws_T;
    std::vector<scalar> ws_alphaT;
    std::vector<scalar> ws_shape_function;
    std::vector<scalar> ws_scv_volume;
    std::vector<scalar> ws_node_volume;

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
        ws_T.resize(nodesPerElement);
        ws_alphaT.resize(nodesPerElement);
        ws_shape_function.resize(numScvIp * nodesPerElement);
        ws_scv_volume.resize(numScvIp);
        ws_node_volume.resize(nodesPerElement);
        const label* scvIpNodeMap = meSCV->ipNodeMap();

        // Create pointers
        scalar* p_coordinates = &ws_coordinates[0];
        scalar* p_disp = &ws_disp[0];
        scalar* p_E = &ws_E[0];
        scalar* p_nu = &ws_nu[0];
        scalar* p_dndx = &ws_dndx[0];
        scalar* p_grad_u = &ws_grad_u[0];
        scalar* p_strain = &ws_strain[0];
        scalar* p_stress = &ws_stress[0];
        scalar* p_T = &ws_T[0];
        scalar* p_alphaT = &ws_alphaT[0];
        scalar* p_shape_function = &ws_shape_function[0];

        // shape functions for the element-mean temperature (thermal stress)
        if (thermalStress)
        {
            meSCV->shape_fcn(p_shape_function);
        }

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

                if (thermalStress)
                {
                    p_T[ni] = *stk::mesh::field_data(*TSTKFieldPtr, node);
                    p_alphaT[ni] =
                        *stk::mesh::field_data(*alphaTSTKFieldPtr, node);
                }
            }

            // Compute gradients
            scalar scs_error = 0.0;
            meSCV->grad_op(1,
                           &ws_coordinates[0],
                           &ws_dndx[0],
                           &ws_deriv[0],
                           &ws_det_j[0],
                           &scs_error);

            // SCV volume of each node: weight of the nodal recovery
            meSCV->determinant(
                1, &ws_coordinates[0], &ws_scv_volume[0], &scs_error);
            std::fill(ws_node_volume.begin(), ws_node_volume.end(), 0.0);
            for (label ip = 0; ip < numScvIp; ++ip)
            {
                ws_node_volume[scvIpNodeMap[ip]] += ws_scv_volume[ip];
            }

            // record elements touching the boundary
            boundaryElement* record = nullptr;
            for (label ni = 0; ni < nodesPerElement; ++ni)
            {
                if (isBoundaryNode[elemNodes[ni].local_offset()])
                {
                    boundaryElements.emplace_back();
                    record = &boundaryElements.back();
                    break;
                }
            }
            if (record)
            {
                record->nodes.assign(elemNodes, elemNodes + nodesPerElement);
                record->coordinates.assign(ws_coordinates.begin(),
                                           ws_coordinates.begin() +
                                               nodesPerElement * SPATIAL_DIM);
                record->centroid.assign(SPATIAL_DIM, 0.0);
                record->gradOp.assign(nodesPerElement * SPATIAL_DIM, 0.0);
                record->Q.assign(nodesPerElement * nQ, 0.0);
                record->weight.assign(ws_node_volume.begin(),
                                      ws_node_volume.end());
                for (label ic = 0; ic < nodesPerElement; ++ic)
                {
                    for (label j = 0; j < SPATIAL_DIM; ++j)
                    {
                        record->centroid[j] +=
                            p_coordinates[ic * SPATIAL_DIM + j] /
                            nodesPerElement;
                        for (label ip = 0; ip < numScvIp; ++ip)
                        {
                            record->gradOp[ic * SPATIAL_DIM + j] +=
                                p_dndx[SPATIAL_DIM * nodesPerElement * ip +
                                       ic * SPATIAL_DIM + j] /
                                numScvIp;
                        }
                    }
                }
            }

            // element-mean temperature, consistent with the mean strain
            scalar TMean = 0.0;
            if (thermalStress)
            {
                for (label ip = 0; ip < numScvIp; ++ip)
                {
                    for (label ic = 0; ic < nodesPerElement; ++ic)
                    {
                        TMean += p_shape_function[ip * nodesPerElement + ic] *
                                 p_T[ic];
                    }
                }
                TMean /= numScvIp;
            }

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

                // strain (and stress), averaged over integration points
                for (label i = 0; i < SPATIAL_DIM * SPATIAL_DIM; ++i)
                {
                    p_strain[i] = 0.0;
                    p_stress[i] = 0.0;
                }

                // plane-strain sigma_zz (F_zz = 1) for von Mises; 0 in 3D
                scalar sigmaZZ = 0.0;

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

                    if (useNeoHookeanStress)
                    {
                        // Deformation gradient: F = I + grad_X(u)
                        scalar F[SPATIAL_DIM * SPATIAL_DIM];
                        for (label i = 0; i < SPATIAL_DIM; ++i)
                        {
                            for (label j = 0; j < SPATIAL_DIM; ++j)
                            {
                                const scalar delta_ij = (i == j) ? 1.0 : 0.0;
                                F[i * SPATIAL_DIM + j] =
                                    delta_ij + p_grad_u[i * SPATIAL_DIM + j];
                            }
                        }

                        // J = det(F)
#if SPATIAL_DIM == 2
                        const scalar J = F[0] * F[3] - F[1] * F[2];
#elif SPATIAL_DIM == 3
                        const scalar J = F[0] * (F[4] * F[8] - F[5] * F[7]) -
                                         F[1] * (F[3] * F[8] - F[5] * F[6]) +
                                         F[2] * (F[3] * F[7] - F[4] * F[6]);
#endif
                        const scalar lnJ = std::log(J);

#if SPATIAL_DIM == 2
                        // plane strain: b_zz = 1
                        sigmaZZ += (lambda / J) * lnJ / numScvIp;
#endif

                        // Cauchy (mu/J)(b-I)+(lambda/J)ln(J)I; E=0.5(F^T F-I)
                        for (label i = 0; i < SPATIAL_DIM; ++i)
                        {
                            for (label j = 0; j < SPATIAL_DIM; ++j)
                            {
                                const scalar delta_ij = (i == j) ? 1.0 : 0.0;

                                scalar b_ij = 0.0;
                                scalar c_ij = 0.0; // (F^T * F)_ij
                                for (label kk = 0; kk < SPATIAL_DIM; ++kk)
                                {
                                    b_ij += F[i * SPATIAL_DIM + kk] *
                                            F[j * SPATIAL_DIM + kk];
                                    c_ij += F[kk * SPATIAL_DIM + i] *
                                            F[kk * SPATIAL_DIM + j];
                                }

                                p_stress[i * SPATIAL_DIM + j] +=
                                    ((mu / J) * (b_ij - delta_ij) +
                                     (lambda / J) * lnJ * delta_ij) /
                                    numScvIp;

                                p_strain[i * SPATIAL_DIM + j] +=
                                    0.5 * (c_ij - delta_ij) / numScvIp;
                            }
                        }
                    }
                    else if (useMooneyRivlinStress)
                    {
                        // Deformation gradient: F = I + grad_X(u)
                        scalar F[SPATIAL_DIM * SPATIAL_DIM];
                        for (label i = 0; i < SPATIAL_DIM; ++i)
                        {
                            for (label j = 0; j < SPATIAL_DIM; ++j)
                            {
                                const scalar delta_ij = (i == j) ? 1.0 : 0.0;
                                F[i * SPATIAL_DIM + j] =
                                    delta_ij + p_grad_u[i * SPATIAL_DIM + j];
                            }
                        }

                        // J = det(F)
#if SPATIAL_DIM == 2
                        const scalar J = F[0] * F[3] - F[1] * F[2];
#elif SPATIAL_DIM == 3
                        const scalar J = F[0] * (F[4] * F[8] - F[5] * F[7]) -
                                         F[1] * (F[3] * F[8] - F[5] * F[6]) +
                                         F[2] * (F[3] * F[7] - F[4] * F[6]);
#endif
                        const scalar lnJ = std::log(J);
                        const scalar Jm23 = std::pow(J, -2.0 / 3.0);
                        const scalar Jm43 = std::pow(J, -4.0 / 3.0);

                        // Left Cauchy-Green tensor b = F * F^T; I1 = trace(b).
                        scalar b[SPATIAL_DIM * SPATIAL_DIM];
                        scalar I1 = 0.0;
                        for (label i = 0; i < SPATIAL_DIM; ++i)
                        {
                            for (label j = 0; j < SPATIAL_DIM; ++j)
                            {
                                scalar b_ij = 0.0;
                                for (label kk = 0; kk < SPATIAL_DIM; ++kk)
                                {
                                    b_ij += F[i * SPATIAL_DIM + kk] *
                                            F[j * SPATIAL_DIM + kk];
                                }
                                b[i * SPATIAL_DIM + j] = b_ij;
                            }
                            I1 += b[i * SPATIAL_DIM + i];
                        }
                        // plane strain: add b_zz = 1 to I1
#if SPATIAL_DIM == 2
                        I1 += 1.0;
#endif
                        // I2 is folded into dev(I1*b - b^2) below

                        // tr(b^2); plane strain adds b_zz^2 = 1
                        scalar b2[SPATIAL_DIM * SPATIAL_DIM];
                        scalar trb2 = 0.0;
                        for (label i = 0; i < SPATIAL_DIM; ++i)
                        {
                            for (label j = 0; j < SPATIAL_DIM; ++j)
                            {
                                scalar b2_ij = 0.0;
                                for (label kk = 0; kk < SPATIAL_DIM; ++kk)
                                {
                                    b2_ij += b[i * SPATIAL_DIM + kk] *
                                             b[kk * SPATIAL_DIM + j];
                                }
                                b2[i * SPATIAL_DIM + j] = b2_ij;
                            }
                            trb2 += b2[i * SPATIAL_DIM + i];
                        }
#if SPATIAL_DIM == 2
                        trb2 += 1.0;
#endif

                        // Flory-split Cauchy stress (physical d = 3)
                        const scalar trX = I1 * I1 - trb2;

#if SPATIAL_DIM == 2
                        // plane strain: b_zz = (b^2)_zz = 1
                        sigmaZZ +=
                            ((2.0 / J) *
                                 (mooneyC1 * Jm23 * (1.0 - I1 / 3.0) +
                                  mooneyC2 * Jm43 * ((I1 - 1.0) - trX / 3.0)) +
                             (mooneyKappa / J) * lnJ) /
                            numScvIp;
#endif
                        for (label i = 0; i < SPATIAL_DIM; ++i)
                        {
                            for (label j = 0; j < SPATIAL_DIM; ++j)
                            {
                                const scalar delta_ij = (i == j) ? 1.0 : 0.0;

                                const scalar devB_ij = b[i * SPATIAL_DIM + j] -
                                                       (I1 / 3.0) * delta_ij;

                                const scalar X_ij =
                                    I1 * b[i * SPATIAL_DIM + j] -
                                    b2[i * SPATIAL_DIM + j];
                                const scalar devX_ij =
                                    X_ij - (trX / 3.0) * delta_ij;

                                p_stress[i * SPATIAL_DIM + j] +=
                                    ((2.0 / J) * (mooneyC1 * Jm23 * devB_ij +
                                                  mooneyC2 * Jm43 * devX_ij) +
                                     (mooneyKappa / J) * lnJ * delta_ij) /
                                    numScvIp;

                                scalar c_ij = 0.0; // (F^T * F)_ij
                                for (label kk = 0; kk < SPATIAL_DIM; ++kk)
                                {
                                    c_ij += F[kk * SPATIAL_DIM + i] *
                                            F[kk * SPATIAL_DIM + j];
                                }
                                p_strain[i * SPATIAL_DIM + j] +=
                                    0.5 * (c_ij - delta_ij) / numScvIp;
                            }
                        }
                    }
                    else if (useSaintVenantKirchhoffStress)
                    {
                        // Deformation gradient: F = I + grad_X(u)
                        scalar F[SPATIAL_DIM * SPATIAL_DIM];
                        for (label i = 0; i < SPATIAL_DIM; ++i)
                        {
                            for (label j = 0; j < SPATIAL_DIM; ++j)
                            {
                                const scalar delta_ij = (i == j) ? 1.0 : 0.0;
                                F[i * SPATIAL_DIM + j] =
                                    delta_ij + p_grad_u[i * SPATIAL_DIM + j];
                            }
                        }

                        // J = det(F)
#if SPATIAL_DIM == 2
                        const scalar J = F[0] * F[3] - F[1] * F[2];
#elif SPATIAL_DIM == 3
                        const scalar J = F[0] * (F[4] * F[8] - F[5] * F[7]) -
                                         F[1] * (F[3] * F[8] - F[5] * F[6]) +
                                         F[2] * (F[3] * F[7] - F[4] * F[6]);
#endif

                        // Green-Lagrange E = 0.5(F^T F - I)
                        scalar Egl[SPATIAL_DIM * SPATIAL_DIM];
                        scalar trE = 0.0;
                        for (label i = 0; i < SPATIAL_DIM; ++i)
                        {
                            for (label j = 0; j < SPATIAL_DIM; ++j)
                            {
                                const scalar delta_ij = (i == j) ? 1.0 : 0.0;
                                scalar c_ij = 0.0;
                                for (label kk = 0; kk < SPATIAL_DIM; ++kk)
                                {
                                    c_ij += F[kk * SPATIAL_DIM + i] *
                                            F[kk * SPATIAL_DIM + j];
                                }
                                Egl[i * SPATIAL_DIM + j] =
                                    0.5 * (c_ij - delta_ij);
                            }
                            trE += Egl[i * SPATIAL_DIM + i];
                        }

#if SPATIAL_DIM == 2
                        // plane strain: E_zz = 0, S_zz = lambda tr(E), F_zz = 1
                        sigmaZZ += lambda * trE / J / numScvIp;
#endif

                        // 2nd PK S = lambda tr(E) I + 2 mu E
                        scalar S[SPATIAL_DIM * SPATIAL_DIM];
                        for (label i = 0; i < SPATIAL_DIM; ++i)
                        {
                            for (label j = 0; j < SPATIAL_DIM; ++j)
                            {
                                const scalar delta_ij = (i == j) ? 1.0 : 0.0;
                                S[i * SPATIAL_DIM + j] =
                                    lambda * trE * delta_ij +
                                    2.0 * mu * Egl[i * SPATIAL_DIM + j];
                            }
                        }

                        // Cauchy sigma = (1/J) F S F^T
                        for (label i = 0; i < SPATIAL_DIM; ++i)
                        {
                            for (label j = 0; j < SPATIAL_DIM; ++j)
                            {
                                scalar sigma_ij = 0.0;
                                for (label k = 0; k < SPATIAL_DIM; ++k)
                                {
                                    for (label l = 0; l < SPATIAL_DIM; ++l)
                                    {
                                        sigma_ij += F[i * SPATIAL_DIM + k] *
                                                    S[k * SPATIAL_DIM + l] *
                                                    F[j * SPATIAL_DIM + l];
                                    }
                                }

                                p_stress[i * SPATIAL_DIM + j] +=
                                    sigma_ij / J / numScvIp;

                                p_strain[i * SPATIAL_DIM + j] +=
                                    Egl[i * SPATIAL_DIM + j] / numScvIp;
                            }
                        }
                    }
                    else
                    {
                        // infinitesimal strain: eps = 0.5*(grad u + grad u^T)
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
                }

                if (!useNeoHookeanStress && !useMooneyRivlinStress &&
                    !useSaintVenantKirchhoffStress)
                {
                    // linear elastic: sigma = lambda*tr(eps)*I + 2*mu*eps
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

#if SPATIAL_DIM == 2
                    // plane strain: eps_zz = 0
                    sigmaZZ = lambda * trace_strain;
#endif

                    if (thermalStress)
                    {
                        const scalar betaDT =
                            thermalStressCoefficient(
                                E, nu, p_alphaT[ni], planeStress) *
                            (TMean - TRefValue);

                        for (label i = 0; i < SPATIAL_DIM; ++i)
                        {
                            p_stress[i * SPATIAL_DIM + i] -= betaDT;
                        }
#if SPATIAL_DIM == 2
                        if (!planeStress)
                        {
                            sigmaZZ -= betaDT;
                        }
#endif
                    }
                }

                // accumulate the element value, weighted by the SCV volume
                const scalar w = ws_node_volume[ni];
                for (label i = 0; i < nTensor; ++i)
                {
                    nodeQ(node, i) += w * p_stress[i];
                    nodeQ(node, nTensor + i) += w * p_strain[i];
                }
                nodeQ(node, 2 * nTensor) += w * sigmaZZ;
                nodeWeight[node.local_offset()] += w;

                if (record)
                {
                    scalar* Q = &record->Q[ni * nQ];
                    for (label i = 0; i < nTensor; ++i)
                    {
                        Q[i] = p_stress[i];
                        Q[nTensor + i] = p_strain[i];
                    }
                    Q[2 * nTensor] = sigmaZZ;
                }
            }
        }
    }

    // interior nodes: SCV-weighted average of the element values
    for (const stk::mesh::Bucket* nodeBucket : nodeBuckets)
    {
        for (stk::mesh::Entity node : *nodeBucket)
        {
            const scalar w = nodeWeight[node.local_offset()];
            if (w > 0.0)
            {
                for (label q = 0; q < nQ; ++q)
                {
                    nodeQ(node, q) /= w;
                }
            }
        }
    }

    const std::vector<const stk::mesh::FieldBase*> recoveredFields = {
        sigmaSTKFieldPtr, epsilonSTKFieldPtr, sigmaVonMisesSTKFieldPtr};
    if (messager::parallel())
    {
        stk::mesh::communicate_field_data(bulkData, recoveredFields);
    }

    // boundary nodes: mean of sigma_e + grad(sigma).(x - x_e) over their
    // elements, grad from the recovered field, iterated to a fixed point
    constexpr label maxBoundaryIterations = 100;
    constexpr scalar boundaryTolerance = 1.0e-10;

    // compact index of the owned/shared boundary nodes (aura nodes are
    // filled by the owner through the communication)
    std::vector<label> boundaryIndex(bulkData.get_size_of_entity_index_space(),
                                     -1);
    std::vector<stk::mesh::Entity> boundaryNodes;
    for (const stk::mesh::Bucket* nodeBucket : nodeBuckets)
    {
        if (!nodeBucket->owned() && !nodeBucket->shared())
        {
            continue;
        }
        for (stk::mesh::Entity node : *nodeBucket)
        {
            if (isBoundaryNode[node.local_offset()])
            {
                boundaryIndex[node.local_offset()] =
                    static_cast<label>(boundaryNodes.size());
                boundaryNodes.push_back(node);
            }
        }
    }

    std::vector<scalar> boundarySum(boundaryNodes.size() * nQ);
    std::vector<scalar> boundaryWeight(boundaryNodes.size());
    std::vector<scalar> gradQ(SPATIAL_DIM * nQ);
    for (label iter = 0; iter < maxBoundaryIterations; ++iter)
    {
        std::fill(boundarySum.begin(), boundarySum.end(), 0.0);
        std::fill(boundaryWeight.begin(), boundaryWeight.end(), 0.0);

        for (const boundaryElement& element : boundaryElements)
        {
            const label nNodes = static_cast<label>(element.nodes.size());

            std::fill(gradQ.begin(), gradQ.end(), 0.0);
            for (label ic = 0; ic < nNodes; ++ic)
            {
                for (label q = 0; q < nQ; ++q)
                {
                    const scalar value = nodeQ(element.nodes[ic], q);
                    for (label j = 0; j < SPATIAL_DIM; ++j)
                    {
                        gradQ[j * nQ + q] +=
                            element.gradOp[ic * SPATIAL_DIM + j] * value;
                    }
                }
            }

            for (label ic = 0; ic < nNodes; ++ic)
            {
                const label ib =
                    boundaryIndex[element.nodes[ic].local_offset()];
                if (ib < 0)
                {
                    continue;
                }

                const scalar w = element.weight[ic];
                for (label q = 0; q < nQ; ++q)
                {
                    scalar value = element.Q[ic * nQ + q];
                    for (label j = 0; j < SPATIAL_DIM; ++j)
                    {
                        value += gradQ[j * nQ + q] *
                                 (element.coordinates[ic * SPATIAL_DIM + j] -
                                  element.centroid[j]);
                    }
                    boundarySum[ib * nQ + q] += w * value;
                }
                boundaryWeight[ib] += w;
            }
        }

        // update; change relative to the largest stress component
        scalar maxChange = 0.0;
        scalar maxStress = 0.0;
        for (std::size_t ib = 0; ib < boundaryNodes.size(); ++ib)
        {
            const stk::mesh::Entity node = boundaryNodes[ib];
            for (label q = 0; q < nQ; ++q)
            {
                const scalar value =
                    boundarySum[ib * nQ + q] / boundaryWeight[ib];
                if (q < nTensor)
                {
                    maxChange =
                        std::max(maxChange, std::abs(value - nodeQ(node, q)));
                    maxStress = std::max(maxStress, std::abs(value));
                }
                nodeQ(node, q) = value;
            }
        }

        if (messager::parallel())
        {
            stk::mesh::communicate_field_data(bulkData, recoveredFields);
        }

        messager::maxReduce(maxChange);
        messager::maxReduce(maxStress);
        if (maxChange <= boundaryTolerance * std::max(maxStress, SMALL))
        {
            break;
        }
    }

    // von Mises sqrt(3/2 dev:dev) of the recovered 3D stress (s4f sigmaEq)
    for (const stk::mesh::Bucket* nodeBucket : nodeBuckets)
    {
        for (stk::mesh::Entity node : *nodeBucket)
        {
            const scalar* stress =
                stk::mesh::field_data(*sigmaSTKFieldPtr, node);
#if SPATIAL_DIM == 2
            const scalar sigmaZZ = nodeQ(node, 2 * nTensor);
#else
            const scalar sigmaZZ = 0.0;
#endif

            scalar trSigma = sigmaZZ;
            for (label i = 0; i < SPATIAL_DIM; ++i)
            {
                trSigma += stress[i * SPATIAL_DIM + i];
            }
            const scalar sigmaMean = trSigma / 3.0;

            scalar devSqr = 0.0;
            for (label i = 0; i < SPATIAL_DIM; ++i)
            {
                for (label j = 0; j < SPATIAL_DIM; ++j)
                {
                    const scalar dev_ij = stress[i * SPATIAL_DIM + j] -
                                          ((i == j) ? sigmaMean : 0.0);
                    devSqr += dev_ij * dev_ij;
                }
            }
#if SPATIAL_DIM == 2
            devSqr += (sigmaZZ - sigmaMean) * (sigmaZZ - sigmaMean);
#endif

            nodeQ(node, 2 * nTensor) = std::sqrt(1.5 * devSqr);
        }
    }
}

} /* namespace accel */
