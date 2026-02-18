// File : meshIO.cpp
// Created : Fri Aug 25 2023 12:55:24 (+0100)
// Author : Mhamad Mahdi Alloush
// Description:
// Copyright (c) 2023 CCFNUM, Lucerne University of Applied Sciences and Arts.
// SPDX-License-Identifier: BSD-3-Clause

// code
#include "boundary.h"
#include "controls.h"
#include "messager.h"
#include "types.h"
#include "zone.h"

namespace accel
{

// Operations

void mesh::read(const YAML::Node& inputNode)
{
    if (messager::master())
    {
        std::cout << "Reading mesh .." << std::endl;
    }

    // Create mesh from exodus file provided by yaml input file
    if (inputNode["mesh"])
    {
        const YAML::Node& mesh_node = inputNode["mesh"];
        fs::path meshFilePath =
            mesh_node["file_path"].template as<std::string>();

        // check for restart
        const auto& restart_ctrl =
            this->controlsRef().solverRef().restartControl_;
        if (restart_ctrl.isRestart_)
        {
            meshFilePath = restart_ctrl.inputFilePath_;
        }

        // automatic domain decomposition
        std::string auto_decomp("None");
        if (mesh_node["automatic_decomposition_type"])
        {
            auto_decomp = mesh_node["automatic_decomposition_type"]
                              .template as<std::string>();
        }
        if (messager::master() && messager::nProcs() > 1 &&
            "None" != auto_decomp)
        {
            std::cout << "Automatic domain decomposition: input Exodus file "
                         "must be a serial file\n";
        }

        stk::ParallelMachine pm = MPI_COMM_WORLD;

        // news for mesh constructs
        stk::mesh::MeshBuilder builder(pm);
        builder.set_aura_option(stk::mesh::BulkData::AUTO_AURA);

        // bulk data and io broker; meta data obtained from bulk data. We use
        // simple fields only that is flat arrays for multidimensional fields
        bulkDataPtr_ = builder.create();

        // Enable simple fields (flat storage)
        bulkDataPtr_->mesh_meta_data().use_simple_fields();

        // Enable declaration of late fields (after data population).
        // This is required for all side fields, which are declared
        // at initialization stage (which is after data population).
        bulkDataPtr_->mesh_meta_data().enable_late_fields();

        ioBrokerPtr_ = new stk::io::StkMeshIoBroker(pm);
        ioBrokerPtr_->set_bulk_data(bulkDataPtr_);

        if ("None" != auto_decomp)
        {
            ioBrokerPtr_->property_add(
                Ioss::Property("DECOMPOSITION_METHOD", auto_decomp));
        }

        // Initialize meta data (from exodus file)
        ioBrokerPtr_->add_mesh_database(meshFilePath,
                                        restart_ctrl.isRestart_
                                            ? stk::io::READ_RESTART
                                            : stk::io::READ_MESH);
        ioBrokerPtr_->create_input_mesh();

        if (restart_ctrl.isRestart_)
        {
            this->controlsRef().deserializeRestartParam(*ioBrokerPtr_);
        }

        // ensure mesh dimensions is consistent with SPATIAL_DIM
        assert(bulkDataPtr_->mesh_meta_data_ptr()->spatial_dimension() ==
               SPATIAL_DIM);

        // Validate YAML input against Exodus file parts
        validateYamlAgainstExodus_(inputNode);

        // check transformation details
        if (mesh_node["transformation"])
        {
            const auto& transformationNode = mesh_node["transformation"];

            if (transformationNode["scale"])
            {
                // set flag to ON
                scaleMesh_ = true;

                // get non-uniform scale
                const auto& scaleNode = transformationNode["scale"];

                // get scale vector
#if SPATIAL_DIM == 3
                scaleVector_[0] = scaleNode["Sx"].template as<scalar>();
                scaleVector_[1] = scaleNode["Sy"].template as<scalar>();
                scaleVector_[2] = scaleNode["Sz"].template as<scalar>();
#elif SPATIAL_DIM == 2
                scaleVector_[0] = scaleNode["Sx"].template as<scalar>();
                scaleVector_[1] = scaleNode["Sy"].template as<scalar>();
#endif
                // scale origin
                scaleOrigin_ =
                    scaleNode["scale_origin"]
                        .template as<std::array<scalar, SPATIAL_DIM>>();

                // export the database with original mesh
                if (scaleNode["export_original"])
                {
                    exportOriginal_ =
                        scaleNode["export_original"].template as<bool>();
                }
            }
        }

        // Check for element validation flag
        if (mesh_node["check_mesh"])
        {
            checkMesh_ = mesh_node["check_mesh"].template as<bool>();
        }

        // Check for element correction flag
        if (mesh_node["enable_correction"])
        {
            enableCorrection_ =
                mesh_node["enable_correction"].template as<bool>();
        }
    }
    else
    {
        errorMsg(
            "path to mesh not provided, i.e. file_path: <path-to-exodus-file>");
    }

    // Read and register zones
    if (inputNode["simulation"])
    {
        if (inputNode["simulation"]["physical_analysis"]["domains"])
        {
            const auto& domainsBlock =
                inputNode["simulation"]["physical_analysis"]["domains"];

            if (messager::master())
            {
                std::cout << "Registering zones" << std::endl;
            }

            for (label iDomain = 0;
                 iDomain < static_cast<label>(domainsBlock.size());
                 iDomain++)
            {
                // Create zone
                std::unique_ptr<zone> zone_ptr = std::make_unique<zone>(
                    this,
                    iDomain,
                    domainsBlock[iDomain]["name"].template as<std::string>());

                // Read zone details from yaml node
                zone_ptr->read(domainsBlock[iDomain]);

                // Turn boolean true for an existing moving zone frame
                if (!this->anyZoneFrameRotating() && zone_ptr->frameRotating())
                {
                    this->setAnyZoneFrameRotating(true);
                }

                // Turn boolean true for an existing moving zone mesh
                if (!this->anyZoneMeshTransforming() &&
                    zone_ptr->meshTransforming())
                {
                    this->setAnyZoneMeshTransforming(true);
                }

                // Turn boolean true for an existing deforming zone mesh
                if (!this->anyZoneMeshDeforming() && zone_ptr->meshDeforming())
                {
                    this->setAnyZoneMeshDeforming(true);
                }

                // Add to the global zones vector
                zoneVector_.push_back(std::move(zone_ptr));
            }

            if (messager::master())
            {
                std::cout << "Finished registering zones" << std::endl;
            }
        }
        else
        {
            errorMsg("simulation->physical_analysis->domains->sequence block "
                     "is not provided in the yaml input file");
        }
    }
    else
    {
        errorMsg("simulation block is not provided in the yaml input file");
    }

    // Collect interior parts: this includes all interior parts found in the
    // mesh database. The code currently assumes that the node-graph used in
    // assembly is based on all the mesh parts. Another use of this variable
    // is for defining fields that are always valid on the whole mesh like
    // volume field.
    // Collect also all boundary parts: this includes all sidesets
    // available in the mesh database.
    if (messager::master())
    {
        std::cout << "Registering mesh parts" << std::endl;
    }

    if (inputNode["simulation"])
    {
        if (inputNode["simulation"]["physical_analysis"])
        {
            // collect interior mesh parts from domains and the relevant
            // boundary mesh parts
            if (inputNode["simulation"]["physical_analysis"]["domains"])
            {
                const auto& domainsBlock =
                    inputNode["simulation"]["physical_analysis"]["domains"];

                for (const auto& domainBlock : domainsBlock)
                {
                    for (std::string elementBlockPartName :
                         domainBlock["location"]
                             .template as<std::vector<std::string>>())
                    {
                        // get part
                        const stk::mesh::Part* part =
                            metaDataRef().get_part(elementBlockPartName);

                        if (part->primary_entity_rank() ==
                            stk::mesh::EntityRank::ELEMENT_RANK)
                        {
                            // add part to vector of interior parts
                            interiorActiveParts_.push_back(part);

                            const auto& boundariesBlock =
                                domainBlock["boundaries"];

                            for (const auto& boundaryBlock : boundariesBlock)
                            {
                                for (std::string boundaryPartName :
                                     boundaryBlock["location"]
                                         .template as<
                                             std::vector<std::string>>())
                                {
                                    const stk::mesh::Part* boundaryPart =
                                        metaDataRef().get_part(
                                            boundaryPartName);

                                    // check validity
                                    if (boundaryPart->primary_entity_rank() !=
                                            metaDataRef().side_rank() ||
                                        boundaryPart->topology() !=
                                            stk::topology::INVALID_TOPOLOGY)
                                    {
                                        errorMsg("invalid boundary part " +
                                                 boundaryPartName);
                                    }

                                    // only add if not stored yet
                                    const auto it =
                                        std::find(boundaryActiveParts_.begin(),
                                                  boundaryActiveParts_.end(),
                                                  boundaryPart);
                                    if (it == boundaryActiveParts_.end())
                                    {
                                        boundaryActiveParts_.push_back(
                                            boundaryPart);

                                        // add wall parts
                                        if (boundaryBlock["type"]
                                                .template as<std::string>() ==
                                            "wall")
                                        {
                                            wallBoundaryActiveParts_.push_back(
                                                boundaryPart);
                                        }

                                        // add symmetry parts
                                        if (boundaryBlock["type"]
                                                .template as<std::string>() ==
                                            "symmetry")
                                        {
                                            symmetryBoundaryActiveParts_
                                                .push_back(boundaryPart);
                                        }
                                    }
                                }
                            }
                        }
                        else
                        {
                            errorMsg("Provided part " + elementBlockPartName +
                                     " is not an element block");
                        }
                    }
                }

#ifndef NDEBUG
                if (messager::master())
                {
                    const auto& allParts = metaDataRef().get_mesh_parts();
                    for (const auto* part : allParts)
                    {
                        if (part->primary_entity_rank() ==
                            stk::mesh::EntityRank::ELEMENT_RANK)
                        {
                            const auto it =
                                std::find(interiorActiveParts_.begin(),
                                          interiorActiveParts_.end(),
                                          part);
                            if (it == interiorActiveParts_.end())
                            {
                                warningMsg("Element block " + part->name() +
                                           " has not been used");
                            }
                        }
                    }
                }
#endif /* NDEBUG */
            }
        }
        else
        {
            errorMsg("simulation->physical_analysis block is not provided in "
                     "the yaml input file");
        }
    }
    else
    {
        errorMsg("simulation block is not provided in the yaml input file");
    }

    if (messager::master())
    {
        std::cout << "Finished registering mesh parts" << std::endl;
    }

    if (messager::master())
    {
        std::cout << "Finished reading mesh .." << std::endl << std::endl;
    }
}

void mesh::write(size_t resultsFileIndex, scalar writeTime)
{
    // restore original scale?
    if (scaleMesh_ && exportOriginal_)
    {
        std::array<scalar, SPATIAL_DIM> invScaleVector;
        for (label i = 0; i < SPATIAL_DIM; i++)
        {
            invScaleVector[i] = 1.0 / scaleVector_[i];
        }
        scale(invScaleVector);
    }

    this->ioBrokerPtr()->process_output_request(resultsFileIndex, writeTime);

    // scale again
    if (scaleMesh_ && exportOriginal_)
    {
        scale(scaleVector_);
    }
}

void mesh::validateYamlAgainstExodus_(const YAML::Node& inputNode)
{
    // Only master rank needs to perform validation
    if (!messager::master())
    {
        return;
    }

    std::cout << "Validating YAML input against Exodus file" << std::endl;

    // Get all available parts from the Exodus file
    const auto& allParts = metaDataRef().get_mesh_parts();

    // Helper: Check if a part name exists in the Exodus file (including
    // aliases)
    auto partExistsInExodus = [this,
                               &allParts](const std::string& partName) -> bool
    {
        // First check using the standard get_part method which should handle
        // aliases
        const stk::mesh::Part* part = metaDataRef().get_part(partName);
        if (part != nullptr)
        {
            return true;
        }

        // Additional check: look through IO database for alternative names
        // Get the input region from IO broker to check for aliases
        try
        {
            auto inputRegion = ioBrokerPtr_->get_input_ioss_region();
            if (inputRegion != nullptr)
            {
                // Check element blocks
                const auto& elementBlocks = inputRegion->get_element_blocks();
                for (const auto* block : elementBlocks)
                {
                    if (block->name() == partName)
                    {
                        return true;
                    }
                    // Check if this block has aliases
                    if (block->get_optional_property(
                            "original_name", std::string("")) == partName)
                    {
                        return true;
                    }
                }

                // Check sidesets
                const auto& sidesets = inputRegion->get_sidesets();
                for (const auto* sideset : sidesets)
                {
                    if (sideset->name() == partName)
                    {
                        return true;
                    }
                    // Check if this sideset has aliases
                    if (sideset->get_optional_property(
                            "original_name", std::string("")) == partName)
                    {
                        return true;
                    }
                }
            }
        }
        catch (...)
        {
            // If IO region access fails, fall back to basic check
        }

        // Fallback: iterate through all parts and check names
        for (const auto* meshPart : allParts)
        {
            if (meshPart->name() == partName)
            {
                return true;
            }
        }

        return false;
    };

    // Helper: Get all element blocks from Exodus file (including aliases)
    auto getAvailableElementBlocks = [this,
                                      &allParts]() -> std::vector<std::string>
    {
        std::vector<std::string> elementBlocks;

        // Add parts from STK mesh
        for (const auto* part : allParts)
        {
            if (part->primary_entity_rank() ==
                stk::mesh::EntityRank::ELEMENT_RANK)
            {
                elementBlocks.push_back(part->name());
            }
        }

        // Also add names from IO database (may include aliases)
        try
        {
            auto inputRegion = ioBrokerPtr_->get_input_ioss_region();
            if (inputRegion != nullptr)
            {
                const auto& ioElementBlocks = inputRegion->get_element_blocks();
                for (const auto* block : ioElementBlocks)
                {
                    // Add if not already in the list
                    if (std::find(elementBlocks.begin(),
                                  elementBlocks.end(),
                                  block->name()) == elementBlocks.end())
                    {
                        elementBlocks.push_back(block->name());
                    }

                    // Add original name if it exists and is different
                    std::string originalName = block->get_optional_property(
                        "original_name", std::string(""));
                    if (!originalName.empty() &&
                        originalName != block->name() &&
                        std::find(elementBlocks.begin(),
                                  elementBlocks.end(),
                                  originalName) == elementBlocks.end())
                    {
                        elementBlocks.push_back(originalName + " (alias for " +
                                                block->name() + ")");
                    }
                }
            }
        }
        catch (...)
        {
            // If IO region access fails, continue with STK parts only
        }

        return elementBlocks;
    };

    // Helper: Get all sidesets touching any of the specified element blocks
    // (including aliases)
    auto getSidesetsForElementBlocks =
        [this](const std::vector<std::string>& elementBlockNames)
        -> std::vector<std::string>
    {
        std::vector<std::string> sidesets;

        // Get sidesets from STK mesh parts
        for (const std::string& elementBlockName : elementBlockNames)
        {
            const stk::mesh::Part* elementBlock =
                metaDataRef().get_part(elementBlockName);
            if (elementBlock)
            {
                const auto& boundaryParts =
                    metaDataRef().get_surfaces_touched_by_block(elementBlock);
                for (const stk::mesh::Part* boundaryPart : boundaryParts)
                {
                    if (boundaryPart->primary_entity_rank() ==
                            metaDataRef().side_rank() &&
                        boundaryPart->topology() ==
                            stk::topology::INVALID_TOPOLOGY)
                    {
                        // Only add if not already in the list
                        if (std::find(sidesets.begin(),
                                      sidesets.end(),
                                      boundaryPart->name()) == sidesets.end())
                        {
                            sidesets.push_back(boundaryPart->name());
                        }
                    }
                }
            }
        }

        // Also add sideset names from IO database (may include aliases)
        try
        {
            auto inputRegion = ioBrokerPtr_->get_input_ioss_region();
            if (inputRegion != nullptr)
            {
                const auto& ioSidesets = inputRegion->get_sidesets();
                for (const auto* sideset : ioSidesets)
                {
                    // Add if not already in the list
                    if (std::find(sidesets.begin(),
                                  sidesets.end(),
                                  sideset->name()) == sidesets.end())
                    {
                        sidesets.push_back(sideset->name());
                    }

                    // Add original name if it exists and is different
                    std::string originalName = sideset->get_optional_property(
                        "original_name", std::string(""));
                    if (!originalName.empty() &&
                        originalName != sideset->name() &&
                        std::find(sidesets.begin(),
                                  sidesets.end(),
                                  originalName) == sidesets.end())
                    {
                        sidesets.push_back(originalName + " (alias for " +
                                           sideset->name() + ")");
                    }
                }
            }
        }
        catch (...)
        {
            // If IO region access fails, continue with STK parts only
        }

        return sidesets;
    };

    // Helper: Check if a sideset is attached to an element block (including
    // aliases)
    auto isSidesetAttachedToElementBlock =
        [this](const std::string& sidesetName,
               const std::string& elementBlockName) -> bool
    {
        const stk::mesh::Part* elementBlock =
            metaDataRef().get_part(elementBlockName);
        if (!elementBlock)
        {
            return false;
        }

        const auto& boundaryParts =
            metaDataRef().get_surfaces_touched_by_block(elementBlock);

        for (const stk::mesh::Part* boundaryPart : boundaryParts)
        {
            // Check primary name
            if (boundaryPart->name() == sidesetName)
            {
                return true;
            }

            // Check if the requested sideset name can be found via get_part
            // This handles aliases - if get_part can find it and it matches
            // this boundary part
            const stk::mesh::Part* requestedPart =
                metaDataRef().get_part(sidesetName);
            if (requestedPart && requestedPart == boundaryPart)
            {
                return true;
            }
        }

        // Additional check through IO database for case-insensitive or alias
        // matching
        try
        {
            auto inputRegion = ioBrokerPtr_->get_input_ioss_region();
            if (inputRegion != nullptr)
            {
                const auto& ioSidesets = inputRegion->get_sidesets();
                for (const auto* sideset : ioSidesets)
                {
                    // Check if this IO sideset matches the requested name
                    // (case-insensitive or alias)
                    if (sideset->name() == sidesetName ||
                        sideset->get_optional_property(
                            "original_name", std::string("")) == sidesetName)
                    {
                        // Now check if this IO sideset corresponds to any of
                        // the boundary parts
                        for (const stk::mesh::Part* boundaryPart :
                             boundaryParts)
                        {
                            if (boundaryPart->name() == sideset->name())
                            {
                                return true;
                            }
                        }
                    }
                }
            }
        }
        catch (...)
        {
            // If IO region access fails, continue with basic check
        }

        return false;
    };

    // Check domains
    if (inputNode["simulation"]["physical_analysis"]["domains"])
    {
        const auto& domainsBlock =
            inputNode["simulation"]["physical_analysis"]["domains"];

        for (const auto& domainBlock : domainsBlock)
        {
            const std::string domainName =
                domainBlock["name"].template as<std::string>();

            // Check domain locations exist in Exodus file
            for (const std::string& location :
                 domainBlock["location"]
                     .template as<std::vector<std::string>>())
            {
                if (!partExistsInExodus(location))
                {
                    std::string errorMessage =
                        "Domain '" + domainName + "': location '" + location +
                        "' does not exist in the Exodus file";

                    // Check if it's not an element block specifically
                    const stk::mesh::Part* part =
                        metaDataRef().get_part(location);
                    if (!part || part->primary_entity_rank() !=
                                     stk::mesh::EntityRank::ELEMENT_RANK)
                    {
                        auto availableElementBlocks =
                            getAvailableElementBlocks();
                        if (!availableElementBlocks.empty())
                        {
                            errorMessage +=
                                "\nAvailable element blocks in Exodus file:";
                            for (const auto& blockName : availableElementBlocks)
                            {
                                errorMessage += "\n - " + blockName;
                            }
                        }
                    }

                    errorMsg(errorMessage);
                }
            }

            // Collect all domain locations for this domain
            std::vector<std::string> domainLocations =
                domainBlock["location"].template as<std::vector<std::string>>();

            // Check boundaries
            if (domainBlock["boundaries"])
            {
                const auto& boundariesBlock = domainBlock["boundaries"];

                for (const auto& boundaryBlock : boundariesBlock)
                {
                    const std::string boundaryName =
                        boundaryBlock["name"].template as<std::string>();

                    for (const std::string& boundaryLocation :
                         boundaryBlock["location"]
                             .template as<std::vector<std::string>>())
                    {
                        // Check if boundary location exists in Exodus file
                        if (!partExistsInExodus(boundaryLocation))
                        {
                            std::string errorMessage =
                                "Domain '" + domainName + "', boundary '" +
                                boundaryName + "': location '" +
                                boundaryLocation +
                                "' does not exist in the Exodus file";

                            // Get all sidesets touching the domain element
                            // blocks
                            auto availableSidesets =
                                getSidesetsForElementBlocks(domainLocations);
                            if (!availableSidesets.empty())
                            {
                                errorMessage += "\nAvailable sidesets touching "
                                                "domain element blocks:";
                                for (const auto& sidesetName :
                                     availableSidesets)
                                {
                                    errorMessage += "\n - " + sidesetName;
                                }
                            }

                            errorMsg(errorMessage);
                        }

                        // Check if boundary is attached to at least one
                        // domain location
                        bool isAttached = false;
                        for (const std::string& domainLocation :
                             domainLocations)
                        {
                            if (isSidesetAttachedToElementBlock(
                                    boundaryLocation, domainLocation))
                            {
                                isAttached = true;
                                break;
                            }
                        }

                        if (!isAttached)
                        {
                            std::string errorMessage =
                                "Domain '" + domainName + "', boundary '" +
                                boundaryName + "': location '" +
                                boundaryLocation +
                                "' is not attached to any of the domain "
                                "locations";

                            // Get all sidesets touching the domain element
                            // blocks
                            auto availableSidesets =
                                getSidesetsForElementBlocks(domainLocations);
                            if (!availableSidesets.empty())
                            {
                                errorMessage += "\nAvailable sidesets touching "
                                                "domain element blocks:";
                                for (const auto& sidesetName :
                                     availableSidesets)
                                {
                                    errorMessage += "\n - " + sidesetName;
                                }
                            }

                            errorMsg(errorMessage);
                        }
                    }
                }
            }
        }
    }

    std::cout << "Finished validating YAML input" << std::endl;
}

} // namespace accel
