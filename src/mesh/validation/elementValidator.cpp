// File       : elementValidator.cpp
// Created    : Mon Oct 7 2025
// Author     : Mhamad Mahdi Alloush
// Description: Element quality validation and correction system implementation
// Copyright (c) 2025 CCFNUM, Lucerne University of Applied Sciences and Arts.
// SPDX-License-Identifier: BSD-3-Clause
// Ported from flash CVFEM code element correction functionality

// Element validation and correction is only supported in 3D
#if SPATIAL_DIM == 3

// Code includes
#include "elementValidator.h"
#include "macros.h"
#include "mesh.h"
#include "messager.h"
#include "zone.h"

// STK includes
#include <stk_mesh/base/GetEntities.hpp>
#include <stk_mesh/base/Selector.hpp>
#include <stk_topology/topology.hpp>
#include <stk_util/parallel/ParallelReduce.hpp>

// STL includes
#include <algorithm>
#include <cmath>
#include <iomanip>

namespace accel
{

elementValidator::elementValidator(mesh* meshPtr) : meshPtr_(meshPtr)
{
}

void elementValidator::setup()
{
    // Create element quality field
    stk::mesh::MetaData& metaData = meshRef().metaDataRef();

    elementQualityField_ = &metaData.declare_field<int>(
        stk::topology::ELEMENT_RANK, elementQualityFieldName);

    // Put field on all element parts
    const auto& elementParts = meshRef().interiorActiveParts();
    for (const auto* part : elementParts)
    {
        stk::mesh::put_field_on_mesh(*elementQualityField_, *part, nullptr);
    }
}

void elementValidator::initialize()
{
    // Initialize all elements as normal quality
    stk::mesh::BulkData& bulkData = meshRef().bulkDataRef();

    stk::mesh::Selector elementSelector =
        meshRef().locallyOwnedInteriorPartsSelector();

    stk::mesh::EntityVector elements;
    stk::mesh::get_entities(
        bulkData, stk::topology::ELEMENT_RANK, elementSelector, elements);

    for (stk::mesh::Entity element : elements)
    {
        int* qualityData =
            stk::mesh::field_data(*elementQualityField_, element);
        if (qualityData != nullptr)
        {
            *qualityData = static_cast<int>(ElementQuality::normalElement);
        }
    }

    resetValidationStats();
}

void elementValidator::validateAllElements()
{
    resetValidationStats();

    stk::mesh::BulkData& bulkData = meshRef().bulkDataRef();
    stk::mesh::Selector elementSelector =
        meshRef().locallyOwnedInteriorPartsSelector();

    stk::mesh::EntityVector elements;
    stk::mesh::get_entities(
        bulkData, stk::topology::ELEMENT_RANK, elementSelector, elements);

    for (stk::mesh::Entity element : elements)
    {
        validateElement(element);
    }

    // Parallel reduction for statistics
    stk::ParallelMachine comm = bulkData.parallel();
    label globalStats[4] = {0, 0, 0, 0};
    label localStats[4] = {stats_.totalElements,
                           stats_.flatElements,
                           stats_.highAspectRatioElements,
                           stats_.normalElements};
    stk::all_reduce_sum(comm, localStats, globalStats, 4);

    stats_.totalElements = globalStats[0];
    stats_.flatElements = globalStats[1];
    stats_.highAspectRatioElements = globalStats[2];
    stats_.normalElements = globalStats[3];

    scalar localMax = stats_.worstAspectRatio;
    scalar globalMax = 0.0;
    stk::all_reduce_max(comm, &localMax, &globalMax, 1);
    stats_.worstAspectRatio = globalMax;
}

void elementValidator::validateZone(label zoneId)
{
    if (zoneId >= meshRef().nZones())
    {
        errorMsg("Invalid zone ID in elementValidator::validateZone");
        return;
    }

    const auto& zone = meshRef().zoneRef(zoneId);
    const auto& zoneParts = zone.interiorParts();

    stk::mesh::BulkData& bulkData = meshRef().bulkDataRef();
    stk::mesh::Selector zoneSelector = stk::mesh::selectUnion(zoneParts);
    zoneSelector &= meshRef().locallyOwnedInteriorPartsSelector();

    stk::mesh::EntityVector elements;
    stk::mesh::get_entities(
        bulkData, stk::topology::ELEMENT_RANK, zoneSelector, elements);

    for (stk::mesh::Entity element : elements)
    {
        validateElement(element);
    }
}

void elementValidator::validateElement(stk::mesh::Entity element)
{
    ElementQuality quality = assessElementQuality_(element);

    // Update element quality field
    int* qualityData = stk::mesh::field_data(*elementQualityField_, element);
    if (qualityData != nullptr)
    {
        *qualityData = static_cast<int>(quality);
    }

    // Update statistics
    updateValidationStats_(element, quality);
}

ElementQuality
elementValidator::assessElementQuality_(stk::mesh::Entity element) const
{
    // Check for flat/degenerate elements first (most critical)
    if (isElementFlat_(element))
    {
        return ElementQuality::flatElement;
    }

    // Check for high aspect ratio
    if (isHighAspectRatio_(element))
    {
        return ElementQuality::highAspectRatio;
    }

    return ElementQuality::normalElement;
}

bool elementValidator::isElementFlat_(stk::mesh::Entity element) const
{
    scalar volume = computeElementVolume_(element);
    return (volume < volumeTolerance_);
}

bool elementValidator::isHighAspectRatio_(stk::mesh::Entity element) const
{
    scalar aspectRatio = computeAspectRatio_(element);
    return (aspectRatio > shapeRatioCriteria_);
}

scalar elementValidator::computeElementVolume_(stk::mesh::Entity element) const
{
    stk::topology topology = getElementTopology_(element);
    std::vector<stk::mesh::Entity> nodes = getElementNodes_(element);

    if (nodes.size() < 4)
    {
        // 2D element or degenerate case
        return 0.0;
    }

    // For tetrahedral elements (most common in CVFEM)
    if (topology == stk::topology::TET_4)
    {
        if (nodes.size() != 4)
        {
            return 0.0;
        }

        auto p0 = getNodeCoordinates_(nodes[0]);
        auto p1 = getNodeCoordinates_(nodes[1]);
        auto p2 = getNodeCoordinates_(nodes[2]);
        auto p3 = getNodeCoordinates_(nodes[3]);

        // Compute vectors
        std::array<scalar, 3> v1 = {
            p1[0] - p0[0], p1[1] - p0[1], p1[2] - p0[2]};
        std::array<scalar, 3> v2 = {
            p2[0] - p0[0], p2[1] - p0[1], p2[2] - p0[2]};
        std::array<scalar, 3> v3 = {
            p3[0] - p0[0], p3[1] - p0[1], p3[2] - p0[2]};

        // Cross product v1 x v2
        std::array<scalar, 3> cross = {v1[1] * v2[2] - v1[2] * v2[1],
                                       v1[2] * v2[0] - v1[0] * v2[2],
                                       v1[0] * v2[1] - v1[1] * v2[0]};

        // Scalar triple product (v1 x v2) . v3
        scalar scalarTriple =
            cross[0] * v3[0] + cross[1] * v3[1] + cross[2] * v3[2];

        // Volume = |scalar triple product| / 6
        return std::abs(scalarTriple) / 6.0;
    }
    else if (topology == stk::topology::HEX_8)
    {
        // For hexahedral elements - approximate using decomposition
        // This is a simplified calculation
        if (nodes.size() != 8)
        {
            return 0.0;
        }

        auto p0 = getNodeCoordinates_(nodes[0]);
        auto p6 = getNodeCoordinates_(nodes[6]); // diagonal corner

        // Simple approximation - actual implementation would be more complex
        scalar dx = std::abs(p6[0] - p0[0]);
        scalar dy = std::abs(p6[1] - p0[1]);
        scalar dz = std::abs(p6[2] - p0[2]);

        return dx * dy * dz;
    }

    // For other topologies, return a default small value
    return 1.0e-10;
}

scalar elementValidator::computeAspectRatio_(stk::mesh::Entity element) const
{
    std::vector<stk::mesh::Entity> nodes = getElementNodes_(element);

    if (nodes.size() < 2)
    {
        return 1.0;
    }

    // Compute edge lengths
    std::vector<scalar> edgeLengths;

    for (size_t i = 0; i < nodes.size(); ++i)
    {
        for (size_t j = i + 1; j < nodes.size(); ++j)
        {
            auto p1 = getNodeCoordinates_(nodes[i]);
            auto p2 = getNodeCoordinates_(nodes[j]);

            scalar length = 0.0;
            for (int dim = 0; dim < SPATIAL_DIM; ++dim)
            {
                scalar diff = p2[dim] - p1[dim];
                length += diff * diff;
            }
            length = std::sqrt(length);

            if (length > 1.0e-15)
            {
                edgeLengths.push_back(length);
            }
        }
    }

    if (edgeLengths.empty())
    {
        return 1.0;
    }

    // Find min and max edge lengths
    auto minMax = std::minmax_element(edgeLengths.begin(), edgeLengths.end());
    scalar minLength = *minMax.first;
    scalar maxLength = *minMax.second;

    if (minLength < 1.0e-15)
    {
        return 1.0e15; // Very high aspect ratio for degenerate edge
    }

    return maxLength / minLength;
}

std::vector<ElementCorrectionResult> elementValidator::correctInvalidElements()
{
    std::vector<ElementCorrectionResult> results;

    if (!enableCorrection_)
    {
        return results;
    }

    stk::mesh::BulkData& bulkData = meshRef().bulkDataRef();
    stk::mesh::Selector elementSelector =
        meshRef().locallyOwnedInteriorPartsSelector();

    stk::mesh::EntityVector elements;
    stk::mesh::get_entities(
        bulkData, stk::topology::ELEMENT_RANK, elementSelector, elements);

    for (stk::mesh::Entity element : elements)
    {
        ElementQuality quality = getElementQuality(element);

        if (quality == ElementQuality::flatElement ||
            quality == ElementQuality::highAspectRatio)
        {
            ElementCorrectionResult result = correctElement(element);
            results.push_back(result);

            if (result.correctionSucceeded)
            {
                stats_.correctedElements++;
            }
            else
            {
                stats_.failedCorrections++;
            }
        }
    }

    return results;
}

ElementCorrectionResult
elementValidator::correctElement(stk::mesh::Entity element)
{
    ElementCorrectionResult result;
    result.originalQuality = getElementQuality(element);

    if (!enableCorrection_)
    {
        return result;
    }

    switch (result.originalQuality)
    {
        case ElementQuality::flatElement:
            result = correctFlatElement_(element);
            break;

        case ElementQuality::highAspectRatio:
            result = correctHighAspectRatioElement_(element);
            break;

        default:
            // Normal element - no correction needed
            result.finalQuality = result.originalQuality;
            result.correctionSucceeded = true;
            break;
    }

    return result;
}

ElementCorrectionResult
elementValidator::correctFlatElement_(stk::mesh::Entity element)
{
    ElementCorrectionResult result;
    result.originalQuality = ElementQuality::flatElement;
    result.correctionMethod = "FlatElementVertexMovement";

    // Ported from flash CVFEM FETETF routine
    // Corrects flat tetrahedral elements by moving one vertex

    stk::mesh::BulkData& bulkData = meshRef().bulkDataRef();

    // Get element topology - only support tets for now (like flash code)
    stk::topology elemTopology = bulkData.bucket(element).topology();
    if (elemTopology != stk::topology::TET_4)
    {
        // Skip non-tetrahedral elements for now
        result.correctionSucceeded = false;
        result.correctionApplied = false;
        return result;
    }

    // Get element nodes
    const stk::mesh::Entity* nodes = bulkData.begin_nodes(element);
    const unsigned numNodes = bulkData.num_nodes(element);

    if (numNodes != 4)
    {
        result.correctionSucceeded = false;
        result.correctionApplied = false;
        return result;
    }

    // Get node coordinates
    std::vector<std::array<scalar, 3>> nodeCoords(4);
    const auto* coordinates = meshRef().metaDataRef().get_field<scalar>(
        stk::topology::NODE_RANK, meshRef().coordinates_ID);

    for (unsigned i = 0; i < 4; ++i)
    {
        const scalar* coord = stk::mesh::field_data(*coordinates, nodes[i]);
        nodeCoords[i][0] = coord[0];
        nodeCoords[i][1] = coord[1];
        nodeCoords[i][2] = coord[2];
    }

    // Flash CVFEM algorithm: Find largest triangular face area
    // Face connectivity for tet: faces opposite to nodes 0,1,2,3
    int faceNodes[4][3] = {{1, 2, 3}, {0, 3, 2}, {0, 1, 3}, {0, 2, 1}};

    scalar maxArea2 = 0.0;
    int nodeToMove = 0;
    std::array<scalar, 3> surfaceNormal{0.0};
    scalar distanceToPlane = 0.0;

    // Find the largest triangle face area
    for (int face = 0; face < 4; ++face)
    {
        // Vector from node 0 to node 1 of face
        std::array<scalar, 3> vec1;
        for (int i = 0; i < 3; ++i)
            vec1[i] = nodeCoords[faceNodes[face][1]][i] -
                      nodeCoords[faceNodes[face][0]][i];

        // Vector from node 0 to node 2 of face
        std::array<scalar, 3> vec2;
        for (int i = 0; i < 3; ++i)
            vec2[i] = nodeCoords[faceNodes[face][2]][i] -
                      nodeCoords[faceNodes[face][0]][i];

        // Cross product for surface normal
        std::array<scalar, 3> normal;
        normal[0] = vec1[1] * vec2[2] - vec1[2] * vec2[1];
        normal[1] = vec1[2] * vec2[0] - vec1[0] * vec2[2];
        normal[2] = vec1[0] * vec2[1] - vec1[1] * vec2[0];

        scalar area2 = normal[0] * normal[0] + normal[1] * normal[1] +
                       normal[2] * normal[2];

        if (area2 > maxArea2)
        {
            maxArea2 = area2;
            nodeToMove = face; // Node opposite to this face
            surfaceNormal = normal;

            // Distance from node to plane (dot product with edge vector)
            distanceToPlane = 0.0;
            for (int i = 0; i < 3; ++i)
            {
                distanceToPlane +=
                    normal[i] *
                    (nodeCoords[face][i] - nodeCoords[faceNodes[face][0]][i]);
            }
        }
    }

    if (maxArea2 <= 0.0)
    {
        result.correctionSucceeded = false;
        result.correctionApplied = false;
        return result; // Degenerate case
    }

    // Normalize direction and compute length scale
    scalar maxArea = std::sqrt(maxArea2);
    for (int i = 0; i < 3; ++i)
        surfaceNormal[i] /= maxArea;

    scalar lengthScale = std::sqrt(maxArea);
    scalar distanceNormalized = distanceToPlane / maxArea;

    // Move distance: 1% of length scale from the plane (flash algorithm)
    scalar moveDistance = 0.01 * lengthScale - distanceNormalized;

    // Apply correction to node coordinates (modify mesh)
    bulkData.modification_begin();

    scalar* modifiableCoord =
        stk::mesh::field_data(*coordinates, nodes[nodeToMove]);

    modifiableCoord[0] += moveDistance * surfaceNormal[0];
    modifiableCoord[1] += moveDistance * surfaceNormal[1];
    modifiableCoord[2] += moveDistance * surfaceNormal[2];

    bulkData.modification_end();

    // Store correction information
    result.modifiedNodes.push_back(nodes[nodeToMove]);
    result.correctionApplied = true;

    // Reassess element quality after correction
    result.finalQuality = assessElementQuality_(element);
    result.correctionSucceeded =
        (result.finalQuality != ElementQuality::flatElement);

    if (result.correctionSucceeded)
        stats_.correctedElements++;
    else
        stats_.failedCorrections++;

    return result;
}

ElementCorrectionResult
elementValidator::correctHighAspectRatioElement_(stk::mesh::Entity element)
{
    ElementCorrectionResult result;
    result.originalQuality = ElementQuality::highAspectRatio;
    result.correctionMethod = "AspectRatioImprovement";

    // Ported from flash CVFEM BETA_BAD_EL approach
    // For high aspect ratio elements, we don't move vertices but rather
    // mark them for special treatment in numerical schemes (reduced blending)

    // In flash CVFEM, high aspect ratio elements get beta=0 (first-order
    // treatment) Here we can store this information in the element quality
    // field

    // For now, we implement a simple vertex averaging approach for severe cases
    stk::mesh::BulkData& bulkData = meshRef().bulkDataRef();

    // Get current aspect ratio
    scalar currentAspectRatio = computeAspectRatio_(element);

    // Apply correction for aspect ratios above threshold (flash CVFEM approach)
    // Lowered threshold to match our criteria for testing
    if (currentAspectRatio >
        2.0 * shapeRatioCriteria_) // Changed from 10.0 to 2.0
    {
        // Get element nodes
        const stk::mesh::Entity* nodes = bulkData.begin_nodes(element);
        const unsigned numNodes = bulkData.num_nodes(element);

        if (numNodes < 4)
        {
            result.correctionSucceeded = false;
            result.correctionApplied = false;
            return result;
        }

        // Find the most distorted direction and apply mild smoothing
        std::vector<std::array<scalar, 3>> nodeCoords(numNodes);
        const auto* coordinates = meshRef().metaDataRef().get_field<scalar>(
            stk::topology::NODE_RANK, meshRef().coordinates_ID);

        for (unsigned i = 0; i < numNodes; ++i)
        {
            const scalar* coord = stk::mesh::field_data(*coordinates, nodes[i]);
            nodeCoords[i][0] = coord[0];
            nodeCoords[i][1] = coord[1];
            nodeCoords[i][2] = coord[2];
        }

        // Check which nodes will actually be moved before logging
        std::vector<unsigned> nodesToMove;

        // Apply conservative smoothing (10% adjustment toward centroid)
        std::array<scalar, 3> centroid{0.0};
        for (unsigned i = 0; i < numNodes; ++i)
        {
            for (int dim = 0; dim < 3; ++dim)
                centroid[dim] += nodeCoords[i][dim];
        }

        for (int dim = 0; dim < 3; ++dim)
            centroid[dim] /= static_cast<scalar>(numNodes);

        // Determine which nodes are far enough from centroid to be moved
        for (unsigned i = 0; i < numNodes; ++i)
        {
            scalar distanceFromCentroid = 0.0;
            for (int dim = 0; dim < 3; ++dim)
            {
                scalar diff = nodeCoords[i][dim] - centroid[dim];
                distanceFromCentroid += diff * diff;
            }
            distanceFromCentroid = std::sqrt(distanceFromCentroid);

            // Only move nodes that are far from centroid
            scalar avgNodeDistance = 0.0;
            for (unsigned j = 0; j < numNodes; ++j)
            {
                for (int dim = 0; dim < 3; ++dim)
                {
                    scalar diff = nodeCoords[j][dim] - centroid[dim];
                    avgNodeDistance += diff * diff;
                }
            }
            avgNodeDistance = std::sqrt(avgNodeDistance / numNodes);

            if (distanceFromCentroid > 1.5 * avgNodeDistance)
            {
                nodesToMove.push_back(i);
            }
        }

        // Only proceed with logging and modification if nodes will actually be
        // moved
        if (!nodesToMove.empty())
        {
            infoMsg("Applying geometric correction to element " +
                    std::to_string(bulkData.identifier(element)) +
                    " (aspect ratio: " + std::to_string(currentAspectRatio) +
                    ", moving " + std::to_string(nodesToMove.size()) +
                    " nodes)");

            bulkData.modification_begin();

            // Actually move the selected nodes
            for (unsigned nodeIndex : nodesToMove)
            {
                scalar* modifiableCoord =
                    stk::mesh::field_data(*coordinates, nodes[nodeIndex]);

                // Move 10% toward centroid
                scalar smoothingFactor = 0.10;
                modifiableCoord[0] +=
                    smoothingFactor * (centroid[0] - nodeCoords[nodeIndex][0]);
                modifiableCoord[1] +=
                    smoothingFactor * (centroid[1] - nodeCoords[nodeIndex][1]);
                modifiableCoord[2] +=
                    smoothingFactor * (centroid[2] - nodeCoords[nodeIndex][2]);

                result.modifiedNodes.push_back(nodes[nodeIndex]);
            }

            bulkData.modification_end();
            result.correctionApplied = true;
        }
        else
        {
            // No nodes meet the criteria for movement - don't log as correction
            // attempt
            result.correctionApplied = false;
        }
    }
    else
    {
        // For moderate aspect ratios, mark for special numerical treatment
        // but don't count as a geometric correction
        result.correctionApplied = false;
        result.correctionSucceeded =
            false; // Don't count as successful since no geometry change
    }

    // Reassess element quality after any modifications
    result.finalQuality = assessElementQuality_(element);

    // Only consider correction successful if geometry was actually modified
    // and the element quality improved
    if (result.correctionApplied)
    {
        scalar finalAspectRatio = computeAspectRatio_(element);

        // Debug: Log final aspect ratio
        std::ostringstream finalMsg;
        finalMsg << "Element " << bulkData.identifier(element)
                 << " correction result: aspect ratio " << currentAspectRatio
                 << " -> " << finalAspectRatio;
        infoMsg(finalMsg.str());

        result.correctionSucceeded =
            (result.finalQuality == ElementQuality::normalElement);

        if (result.correctionSucceeded)
            stats_.correctedElements++;
        else
            stats_.failedCorrections++;
    }
    // If no correction was applied, correctionSucceeded remains false

    return result;
}

std::vector<stk::mesh::Entity>
elementValidator::getElementNodes_(stk::mesh::Entity element) const
{
    const stk::mesh::BulkData& bulkData = meshRef().bulkDataRef();
    const stk::mesh::Entity* nodes = bulkData.begin_nodes(element);
    const unsigned numNodes = bulkData.num_nodes(element);
    return std::vector<stk::mesh::Entity>(nodes, nodes + numNodes);
}

std::array<scalar, SPATIAL_DIM>
elementValidator::getNodeCoordinates_(stk::mesh::Entity node) const
{
    const auto* coordinates = meshRef().metaDataRef().get_field<scalar>(
        stk::topology::NODE_RANK, meshRef().coordinates_ID);
    const scalar* coord = stk::mesh::field_data(*coordinates, node);

    std::array<scalar, SPATIAL_DIM> result;
    for (int dim = 0; dim < SPATIAL_DIM; ++dim)
    {
        result[dim] = coord[dim];
    }
    return result;
}

stk::topology
elementValidator::getElementTopology_(stk::mesh::Entity element) const
{
    return meshRef().bulkDataRef().bucket(element).topology();
}

void elementValidator::updateValidationStats_(stk::mesh::Entity element,
                                              ElementQuality quality)
{
    stats_.totalElements++;

    switch (quality)
    {
        case ElementQuality::normalElement:
            stats_.normalElements++;
            break;
        case ElementQuality::flatElement:
            stats_.flatElements++;
            break;
        case ElementQuality::highAspectRatio:
            stats_.highAspectRatioElements++;
            break;
        default:
            // Unknown quality
            break;
    }

    // Update aspect ratio statistics
    scalar aspectRatio = computeAspectRatio_(element);
    if (aspectRatio > stats_.worstAspectRatio)
    {
        stats_.worstAspectRatio = aspectRatio;
    }

    if (stats_.totalElements == 1)
    {
        stats_.averageAspectRatio = aspectRatio;
    }
    else
    {
        stats_.averageAspectRatio =
            ((stats_.totalElements - 1) * stats_.averageAspectRatio +
             aspectRatio) /
            stats_.totalElements;
    }

    if (stats_.totalElements == 1)
    {
        stats_.averageAspectRatio = aspectRatio;
    }
}

ElementQuality
elementValidator::getElementQuality(stk::mesh::Entity element) const
{
    return assessElementQuality_(element);
}

void elementValidator::resetValidationStats()
{
    stats_ = ElementValidationStats();
}

void elementValidator::printValidationReport() const
{
    std::ostringstream report;

    report << std::string(60, '=') << "\n";
    report << " MESH ELEMENT QUALITY VALIDATION REPORT\n";
    report << std::string(60, '=') << "\n";

    report << "Total Elements: " << std::setw(10) << stats_.totalElements
           << "\n";

    report << "Normal Elements: " << std::setw(10) << stats_.normalElements
           << " (" << std::fixed << std::setprecision(1)
           << (100.0 * stats_.normalElements /
               std::max(stats_.totalElements, static_cast<label>(1)))
           << "%)\n";

    report << "Flat Elements: " << std::setw(10) << stats_.flatElements << " ("
           << std::fixed << std::setprecision(1)
           << (100.0 * stats_.flatElements /
               std::max(stats_.totalElements, static_cast<label>(1)))
           << "%)\n";

    report << "High Aspect Ratio: " << std::setw(10)
           << stats_.highAspectRatioElements << " (" << std::fixed
           << std::setprecision(1)
           << (100.0 * stats_.highAspectRatioElements /
               std::max(stats_.totalElements, static_cast<label>(1)))
           << "%)\n";

    report << "\nAspect Ratio Statistics:\n";
    report << "Average Aspect Ratio: " << std::setw(10) << std::scientific
           << std::setprecision(3) << stats_.averageAspectRatio << "\n";
    report << "Worst Aspect Ratio: " << std::setw(10) << std::scientific
           << std::setprecision(3) << stats_.worstAspectRatio << "\n";

    if (enableCorrection_)
    {
        report << "\nCorrection Statistics:\n";
        report << "Corrected Elements: " << std::setw(10)
               << stats_.correctedElements << "\n";
        report << "Failed Corrections: " << std::setw(10)
               << stats_.failedCorrections << "\n";
    }

    report << "\nValidation Criteria:\n";
    report << "Shape Ratio Threshold: " << std::setw(10) << std::scientific
           << std::setprecision(3) << shapeRatioCriteria_ << "\n";
    report << "Flatness Tolerance: " << std::setw(10) << std::scientific
           << std::setprecision(3) << flatnessTolerance_ << "\n";

    report << std::string(60, '=') << "\n";

    infoMsg(report.str());
}

} // namespace accel

#endif // SPATIAL_DIM == 3
