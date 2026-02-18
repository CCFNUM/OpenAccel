// File : elementValidator.h
// Created : Mon Oct 7 2025
// Author : Mhamad Mahdi Alloush
// Description: Element quality validation and correction system
// Copyright (c) 2025 CCFNUM, Lucerne University of Applied Sciences and Arts.
// SPDX-License-Identifier: BSD-3-Clause
// Ported from flash CVFEM code element correction functionality

#ifndef ELEMENTVALIDATOR_H
#define ELEMENTVALIDATOR_H

// Element validation and correction is only supported in 3D
#if SPATIAL_DIM == 3

// STL includes
#include <array>
#include <memory>
#include <vector>

// STK includes
#include <stk_mesh/base/BulkData.hpp>
#include <stk_mesh/base/Entity.hpp>
#include <stk_mesh/base/Field.hpp>
#include <stk_mesh/base/MetaData.hpp>

// Code includes
#include "types.h"

namespace accel
{

class mesh;

// Element quality flags (from flash CVFEM code)
enum class ElementQuality : int
{
    flatElement = -1,    // Flat/degenerate element (needs correction)
    highAspectRatio = 0, // High aspect ratio element
    reserved1 = 1,       // Reserved for future use
    reserved2 = 2,       // Reserved for future use
    reserved3 = 3,       // Reserved for future use
    reserved4 = 4,       // Reserved for future use
    reserved5 = 5,       // Reserved for future use
    reserved6 = 6,       // Reserved for future use
    reserved7 = 7,       // Reserved for future use
    reserved8 = 8,       // Reserved for future use
    reserved9 = 9,       // Reserved for future use
    normalElement = 10   // Normal/good element (default)
};

// Element correction results
struct ElementCorrectionResult
{
    bool correctionApplied = false;
    bool correctionSucceeded = false;
    ElementQuality originalQuality = ElementQuality::normalElement;
    ElementQuality finalQuality = ElementQuality::normalElement;
    std::string correctionMethod;
    std::vector<stk::mesh::Entity> modifiedNodes;
};

// Element validation statistics
struct ElementValidationStats
{
    label totalElements = 0;
    label flatElements = 0;
    label highAspectRatioElements = 0;
    label normalElements = 0;
    label correctedElements = 0;
    label failedCorrections = 0;
    scalar averageAspectRatio = 0.0;
    scalar worstAspectRatio = 0.0;
    stk::mesh::Entity worstElement;
};

/**
 * @brief Element quality validation and correction class
 *
 * This class implements element quality checking and correction algorithms
 * ported from the flash CVFEM code. It can identify flat/degenerate elements,
 * high aspect ratio elements, and apply correction strategies.
 */
class elementValidator
{
private:
    // Reference to mesh
    mesh* meshPtr_;

    // Quality assessment parameters
    scalar shapeRatioCriteria_ =
        1.5; // Threshold for high aspect ratio (TESTING: lowered to trigger
    // corrections)
    scalar flatnessTolerance_ = 1.0e-12; // Tolerance for flatness detection
    scalar volumeTolerance_ = 1.0e-15;   // Minimum acceptable element volume

    // Correction parameters
    bool enableCorrection_ = false;
    label maxCorrectionIterations_ = 5;
    scalar correctionFactor_ = 0.1; // Factor for vertex movement

    // Element quality field
    stk::mesh::Field<int>* elementQualityField_ = nullptr;

    // Validation statistics
    ElementValidationStats stats_;

    // Private methods for quality assessment
    ElementQuality assessElementQuality_(stk::mesh::Entity element) const;

    scalar computeAspectRatio_(stk::mesh::Entity element) const;

    scalar computeElementVolume_(stk::mesh::Entity element) const;

    bool isElementFlat_(stk::mesh::Entity element) const;

    bool isHighAspectRatio_(stk::mesh::Entity element) const;

    // Private methods for element correction
    ElementCorrectionResult correctFlatElement_(stk::mesh::Entity element);

    ElementCorrectionResult
    correctHighAspectRatioElement_(stk::mesh::Entity element);

    bool moveVertexForCorrection_(
        stk::mesh::Entity vertex,
        const std::array<scalar, SPATIAL_DIM>& displacement);

    std::array<scalar, SPATIAL_DIM>
    computeCorrectionDisplacement_(stk::mesh::Entity element,
                                   stk::mesh::Entity vertex) const;

    // Utility methods
    std::vector<stk::mesh::Entity>
    getElementNodes_(stk::mesh::Entity element) const;

    std::array<scalar, SPATIAL_DIM>
    getNodeCoordinates_(stk::mesh::Entity node) const;

    void
    setNodeCoordinates_(stk::mesh::Entity node,
                        const std::array<scalar, SPATIAL_DIM>& coordinates);

    stk::topology getElementTopology_(stk::mesh::Entity element) const;

    void updateValidationStats_(stk::mesh::Entity element,
                                ElementQuality quality);

public:
    // Static field identifiers
    static constexpr char elementQualityFieldName[] = "element_quality";

    // Constructor
    explicit elementValidator(mesh* meshPtr);

    // Destructor
    ~elementValidator() = default;

    // Delete copy constructor and assignment operator
    elementValidator(const elementValidator&) = delete;
    elementValidator& operator=(const elementValidator&) = delete;

    // Setup and initialization
    void setup();
    void initialize();

    // Main validation interface
    void validateAllElements();
    void validateZone(label zoneId);
    void validateElement(stk::mesh::Entity element);

    // Element correction interface
    std::vector<ElementCorrectionResult> correctInvalidElements();
    ElementCorrectionResult correctElement(stk::mesh::Entity element);

    // Quality assessment interface
    ElementQuality getElementQuality(stk::mesh::Entity element) const;
    scalar getElementAspectRatio(stk::mesh::Entity element) const;
    scalar getElementVolume(stk::mesh::Entity element) const;

    // Configuration interface
    void setShapeRatioCriteria(scalar criteria)
    {
        shapeRatioCriteria_ = criteria;
    }

    scalar getShapeRatioCriteria() const
    {
        return shapeRatioCriteria_;
    }

    void setFlatnessTolerance(scalar tolerance)
    {
        flatnessTolerance_ = tolerance;
    }

    scalar getFlatnessTolerance() const
    {
        return flatnessTolerance_;
    }

    void setVolumeTolerance(scalar tolerance)
    {
        volumeTolerance_ = tolerance;
    }

    scalar getVolumeTolerance() const
    {
        return volumeTolerance_;
    }

    void enableElementCorrection(bool enable)
    {
        enableCorrection_ = enable;
    }

    bool isElementCorrectionEnabled() const
    {
        return enableCorrection_;
    }

    void setMaxCorrectionIterations(label iterations)
    {
        maxCorrectionIterations_ = iterations;
    }

    label getMaxCorrectionIterations() const
    {
        return maxCorrectionIterations_;
    }

    void setCorrectionFactor(scalar factor)
    {
        correctionFactor_ = factor;
    }

    scalar getCorrectionFactor() const
    {
        return correctionFactor_;
    }

    // Statistics and reporting
    const ElementValidationStats& getValidationStats() const
    {
        return stats_;
    }

    void resetValidationStats();
    void printValidationReport() const;

    // Access to mesh
    mesh& meshRef()
    {
        return *meshPtr_;
    }

    const mesh& meshRef() const
    {
        return *meshPtr_;
    }

    // Field access
    stk::mesh::Field<int>& elementQualityFieldRef()
    {
        return *elementQualityField_;
    }

    const stk::mesh::Field<int>& elementQualityFieldRef() const
    {
        return *elementQualityField_;
    }
};

} // namespace accel

#endif // SPATIAL_DIM == 3

#endif // ELEMENTVALIDATOR_H