// File : nodeField.h
// Created : Fri Aug 25 2023 12:55:24 (+0100)
// Author : Mhamad Mahdi Alloush
// Description: Templated node field with boundary/initial conditions and
// gradient support
// Copyright (c) 2023 CCFNUM, Lucerne University of Applied
// Sciences and Arts.
// SPDX-License-Identifier: BSD-3-Clause

#ifndef NODEFIELD_H
#define NODEFIELD_H

// code
#include "boundary.h"
#include "controls.h"
#ifdef HAS_INTERFACE
#include "dataTransfer.h"
#include "dgInfo.h"
#include "interface.h"
#endif /* HAS_INTERFACE */
#include "dataHandler.h"
#include "nodeSideField.h"
#include "scatterToSurface.h"
#include "sideField.h"

// external
#include "exprtk.hpp"

namespace accel
{

class initialConditionDictionary : public dataHandler
{
private:
    initialConditionOption type_ = initialConditionOption::null;

public:
    // Constructors

    initialConditionDictionary() = default;

    // Access

    initialConditionOption type() const
    {
        return type_;
    };

    // Set

    void setType(initialConditionOption type)
    {
        type_ = type;
    };
};

class boundaryConditionDictionary : public dataHandler
{
private:
    label index_ = -1; // index of boundary in corresponding zone

    std::string name_; // simply name of boundary

    boundaryConditionType type_ = boundaryConditionType::zeroGradient;

public:
    // Constructors

    boundaryConditionDictionary() = default;

    // Access

    label index() const
    {
        return index_;
    };

    std::string name() const
    {
        return name_;
    };

    boundaryConditionType type() const
    {
        return type_;
    };

    // Set

    void setIndex(label index)
    {
        index_ = index;
    };

    void setName(std::string name)
    {
        name_ = name;
    };

    void setType(boundaryConditionType type)
    {
        type_ = type;
    };
};

template <size_t N, size_t M = 1>
class nodeField : public field<scalar, N>
{
protected:
    std::unique_ptr<nodeField<N, M>> prevIterFieldPtr_ = nullptr;

    std::unique_ptr<nodeField<N, M>> prevTimeFieldPtr_ = nullptr;

    std::unique_ptr<nodeField<N>> blendingFactorFieldPtr_ = nullptr;

    std::unique_ptr<nodeField<N>> gradientLimiterFieldPtr_ = nullptr;

    std::unique_ptr<nodeField<N, M>> maxValueFieldPtr_ = nullptr;

    std::unique_ptr<nodeField<N, M>> minValueFieldPtr_ = nullptr;

    std::unique_ptr<nodeField<M>> gradFieldPtr_ = nullptr;

    std::unique_ptr<nodeSideField<scalar, N>> nodeSideFieldPtr_ = nullptr;

    std::unique_ptr<sideField<scalar, N>> sideFieldPtr_ = nullptr;

    std::unique_ptr<sideField<scalar, N>> sideFluxFieldPtr_ = nullptr;

#ifdef HAS_INTERFACE
    std::vector<std::unique_ptr<dataTransfer>> dataTransferVector_;
#endif /* HAS_INTERFACE */

    std::vector<std::vector<boundaryConditionDictionary>>
        boundaryConditionsDictionaryArray_;

    std::vector<initialConditionDictionary> initialConditionsDictionaryArray_;

    // advection scheme (upwind or high-resolution)
    advectionSchemeType advectionScheme_ = advectionSchemeType::upwind;

    // interpolation scheme: shifted ip = linearLinear
    interpolationSchemeType interpolationScheme_ =
        interpolationSchemeType::linearLinear;

    // gradient interpolation scheme: shifted ip = linearLinear
    interpolationSchemeType gradientInterpolationScheme_ =
        interpolationSchemeType::linearLinear;

    scalar blendingFactorMax_ = 1.0; // maximum blending factor

    scalar gradURF_ = 0.5625; // gradient under-relaxation

    // values of nodes at boundaries are corrected or conservative?
    bool correctedBoundaryNodeValues_ = false;

    // flag to mark if the field operates in all mediums or a single type of
    // them, for instance, velocity and pressure may only apply to fluid domains
    bool mediumIndependent_ = true;

    // flag to limit the gradient to prevent new extrema
    bool limitGradient_ = false;

    // flag to correct gradient at symmetry planes
    bool correctGradient_ = false;

    // flag to calculate gradient in terms of the incremental change relative to
    // the node: more accurate gradient
    bool incrementalGradientChange_ = false;

    // flag to mark initialzed field at a zone
    std::vector<bool> isInitialized_;

    void putFieldOnRegisteredParts_();

    // set the values given an array (scalar* arr): internal use only

    void setToValue_(const scalar* val);

    void setToValue_(const scalar* val, const stk::mesh::Part& part);

    void setToValue_(const scalar* val, const stk::mesh::PartVector parts);

    virtual void limitGradientField_(label iZone);

    virtual void correctGradientField_(label iZone) {};

    virtual std::string getCoordinatesID_(label iZone) const
    {
        return mesh::coordinates_ID;
    }

    virtual std::string getDualNodalVolumeID_(label iZone) const
    {
        return mesh::dual_nodal_volume_ID;
    }

    virtual std::string getExposedAreaVectorID_(label iZone) const
    {
        return mesh::exposed_area_vector_ID;
    }

public:
    // Constructors

    nodeField(mesh* meshPtr,
              std::string name,
              unsigned numberOfStates,
              bool prevIter,
              bool highResolution = false,
              bool computeGradient = false,
              bool correctedBoundaryNodeValues = false);

    nodeField(mesh* meshPtr, STKScalarField* stkField_ptr);

    // Operations

    void setZone(label iZone) override;

    void putFieldOnPart(const stk::mesh::Part& part);

    bool definedOn(const stk::mesh::Part& part) const;

    bool definedOn(const stk::mesh::PartVector parts) const;

    // set the values given various inputs

    void setToValue(std::initializer_list<scalar> val);

    void setToValue(std::initializer_list<scalar> val,
                    const stk::mesh::Part& part);

    void setToValue(std::initializer_list<scalar> val,
                    const stk::mesh::PartVector parts);

    void setToValue(std::vector<scalar> val);

    void setToValue(std::vector<scalar> val, const stk::mesh::Part& part);

    void setToValue(std::vector<scalar> val, const stk::mesh::PartVector parts);

    void setToValue(std::array<scalar, N> val);

    void setToValue(std::array<scalar, N> val, const stk::mesh::Part& part);

    void setToValue(std::array<scalar, N> val,
                    const stk::mesh::PartVector parts);

    void registerSideFields(label iZone, label iBoundary);

    void registerNodeSideField(label iZone, label iBoundary);

    void registerSideField(label iZone, label iBoundary);

    void registerSideFluxField(label iZone, label iBoundary);

#ifdef HAS_INTERFACE
    void registerSideFieldsForInterfaceSide(label iInterface,
                                            bool master,
                                            bool onlyIfNonoverlap = false);

    void registerSideFluxFieldForInterfaceSide(label iInterface,
                                               bool master,
                                               bool onlyIfNonoverlap = false);
#endif /* HAS_INTERFACE */

    nodeField& operator=(const nodeField& fld);

    virtual void setupGradientField();

    virtual void setupGradientLimiterField();

    virtual void setupBlendingFactorField(bool dummy = false);

    virtual void setupMinMaxFields();

    // Initialize

    virtual void initialize(label iZone, bool force = false);

    virtual void initializeField(label iZone);

    virtual void initializeSideFields(label iZone);

    // Update

    virtual void update(label iZone);

    virtual void updateField(label iZone);

    virtual void updateSideFields(label iZone);

#ifdef HAS_INTERFACE
    virtual void updateInterfaceSideField(label iInterface, bool master);
#endif /* HAS_INTERFACE */

    virtual void updateBoundarySideField(label iZone, label iBoundary);

    virtual void updateBoundarySideFieldSpecifiedValue(label iZone,
                                                       label iBoundary);

    virtual void updateSideFluxField(label iZone, label iBoundary);

    // Prev field updates

    void updatePrevIterField(label iZone);

    void updatePrevTimeField(label iZone);

    void restorePrevIterField(label iZone);

    void restorePrevTimeField(label iZone);

    // Sub-field updates

    void updateBlendingFactorField(label iZone);

    void updateMinMaxFields(label iZone,
                            bool applyExtremaExpansion = false,
                            scalar extremaCoeff = 0.0,
                            scalar smoothLimiter = 0.0);

    void updateGradientField(label iZone);

    // Synchronize

    void synchronize(label iZone);

#ifdef HAS_INTERFACE
    // Data transfer across interface

    void transfer(label iInterface,
                  dataTransferType type = dataTransferType::copy,
                  bool reverse = false);
#endif /* HAS_INTERFACE */

    // Stats

    void updateScale() override;

    virtual scalar max();

    virtual scalar min();

    // Other operations

    void correctBoundaryNodes(label iZone, label iBoundary);

#ifdef HAS_INTERFACE
    void correctInterfaceNodes(label iInterface, bool master);
#endif /* HAS_INTERFACE */

    void relax(label iZone, const scalar urf);

    // Access

    nodeField<N, M>& prevIterRef();

    const nodeField<N, M>& prevIterRef() const;

    nodeField<N, M>& prevTimeRef();

    const nodeField<N, M>& prevTimeRef() const;

    nodeField<N>& blendingFactorRef();

    const nodeField<N>& blendingFactorRef() const;

    nodeField<N>* blendingFactorPtr();

    const nodeField<N>* blendingFactorPtr() const;

    nodeField<N>& gradientLimiterRef();

    const nodeField<N>& gradientLimiterRef() const;

    nodeField<N>* gradientLimiterPtr();

    const nodeField<N>* gradientLimiterPtr() const;

    nodeField<N, M>& maxValueRef();

    const nodeField<N, M>& maxValueRef() const;

    nodeField<N, M>& minValueRef();

    const nodeField<N, M>& minValueRef() const;

    nodeField<M>& gradRef();

    const nodeField<M>& gradRef() const;

    nodeSideField<scalar, N>& nodeSideFieldRef();

    const nodeSideField<scalar, N>& nodeSideFieldRef() const;

    sideField<scalar, N>& sideFieldRef();

    const sideField<scalar, N>& sideFieldRef() const;

    sideField<scalar, N>& sideFluxFieldRef();

    const sideField<scalar, N>& sideFluxFieldRef() const;

    bool isShifted() const
    {
        return interpolationScheme_ == interpolationSchemeType::linearLinear;
    };

    bool isGradientShifted() const
    {
        return gradientInterpolationScheme_ ==
               interpolationSchemeType::linearLinear;
    };

    bool isHighResolution() const
    {
        return advectionScheme_ != advectionSchemeType::upwind;
    };

    bool correctedBoundaryNodeValues() const
    {
        return correctedBoundaryNodeValues_;
    };

    // Boundary-related access

    label nBoundaries() const
    {
        return boundaryConditionsDictionaryArray_.size();
    };

    boundaryConditionDictionary& boundaryConditionRef(label iZone,
                                                      label iBoundary)
    {
        return boundaryConditionsDictionaryArray_[iZone][iBoundary];
    };

    const boundaryConditionDictionary&
    boundaryConditionRef(label iZone, label iBoundary) const
    {
        return boundaryConditionsDictionaryArray_[iZone][iBoundary];
    };

    // Initialization-related access

    initialConditionDictionary& initialConditionRef(label iZone)
    {
        return initialConditionsDictionaryArray_[iZone];
    };

    const initialConditionDictionary& initialConditionRef(label iZone) const
    {
        return initialConditionsDictionaryArray_[iZone];
    };

    // Getters and setters

    interpolationSchemeType interpolationScheme() const
    {
        return interpolationScheme_;
    }

    void setInterpolationScheme(interpolationSchemeType interpolationScheme)
    {
        interpolationScheme_ = interpolationScheme;
    }

    interpolationSchemeType gradientInterpolationScheme() const
    {
        return gradientInterpolationScheme_;
    }

    void setGradientInterpolationScheme(
        interpolationSchemeType gradientInterpolationScheme)
    {
        gradientInterpolationScheme_ = gradientInterpolationScheme;
    }

    bool isInitialized(label iZone) const
    {
        assert(isInitialized_.size() > iZone);
        assert(this->isZoneSet(iZone));
        return isInitialized_[iZone];
    }

    void setIsInitialized(label iZone, bool state = true)
    {
        assert(isInitialized_.size() > iZone);
        assert(this->isZoneSet(iZone));
        isInitialized_[iZone] = state;
    }

    bool mediumIndependent() const
    {
        return mediumIndependent_;
    };

    void setMediumIndependent(bool state)
    {
        mediumIndependent_ = state;
    };

    scalar gradURF() const
    {
        return gradURF_;
    }

    void setGradURF(scalar urf)
    {
        gradURF_ = urf;
    }
};

// Out-of-line definitions

template <size_t N, size_t M>
std::ostream& operator<<(std::ostream& os, const nodeField<N, M>& field);

template <size_t N, size_t M>
std::ostream& operator<<(std::ostream& os,
                         const boundaryConditionDictionary& dict);

} // namespace accel

#include "nodeField.hpp"

#endif // NODEFIELD_H
