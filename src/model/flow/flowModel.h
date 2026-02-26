// File       : flowModel.h
// Created    : Fri Mar 15 2024 15:06:38 (+0100)
// Author     : Fabian Wermelinger
// Description: Abstract base model for Navier-Stokes equations. Implements
// common model code.
// Copyright (c) 2024 CCFNUM, Lucerne University of Applied Sciences and Arts.
// SPDX-License-Identifier: BSD-3-Clause
#ifndef FLOWMODEL_H
#define FLOWMODEL_H

// code
#include "model.h"

namespace accel
{

class flowModel : public model
{
public:
    static constexpr char du_ID[] = "du";
    static constexpr char duTilde_ID[] = "du_tilde";
    static constexpr char F_ID[] = "body_forces";
    static constexpr char FOriginal_ID[] = "body_forces_original";

    flowModel(realm* realm);
    ~flowModel();

    using fieldBroker::betaRef;
    using fieldBroker::CoRef;
    using fieldBroker::cpRef;
    using fieldBroker::muEffRef;
    using fieldBroker::muRef;
    using fieldBroker::p0Ref;
    using fieldBroker::pRef;
    using fieldBroker::psiRef;
    using fieldBroker::TRef;
    using fieldBroker::uWallCoeffsRef;
    using fieldBroker::wallShearStressRef;

    // methods

    void updatePressureScale();

    void setupBodyForces(const std::shared_ptr<domain> domain);

    virtual void computeBodyForces(const std::shared_ptr<domain> domain);

    virtual void redistributeBodyForces(const std::shared_ptr<domain> domain);

    void transformMassFlowRateToRelative(const std::shared_ptr<domain> domain);

    void transformMassFlowRateToAbsolute(const std::shared_ptr<domain> domain);

    void updateFlowReversalFlag(const std::shared_ptr<domain> domain);

    void updateMassDivergenceField(const std::shared_ptr<domain> domain);

    void updateSideMassFlowRateFraction(const std::shared_ptr<domain> domain);

    // initialize

    virtual void
    initializePressure(const std::shared_ptr<domain> domain) override;

    virtual void
    initializeVelocity(const std::shared_ptr<domain> domain) override;

    virtual void
    initializeDensity(const std::shared_ptr<domain> domain) override;

    virtual void
    initializeDynamicViscosity(const std::shared_ptr<domain> domain) override;

    // update

    virtual void updatePressure(const std::shared_ptr<domain> domain) override;

    virtual void updateVelocity(const std::shared_ptr<domain> domain) override;

    virtual void updateDensity(const std::shared_ptr<domain> domain) override;

    virtual void
    updateDynamicViscosity(const std::shared_ptr<domain> domain) override;

    virtual void
    updateEffectiveDynamicViscosity(const std::shared_ptr<domain> domain);

    void updateWallShearStress(const std::shared_ptr<domain> domain);

    void resetCourantNumber();

#ifndef NDEBUG
    void zeroSCLCheckField(const std::shared_ptr<domain> domain);
    void syncSCLCheckField(const std::shared_ptr<domain> domain);
#endif /* NDEBUG */

protected:
    STKScalarField* FSTKFieldPtr_;
    STKScalarField* FOriginalSTKFieldPtr_; // workspace for balanced force

    static constexpr char COMMENT[] = "# ";
    std::unique_ptr<std::ostream> mDot_stream_; // file stream for mass flux

    void reportFlowData_();

    // generic mass flux initialization

    virtual void initializeMassFlowRateInterior_(
        const std::shared_ptr<domain> domain) override;

#ifdef HAS_INTERFACE
    virtual void initializeMassFlowRateInterfaceSideField_(
        const std::shared_ptr<domain> domain,
        const interfaceSideInfo* interfaceSideInfoPtr) override;
#endif /* HAS_INTERFACE */

    virtual void
    initializeMassFlowRateBoundaryField_(const std::shared_ptr<domain> domain,
                                         const boundary* boundary) override;

    virtual void
    initializeMassFlowRateInterior_(const std::shared_ptr<domain> domain,
                                    elementField<scalar, 1>& mDotField,
                                    const nodeField<1, SPATIAL_DIM>& rhoField);

#ifdef HAS_INTERFACE
    virtual void initializeMassFlowRateInterfaceSideField_(
        const std::shared_ptr<domain> domain,
        const interfaceSideInfo* interfaceSideInfoPtr,
        sideField<scalar, 1>& mDotSideField,
        const nodeField<1, SPATIAL_DIM>& rhoField);
#endif /* HAS_INTERFACE */

    virtual void initializeMassFlowRateBoundaryField_(
        const std::shared_ptr<domain> domain,
        const boundary* boundary,
        sideField<scalar, 1>& mDotSideField,
        const nodeField<1, SPATIAL_DIM>& rhoField);

    // generic mass flux update

    virtual void
    updateMassFlowRateInterior_(const std::shared_ptr<domain> domain) override;

#ifdef HAS_INTERFACE
    virtual void updateMassFlowRateInterfaceSideField_(
        const std::shared_ptr<domain> domain,
        const interfaceSideInfo* interfaceSideInfoPtr) override;
#endif /* HAS_INTERFACE */

    virtual void
    updateMassFlowRateBoundaryField_(const std::shared_ptr<domain> domain,
                                     const boundary* boundary) override;

    virtual void
    updateMassFlowRateInterior_(const std::shared_ptr<domain> domain,
                                elementField<scalar, 1>& mDotField,
                                const nodeField<1, SPATIAL_DIM>& rhoField);

#ifdef HAS_INTERFACE
    virtual void updateMassFlowRateInterfaceSideField_(
        const std::shared_ptr<domain> domain,
        const interfaceSideInfo* interfaceSideInfoPtr,
        sideField<scalar, 1>& mDotSideField,
        const nodeField<1, SPATIAL_DIM>& rhoField);
#endif /* HAS_INTERFACE */

    virtual void
    updateMassFlowRateBoundaryField_(const std::shared_ptr<domain> domain,
                                     const boundary* boundary,
                                     sideField<scalar, 1>& mDotSideField,
                                     const nodeField<1, SPATIAL_DIM>& rhoField);

    virtual void updateMassFlowRateBoundaryFieldInletSpecifiedVelocity_(
        const std::shared_ptr<domain> domain,
        const boundary* boundary,
        sideField<scalar, 1>& mDotSideField,
        const nodeField<1, SPATIAL_DIM>& rhoField);

    virtual void updateMassFlowRateBoundaryFieldInletSpecifiedPressure_(
        const std::shared_ptr<domain> domain,
        const boundary* boundary,
        sideField<scalar, 1>& mDotSideField,
        const nodeField<1, SPATIAL_DIM>& rhoField);

    virtual void
    updateMassFlowRateBoundaryFieldInletSpecifiedVelocityAndPressure_(
        const std::shared_ptr<domain> domain,
        const boundary* boundary,
        sideField<scalar, 1>& mDotSideField,
        const nodeField<1, SPATIAL_DIM>& rhoField);

    virtual void updateMassFlowRateBoundaryFieldInletSpecifiedMassFlowRate_(
        const std::shared_ptr<domain> domain,
        const boundary* boundary,
        sideField<scalar, 1>& mDotSideField,
        const nodeField<1, SPATIAL_DIM>& rhoField);

    virtual void updateMassFlowRateBoundaryFieldOutletSpecifiedPressure_(
        const std::shared_ptr<domain> domain,
        const boundary* boundary,
        sideField<scalar, 1>& mDotSideField,
        const nodeField<1, SPATIAL_DIM>& rhoField);

    virtual void updateMassFlowRateBoundaryFieldOutletSpecifiedMassFlowRate_(
        const std::shared_ptr<domain> domain,
        const boundary* boundary,
        sideField<scalar, 1>& mDotSideField,
        const nodeField<1, SPATIAL_DIM>& rhoField,
        scalar F);

    virtual void updateMassFlowRateBoundaryFieldOutletOutflow_(
        const std::shared_ptr<domain> domain,
        const boundary* boundary,
        sideField<scalar, 1>& mDotSideField,
        const nodeField<1, SPATIAL_DIM>& rhoField);

    virtual void updateMassFlowRateBoundaryFieldOpeningPressure_(
        const std::shared_ptr<domain> domain,
        const boundary* boundary,
        sideField<scalar, 1>& mDotSideField,
        const nodeField<1, SPATIAL_DIM>& rhoField);

    void
    transformMassFlowRateToRelative_(const std::shared_ptr<domain> domain,
                                     elementField<scalar, 1>& mDotField,
                                     const nodeField<1, SPATIAL_DIM>& rhoField);
    void
    transformMassFlowRateToAbsolute_(const std::shared_ptr<domain> domain,
                                     elementField<scalar, 1>& mDotField,
                                     const nodeField<1, SPATIAL_DIM>& rhoField);

    void updateMassDivergenceField_(const std::shared_ptr<domain> domain,
                                    const elementField<scalar, 1>& mDotField,
                                    const nodeField<1, SPATIAL_DIM>& rhoField,
                                    nodeField<1>& divField);

    void updateFlowReversalFlag_(const std::shared_ptr<domain> domain,
                                 sideField<scalar, 1>& mDotSideField,
                                 bool ignoreFlagUpdate = false);

    void
    updateSideMassFlowRateFraction_(const std::shared_ptr<domain> domain,
                                    const nodeField<1, SPATIAL_DIM>& rhoField,
                                    massFlowRate& mDotField);

    void updateTotalPressureField_(const std::shared_ptr<domain> domain);

    void updateCourantNumberField_(const std::shared_ptr<domain> domain);

    void updateMachNumberField_(const std::shared_ptr<domain> domain);

    void updateRelativeVelocityField_(const std::shared_ptr<domain> domain);

    void updateMassImbalance_(const std::shared_ptr<domain> domain);

    void updateInterfaceMassImbalance_(const std::shared_ptr<domain> domain);

private:
    // flow data for reporting
    struct FlowData
    {
        scalar maxCo = -1;
    };

    FlowData flowData_;

    // boundary patch data for reporting
    MPI_Datatype MPIFlowBoundaryData;
    MPI_Op MPIFlowBoundaryData_SUM;

    struct FlowBoundaryData
    {
        scalar inflow = 0.0;
        scalar outflow = 0.0;
        scalar inflow_area = 0.0;
        scalar outflow_area = 0.0;
        scalar total_area = 0.0;

        unsigned long int blocked_ip_count = 0;
        unsigned long int total_ip_count = 0;

        char name[32] = "NA";
        char type[32] = "NA";
    };

    static void sumFlowBoundaryData_(void* a, void* b, int* n, MPI_Datatype*)
    {
        const FlowBoundaryData* a_ = static_cast<const FlowBoundaryData*>(a);
        FlowBoundaryData* b_ = static_cast<FlowBoundaryData*>(b);
        for (int i = 0; i < *n; i++)
        {
            b_[i].inflow += a_[i].inflow;
            b_[i].outflow += a_[i].outflow;
            b_[i].inflow_area += a_[i].inflow_area;
            b_[i].outflow_area += a_[i].outflow_area;
            b_[i].total_area += a_[i].total_area;

            b_[i].blocked_ip_count += a_[i].blocked_ip_count;
            b_[i].total_ip_count += a_[i].total_ip_count;

            for (label j = 0; j < 32; j++)
            {
                b_[i].name[j] = a_[i].name[j];
                b_[i].type[j] = a_[i].type[j];
            }
        }
    }

    std::vector<FlowBoundaryData> flowBoundaryDataVector_;
    std::vector<FlowBoundaryData> flowInterfaceDataVector_;

    void updateMassFlowRateBoundaryFieldInletSpecifiedVelocity_(
        const std::shared_ptr<domain> domain,
        const boundary* boundary);

    void updateMassFlowRateBoundaryFieldInletSpecifiedPressure_(
        const std::shared_ptr<domain> domain,
        const boundary* boundary);

    void updateMassFlowRateBoundaryFieldInletSpecifiedVelocityAndPressure_(
        const std::shared_ptr<domain> domain,
        const boundary* boundary);

    void updateMassFlowRateBoundaryFieldInletSpecifiedMassFlowRate_(
        const std::shared_ptr<domain> domain,
        const boundary* boundary);

    void updateMassFlowRateBoundaryFieldOutletSpecifiedPressure_(
        const std::shared_ptr<domain> domain,
        const boundary* boundary);

    void updateMassFlowRateBoundaryFieldOutletSpecifiedMassFlowRate_(
        const std::shared_ptr<domain> domain,
        const boundary* boundary);

    void updateMassFlowRateBoundaryFieldOutletOutflow_(
        const std::shared_ptr<domain> domain,
        const boundary* boundary);

    void updateMassFlowRateBoundaryFieldOpeningPressure_(
        const std::shared_ptr<domain> domain,
        const boundary* boundary);

    void updateVelocitySideFields_(const std::shared_ptr<domain> domain);

    void updatePressureSideFields_(const std::shared_ptr<domain> domain);

    void updateDensitySideFields_(const std::shared_ptr<domain> domain);

    void updateVelocityBoundarySideFieldMassFlowRate_(
        const std::shared_ptr<domain> domain,
        const boundary* boundary);

    void updatePressureBoundarySideFieldAverageStaticPressure_(
        const std::shared_ptr<domain> domain,
        const boundary* boundary);

    void updatePressureBoundarySideFieldInletTotalPressure_(
        const std::shared_ptr<domain> domain,
        const boundary* boundary);

    void updatePressureBoundarySideFieldOpening_(
        const std::shared_ptr<domain> domain,
        const boundary* boundary);
};

} /* namespace accel */

#endif // FLOWMODEL_H
