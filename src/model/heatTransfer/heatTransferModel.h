// File       : heatTransferModel.h
// Created    : Thu Feb 01 2024 16:41:11 (+0100)
// Author     : Fabian Wermelinger
// Description: Base class for heat transfer problems
// Copyright (c) 2024 CCFNUM, Lucerne University of Applied Sciences and Arts.
// SPDX-License-Identifier: BSD-3-Clause

#ifndef HEATTRANSFERMODEL_H
#define HEATTRANSFERMODEL_H

// code
#include "model.h"

namespace accel
{

class heatTransferModel : public model
{
protected:
    void updateTotalTemperatureField_(const std::shared_ptr<domain> domain);

    void reportHeatData_();

    void updateHeatImbalance_(const std::shared_ptr<domain> domain);

#ifdef HAS_INTERFACE
    void updateInterfaceHeatImbalance_(const std::shared_ptr<domain> domain);
#endif /* HAS_INTERFACE */

public:
    heatTransferModel(realm* realm);

    using fieldBroker::betaRef;
    using fieldBroker::cpRef;
    using fieldBroker::gammaRef;
    using fieldBroker::h0Ref;
    using fieldBroker::hRef;
    using fieldBroker::lambdaEffRef;
    using fieldBroker::lambdaRef;
    using fieldBroker::MaRef;
    using fieldBroker::muEffRef;
    using fieldBroker::pRef;
    using fieldBroker::psiRef;
    using fieldBroker::qDotRef;
    using fieldBroker::rhoRef;
    using fieldBroker::T0Ref;
    using fieldBroker::TRef;
    using fieldBroker::TWallCoeffsRef;

    // TODO: Implement utility methods for use in
    // equations that inherit from this model, e.g., computation of enthalpy
    // from temperature, post processing routines, etc.

    // setup

    virtual void
    setupSpecificHeatCapacity(const std::shared_ptr<domain> domain) override;

    virtual void
    setupThermalConductivity(const std::shared_ptr<domain> domain) override;

    // initialize

    virtual void
    initializeSpecificEnthalpy(const std::shared_ptr<domain> domain) override;

    virtual void initializeSpecificTotalEnthalpy(
        const std::shared_ptr<domain> domain) override;

    virtual void initializeSpecificHeatCapacity(
        const std::shared_ptr<domain> domain) override;

    virtual void initializeThermalConductivity(
        const std::shared_ptr<domain> domain) override;

    virtual void
    initializeCompressibility(const std::shared_ptr<domain> domain) override;

    // update

    virtual void
    updateTemperature(const std::shared_ptr<domain> domain) override;

    virtual void
    updateSpecificEnthalpy(const std::shared_ptr<domain> domain) override;

    virtual void
    updateSpecificTotalEnthalpy(const std::shared_ptr<domain> domain) override;

    virtual void
    updateSpecificHeatCapacity(const std::shared_ptr<domain> domain) override;

    virtual void
    updateThermalConductivity(const std::shared_ptr<domain> domain) override;

    virtual void
    updateCompressibility(const std::shared_ptr<domain> domain) override;

    virtual void
    updateEffectiveThermalConductivity(const std::shared_ptr<domain> domain);

    virtual void updateSpecificHeatCapacityGradientField(
        const std::shared_ptr<domain> domain) override;

private:
    void updateTemperatureSideFields_(const std::shared_ptr<domain> domain);

    void
    updateSpecificEnthalpySideFields_(const std::shared_ptr<domain> domain);

    void updateSpecificTotalEnthalpySideFields_(
        const std::shared_ptr<domain> domain);

    void updateTemperatureBoundarySideFieldInletTotalTemperature_(
        const std::shared_ptr<domain> domain,
        const boundary* boundary);

    void updateTemperatureBoundarySideFieldOpening_(
        const std::shared_ptr<domain> domain,
        const boundary* boundary);

    void updateSpecificEnthalpyBoundarySideFieldWallSpecifiedValue_(
        const std::shared_ptr<domain> domain,
        const boundary* boundary);

    void updateSpecificEnthalpyBoundarySideFieldInletSpecifiedValue_(
        const std::shared_ptr<domain> domain,
        const boundary* boundary);

    void updateSpecificEnthalpyBoundarySideFieldOpening_(
        const std::shared_ptr<domain> domain,
        const boundary* boundary);

    void updateSpecificTotalEnthalpyBoundarySideFieldWallSpecifiedValue_(
        const std::shared_ptr<domain> domain,
        const boundary* boundary);

    void updateSpecificTotalEnthalpyBoundarySideFieldInletSpecifiedValue_(
        const std::shared_ptr<domain> domain,
        const boundary* boundary);

    void updateSpecificTotalEnthalpyBoundarySideFieldOpening_(
        const std::shared_ptr<domain> domain,
        const boundary* boundary);

    // boundary patch data for reporting
    MPI_Datatype MPIHeatBoundaryData;
    MPI_Op MPIHeatBoundaryData_SUM;

    struct HeatBoundaryData
    {
        scalar in = 0.0;
        scalar out = 0.0;
        scalar in_area = 0.0;
        scalar out_area = 0.0;
        scalar total_area = 0.0;

        char name[32] = "NA";
        char type[32] = "NA";
    };

    static void sumHeatBoundaryData_(void* a, void* b, int* n, MPI_Datatype*)
    {
        const HeatBoundaryData* a_ = static_cast<const HeatBoundaryData*>(a);
        HeatBoundaryData* b_ = static_cast<HeatBoundaryData*>(b);
        for (int i = 0; i < *n; i++)
        {
            b_[i].in += a_[i].in;
            b_[i].out += a_[i].out;
            b_[i].in_area += a_[i].in_area;
            b_[i].out_area += a_[i].out_area;
            b_[i].total_area += a_[i].total_area;

            for (label j = 0; j < 32; j++)
            {
                b_[i].name[j] = a_[i].name[j];
                b_[i].type[j] = a_[i].type[j];
            }
        }
    }

    std::vector<HeatBoundaryData> heatBoundaryDataVector_;
    std::vector<HeatBoundaryData> heatInterfaceDataVector_;
};

} /* namespace accel */

#endif // HEATTRANSFERMODEL_H
