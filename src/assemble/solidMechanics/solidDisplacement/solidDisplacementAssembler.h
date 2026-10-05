// File       : solidDisplacementAssembler.h
// Created    : Thu Dec 04 2025 08:42:10 (+0100)
// Author     : Mhamad Mahdi Alloush
// Description: Assembler for the solid displacement equation in structural
//              mechanics
// Copyright 2026 CCFNUM HSLU T&A. All Rights Reserved.

#ifndef SOLIDDISPLACEMENTASSEMBLER_H
#define SOLIDDISPLACEMENTASSEMBLER_H

#include "phiAssembler.h"
#include "solidMechanicsModel.h"

namespace accel
{

class solidDisplacementAssembler : public phiAssembler<SPATIAL_DIM>
{
private:
    solidMechanicsModel* model_;

public:
    using Base = phiAssembler<SPATIAL_DIM>;

protected:
    using Context = typename Base::Context;
    using Matrix = typename Base::Matrix;
    using Vector = typename Base::Vector;

    std::string getCoordinatesID_(const domain* domain) const override
    {
        bool useOrig = false;
        if (domain->solidMechanics_.formulation_ ==
                kinematicFormulationType::totalLagrangian &&
            field_broker_->controlsRef().isTransient())
        {
            useOrig = true;
        }
        return useOrig ? mesh::original_coordinates_ID : mesh::coordinates_ID;
    }

    std::string getDualNodalVolumeID_(const domain* domain) const override
    {
        bool useOrig = false;
        if (domain->solidMechanics_.formulation_ ==
                kinematicFormulationType::totalLagrangian &&
            field_broker_->controlsRef().isTransient())
        {
            useOrig = true;
        }
        return useOrig ? mesh::original_dual_nodal_volume_ID
                       : mesh::dual_nodal_volume_ID;
    }

    std::string getExposedAreaVectorID_(const domain* domain) const override
    {
        bool useOrig = false;
        if (domain->solidMechanics_.formulation_ ==
                kinematicFormulationType::totalLagrangian &&
            field_broker_->controlsRef().isTransient())
        {
            useOrig = true;
        }
        return useOrig ? mesh::original_exposed_area_vector_ID
                       : mesh::exposed_area_vector_ID;
    }

public:
    solidDisplacementAssembler(solidMechanicsModel* model)
        : Base(model), model_(model)
    {
    }

    void postAssemble(const domain* domain, Context* ctx) override;

protected:
    void applySymmetryConditions_(const domain* domain, Context* ctx) override;

private:
    void assembleNodeTermsFusedSteady_(const domain* domain,
                                       Context* ctx) override;

    void assembleNodeTermsFusedFirstOrderUnsteady_(const domain* domain,
                                                   Context* ctx) override;

    void assembleNodeTermsFusedSecondOrderUnsteady_(const domain* domain,
                                                    Context* ctx) override;

    void assembleElemTermsInterior_(const domain* domain,
                                    Context* ctx) override;
    void assembleElemTermsInterfaceSide_(
        const domain* domain,
        const interfaceSideInfo* interfaceSideInfoPtr,
        Context* ctx) override;

    void assembleElemTermsBoundaryWallSpecifiedFlux_(const domain* domain,
                                                     const boundary* boundary,
                                                     Context* ctx) override;
};

} /* namespace accel */

#endif // SOLIDDISPLACEMENTASSEMBLER_H
