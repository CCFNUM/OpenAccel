// File : solidDisplacementAssembler.h
// Created : Thu Dec 04 2025 08:42:10 (+0100)
// Author : Mhamad Mahdi Alloush
// Description: Assembler for the solid displacement equation in structural
// mechanics
// Copyright (c) 2025 CCFNUM, Lucerne University of Applied Sciences
// and Arts.
// SPDX-License-Identifier: BSD-3-Clause

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
        return domain->solidMechanics_.formulation_ ==
                       kinematicFormulationType::updatedLagrangian
                   ? mesh::coordinates_ID
                   : mesh::original_coordinates_ID;
    }

    std::string getDualNodalVolumeID_(const domain* domain) const override
    {
        return domain->solidMechanics_.formulation_ ==
                       kinematicFormulationType::updatedLagrangian
                   ? mesh::dual_nodal_volume_ID
                   : mesh::original_dual_nodal_volume_ID;
    }

    std::string getExposedAreaVectorID_(const domain* domain) const override
    {
        return domain->solidMechanics_.formulation_ ==
                       kinematicFormulationType::updatedLagrangian
                   ? mesh::exposed_area_vector_ID
                   : mesh::original_exposed_area_vector_ID;
    }

public:
    solidDisplacementAssembler(solidMechanicsModel* model)
        : Base(model), model_(model)
    {
    }

protected:
    void postAssemble_(const domain* domain, Context* ctx) override;

    void applySymmetryConditions_(const domain* domain, Vector& b);

private:
    void assembleNodeTermsFusedSteady_(const domain* domain,
                                       Context* ctx) override;

    void assembleNodeTermsFusedFirstOrderUnsteady_(const domain* domain,
                                                   Context* ctx) override;

    void assembleNodeTermsFusedSecondOrderUnsteady_(const domain* domain,
                                                    Context* ctx) override;

    void assembleElemTermsInterior_(const domain* domain,
                                    Context* ctx) override;
#ifdef HAS_INTERFACE
    void assembleElemTermsInterfaceSide_(
        const domain* domain,
        const interfaceSideInfo* interfaceSideInfoPtr,
        Context* ctx) override;
#endif /* HAS_INTERFACE */

    void assembleElemTermsBoundary_(const domain* domain,
                                    Context* ctx) override;
};

} /* namespace accel */

#endif // SOLIDDISPLACEMENTASSEMBLER_H
