// File       : displacementDiffusionAssembler.h
// Created    : Tue Nov 26 2024
// Author     : Mhamad Mahdi Alloush
// Description: Assembly of the displacement diffusion linear system
// Copyright 2024 CCFNUM HSLU T&A. All Rights Reserved.

#ifndef DISPLACEMENTDIFFUSIONASSEMBLER_H
#define DISPLACEMENTDIFFUSIONASSEMBLER_H

#include "phiAssembler.h"

namespace accel
{

class displacementDiffusionAssembler : public phiAssembler<SPATIAL_DIM>
{
public:
    using Base = phiAssembler<SPATIAL_DIM>;

    using Base::phiAssembler;

    void postAssemble(const domain* domain, Context* ctx) override;

protected:
    void applySymmetryConditions_(const domain* domain, Context* ctx) override;

protected:
    // kernel drivers
    void assembleNodeTermsFused_(const domain* domain, Context* ctx) override
    {
    }

    void assembleElemTermsInterfaces_(const domain* domain,
                                      Context* ctx) override;

    void assembleElemTermsInterfaceSideSpecifiedValue_(
        const domain* domain,
        const interfaceSideInfo* interfaceSideInfoPtr,
        Context* ctx);

    // Boundary conditions
    void assembleElemTermsBoundary_(const domain* domain,
                                    Context* ctx) override;

    std::string getCoordinatesID_(const domain* domain) const override
    {
        bool relDisp = true;
        if (!domain->zonePtr()
                 ->deformationRef()
                 .displacementRelativeToPreviousMesh())
        {
            relDisp = false;
        }

        return relDisp ? mesh::coordinates_ID : mesh::original_coordinates_ID;
    }

    std::string getDualNodalVolumeID_(const domain* domain) const override
    {
        bool relDisp = true;
        if (!domain->zonePtr()
                 ->deformationRef()
                 .displacementRelativeToPreviousMesh())
        {
            relDisp = false;
        }

        return relDisp ? mesh::dual_nodal_volume_ID
                       : mesh::original_dual_nodal_volume_ID;
    }

    std::string getExposedAreaVectorID_(const domain* domain) const override
    {
        bool relDisp = true;
        if (!domain->zonePtr()
                 ->deformationRef()
                 .displacementRelativeToPreviousMesh())
        {
            relDisp = false;
        }

        return relDisp ? mesh::exposed_area_vector_ID
                       : mesh::original_exposed_area_vector_ID;
    }
};

} /* namespace accel */

#endif // DISPLACEMENTDIFFUSIONASSEMBLER_H
