// File       : displacementDiffusionAssemblerElemInterfaceConditions.cpp
// Created    : Tue Nov 26 2024
// Author     : Mhamad Mahdi Alloush
// Description:
// Copyright 2024 CCFNUM HSLU T&A. All Rights Reserved.

#ifdef HAS_INTERFACE
#include "displacementDiffusionAssembler.h"

namespace accel
{

void displacementDiffusionAssembler::assembleElemTermsInterfaces_(
    const domain* domain,
    Context* ctx)
{
    for (const interface* interf : domain->interfacesRef())
    {
        if (interf->isInternal())
        {
            assembleElemTermsInterfaceSide_(
                domain, interf->masterInfoPtr(), ctx);
            assembleElemTermsInterfaceSide_(
                domain, interf->slaveInfoPtr(), ctx);
        }
        else
        {
            if (interf->isFluidSolidType())
            {
                // strong dirichlet treatment
            }
            else
            {
                assembleElemTermsInterfaceSide_(
                    domain, interf->interfaceSideInfoPtr(domain->index()), ctx);
            }
        }
    }
}

} /* namespace accel */
#endif /* HAS_INTERFACE */
