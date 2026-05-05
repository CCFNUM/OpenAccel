// File       : volumeFractionAssemblerElemInterfaceConditions.cpp
// Created    : Thu Apr 17 2025 22:07:40 (+0100)
// Author     : Mhamad Mahdi Alloush
// Description:
// Copyright 2025 CCFNUM HSLU T&A. All Rights Reserved.

#ifdef HAS_INTERFACE

#include "interface.h"
#include "volumeFractionAssembler.h"

namespace accel
{

void volumeFractionAssembler::assembleElemTermsInterfaceSide_(
    const domain* domain,
    const interfaceSideInfo* interfaceSideInfoPtr,
    Context* ctx)
{
    if (interfaceSideInfoPtr->interfPtr()->isFluidSolidType())
        return;

    phiAssembler<1>::assembleElemTermsInterfaceSide_(
        domain, interfaceSideInfoPtr, ctx);
}

} // namespace accel

#endif /* HAS_INTERFACE */
