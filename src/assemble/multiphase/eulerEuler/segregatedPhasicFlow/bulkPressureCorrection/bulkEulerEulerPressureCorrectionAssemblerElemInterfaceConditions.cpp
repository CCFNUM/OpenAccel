// File       : bulkEulerEulerPressureCorrectionAssemblerElemInterfaceConditions.cpp
// Description: Interface flux terms for common Euler-Euler pressure correction

#include "bulkEulerEulerPressureCorrectionAssembler.h"

namespace accel
{

void bulkEulerEulerPressureCorrectionAssembler::assembleElemTermsInterfaces_(
    const domain* domain,
    Context* ctx)
{
    for (const interface* interf : domain->interfacesRef())
    {
        if (interf->isFluidSolidType())
        {
            continue;
        }
        if (interf->isInternal())
        {
            assembleSideFlux_(
                domain, interf->masterInfoRef().currentPartVec_, ctx);
            assembleSideFlux_(
                domain, interf->slaveInfoRef().currentPartVec_, ctx);
        }
        else
        {
            assembleSideFlux_(
                domain,
                interf->interfaceSideInfoPtr(domain->index())
                    ->currentPartVec_,
                ctx);
        }
    }
}

} // namespace accel
