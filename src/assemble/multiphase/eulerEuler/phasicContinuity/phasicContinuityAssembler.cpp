// File       : phasicContinuityAssembler.cpp
// Description: Assembler for conservative Euler-Euler phasic continuity

#include "phasicContinuityAssembler.h"

namespace accel
{

phasicContinuityAssembler::phasicContinuityAssembler(
    eulerEulerModel* model,
    label phaseIndex)
    : phiAssembler<1>(model), model_(model), phaseIndex_(phaseIndex)
{
}

} // namespace accel
