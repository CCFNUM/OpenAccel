// File       : phasicContinuityEquation.h
// Description: Conservative Euler-Euler phase-fraction transport

#ifndef PHASICCONTINUITYEQUATION_H
#define PHASICCONTINUITYEQUATION_H

#include "eulerEulerModel.h"
#include "equation.h"
#include "linearSystem.h"
#include "phasicContinuityAssembler.h"

namespace accel
{

class phasicContinuityEquation : public equation, public linearSystem<1>
{
public:
    phasicContinuityEquation(realm* realm,
                             eulerEulerModel* model,
                             label phaseIndex);

    void checkDomain(const std::shared_ptr<domain> domain) override;
    bool isConverged() const override;
    void setup() override;
    void initialize() override;
    void postInitialize() override;
    void preSolve() override;
    void solve() override;
    void preTimeStep() override;
    void printScales() override;

    equationID getID() override
    {
        return equationID::volumeFraction;
    }

protected:
    void setResidualScales_() override;

private:
    eulerEulerModel* model_;
    label phaseIndex_;
    std::unique_ptr<phasicContinuityAssembler> assembler_;
};

} // namespace accel

#endif // PHASICCONTINUITYEQUATION_H
