// File       : phasicNavierStokesEquation.h
// Description: One conservative momentum equation for an Euler-Euler phase

#ifndef PHASICNAVIERSTOKESEQUATION_H
#define PHASICNAVIERSTOKESEQUATION_H

#include "eulerEulerModel.h"
#include "equation.h"
#include "linearSystem.h"
#include "phasicNavierStokesAssembler.h"

namespace accel
{

class phasicNavierStokesEquation : public equation,
                                   public linearSystem<SPATIAL_DIM>
{
public:
    phasicNavierStokesEquation(realm* realm,
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
        return equationID::coupledNavierStokes;
    }

protected:
    void setResidualScales_() override;

private:

    eulerEulerModel* model_;
    label phaseIndex_;
    std::unique_ptr<phasicNavierStokesAssembler> assembler_;
};

} // namespace accel

#endif // PHASICNAVIERSTOKESEQUATION_H
