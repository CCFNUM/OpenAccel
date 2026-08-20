// File       : segregatedEulerEulerFlowEquations.h
// Description: Euler-Euler equation-system registration and field lifecycle

#ifndef SEGREGATEDEULEREULERFLOWEQUATIONS_H
#define SEGREGATEDEULEREULERFLOWEQUATIONS_H

#include "equation.h"
#include "eulerEulerModel.h"
#include "bulkEulerEulerPressureCorrectionEquation.h"
#include "phasicContinuityEquation.h"
#include "phasicNavierStokesEquation.h"

namespace accel
{

class segregatedEulerEulerFlowEquations : public equation,
                                           public eulerEulerModel
{
public:
    static constexpr equationID ID = equationID::segregatedEulerEulerFlow;

    explicit segregatedEulerEulerFlowEquations(realm* realm);

    void addDomain(std::shared_ptr<domain> domain) override;
    bool isConverged() const override;
    void setup() override;
    void initialize() override;
    void postInitialize() override;
    void preTimeStep() override;
    void preSolve() override;
    void solve() override;
    void printScales() override;

    equationID getID() override
    {
        return ID;
    }

private:
    std::vector<std::unique_ptr<phasicNavierStokesEquation>> momentumEquations_;
    std::vector<std::unique_ptr<phasicContinuityEquation>> continuityEquations_;
    std::unique_ptr<bulkEulerEulerPressureCorrectionEquation>
        pressureCorrectionEquation_;
};

} // namespace accel

#endif // SEGREGATEDEULEREULERFLOWEQUATIONS_H
