// File       : bulkEulerEulerPressureCorrectionEquation.h
// Description: Common pressure correction assembled from all phase fluxes

#ifndef BULKEULEREULERPRESSURECORRECTIONEQUATION_H
#define BULKEULEREULERPRESSURECORRECTIONEQUATION_H

#include "bulkEulerEulerPressureCorrectionAssembler.h"
#include "eulerEulerModel.h"
#include "equation.h"
#include "linearSystem.h"

namespace accel
{

class bulkEulerEulerPressureCorrectionEquation : public equation,
                                                 public linearSystem<1>
{
public:
    bulkEulerEulerPressureCorrectionEquation(realm* realm,
                                              eulerEulerModel* model);

    void checkDomain(const std::shared_ptr<domain> domain) override;
    bool isConverged() const override;
    void setup() override;
    void initialize() override;
    void postInitialize() override;
    void preSolve() override;
    void solve() override;
    void postSolve() override;
    void preTimeStep() override;
    void printScales() override;

    equationID getID() override
    {
        return equationID::pressureCorrection;
    }

protected:
    void setResidualScales_() override;

private:
    void storePressureCorrection_(
        const domain* domain,
        ::linearSolver::coefficients<linearSystem::BLOCKSIZE>* coefficients);

    eulerEulerModel* model_;
    std::unique_ptr<bulkEulerEulerPressureCorrectionAssembler> assembler_;
};

} // namespace accel

#endif // BULKEULEREULERPRESSURECORRECTIONEQUATION_H
