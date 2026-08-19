#ifndef SMAGORINSKYLESEQUATIONS_H
#define SMAGORINSKYLESEQUATIONS_H

#include "SmagorinskyLESModel.h"
#include "equation.h"

namespace accel
{

class smagorinskyLESEquations : public equation,
                                public SmagorinskyLESModel
{
public:
    static constexpr equationID ID = equationID::smagorinskyLES;

    explicit smagorinskyLESEquations(realm* realm);

    bool isConverged() const override;
    void setup() override;
    void initialize() override;
    void postInitialize() override;
    void preTimeStep() override;
    void preSolve() override;
    void solve() override;

    equationID getID() override
    {
        return ID;
    }
};

} // namespace accel

#endif // SMAGORINSKYLESEQUATIONS_H
