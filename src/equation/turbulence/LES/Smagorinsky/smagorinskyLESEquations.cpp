#include "smagorinskyLESEquations.h"

namespace accel
{

smagorinskyLESEquations::smagorinskyLESEquations(realm* realm)
    : equation("Smagorinsky LES"), SmagorinskyLESModel(realm)
{
}

bool smagorinskyLESEquations::isConverged() const
{
    return true;
}

void smagorinskyLESEquations::setup()
{
    equation::isCreated_ = true;
}

void smagorinskyLESEquations::initialize()
{
    equation::isInitialized_ = true;
}

void smagorinskyLESEquations::postInitialize()
{
}

void smagorinskyLESEquations::preTimeStep()
{
}

void smagorinskyLESEquations::preSolve()
{
    FOREACH_DOMAIN(updateTurbulentDynamicViscosity);
    FOREACH_DOMAIN(updateEffectiveDynamicViscosity);
}

void smagorinskyLESEquations::solve()
{
}

} // namespace accel
