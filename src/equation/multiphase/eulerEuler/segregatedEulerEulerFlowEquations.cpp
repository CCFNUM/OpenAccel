// File       : segregatedEulerEulerFlowEquations.cpp
// Description: Euler-Euler equation-system registration and field lifecycle

#include "segregatedEulerEulerFlowEquations.h"

namespace accel
{

segregatedEulerEulerFlowEquations::segregatedEulerEulerFlowEquations(
    realm* realm)
    : equation("Segregated Euler Euler Flow"), eulerEulerModel(realm)
{
    pressureCorrectionEquation_ =
        std::make_unique<bulkEulerEulerPressureCorrectionEquation>(realm, this);
    pressureCorrectionEquation_->setFallbackName(
        canonicalEquationIDName(equation::name_));

    momentumEquations_.resize(nPhases());
    continuityEquations_.resize(nPhases());
    for (label phase = 0; phase < nPhases(); ++phase)
    {
        const label phaseIndex = this->phaseIndex(phase);
        momentumEquations_[phase] =
            std::make_unique<phasicNavierStokesEquation>(
                realm, this, phaseIndex);
        momentumEquations_[phase]->setFallbackName(
            canonicalEquationIDName(equation::name_));

        bool secondarySomewhere = false;
        for (const auto& domain : realm->simulationRef().domainVector())
        {
            secondarySomewhere =
                secondarySomewhere ||
                (isPhaseActive(domain, phaseIndex) &&
                 !isPrimaryPhase(domain, phaseIndex));
        }
        if (secondarySomewhere)
        {
            continuityEquations_[phase] =
                std::make_unique<phasicContinuityEquation>(
                    realm, this, phaseIndex);
            continuityEquations_[phase]->setFallbackName(
                canonicalEquationIDName(equation::name_));
        }
    }

    subIters_ = controlsRef()
                    .solverRef()
                    .solverControl_.advancedOptions_.equationControls_
                    .subIterations_.segregatedFlow_;
}

void segregatedEulerEulerFlowEquations::addDomain(
    std::shared_ptr<domain> domain)
{
    equation::addDomain(domain);
    pressureCorrectionEquation_->addDomain(domain);
    for (label phase = 0; phase < nPhases(); ++phase)
    {
        const label phaseIndex = this->phaseIndex(phase);
        if (!isPhaseActive(domain, phaseIndex))
        {
            continue;
        }
        momentumEquations_[phase]->addDomain(domain);
        if (!isPrimaryPhase(domain, phaseIndex))
        {
            assert(continuityEquations_[phase]);
            continuityEquations_[phase]->addDomain(domain);
        }
    }
}

bool segregatedEulerEulerFlowEquations::isConverged() const
{
    bool converged = pressureCorrectionEquation_->isConverged();
    for (label phase = 0; phase < nPhases(); ++phase)
    {
        converged =
            converged && momentumEquations_[phase]->isConverged();
        if (continuityEquations_[phase])
        {
            converged =
                converged && continuityEquations_[phase]->isConverged();
        }
    }
    return converged;
}

void segregatedEulerEulerFlowEquations::setup()
{
    for (const auto& domain : domainVector_)
    {
        setupPhaseFields(domain);
        // inherited `body_forces` scratch field used to stage each phase
        // through flowModel::redistributeBodyForces
        setupBodyForces(domain);
    }
    for (auto& momentumEquation : momentumEquations_)
    {
        momentumEquation->setup();
    }
    for (auto& continuityEquation : continuityEquations_)
    {
        if (continuityEquation)
        {
            continuityEquation->setup();
        }
    }
    pressureCorrectionEquation_->setup();

    isCreated_ = pressureCorrectionEquation_->isCreated();
    for (label phase = 0; phase < nPhases(); ++phase)
    {
        isCreated_ = isCreated_ && momentumEquations_[phase]->isCreated();
        if (continuityEquations_[phase])
        {
            isCreated_ =
                isCreated_ && continuityEquations_[phase]->isCreated();
        }
    }
}

void segregatedEulerEulerFlowEquations::initialize()
{
    pressureCorrectionEquation_->initialize();
    for (const auto& domain : domainVector_)
    {
        initializePhaseFields(domain);
    }
    for (auto& momentumEquation : momentumEquations_)
    {
        momentumEquation->initialize();
    }
    for (auto& continuityEquation : continuityEquations_)
    {
        if (continuityEquation)
        {
            continuityEquation->initialize();
        }
    }
}

void segregatedEulerEulerFlowEquations::postInitialize()
{
    pressureCorrectionEquation_->postInitialize();
    for (auto& momentumEquation : momentumEquations_)
    {
        momentumEquation->postInitialize();
    }
    for (auto& continuityEquation : continuityEquations_)
    {
        if (continuityEquation)
        {
            continuityEquation->postInitialize();
        }
    }
    isInitialized_ = pressureCorrectionEquation_->isInitialized();
    for (label phase = 0; phase < nPhases(); ++phase)
    {
        isInitialized_ =
            isInitialized_ && momentumEquations_[phase]->isInitialized();
        if (continuityEquations_[phase])
        {
            isInitialized_ =
                isInitialized_ && continuityEquations_[phase]->isInitialized();
        }
    }
}

void segregatedEulerEulerFlowEquations::preTimeStep()
{
    pressureCorrectionEquation_->preTimeStep();
    for (const auto& domain : domainVector_)
    {
        for (const label phaseIndex : activePhaseIndices(domain))
        {
            updatePhaseMassCoefficientPrevTime(domain, phaseIndex);
            updateVelocityPrevTimeField(domain, phaseIndex);
            updateVolumeFractionPrevTimeField(domain, phaseIndex);
            updateDensityPrevTimeField(domain, phaseIndex);
            updateDynamicViscosityPrevTimeField(domain, phaseIndex);
        }
    }
}

void segregatedEulerEulerFlowEquations::preSolve()
{
}

void segregatedEulerEulerFlowEquations::solve()
{
    for (const auto& domain : domainVector_)
    {
        for (const label phaseIndex : activePhaseIndices(domain))
        {
            const label localPhaseIndex =
                domain->globalToLocalMaterialIndex(phaseIndex);
            if (domain->isMaterialCompressible(localPhaseIndex))
            {
                fieldBroker::updateDensity(domain, phaseIndex);
                fieldBroker::updateCompressibility(domain, phaseIndex);
            }
            fieldBroker::updateDynamicViscosity(domain, phaseIndex);
            updatePhaseMassCoefficient(domain, phaseIndex);
        }
        // Must precede the drag: it supplies the volume fraction the drag
        // law shares with the body force.
        updatePhaseSmoothedVolumeFraction(domain);
        updateInterphaseMomentumSources(domain);
        updatePhaseBodyForces(domain);
    }

    if (!controlsRef().solverRef().solverControl_.expertParameters_.freezeFlow_)
    {
        for (auto& momentumEquation : momentumEquations_)
        {
            momentumEquation->preSolve();
            momentumEquation->solve();
            momentumEquation->postSolve();
        }

        for (const auto& domain : domainVector_)
        {
            for (const label phaseIndex : activePhaseIndices(domain))
            {
                updatePhaseMassFlux(domain, phaseIndex, true);
                updatePhaseFluxDivergence(domain, phaseIndex);
            }
        }

        for (label subIteration = 0;
             subIteration < pressureCorrectionEquation_->subIters();
             ++subIteration)
        {
            pressureCorrectionEquation_->preSolve();
            pressureCorrectionEquation_->solve();
            pressureCorrectionEquation_->postSolve();

            for (const auto& domain : domainVector_)
            {
                for (const label phaseIndex : activePhaseIndices(domain))
                {
                    updatePhaseMassFlux(domain, phaseIndex, true);
                    updatePhaseFluxDivergence(domain, phaseIndex);
                    correctPhaseVelocity(domain, phaseIndex);
                }
                updatePressureGradientField(domain);
            }
            for (label phase = 0; phase < nPhases(); ++phase)
            {
                updateVelocityScale(this->phaseIndex(phase));
            }
            if (pressureCorrectionEquation_->isConverged())
            {
                break;
            }
        }
    }

    for (auto& continuityEquation : continuityEquations_)
    {
        if (!continuityEquation)
        {
            continue;
        }
        continuityEquation->preSolve();
        continuityEquation->solve();
        continuityEquation->postSolve();
    }

    for (const auto& domain : domainVector_)
    {
        applyVolumeFractionClosure(domain);
        for (const label phaseIndex : activePhaseIndices(domain))
        {
            updateVolumeFractionGradientField(domain, phaseIndex);
            updatePhaseMassCoefficient(domain, phaseIndex);
            updatePhaseMassFlux(domain, phaseIndex, true);
            updatePhaseFluxDivergence(domain, phaseIndex);
            updateVelocityGradientField(domain, phaseIndex);
            updateVelocityBlendingFactorField(domain, phaseIndex);
        }
        updatePhaseSmoothedVolumeFraction(domain);
        updateInterphaseMomentumSources(domain);
    }
}

void segregatedEulerEulerFlowEquations::printScales()
{
    pressureCorrectionEquation_->printScales();
    for (auto& momentumEquation : momentumEquations_)
    {
        momentumEquation->printScales();
    }
    for (auto& continuityEquation : continuityEquations_)
    {
        if (continuityEquation)
        {
            continuityEquation->printScales();
        }
    }
}

} // namespace accel
