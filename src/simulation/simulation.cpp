// File       : simulation.cpp
// Created    : Tue Jun 11 2024 15:06:38 (+0100)
// Author     : Mhamad Mahdi Alloush
// Description:
// Copyright 2024 CCFNUM HSLU T&A. All Rights Reserved.

// code
#include "simulation.h"
#include "domain.h"
#include "elementField.h"
#include "fsiConvergence.h"
#include "interface.h"
#include "mesh.h"
#include "realm.h"
#include "version.h"

#include <csignal>

extern volatile sig_atomic_t g_signalSent;

namespace accel
{

// Constructors

simulation::simulation(const int argc, const char* argv[]) : verbose_(0)
{
    // Argument parsing
    std::string inputFileName;
    stk::OptionsSpecification desc("Accel-stk Supported Options:");
    desc.add_options()("input,i",
                       "Analysis input file",
                       stk::DefaultValue<std::string>("input.i"),
                       stk::TargetPointer<std::string>(&inputFileName))(
        "residuals,r", "plot residuals")("scales,s", "print scales");

    stk::parse_command_line_args(argc, argv, desc, args_);

    printSolverHeader();

    inputFilePath_ = fs::path(inputFileName);
    inputNode_ = YAML::LoadFile(inputFilePath_.c_str());

    if (getYAMLSimulationNode()["verbose"])
    {
        verbose_ = getYAMLSimulationNode()["verbose"].template as<int>();
    }

    if (messager::master())
    {
        if (verbose_ > 0)
        {
            std::cout << "Git hash: " << accel::GIT_HASH << '\n';
            std::cout << "Git describe: " << accel::GIT_DESCRIBE << '\n';
        }
        std::cout << "Command line:";
        for (int i = 0; i < argc; i++)
        {
            std::cout << " " << argv[i];
        }
        std::cout << std::endl;
    }

    // 1) Create components: read/setup
    create_();

    // 2) Initialize components: populate/fill
    initialize_();
}

// Destructor

simulation::~simulation()
{
    if (gpPtr_)
    {
        gpPtr_->sendcommand("quit");
    }
}

// Methods
void simulation::initialize_()
{
    // Initialize mesh
    meshPtr_->initialize();

    // Setup domains
    for (auto& domain : domainVector_)
    {
        domain->setup();
    }

    // Setup equations
    assert(domainVector_.size() > 0); // at least one domain is required
    for (auto& equation : equationVector_)
    {
        if (messager::master())
        {
            std::cout << "Setup equation `" + equation->name() +
                             "` on realm `" + realmPtr_->name() + "`\n";
        }
        equation->setup();
    }

    // FSI outputs are written by fsiConvergence, which only runs with
    // fsi_convergence enabled: tell the user if they were requested without
    {
        const auto& fsiConvCfg = controlsRef()
                                     .solverRef()
                                     .solverControl_.basicSettings_
                                     .convergenceCriteria_.fsiConvergence_;
        bool forceBalanceRequested = false;
        for (const interface* interf : meshRef().interfaceVector())
        {
            if (interf->isFluidSolidType() && interf->writeForceBalance())
            {
                forceBalanceRequested = true;
            }
        }
        if (!fsiConvCfg.enabled_ &&
            (forceBalanceRequested || fsiConvCfg.diagnosticsEnabled_) &&
            messager::master())
        {
            warningMsg("write_force_balance / fsi_convergence diagnostics "
                       "have no effect: fsi_convergence is not enabled");
        }
    }

    // setup dynamic mesh components in case of deforming/moving mesh
    if (meshMotionPtr_)
    {
        meshMotionPtr_->setup();
    }

    // setup wall distance calculator if required
    if (wallDistancePtr_)
    {
        wallDistancePtr_->setup();
    }

    // initialize dynamic mesh components
    if (meshMotionPtr_)
    {
        meshMotionPtr_->initialize();
    }

    // initialize wall distance calculator if required
    if (wallDistancePtr_)
    {
        wallDistancePtr_->initialize();
    }

    // Initialize equations
    for (auto& equation : equationVector_)
    {
        if (messager::master())
        {
            std::cout << std::endl
                      << "Initializing equation `" + equation->name() +
                             "` on realm `" + realmPtr_->name() + "`\n\n";
        }
        equation->initialize();
    }

    for (auto& equation : equationVector_)
    {
        equation->postInitialize();
    }

    // Initialize the realm
    realmPtr_->initialize();

    // initialize output
    initializeOutput_();
}

void simulation::run()
{
    if (controlsRef().isTransient())
    {
        runTransient();
    }
    else
    {
        runSteadyState();
    }

    // ensure results are written if `match_final_time` is enabled
    writeResults(controlsRef().solverRef().outputControl_.matchFinalTime_);

    // force a restart snapshot for final state
    writeRestart(true);

    printSolverFooter();
}

void simulation::runSteadyState()
{
    auto& io = controlsRef().solverRef().outputControl_;
    const auto& restart = controlsRef().solverRef().restartControl_;

    label minIterations =
        controlsRef()
            .solverRef()
            .solverControl_.basicSettings_.convergenceControl_.minIterations_;
    label maxIterations =
        controlsRef()
            .solverRef()
            .solverControl_.basicSettings_.convergenceControl_.maxIterations_;

    for (controlsRef().iter = 1; controlsRef().iter <= maxIterations;
         controlsRef().iter++)
    {
        io.writeCounter_ = ++controlsRef().globalIter;
        io.writeTime_ = io.writeCounter_;

        if (messager::master())
        {
            std::cout << std::endl << "Iter = " << controlsRef().iter;
            if (restart.isRestart_)
            {
                std::cout << " (global iter = " << controlsRef().globalIter
                          << ")";
            }
            std::cout << std::endl << std::endl;
        }

        auto elapsedIterTimeStart = std::chrono::high_resolution_clock::now();

        // pre work: nothing really ..
        preWork();

        // assemble and solve systems now
        assembleAndSolveSystems();

        auto elapsedIterTimeEnd = std::chrono::high_resolution_clock::now();
        std::chrono::duration<scalar> elapsedSeconds =
            elapsedIterTimeEnd - elapsedIterTimeStart;

        if (messager::master())
        {
            std::cout << std::endl
                      << " Iter CPU Time: " << std::scientific
                      << elapsedSeconds.count() << "[s]" << std::endl
                      << std::endl;
        }

        if (verbose_ > 3)
        {
            getProfiler().printReport("Simulation profile:\n");
        }

        plotResiduals();

        if (checkConvergence() && controlsRef().iter >= minIterations)
        {
            // force post work: post-process and write
            postWork();

            if (messager::master())
            {
                std::cout << "\nConverged.\n";
            }

            break;
        }

        // post work: post-process and write
        postWork();

        // break iteration loop after signal sent once
        if (g_signalSent)
        {
            MPI_Barrier(messager::comm());
            break;
        }
    }
}

void simulation::runTransient()
{
    auto& io = controlsRef().solverRef().outputControl_;

    while (controlsRef().getTotalTime() - controlsRef().time >
           50 * ::accel::SMALL)
    {
        controlsRef().advanceAndSetTimestep();
        controlsRef().time += controlsRef().getTimestep();
        io.writeCounter_ = controlsRef().getTimeStepCount();
        io.writeTime_ = controlsRef().time;

        label minCoeffIterations = controlsRef()
                                       .solverRef()
                                       .solverControl_.basicSettings_
                                       .convergenceControl_.minIterations_;
        label maxCoeffIterations = controlsRef()
                                       .solverRef()
                                       .solverControl_.basicSettings_
                                       .convergenceControl_.maxIterations_;

        // optional outer-iteration cap until all interfaces are coupled
        {
            const label maxUncoupled =
                controlsRef()
                    .solverRef()
                    .solverControl_.basicSettings_.convergenceControl_
                    .maxIterationsUncoupled_;
            if (maxUncoupled > 0 && hasFsiInterface() && !allFsiCoupled())
            {
                maxCoeffIterations = maxUncoupled;
            }
        }

        if (messager::master())
            std::cout << std::endl
                      << "Time = " << std::scientific << std::setprecision(2)
                      << controlsRef().time << "[s]" << std::endl
                      << std::endl;

        // pre work: update time states
        preWork();

        for (controlsRef().iter = 1; controlsRef().iter <= maxCoeffIterations;
             controlsRef().iter++)
        {
            if (messager::master())
            {
                std::cout << std::endl
                          << "Iter = "
                          << controlsRef().iter
                          // show global iterations unconditionally for
                          // transient runs
                          << " (global iter = " << controlsRef().globalIter + 1
                          << ")" << std::endl
                          << std::endl;
            }

            // per-iteration mesh motion (not frozen per timestep)
            const bool updateMeshPerIteration =
                (this->meshRef().anyZoneMeshDeforming() ||
                 this->meshRef().anyZoneMeshTransforming()) &&
                !controlsRef()
                     .solverRef()
                     .solverControl_.advancedOptions_.equationControls_
                     .meshMotion_.freezePerTimestep_;

            const bool isFsiRun = hasFsiInterface();

            // non-FSI run: update mesh in case of deforming/moving mesh,
            // before assembling the systems
            if (updateMeshPerIteration && !isFsiRun)
            {
                meshMotionPtr_->update();
            }

            controlsRef().globalIter++;

            auto elapsedIterTimeStart =
                std::chrono::high_resolution_clock::now();

            // assemble and solve systems now
            assembleAndSolveSystems();

            // FSI: mesh motion after the solid solve (Tukovic Alg. 1 step 7),
            // once coupled
            if (updateMeshPerIteration && isFsiRun && allFsiCoupled())
            {
                meshMotionPtr_->update();
            }

            auto elapsedIterTimeEnd = std::chrono::high_resolution_clock::now();
            std::chrono::duration<scalar> elapsedSeconds =
                elapsedIterTimeEnd - elapsedIterTimeStart;

            if (messager::master())
            {
                std::cout << std::endl
                          << " Iter CPU Time: " << std::scientific
                          << elapsedSeconds.count() << "[s]" << std::endl
                          << std::endl;
            }

            if (verbose_ > 3)
            {
                getProfiler().printReport("Simulation profile:\n");
            }

            plotResiduals();

            if (checkConvergence() && controlsRef().iter >= minCoeffIterations)
            {
                if (messager::master())
                {
                    std::cout << "\nConverged.\n";
                }

                break;
            }

            if (g_signalSent)
            {
                break;
            }
        }

        // post work: post-process and write
        postWork();

        // break iteration loop after signal sent once
        if (g_signalSent)
        {
            MPI_Barrier(messager::comm());
            break;
        }
    }
}

void simulation::assembleAndSolveSystems()
{
    if (messager::master() && verbose_ > 1)
    {
        std::cout << "Solving realm: " << realmPtr_->name() << "\n";
    }

    // iterate realm equations:
    // 1.) pre-solve tasks
    // 2.) solve equation
    // 3.) post-solve tasks

    for (auto& equation : equationVector_)
    {
        // freeze the solid until all interfaces are coupled (s4f coupled_)
        if (equation->getID() == equationID::solidDisplacement)
        {
            if (!allFsiCoupled())
            {
                // once per uncoupled step (first outer iteration)
                if (messager::master() && controlsRef().iter == 1)
                {
                    std::cout << "\n[FSI] t = " << controlsRef().time
                              << " < coupling_start_time: solid frozen, "
                              << "skipping " << equation->name() << "\n";
                }
                continue;
            }
        }

        if (messager::master() && equation->subIters() > 1)
        {
            std::cout << std::endl
                      << "Sub-iterating " + equation->name() << ":\n";
        }

        for (label subIter = 1; subIter <= equation->subIters(); subIter++)
        {
            equation->setCurrentSubIter(subIter);

            // only print sub-iter if it is active (> 1)
            if (messager::master() && equation->subIters() > 1)
            {
                std::cout << std::endl << " sub-iter: " << subIter << "\n";
            }

            // 1.)
            if (messager::master() && verbose_ > 1)
            {
                std::cout << "\t pre-solve equation: " << equation->name()
                          << "\n";
            }
            equation->preSolve();

            // 2.)
            if (messager::master() && verbose_ > 1)
            {
                std::cout << "\t solve equation: " << equation->name() << "\n";
            }
            // FIXME: This check may not apply
            // for segregated equation queues(?):
            // if (!equation->isConverged()) {
            // equation->solve();
            // }
            equation->solve(); // conservative choice

            // 3.)
            if (messager::master() && verbose_ > 1)
            {
                std::cout << "\t post-solve equation: " << equation->name()
                          << "\n";
            }
            equation->postSolve();

            // if converged .. break sub-iter loop
            if (equation->isConverged())
            {
                break;
            }
        }
    }
}

void simulation::preWork()
{
    if (controlsRef().isTransient())
    {
        if (fsiConvergence* fsiConv = getFsiConvergence_())
        {
            fsiConv->resetForTimeStep();
        }

        // equation dependent pre-timestep tasks for transient simulations
        // (typically updating temporal fields)
        for (auto& equation : equationVector_)
        {
            equation->preTimeStep();
        }

        // reset/update mesh motion quantities
        if (this->meshRef().anyZoneMeshDeforming() ||
            this->meshRef().anyZoneMeshTransforming())
        {
            meshMotionPtr_->reset();

            if (controlsRef()
                    .solverRef()
                    .solverControl_.advancedOptions_.equationControls_
                    .meshMotion_.freezePerTimestep_)
            {
                meshMotionPtr_->update();
            }
        }

        // update wall distance if required
        if (wallDistancePtr_)
        {
            wallDistancePtr_->update();
        }
    }
    else
    {
        // update wall distance if required
        if (wallDistancePtr_)
        {
            wallDistancePtr_->update();
        }
    }
}

void simulation::postWork()
{
    // 1. always required for steady and transient simulations
    postProcessPtr_->update();

    // IO
    {
        const auto& outputCtrl = controlsRef().solverRef().outputControl_;
        const auto& outFreq = outputCtrl.outputFrequency_;
        const bool write_now =
            (controlsRef().isTransient() &&
             outFreq.option_ == outputFrequencyType::timeInterval)
                ? writeNowTime(outFreq.timeInterval_)
                : writeNow(outFreq.timestepInterval_);
        writeResults(write_now);
    }
    writeRestart(
        writeNow(controlsRef().solverRef().outputControl_.restartFrequency_));

    // 2. equation dependent post-timestep tasks for transient simulations only
    if (controlsRef().isTransient())
    {
        for (auto& equation : equationVector_)
        {
            equation->postTimeStep();
        }
    }

    // 3. print scales
    if (printScales_)
    {
        for (auto& equation : equationVector_)
        {
            equation->printScales();
        }
    }
}

bool simulation::hasFsiInterface() const
{
    for (const interface* interf : meshRef().interfaceVector())
    {
        if (interf->isFluidSolidType())
        {
            return true;
        }
    }
    return false;
}

bool simulation::allFsiCoupled() const
{
    const scalar tNow = controlsRef().time;
    for (const interface* interf : meshRef().interfaceVector())
    {
        if (interf->isFluidSolidType() &&
            tNow < (interf->couplingStartTime() - ::accel::SMALL))
        {
            return false;
        }
    }
    return true;
}

bool simulation::checkConvergence()
{
    fsiConvergence* fsiConv = getFsiConvergence_();
    // FSI criterion only once all interfaces are coupled
    if (fsiConv && allFsiCoupled())
    {
        fsiConv->update();
        if (fsiConv->enabled() && fsiConv->hasCriteria())
        {
            return fsiConv->isConverged();
        }
    }

    bool convergence = true;
    for (auto& equation : equationVector_)
    {
        convergence = convergence && equation->isConverged();
    }

    return convergence;
}

fsiConvergence* simulation::getFsiConvergence_()
{
    const auto& fsiConv =
        controlsRef()
            .solverRef()
            .solverControl_.basicSettings_.convergenceCriteria_.fsiConvergence_;

    // created when fsi_convergence is enabled, also without criteria
    if (!fsiConv.enabled_)
    {
        return nullptr;
    }

    if (!fsiConvergencePtr_)
    {
        fsiConvergencePtr_ = std::make_unique<fsiConvergence>(*this);
    }
    return fsiConvergencePtr_.get();
}

// Access

label simulation::nDomains() const
{
    return domainVector_.size();
}

equation* simulation::equationPtr(equationID id)
{
    for (const auto& eq : equationVector_)
    {
        if (eq->getID() == id)
        {
            return eq.get();
        }
    }
    return nullptr;
}

domain& simulation::domainRef(label iDomain)
{
    return *domainVector_[iDomain].get();
}

const domain& simulation::domainRef(label iDomain) const
{
    return *domainVector_[iDomain].get();
}

domain& simulation::domainRef(std::string name)
{
    label idx = -1;
    for (auto dm : domainVector_)
    {
        if (dm->name() == name)
        {
            idx = dm->index();
        }
    }
    assert(idx != -1);
    return *domainVector_[idx].get();
}

const domain& simulation::domainRef(std::string name) const
{
    label idx = -1;
    for (auto dm : domainVector_)
    {
        if (dm->name() == name)
        {
            idx = dm->index();
        }
    }
    assert(idx != -1);
    return *domainVector_[idx].get();
}

controls* simulation::controlsPtr()
{
    return controlsPtr_.get();
}

const controls* simulation::controlsPtr() const
{
    return controlsPtr_.get();
}

controls& simulation::controlsRef()
{
    return *controlsPtr_.get();
}

const controls& simulation::controlsRef() const
{
    return *controlsPtr_.get();
}

mesh* simulation::meshPtr()
{
    return meshPtr_.get();
}

const mesh* simulation::meshPtr() const
{
    return meshPtr_.get();
}

mesh& simulation::meshRef()
{
    return *meshPtr_.get();
}

const mesh& simulation::meshRef() const
{
    return *meshPtr_.get();
}

meshMotion* simulation::meshMotionPtr()
{
    return meshMotionPtr_.get();
}

const meshMotion* simulation::meshMotionPtr() const
{
    return meshMotionPtr_.get();
}

meshMotion& simulation::meshMotionRef()
{
    return *meshMotionPtr_.get();
}

const meshMotion& simulation::meshMotionRef() const
{
    return *meshMotionPtr_.get();
}

overrides* simulation::overridesPtr()
{
    return overridesPtr_.get();
}

const overrides* simulation::overridesPtr() const
{
    return overridesPtr_.get();
}

overrides& simulation::overridesRef()
{
    return *overridesPtr_.get();
}

const overrides& simulation::overridesRef() const
{
    return *overridesPtr_.get();
}

} // namespace accel
