// File : simulation.cpp
// Created : Tue Jun 11 2024 15:06:38 (+0100)
// Author : Mhamad Mahdi Alloush
// Description:
// Copyright (c) 2024 CCFNUM, Lucerne University of Applied Sciences and Arts.
// SPDX-License-Identifier: BSD-3-Clause

// code
#include "simulation.h"
#include "domain.h"
#include "elementField.h"
#include "mesh.h"
#include "realm.h"

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

    printSolverHeader(argc, argv);

    inputFilePath_ = fs::path(inputFileName);
    inputNode_ = YAML::LoadFile(inputFilePath_.c_str());

    if (getYAMLSimulationNode()["verbose"])
    {
        verbose_ = getYAMLSimulationNode()["verbose"].template as<int>();
    }

    // 1) Create components: read/setup
    create_();

    // 2) Initialize components: populate/fill
    initialize_();
}

// Destructor

simulation::~simulation()
{
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
        io_write_counter_ = ++controlsRef().globalIter;
        io_write_time_ = io_write_counter_;

        if (messager::master())
            std::cout << std::endl
                      << "Iter = " << controlsRef().iter << std::endl
                      << std::endl;

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
    }
}

void simulation::runTransient()
{
    while (controlsRef().getTotalTime() - controlsRef().time >
           50 * ::accel::SMALL)
    {
        controlsRef().advanceAndSetTimestep();
        controlsRef().time += controlsRef().getTimestep();
        io_write_counter_ = controlsRef().getTimeStepCount();
        io_write_time_ = controlsRef().time;

        label minCoeffIterations = controlsRef()
                                       .solverRef()
                                       .solverControl_.basicSettings_
                                       .convergenceControl_.minIterations_;
        label maxCoeffIterations = controlsRef()
                                       .solverRef()
                                       .solverControl_.basicSettings_
                                       .convergenceControl_.maxIterations_;

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
                std::cout << std::endl
                          << "Iter = " << controlsRef().iter << std::endl
                          << std::endl;

            // update mesh in case of deforming/moving mesh
            if (this->meshRef().anyZoneMeshDeforming() ||
                this->meshRef().anyZoneMeshTransforming())
            {
                if (!controlsRef()
                         .solverRef()
                         .solverControl_.advancedOptions_.equationControls_
                         .meshMotion_.freezePerTimestep_)
                {
                    meshMotionPtr_->update();
                }
            }

            controlsRef().globalIter++;

            auto elapsedIterTimeStart =
                std::chrono::high_resolution_clock::now();

            // assemble and solve systems now
            assembleAndSolveSystems();

            auto elapsedIterTimeEnd = std::chrono::high_resolution_clock::now();
            std::chrono::duration<scalar> elapsedSeconds =
                elapsedIterTimeEnd - elapsedIterTimeStart;

            if (messager::master())
            {
                std::cout << std::endl
                          << " Iter CPU Time: " << std::scientific
                          << elapsedSeconds.count() << std::endl
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
        }

        // post work: post-process and write
        postWork();
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
        if (messager::master() && equation->subIters() > 1)
        {
            std::cout << std::endl
                      << "Sub-iterating " + equation->name() << ":\n";
        }

        for (label subIter = 1; subIter <= equation->subIters(); subIter++)
        {
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
            // FIXME: [2024-03-07] This check may not apply
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

bool simulation::checkConvergence()
{
    bool convergence = true;
    for (auto& equation : equationVector_)
    {
        convergence = convergence && equation->isConverged();
    }

    return convergence;
}

// Access

label simulation::nDomains() const
{
    return domainVector_.size();
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
