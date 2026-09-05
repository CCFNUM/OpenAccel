// File       : controls.h
// Created    : Fri Aug 25 2023 12:55:24 (+0100)
// Author     : Mhamad Mahdi Alloush
// Description: Simulation controls for time stepping, solver settings, and
//              output
// Copyright 2023 CCFNUM HSLU T&A. All Rights Reserved.

#ifndef CONTROLS_H
#define CONTROLS_H

// code
#include "Profiler.h"
#include "types.h"

namespace accel
{

struct analysisTypeDictionary
{
    // number of entries in timestep ring buffer (sufficient for 4th order
    // multi-step methods, increase DT_ENTRIES for a larger history should it
    // ever be required)
    static constexpr int DT_ENTRIES = 4;

    bool transient_{false};
    label timeStepCount_{0};      // timestep ID
    scalar totalTime_{0};         // final time
    scalar initialTimestep_{0.0}; // what is specified in YAML input
    std::array<scalar, DT_ENTRIES> timestep_{0}; // buffer of timesteps

    struct timeStepsDictionary
    {
        struct timestepAdaptationDictionary
        {
            timestepAdaptationType option_ = timestepAdaptationType::maxCourant;
            scalar courantNumber_{5.0};
            scalar minTimestep_{0.0};
            scalar maxTimestep_{VBIG};
            scalar timestepDecreaseFactor_{0.8};
            scalar timestepIncreaseFactor_{1.06};
        };

        timestepMode mode_{timestepMode::constant};
        std::list<scalar> startTime_{0};
        std::list<scalar> intervalLength_{0};
        std::list<scalar> timestepInterval_{0};
        scalar period_{VBIG};
        label timestepUpdateFrequency_{1};
        timestepAdaptationDictionary timestepAdaptation_;
    };

    timeStepsDictionary timeSteps_;
};

struct solverDictionary
{
    struct solverControlDictionary
    {
        struct basicSettingsDictionary
        {
            advectionSchemeType advectionScheme_ = advectionSchemeType::upwind;
            advectionSchemeType turbulenceNumerics_ =
                advectionSchemeType::upwind;
            transientSchemeType transientScheme_ =
                transientSchemeType::firstOrderBackwardEuler;

            struct relaxationParametersDictionary
            {
                scalar relaxMass_ = 1.0;
                scalar wallScaleRelaxationFactor_ = 1.0;
                scalar energyRelaxationFactor_ = 1.0;
                scalar velocityRelaxationFactor_ = 1.0;
                scalar pressureRelaxationFactor_ = 1.0;
                scalar turbulenceRelaxationFactor_ = 1.0;
                scalar solidDisplacementRelaxationFactor_ = 1.0;
            };

            struct convergenceControlDictionary
            {
                label minIterations_ = 1;
                label maxIterations_ = 100;
                scalar physicalTimescale_ = 1.0;

                relaxationParametersDictionary relaxationParameters_;
            };

            struct interpolationSchemeTypeDictionary
            {
                interpolationSchemeType velocityInterpolationType_ =
                    interpolationSchemeType::trilinear;
                interpolationSchemeType pressureInterpolationType_ =
                    interpolationSchemeType::linearLinear;
                interpolationSchemeType temperatureInterpolationType_ =
                    interpolationSchemeType::linearLinear;
                interpolationSchemeType
                    turbulentKineticEnergyInterpolationType_ =
                        interpolationSchemeType::linearLinear;
                interpolationSchemeType
                    turbulentEddyFrequencyInterpolationType_ =
                        interpolationSchemeType::linearLinear;
                interpolationSchemeType
                    turbulentDissipationRateInterpolationType_ =
                        interpolationSchemeType::linearLinear;
                interpolationSchemeType
                    transitionOnsetReynoldsNumberInterpolationType_ =
                        interpolationSchemeType::linearLinear;
                interpolationSchemeType
                    turbulentIntermittencyInterpolationType_ =
                        interpolationSchemeType::linearLinear;
                interpolationSchemeType wallScaleInterpolationType_ =
                    interpolationSchemeType::linearLinear;
                interpolationSchemeType volumeFractionInterpolationType_ =
                    interpolationSchemeType::linearLinear;
                interpolationSchemeType displacementInterpolationType_ =
                    interpolationSchemeType::linearLinear;

                // Gradient interpolation scheme types
                interpolationSchemeType velocityGradientInterpolationType_ =
                    interpolationSchemeType::linearLinear;
                interpolationSchemeType pressureGradientInterpolationType_ =
                    interpolationSchemeType::linearLinear;
                interpolationSchemeType temperatureGradientInterpolationType_ =
                    interpolationSchemeType::linearLinear;
                interpolationSchemeType
                    turbulentKineticEnergyGradientInterpolationType_ =
                        interpolationSchemeType::linearLinear;
                interpolationSchemeType
                    turbulentEddyFrequencyGradientInterpolationType_ =
                        interpolationSchemeType::linearLinear;
                interpolationSchemeType
                    turbulentDissipationRateGradientInterpolationType_ =
                        interpolationSchemeType::linearLinear;
                interpolationSchemeType
                    transitionOnsetReynoldsNumberGradientInterpolationType_ =
                        interpolationSchemeType::linearLinear;
                interpolationSchemeType
                    turbulentIntermittencyGradientInterpolationType_ =
                        interpolationSchemeType::linearLinear;
                interpolationSchemeType wallScaleGradientInterpolationType_ =
                    interpolationSchemeType::linearLinear;
                interpolationSchemeType
                    volumeFractionGradientInterpolationType_ =
                        interpolationSchemeType::linearLinear;
                interpolationSchemeType displacementGradientInterpolationType_ =
                    interpolationSchemeType::linearLinear;
            };

            struct convergenceCriteriaDictionary
            {
                residualType residualType_ = residualType::RMS;
                scalar residualTarget_ = 1e-4;

                struct physicsConvergenceDictionary
                {
                    bool enabled_ = false;
                    bool writeResiduals_ = false;
                    std::vector<physicsConvergenceType> criteria_;
                    scalar fsiInterfaceResidualTarget_ = 1e-3;
                    scalar fsiForceResidualTarget_ = 1e-3;
                };

                physicsConvergenceDictionary physicsConvergence_;
            };

            convergenceControlDictionary convergenceControl_;
            interpolationSchemeTypeDictionary interpolationSchemeType_;
            convergenceCriteriaDictionary convergenceCriteria_;
        };

        struct advancedOptionsDictionary
        {
            struct pressureLevelInformationDictionary
            {
                pressureLevelInformationSpecification option_ =
                    pressureLevelInformationSpecification::automatic;
                std::vector<scalar> cartesianCoordinates_ =
                    std::vector<scalar>(SPATIAL_DIM, 0.0);
                scalar relativePressureLevel_ = 0.0;
            };

            struct interfaceTransferDictionary
            {
                scalar searchTolerance_ = 1e-4;
                scalar searchExpansionFactor_ = 1.5;
                bool forceResearch_ = false;
                bool conservativeFluxTransfer_ = false;
                label verbose_ = 0;
            };

            struct solidMechanicsDictionary
            {
                label verbose_ = 0;
            };

            struct equationControlsDictionary
            {
                struct subIterationsDictionary
                {
                    label pressureCorrection_ = 1;
                    label solidDisplacement_ = 1;
                    label segregatedFlow_ = 1;
                    label volumeFraction_ = 1;
                };

                struct volumeFractionSmoothingDictionary
                {
                    bool smoothVolumeFraction_ = false;
                    label smoothingIterations_ = 3;
                    scalar fourierNumber_ = 0.25;
                    // curvature-extension smoother: false = box_average
                    // (default), true = laplacian (reuses the VOF diffusion
                    // smoother)
                    bool curvatureSmootherLaplacian_ = false;
                    label curvatureSmoothingIterations_ = 40;
                };

                struct meshMotionDictionary
                {
                    bool freezePerTimestep_ = true;
                    label maxSmoothingIterations_ = 5;
                };

                subIterationsDictionary subIterations_;
                volumeFractionSmoothingDictionary volumeFractionSmoothing_;
                meshMotionDictionary meshMotion_;
            };

            struct domainDecompositionDictionary
            {
                std::string method_{""};
                std::map<std::string, std::string> properties_;
            };

            pressureLevelInformationDictionary pressureLevelInformation_;
            interfaceTransferDictionary interfaceTransfer_;
            solidMechanicsDictionary solidMechanics_;
            equationControlsDictionary equationControls_;
            domainDecompositionDictionary domainDecomposition_;
        };

        struct expertParametersDictionary
        {
            bool disableMomentumPredictor_ = false;
            bool printMomentumInterfaceImbalance_ = false;
            bool consistent_ = false;
            bool limitGradients_ = false;
            bool correctGradients_ = false;
            bool incrementalGradientChange_ = true;
            bool relaxGradients_ = true;
            bool falseMassAccumulation_ = true;
            bool fractionalStepMethod_ = false;
            bool coriolisProductionTurbulence_ = false;
            bool bodyForceRedistribution_ = true;
            wallDistanceMethod wallDistanceMethod_ =
                wallDistanceMethod::poisson;
            bool strongDirichletWallScale_ = false;
            scalar volumeFractionBlendingFactorMax_ = 2.0;
            bool bandwidthReduction_ = true;
            bool nodeReordering_ = false;
            bool forceWallDistanceCalculation_ = false;
            bool disablePhysics_ = false;
            bool freezeFlow_ = false;
            bool freezePressure_ = false;
            bool freezeEnergy_ = false;
            // disable subset node graphs (default): subsets are opt-in
            bool forceFullNodeGraph_ = true;
            bool nso_ = false;
            scalar nsoFourthOrderFac_ = 1.0;
            bool highSpeedBlendDamping_ = false;
            // Cap on the high-resolution blend beta in
            // phi_ip = phi_upwind + beta * (grad phi . dx): 0 upwind, 1 linear
            scalar blendFactorMax_ = 1.0;
            nonconformalMethod nonconformalMethod_ =
                nonconformalMethod::discontinuousGalerkin;
            gradientAveragingType cvpgType_ = gradientAveragingType::arithAver;
            solidAssemblerType solidAssemblerType_ = solidAssemblerType::cvfem;
        };

        basicSettingsDictionary basicSettings_;
        advancedOptionsDictionary advancedOptions_;
        expertParametersDictionary expertParameters_;
    };

    struct outputControlDictionary
    {
        struct outputFrequencyDictionary
        {
            outputFrequencyType option_{outputFrequencyType::timestepInterval};
            scalar timeInterval_{1.0};
            label timestepInterval_{1};
        };

        std::string filePath_;
        std::string restartFileName_{"restart.bin"};
        outputFrequencyDictionary outputFrequency_;
        label restartFrequency_{25};
        bool matchFinalTime_{false};
        bool correctedBoundaryValues_ = false;
        bool writeTimestepInfo_{false};
        std::vector<std::string> outputFields_ = {};

        // STK specific (if restart this becomes stk::io::APPEND_RESULTS)
        stk::io::DatabasePurpose writeMode_{stk::io::WRITE_RESULTS};
        std::unique_ptr<Ioss::PropertyManager> propertyManagerPtr_{nullptr};

        // I/O control
        size_t fileIndex_{0};
        label lastRestart_{-1};
        label lastResults_{-1};
        label writeCounter_{0};
        scalar writeTime_{0.0};
        scalar lastResultsTime_{-1.0};
    };

    struct restartControlDictionary
    {
        std::string inputFilePath_;
        bool isRestart_{false};
        bool writeInitial_{false};
        scalar restartTime_{0.0};
        label keepNRestartSnapshots_{4};

        // STK specific (see controlsIO.cpp)
        stk::io::MeshField::TimeMatchOption timeMatchOption_{
            stk::io::MeshField::CLOSEST};
        std::unique_ptr<Ioss::PropertyManager> propertyManagerPtr_{nullptr};

        // state
        size_t fileIndex_{0};
        std::set<std::string> fields_; // registered restart fields
    };

    solverControlDictionary solverControl_;
    outputControlDictionary outputControl_;
    restartControlDictionary restartControl_;
};

class controls
{
public:
    // Constructors

    controls(fs::path working_directory = "./");

    // Destructor

    ~controls();

    // IO

    void read(YAML::Node inputNode);

    // Access

    const solverDictionary& solverRef() const
    {
        return solver_;
    };

    solverDictionary& solverRef()
    {
        return solver_;
    };

    // Matrix-stencil reduction and node renumbering are independent controls.
    bool isReducedStencil() const;

    bool isRenumbered() const;

    bool isCvfemSolidMechanics() const;

    bool isTransient() const;

    bool isHighResolution() const;

    bool isHighResolutionTurbulenceNumerics() const;

    bool isNSO() const;

    bool useAutomaticDomainDecomposition() const;

    label getNumberOfStates() const;

    Profiler& getProfiler()
    {
        return profiler_;
    }

    stk::util::ParameterList& getRestartParam();

    void setRestartParam();

    void deserializeRestartParam(const stk::io::StkMeshIoBroker& io_broker);

    void registerRestartField(const std::string& fieldName);

    scalar getTotalTime() const;

    void setTimestep(const scalar dt);

    scalar getTimestep(const int i = 0) const;

    scalar getPhysicalTimescale() const;

    label getTimeStepCount() const;

    void advanceAndSetTimestep();

    void updateMaxCourant(const scalar maxCourant);

    scalar getMaxCourant() const;

    void resetMaxCourant();

    // IO paths
    fs::path getWorkingDirectory() const;

    fs::path getPostProcessingDirectory() const;

    fs::path getResidualDirectory() const;

    fs::path getRestartDirectory() const;

    fs::path getResultsDirectory() const;

    fs::path getAdaptiveTimesteppingDirectory() const;

    // Run loop-related public variables

    label iter{0};

    label globalIter{0};

    scalar time{0};

private:
    fs::path workingDirectory_;

    analysisTypeDictionary analysisType_;

    solverDictionary solver_;

    Profiler profiler_;

    stk::util::ParameterList restartParameter_;

    scalar maxCourant_{-1};

    std::unique_ptr<std::ofstream> dtLogFile_;

    struct LogFileInfo
    {
        LogFileInfo(timestepMode m)
            : mode(m), period(0), dt(0), maxCourant(0), intervalStart(0),
              intervalLength(0), action("unchanged")
        {
        }

        const timestepMode mode;
        size_t period;
        scalar dt, maxCourant, intervalStart, intervalLength;
        std::string action{"unchanged"};
    };

    void defineRestartParam_();

    int timestepPosition_(const int i) const;

    void writeTimestepLog_(const LogFileInfo& info);

    void periodicIntervalTimestep_(LogFileInfo& info) const;

    void specifiedIntervalTimestep_(LogFileInfo& info);

    void courantAdaptiveTimestep_(LogFileInfo& info) const;
};

} // namespace accel

#endif // CONTROLS_H
