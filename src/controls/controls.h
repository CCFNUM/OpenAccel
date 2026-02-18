// File : controls.h
// Created : Fri Aug 25 2023 12:55:24 (+0100)
// Author : Mhamad Mahdi Alloush
// Description: Simulation controls for time stepping, solver settings, and
// output
// Copyright (c) 2023 CCFNUM, Lucerne University of Applied Sciences and
// Arts. SPDX-License-Identifier: BSD-3-Clause

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
            bool reducedStencil_;

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

            struct equationControlsDictionary
            {
                struct subIterationsDictionary
                {
                    label pressureCorrection_ = 1;
                    label solidDisplacement_ = 1;
                    label segregatedFlow_ = 1;
                    label volumeFraction_ = 1;
                };

                struct accelerationDictionary
                {
                    struct solidDisplacementAccelerationDictionary
                    {
                        bool aitkenEnabled_ = false;
                        scalar aitkenInitialOmega_ = 1.0;
                        scalar aitkenOmegaMin_ = 0.1;
                        scalar aitkenOmegaMax_ = 1.0;
                    };

                    solidDisplacementAccelerationDictionary solidDisplacement_;
                };

                struct volumeFractionSmoothingDictionary
                {
                    bool smoothVolumeFraction_ = false;
                    label smoothingIterations_ = 3;
                    scalar fourierNumber_ = 0.25;
                };

                struct meshMotionDictionary
                {
                    bool freezePerTimestep_ = true;
                    label maxSmoothingIterations_ = 5;
                };

                subIterationsDictionary subIterations_;
                accelerationDictionary acceleration_;
                volumeFractionSmoothingDictionary volumeFractionSmoothing_;
                meshMotionDictionary meshMotion_;
            };

            pressureLevelInformationDictionary pressureLevelInformation_;
            equationControlsDictionary equationControls_;
        };

        struct expertParametersDictionary
        {
            bool disableMomentumPredictor_ = false;
            bool consistent_ = false;
            bool limitGradients_ = false;
            bool correctGradients_ = false;
            bool incrementalGradientChange_ = true;
            bool falseMassAccumulation_ = true;
            bool fractionalStepMethod_ = false;
            bool coriolisProductionTurbulence_ = false;
            bool bodyForceRedistribution_ = true;
            bool geometricWallDistanceCalculation_ = false;
            bool strongDirichletWallScale_ = false;
            scalar volumeFractionBlendingFactorMax_ = 2.0;
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
    };

    solverControlDictionary solverControl_;
    outputControlDictionary outputControl_;
    restartControlDictionary restartControl_;
};

class controls
{
public:
    // Constructors

    controls();

    // Destructor

    ~controls();

    // IO

    void read(YAML::Node inputNode);

    // Access

    const solverDictionary& solverRef() const
    {
        return solver_;
    };

    bool isReducedStencil() const;

    bool isTransient() const;

    bool isHighResolution() const;

    bool isHighResolutionTurbulenceNumerics() const;

    label getNumberOfStates() const;

    Profiler& getProfiler()
    {
        return profiler_;
    }

    stk::util::ParameterList& getRestartParam();

    void setRestartParam();

    void deserializeRestartParam(const stk::io::StkMeshIoBroker& io_broker);

    scalar getTotalTime() const;

    void setTimestep(const scalar dt);

    scalar getTimestep(const int i = 0) const;

    scalar getPhysicalTimescale() const;

    label getTimeStepCount() const;

    void advanceAndSetTimestep();

    void updateMaxCourant(const scalar maxCourant);

    scalar getMaxCourant() const;

    void resetMaxCourant();

    // Run loop-related public variables

    label iter{0};

    label globalIter{0};

    scalar time{0};

private:
    analysisTypeDictionary analysisType_;

    solverDictionary solver_;

    Profiler profiler_;

    stk::util::ParameterList restartParameter_;

    scalar maxCourant_{-1};

    int timestepPosition_(const int i) const;

    scalar periodicIntervalTimestep_() const;

    scalar specifiedIntervalTimestep_();
};

} // namespace accel

#endif // CONTROLS_H
