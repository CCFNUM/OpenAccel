// File       : simulation.h
// Created    : Tue Jun 11 2024 15:06:38 (+0100)
// Author     : Mhamad Mahdi Alloush
// Description: This is the universe in which a full simulation can be made
// Copyright (c) 2024 CCFNUM, Lucerne University of Applied Sciences and Arts.
// SPDX-License-Identifier: BSD-3-Clause

#ifndef SIMULATION_H
#define SIMULATION_H

// code
#include "controls.h"
#include "equation.h"
#include "meshMotion.h"
#include "overrides.h"
#include "postProcess.h"
#include "realm.h"
#include "rigidBody.h"
#include "types.h"
#include "wallDistance.h"

// external
#include "gplot++.h"

namespace accel
{

class simulation
{
public:
    struct residualPlotItem
    {
        std::string data_file;
        std::string legend_name;
        unsigned int xdata_idx; // 0-based data column index
        unsigned int ydata_idx; // 0-based data column index
    };

private:
    friend class domain;

    // Primary members

    // analysis type, solver settings, etc.
    std::unique_ptr<controls> controlsPtr_ = nullptr;

    // stk mesh containers (IO borker, bulk data, meta data, etc), and geometric
    // fields
    std::unique_ptr<mesh> meshPtr_ = nullptr;

    // mesh motion for managing mesh deformation or rigid motion
    std::unique_ptr<meshMotion> meshMotionPtr_ = nullptr;

    // rigid body motion
    std::vector<std::unique_ptr<rigidBody>> rigidBodyVector_;

    // wall distance calcualtor
    std::unique_ptr<wallDistance> wallDistancePtr_ = nullptr;

    // the realm maintains a set of equations which are solved on the mesh that
    // it references
    std::unique_ptr<realm> realmPtr_;

    // domain: contains region specific properties. Each domain has a unique
    // mesh zone. The mesh zone is owned by the mesh.
    std::vector<std::shared_ptr<domain>> domainVector_;

    // equations that are solved on the realm
    std::vector<std::unique_ptr<equation>> equationVector_;

    std::unique_ptr<postProcess> postProcessPtr_ = nullptr;

    std::unordered_map<std::string, label> registeredMaterialsMap_;

    std::unique_ptr<overrides> overridesPtr_ = nullptr;

    // IO
    std::set<std::string> restartFields_;
    label io_last_restart_;
    label io_last_results_;
    label io_write_counter_;
    scalar io_write_time_;
    scalar io_last_results_time_{-1.0};

    // path to input file: file contains all information of simulation
    fs::path inputFilePath_;

    YAML::Node inputNode_;

    // manage the output exodus file
    std::unique_ptr<Ioss::PropertyManager> resultsPropertyManagerPtr_ = nullptr;
    std::unique_ptr<Ioss::PropertyManager> restartPropertyManagerPtr_;

    // index of the results file
    size_t resultsFileIndex_;
    size_t restartFileIndex_;

    bool plotRes_ = false;
    bool plotResInitialized_ = false;

    bool printScales_ = false;

    int verbose_;

    std::unique_ptr<Gnuplot> gp_ptr_ = nullptr;

    std::vector<residualPlotItem> plot_items_;

    stk::ParsedOptions args_;

    // Scaled residuals for equations: used here to check convergence

    std::unordered_map<std::string, scalar> scalarResiduals_;

    // Private methods

    // read and create settings/directories
    void createControls_();

    void createDirectories_();

    // read and create mesh
    void createMesh_();

    // create realm
    void createRealm_();

    // read and create equations
    void createEquations_();

    // create add-on equations or methods
    void createAddons_();

    void collectEquations_();

    bool findEquation_(equationID);

    // read and create domains
    void createDomains_();

    // post process objects
    void createPostProcess_();

    // physical overrides
    void createOverrides_();

    // correct boundary values for output fields
    void applyArbitraryFieldCorrections_();

    // correct boundary values for output fields
    void correctBoundaryValues_();

    // restore conservative boundary values for output fields
    void restoreBoundaryValues_();

    // Create components
    void create_();

    // Initialize components
    void initialize_();

    // Initialize output structure and detect output fields to be written
    void initializeOutput_();

public:
    // Constructors

    simulation(const int argc, const char* argv[]);

    // Destructor

    ~simulation();

    // Operations

    void run();

    void runSteadyState();

    void runTransient();

    void assembleAndSolveSystems();

    void preWork();

    void postWork();

    bool writeNow(const label freq)
    {
        return (freq > 0) && (io_write_counter_ % freq == 0);
    }

    bool writeNowTime(const scalar interval)
    {
        // Use small tolerance to handle floating-point accumulation errors
        const scalar tolerance = interval * 1.0e-6;
        return (interval > 0.0) && (io_last_results_time_ < 0.0 ||
                                    (io_write_time_ - io_last_results_time_) >=
                                        (interval - tolerance));
    }

    void writeResults(const bool write_condition);

    void writeRestart(const bool write_condition);

    void printSolverHeader(const int argc, const char* argv[]);

    void printSolverFooter();

    bool checkConvergence();

    fs::path getSimulationDirectory() const;

    fs::path getPostProcessingDirectory() const;

    fs::path getResidualDirectory() const;

    void initializeResidualPlot();

    void updateResidualPlot();

    std::string residualPlotCommand_() const;

    void plotResiduals();

    void addPlotItem(const residualPlotItem& item);

    Profiler& getProfiler()
    {
        return controlsRef().getProfiler();
    }

    label getIterationCount() const
    {
        return controlsRef().iter;
    }

    label getGlobalIterationCount() const
    {
        return controlsRef().globalIter;
    }

    scalar getSimulationTime() const
    {
        return controlsRef().time;
    }

    YAML::Node getYAMLSimulationNode() const
    {
        return inputNode_["simulation"];
    }

    YAML::Node getYAMLPhysicalAnalysisNode() const
    {
        return inputNode_["simulation"]["physical_analysis"];
    }

    YAML::Node getYAMLSolverControlNode() const
    {
        return inputNode_["simulation"]["solver"]["solver_control"];
    }

    YAML::Node getYAMLMaterialLibrary() const
    {
        return inputNode_["simulation"]["material_library"];
    }

    void registerRestartField(const std::string& fieldName)
    {
        restartFields_.insert(fieldName);
    }

    void registerMaterial(const std::string& materialName)
    {
        if (!registeredMaterialsMap_.contains(materialName))
        {
            registeredMaterialsMap_[materialName] =
                registeredMaterialsMap_.size();
        }
    }

    label materialIndex(const std::string& materialName) const
    {
        auto it = registeredMaterialsMap_.find(materialName);
        if (it != registeredMaterialsMap_.end())
        {
            return it->second;
        }
        else
        {
            errorMsg("Material not registered");
        }

        return -1;
    }

    std::string materialName(label materialIndex) const
    {
        for (const auto& pair : registeredMaterialsMap_)
        {
            if (pair.second == materialIndex)
            {
                return pair.first; // Return the key if the value matches
            }
        }

        errorMsg("Material of index " + std::to_string(materialIndex) +
                 " not registered");

        return "";
    }

    label nRegisteredMaterials() const
    {
        return registeredMaterialsMap_.size();
    }

    // Access

    label nDomains() const;

    domain& domainRef(label iZone);

    const domain& domainRef(label iZone) const;

    domain& domainRef(std::string name);

    const domain& domainRef(std::string name) const;

    controls* controlsPtr();

    const controls* controlsPtr() const;

    controls& controlsRef();

    const controls& controlsRef() const;

    mesh* meshPtr();

    const mesh* meshPtr() const;

    mesh& meshRef();

    const mesh& meshRef() const;

    const std::vector<std::shared_ptr<domain>>& domainVector() const
    {
        return domainVector_;
    }

    meshMotion* meshMotionPtr();

    const meshMotion* meshMotionPtr() const;

    meshMotion& meshMotionRef();

    const meshMotion& meshMotionRef() const;

    overrides* overridesPtr();

    const overrides* overridesPtr() const;

    overrides& overridesRef();

    const overrides& overridesRef() const;

    label nRigidBodies() const
    {
        return rigidBodyVector_.size();
    }

    rigidBody& rigidBodyRef(label iRigidBody)
    {
        return *rigidBodyVector_[iRigidBody].get();
    }

    const rigidBody& rigidBodyRef(label iRigidBody) const
    {
        return *rigidBodyVector_[iRigidBody].get();
    }

    rigidBody* rigidBodyPtr(std::string name)
    {
        for (label iRigidBody = 0; iRigidBody < nRigidBodies(); iRigidBody++)
        {
            if (rigidBodyRef(iRigidBody).name() == name)
            {
                return rigidBodyVector_[iRigidBody].get();
            }
        }
        errorMsg("rigid body not available");
        return nullptr;
    }

    const rigidBody* rigidBodyPtr(std::string name) const
    {
        for (label iRigidBody = 0; iRigidBody < nRigidBodies(); iRigidBody++)
        {
            if (rigidBodyRef(iRigidBody).name() == name)
            {
                return rigidBodyVector_[iRigidBody].get();
            }
        }
        errorMsg("rigid body not available");
        return nullptr;
    }
};

} // namespace accel

#endif // SIMULATION_H
