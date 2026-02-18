// File : simulationIO.cpp
// Created : Tue Jun 11 2024 15:06:38 (+0100)
// Author : Mhamad Mahdi Alloush
// Description:
// Copyright (c) 2024 CCFNUM, Lucerne University of Applied Sciences and Arts.
// SPDX-License-Identifier: BSD-3-Clause

// code
#include "domain.h"
#include "elementField.h"
#include "equation.h"
#include "fluidEquationsIncludes.h"
#include "git_revision.h"
#include "mesh.h"
#include "messager.h"
#include "simulation.h"
#include "types.h"

#ifdef HAS_GNUPLOT
#include <algorithm>
#include <sstream>
#endif

namespace accel
{

#ifdef HAS_GNUPLOT
std::string simulation::residualPlotCommand_() const
{
    std::ostringstream cmd;
    cmd << "plot ";

    for (size_t i = 0; i < plot_items_.size(); ++i)
    {
        const auto& item = plot_items_[i];

        std::string legend = item.legend_name;
        std::replace(legend.begin(), legend.end(), '\'', ' ');

        cmd << "'" << item.data_file << "' using " << (item.xdata_idx + 1)
            << ":" << (item.ydata_idx + 1) << " with linespoints title '"
            << legend << "'";

        if (i + 1 < plot_items_.size())
        {
            cmd << ", ";
        }
    }

    return cmd.str();
}
#endif

// Read and register information

void simulation::create_()
{
    createDirectories_();
    createControls_(); // must call before createMesh_ for correct restarts!
    createMesh_();
    createRealm_();
    createDomains_();
    createEquations_();
    createAddons_();
    createPostProcess_();
    createOverrides_();
}

void simulation::createControls_()
{
    controlsPtr_ = std::make_unique<controls>();
    controlsPtr_->read(getYAMLSimulationNode());
}

void simulation::createMesh_()
{
    meshPtr_ = std::make_unique<mesh>(controlsPtr_.get());
    meshPtr_->read(inputNode_);
    meshPtr_->setup();
}

void simulation::createRealm_()
{
    // Allocate the realm (currently one default realm). A realm is a container
    // in which equations and domains are manipulated.
    realmPtr_ = std::make_unique<::accel::realm>(this, "default");
    if (messager::master())
    {
        std::cout << "Created realm: `" + realmPtr_->name() + "`:\n";
    }
}

void simulation::createDomains_()
{
    // steps:
    // 1. Create domains (in same order as mesh zones were created before this
    // call). Domains will reference equation(s) and material they act on.
    // 2. Assign interfaces to domains. In case of two domains separated by an
    // interface, the domain with the lower index will be assigned the
    // interface.
    const YAML::Node sim_root = getYAMLSimulationNode();

    // 1.) create domains (at this point mesh did setup zones based on
    // sim_root["physical_analysis"]["domains"] YAML sequence, so
    // this node is guaranteed to exist, not checking)
    if (messager::master())
    {
        std::cout << "Setting up simulation domains:\n";
    }

    std::map<std::string, label> name2index; // auxiliary data structure for
    // interface assignment

    const YAML::Node domains = sim_root["physical_analysis"]["domains"];
    for (YAML::const_iterator domain = domains.begin(); domain != domains.end();
         ++domain)
    {
        const label index = static_cast<label>(domainVector_.size());

        std::string name;
        if ((*domain)["name"])
        {
            name = (*domain)["name"].template as<std::string>();
        }
        else
        {
            errorMsg("simulation: `name` field not set for domain index " +
                     std::to_string(index));
        }

        if (messager::master())
        {
            std::cout << "\t adding domain `" + name +
                             "` (index = " + std::to_string(index) + ")\n";
        }

        const auto it = name2index.find(name);
        if (it != name2index.end())
        {
            // domain name must be unique
            errorMsg("simulation: ambiguous name `" + name +
                     "` for domain index " + std::to_string(index));
        }
        name2index[name] = index;

        // create domain
        domainVector_.push_back(std::make_shared<::accel::domain>(
            this, meshRef().zonePtr(index), *domain));
    }

    // check consistency
    const label n_domains = static_cast<label>(domainVector_.size());
    if (meshRef().nZones() != n_domains)
    {
        errorMsg("simulation: domain and mesh zone mismatch during domain "
                 "creation (nZones = " +
                 std::to_string(meshRef().nZones()) +
                 "; nDomains = " + std::to_string(n_domains) + ")");
    }

    // check if there is a dupilcate domain name
    for (const auto& domain1 : domainVector_)
    {
        for (const auto& domain2 : domainVector_)
        {
            if (domain1->index() != domain2->index() &&
                domain1->name() == domain2->name())
            {
                errorMsg("detected duplicate domain names");
            }
        }
    }
}

void simulation::createEquations_()
{
    // steps:
    // 1. Create equations (equations are solved on the mesh that is referenced
    // by the realm)
    // 2. Link domains to equations (equations will only be assembled for linked
    // domains)
    // 3. Physics-sanity check: every fluid-fluid interface must be pairing two
    // domains with exactly same physics

    // 1.) create equations
    const YAML::Node sim_root = getYAMLSimulationNode();
    collectEquations_();

    // 2.) link domains to equations
    for (auto& equation : equationVector_)
    {
        const equationID id = equation->getID();

        // equations in queue must have a valid ID
        assert(id != equationID::noID);
        for (auto domain : domainVector_)
        {
            if (domain->hasEquation(id))
            {
                // register domain for equation
                equation->addDomain(domain);
            }
        }
    }

    if (messager::master())
    {
        std::cout << "\tEquation queue = [\n";
        for (const auto& equation : equationVector_)
        {
            std::cout << "\t\t" << equation->name() << '\n';
        }
        std::cout << "\t]\n";
    }
}

void simulation::createAddons_()
{
    // mesh motion
    if (this->meshRef().anyZoneMeshDeforming() ||
        this->meshRef().anyZoneMeshTransforming())
    {
        meshMotionPtr_ = std::make_unique<meshMotion>(realmPtr_.get());
    }

    // rigid bodies
    if (getYAMLSimulationNode()["physical_analysis"]["rigid_bodies"])
    {
        const YAML::Node rigidBodiesBlockArray =
            getYAMLSimulationNode()["physical_analysis"]["rigid_bodies"];

        for (const auto& rigidBodyBlock : rigidBodiesBlockArray)
        {
            std::unique_ptr<rigidBody> rigidBodyPtr =
                std::make_unique<rigidBody>(realmPtr_.get());
            rigidBodyPtr->read(rigidBodyBlock);
            rigidBodyVector_.push_back(std::move(rigidBodyPtr));
        }
    }

    // wall distance
    for (const auto& domain : domainVector_)
    {
        if (domain->isWallDistanceRequired())
        {
            wallDistancePtr_ = std::make_unique<wallDistance>(realmPtr_.get());
            break;
        }
    }
}

void simulation::collectEquations_()
{
    for (const auto& domain : domainVector_)
    {
        // segregated NS-pcorr equations
        if (domain->hasEquation(equationID::segregatedFlow))
        {
            if (!findEquation_(equationID::segregatedFlow))
            {
                equationVector_.push_back(
                    std::make_unique<segregatedFlowEquations>(realmPtr_.get()));
            }
        }

        // segregated free surface flow equations
        if (domain->hasEquation(equationID::segregatedFreeSurface))
        {
            if (!findEquation_(equationID::segregatedFreeSurface))
            {
                equationVector_.push_back(
                    std::make_unique<segregatedFreeSurfaceFlowEquations>(
                        realmPtr_.get()));
            }
        }

        // segregated shear stress transport equations
        if (domain->hasEquation(equationID::segregatedShearStressTransport))
        {
            if (!findEquation_(equationID::segregatedShearStressTransport))
            {
                equationVector_.push_back(
                    std::make_unique<segregatedShearStressTransportEquations>(
                        realmPtr_.get()));
            }
        }

        // segregated transition shear stress transport equations
        if (domain->hasEquation(
                equationID::segregatedTransitionShearStressTransport))
        {
            if (!findEquation_(
                    equationID::segregatedTransitionShearStressTransport))
            {
                equationVector_.push_back(
                    std::make_unique<
                        segregatedTransitionShearStressTransportEquations>(
                        realmPtr_.get()));
            }
        }

        // segregated correlation-based transition shear stress transport
        // equations (Menter 2015)
        if (domain->hasEquation(
                equationID::
                    segregatedCorrelationTransitionShearStressTransport))
        {
            if (!findEquation_(
                    equationID::
                        segregatedCorrelationTransitionShearStressTransport))
            {
                equationVector_.push_back(
                    std::make_unique<
                        segregatedCorrelationTransitionShearStressTransportEquations>(
                        realmPtr_.get()));
            }
        }

        // segregated k-epsilon equations
        if (domain->hasEquation(equationID::segregatedKEpsilon))
        {
            if (!findEquation_(equationID::segregatedKEpsilon))
            {
                equationVector_.push_back(
                    std::make_unique<segregatedKEpsilonEquations>(
                        realmPtr_.get()));
            }
        }

        // thermal energy equation
        if (domain->hasEquation(equationID::thermalEnergy))
        {
            if (!findEquation_(equationID::thermalEnergy))
            {
#ifdef WITH_THERMAL_TEMPERATURE
                equationVector_.push_back(
                    std::make_unique<thermalTemperatureEquation>(
                        realmPtr_.get()));
#else
                equationVector_.push_back(
                    std::make_unique<thermalEnergyEquation>(realmPtr_.get()));
#endif
            }
        }

        // total energy equation
        if (domain->hasEquation(equationID::totalEnergy))
        {
            if (!findEquation_(equationID::totalEnergy))
            {
                equationVector_.push_back(
                    std::make_unique<totalEnergyEquation>(realmPtr_.get()));
            }
        }

        // solid displacement equation
        if (domain->hasEquation(equationID::solidDisplacement))
        {
            if (!findEquation_(equationID::solidDisplacement))
            {
                equationVector_.push_back(
                    std::make_unique<solidDisplacementEquation>(
                        realmPtr_.get()));
            }
        }
    }
}

bool simulation::findEquation_(equationID id)
{
    for (const auto& eq : equationVector_)
    {
        if (eq->getID() == id)
        {
            return true;
        }
    }
    return false;
}

void simulation::createPostProcess_()
{
    postProcessPtr_ = std::make_unique<postProcess>(
        realmPtr_.get(), this->getPostProcessingDirectory());
    postProcessPtr_->read(getYAMLSimulationNode());
}

void simulation::createOverrides_()
{
    overridesPtr_ = std::make_unique<overrides>(realmPtr_.get());
    overridesPtr_->read(getYAMLSimulationNode());
}

void simulation::initializeOutput_()
{
    const bool is_restart =
        controlsRef().solverRef().restartControl_.isRestart_;

    // results
    resultsPropertyManagerPtr_ = std::make_unique<Ioss::PropertyManager>();
    resultsFileIndex_ = meshRef().ioBrokerPtr()->create_output_mesh(
        controlsRef().solverRef().outputControl_.filePath_,
        controlsRef().solverRef().outputControl_.writeMode_,
        *resultsPropertyManagerPtr_.get());
    meshRef().ioBrokerPtr()->use_nodeset_for_part_nodes_fields(
        resultsFileIndex_, false);

    // Add result fields
    const std::vector<std::string>& outputFieldNames =
        controlsRef().solverRef().outputControl_.outputFields_;
    const auto& metaData = meshRef().metaDataRef();

    const stk::mesh::FieldVector& fields = metaData.get_fields();
    std::vector<std::string> registeredFieldNames;
    for (const auto& field : fields)
    {
        auto iter = std::find(
            outputFieldNames.begin(), outputFieldNames.end(), field->name());
        if (iter != outputFieldNames.end())
        {
            meshRef().ioBrokerPtr()->add_field(
                resultsFileIndex_, *field, field->name());
            registeredFieldNames.push_back(field->name());
        }
    }

    // check queried fields that are not available
    if (messager::master())
    {
        for (const auto& outputFieldName : outputFieldNames)
        {
            auto iter = std::find(registeredFieldNames.begin(),
                                  registeredFieldNames.end(),
                                  outputFieldName);
            if (iter == registeredFieldNames.end())
            {
                warningMsg(
                    "simulation: output field `" + outputFieldName +
                    "` is not defined in mesh database: ignored for output");
            }
        }
    }

#ifndef NDEBUG
    // unconditional write of rank ID cell field
    meshRef().ioBrokerPtr()->add_field(
        resultsFileIndex_,
        *stk::mesh::get_field_by_name<int>(::accel::mesh::rank_ID, metaData),
        ::accel::mesh::rank_ID);
#endif /* NDEBUG */

    // Write result fields at time = t_0
    io_write_counter_ = 0;
    io_last_results_ = -1;
    io_last_results_time_ = -1.0;
    io_write_time_ = controlsRef().time;
    const auto& outputCtrl = controlsRef().solverRef().outputControl_;
    const auto& outFreq = outputCtrl.outputFrequency_;
    const bool results_freq =
        (outFreq.option_ == outputFrequencyType::timeInterval)
            ? (outFreq.timeInterval_ > 0.0)
            : (outFreq.timestepInterval_ > 0);
    const bool restart_write_initial =
        controlsRef().solverRef().restartControl_.writeInitial_;
    writeResults((results_freq && !is_restart) || restart_write_initial);

    // restarts
    io_last_restart_ = -1;

    // compose restart path name
    auto restart_path = this->getSimulationDirectory();
    std::string restart_dir_name("restart."); // default restart dir name
    const size_t len = restart_dir_name.size();
    int nrestart_dirs = 0;
    for (const auto& entry : fs::directory_iterator{restart_path})
    {
        if (restart_dir_name.compare(0, len, entry.path().filename(), 0, len) ==
            0)
        {
            ++nrestart_dirs;
        }
    }
    std::ostringstream oss;
    oss << restart_dir_name << std::setfill('0') << std::setw(3)
        << nrestart_dirs;
    restart_path /= oss.str();
    restart_path /= controlsRef().solverRef().outputControl_.restartFileName_;

    restartPropertyManagerPtr_ = std::make_unique<Ioss::PropertyManager>();
    restartFileIndex_ = meshRef().ioBrokerPtr()->create_output_mesh(
        restart_path,
        stk::io::WRITE_RESTART,
        *restartPropertyManagerPtr_.get());
    for (const auto& field_name : restartFields_)
    {
        stk::mesh::FieldBase* field =
            stk::mesh::get_field_by_name(field_name, metaData);
        assert(field);
        meshRef().ioBrokerPtr()->add_field(
            restartFileIndex_, *field, field_name);
    }

    for (const auto& param : controlsRef().getRestartParam())
    {
        assert(param.second.toRestartFile);
        meshRef().ioBrokerPtr()->add_global(
            restartFileIndex_, param.first, param.second);
    }

    meshRef()
        .ioBrokerPtr()
        ->get_output_ioss_region(restartFileIndex_)
        ->get_database()
        ->set_cycle_count(
            controlsRef().solverRef().restartControl_.keepNRestartSnapshots_);
}

void simulation::writeResults(const bool write_condition)
{
    if (write_condition && io_write_counter_ != io_last_results_)
    {
        // arbitrary field corrections/post-process: done with care
        applyArbitraryFieldCorrections_();

        // overwrite boundary node values with specified ones (hybrid)
        correctBoundaryValues_();

        // write solution to database file
        meshRef().write(resultsFileIndex_, io_write_time_);

        // restore conservative boundary node values
        restoreBoundaryValues_();

        //
        io_last_results_ = io_write_counter_;
        io_last_results_time_ = io_write_time_;
    }
}

void simulation::writeRestart(const bool write_condition)
{
    if (write_condition && io_write_counter_ != io_last_restart_)
    {
        controlsRef().setRestartParam();

        meshRef().ioBrokerPtr()->begin_output_step(restartFileIndex_,
                                                   io_write_time_);
        meshRef().ioBrokerPtr()->write_defined_output_fields(restartFileIndex_);
        for (const auto& param : controlsRef().getRestartParam())
        {
            assert(param.second.toRestartFile);
            meshRef().ioBrokerPtr()->write_global(
                restartFileIndex_, param.first, param.second);
        }
        meshRef().ioBrokerPtr()->end_output_step(restartFileIndex_);
        io_last_restart_ = io_write_counter_;
    }
}

void simulation::applyArbitraryFieldCorrections_()
{
}

void simulation::correctBoundaryValues_()
{
    if (!controlsRef().solverRef().outputControl_.correctedBoundaryValues_)
        return;

    const std::vector<std::string>& outputFields =
        controlsRef().solverRef().outputControl_.outputFields_;
    const auto& meta_data = meshRef().metaDataRef();
    for (label iField = 0; iField < static_cast<label>(outputFields.size());
         iField++)
    {
        // detect if a node side field is available
        const STKScalarField* nodeSideSTKFieldPtr =
            stk::mesh::get_field_by_name<scalar>(
                outputFields[iField] + "_node_side", meta_data);

        if (nodeSideSTKFieldPtr)
        {
            STKScalarField* stkFieldPtr = stk::mesh::get_field_by_name<scalar>(
                outputFields[iField], meta_data);

            // prev field will be used to store conservative boundary data to be
            // able to restore them at a later stage
            STKScalarField* stkFieldPtrPrev =
                stk::mesh::get_field_by_name<scalar>(
                    outputFields[iField] + "_prev_iter", meta_data);

            label fieldDim = stkFieldPtr->max_size();

            const stk::mesh::Selector selAllNodes =
                meshPtr_->metaDataRef().universal_part() &
                stk::mesh::selectField(*nodeSideSTKFieldPtr);
            const stk::mesh::BucketVector& nodeBuckets =
                meshPtr_->bulkDataRef().get_buckets(stk::topology::NODE_RANK,
                                                    selAllNodes);

            for (stk::mesh::Bucket::size_type ib = 0; ib < nodeBuckets.size();
                 ++ib)
            {
                const stk::mesh::Bucket& bucket = *nodeBuckets[ib];
                const stk::mesh::Bucket::size_type n_entities = bucket.size();

                for (stk::mesh::Bucket::size_type i = 0; i < n_entities; ++i)
                {
                    const stk::mesh::Entity node = bucket[i];
                    const scalar* nodeSideVal =
                        stk::mesh::field_data(*nodeSideSTKFieldPtr, node);
                    scalar* val = stk::mesh::field_data(*stkFieldPtr, node);
                    scalar* valPrev =
                        stk::mesh::field_data(*stkFieldPtrPrev, node);

                    for (label i = 0; i < fieldDim; i++)
                    {
                        valPrev[i] = val[i];
                        val[i] = nodeSideVal[i];
                    }
                }
            }
        }
    }
}

void simulation::restoreBoundaryValues_()
{
    if (!controlsRef().solverRef().outputControl_.correctedBoundaryValues_)
        return;

    const std::vector<std::string>& outputFields =
        controlsRef().solverRef().outputControl_.outputFields_;
    const auto& meta_data = meshRef().metaDataRef();
    for (label iField = 0; iField < static_cast<label>(outputFields.size());
         iField++)
    {
        // detect if a node side field is available
        const STKScalarField* nodeSideSTKFieldPtr =
            stk::mesh::get_field_by_name<scalar>(
                outputFields[iField] + "_node_side", meta_data);

        if (nodeSideSTKFieldPtr)
        {
            STKScalarField* stkFieldPtr = stk::mesh::get_field_by_name<scalar>(
                outputFields[iField], meta_data);

            const STKScalarField* stkFieldPtrPrev =
                stk::mesh::get_field_by_name<scalar>(
                    outputFields[iField] + "_prev_iter", meta_data);

            label fieldDim = stkFieldPtr->max_size();

            const stk::mesh::Selector selAllNodes =
                meshPtr_->metaDataRef().universal_part() &
                stk::mesh::selectField(*nodeSideSTKFieldPtr);
            const stk::mesh::BucketVector& nodeBuckets =
                meshPtr_->bulkDataRef().get_buckets(stk::topology::NODE_RANK,
                                                    selAllNodes);

            for (stk::mesh::Bucket::size_type ib = 0; ib < nodeBuckets.size();
                 ++ib)
            {
                const stk::mesh::Bucket& bucket = *nodeBuckets[ib];
                const stk::mesh::Bucket::size_type n_entities = bucket.size();

                for (stk::mesh::Bucket::size_type i = 0; i < n_entities; ++i)
                {
                    const stk::mesh::Entity node = bucket[i];
                    scalar* val = stk::mesh::field_data(*stkFieldPtr, node);
                    const scalar* valPrev =
                        stk::mesh::field_data(*stkFieldPtrPrev, node);

                    for (label i = 0; i < fieldDim; i++)
                    {
                        val[i] = valPrev[i];
                    }
                }
            }
        }
    }
}

void simulation::createDirectories_()
{
    plotResInitialized_ = false;
    if (messager::master())
    {
        // Set booleans
        plotRes_ = args_.count("residuals");
        printScales_ = args_.count("scales");

        // Post processing directory
        fs::create_directory(this->getPostProcessingDirectory());

        // Residuals directory and file
        fs::create_directories(this->getResidualDirectory());
    }
}

void simulation::initializeResidualPlot()
{
#ifdef HAS_GNUPLOT
    try
    {
        gp_ptr_ = std::make_unique<Gnuplot>();
        gp_ptr_->sendcommand("set title 'Data Monitoring'");
        gp_ptr_->sendcommand("set xlabel 'Iter'");
        gp_ptr_->sendcommand("set ylabel 'scaled residuals'");
        gp_ptr_->sendcommand("set grid");
        gp_ptr_->sendcommand("set logscale y");
        gp_ptr_->sendcommand("set format y \"1e{%T}\"");
        gp_ptr_->sendcommand("set yrange [1e-10:1]");

        if (!plot_items_.empty())
        {
            gp_ptr_->sendcommand(residualPlotCommand_());
            gp_ptr_->show();
        }
    }
    catch (const std::exception& ex)
    {
        std::cout << ex.what() << std::endl;
        return;
    }
#else
    std::cout << "Residual plotting requested but gnuplot support is "
                 "disabled at compile time.\n";
#endif
}

void simulation::updateResidualPlot()
{
#ifdef HAS_GNUPLOT
    if (gp_ptr_ && !plot_items_.empty())
    {
        gp_ptr_->sendcommand(residualPlotCommand_());
        gp_ptr_->show();
    }
#endif
}

fs::path simulation::getSimulationDirectory() const
{
    return fs::absolute(inputFilePath_).parent_path();
}

fs::path simulation::getPostProcessingDirectory() const
{
    fs::path pp_path(this->getSimulationDirectory() / "postProcessing");
    return pp_path;
}

fs::path simulation::getResidualDirectory() const
{
    fs::path res_path(this->getPostProcessingDirectory() / "0" / "residuals");
    return res_path;
}

void simulation::plotResiduals()
{
    if (messager::master())
    {
        if (plotRes_)
        {
#ifdef HAS_GNUPLOT
            // Initialize window first
            if (!plotResInitialized_)
            {
                initializeResidualPlot();
                plotResInitialized_ = true;
            }

            // Plot residuals
            updateResidualPlot();
#endif
        }
    }
}

void simulation::addPlotItem(const residualPlotItem& item)
{
#ifdef HAS_GNUPLOT
    plot_items_.push_back(item);
#endif
}

// clang-format off

void simulation::printSolverHeader(const int argc, const char* argv[])
{
 if (messager::master())
 {
 std::cout << "╔══════════════════════════════════════════════════════════════════════╗" << std::endl;
 std::cout << "║                          OpenAccel " << SPATIAL_DIM
 << "D                                ║" << std::endl;
 std::cout << "║          Parallel fluid flow CFD package based on CVFEM              ║" << std::endl;
 std::cout << "║                    |<| Powered by Trilinos |>|                       ║" << std::endl;
 std::cout << "╚══════════════════════════════════════════════════════════════════════╝" << std::endl;
 std::cout << "Revision: " << accel::git_revision << std::endl;
 std::cout << "Command line:";
 for (int i = 0; i < argc; i++) {
 std::cout << " " << argv[i];
 }
 std::cout << '\n';
 std::cout << std::endl;
 }
}

void simulation::printSolverFooter()
{
 if (messager::master())
 {
 std::cout << "╔══════════════════════════════════════════════════════════════════════╗" << std::endl;
 std::cout << "║                       Simulation is complete                         ║" << std::endl;
 std::cout << "╚══════════════════════════════════════════════════════════════════════╝" << std::endl;
 std::cout << std::endl;
 }
}

// clang-format on

} // namespace accel
