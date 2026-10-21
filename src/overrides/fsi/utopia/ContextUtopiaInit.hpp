// File       : ContextUtopiaInit.hpp
// Created    : Mon Oct 14 2024 10:05:22 (+0200)
// Author     : Fabian Wermelinger
// Description: Utopia FSI initialization
// Copyright 2024 CCFNUM HSLU T&A. All Rights Reserved.

// no include guard: must only be included in ContextUtopia.h

#include <cassert>
#include <vector>

#include "utopia_stk_DofMap.hpp"
#ifdef UTOPIA_ENABLE_SFEM
#include "utopia_sfem_stk_SDFObstacle.hpp"
#endif // UTOPIA_ENABLE_SFEM

namespace accel
{
namespace fsi
{

template <size_t N>
void ContextUtopia<N>::init(::stk::mesh::BulkData& data,
                            const YAML::Node& fsi_node,
                            simulation& sim)
{
    UTOPIA_TRACE_REGION_BEGIN("ContextUtopia::init");

    const controls& ctrl = sim.controlsRef();
    assert(ctrl.isTransient()); // Utopia currently only supports
                                // transient simulations

    stime_ = std::make_shared<UtopiaTime>();
    stime_->restart(ctrl.time);
    stime_->delta() = ctrl.getTimestep();

    this->restart_enabled_ = ctrl.solverRef().restartControl_.isRestart_;

    // structure restarts filename (hardcoded here)
    const std::string restart_filename("restart_structure.bin");

    // setup norm_increment_ output stream
    sim_ = &sim; // keep a const reference for iteration and time
    const fs::path fname =
        sim.getPostProcessingDirectory() / "utopia_norm_increment.out";
    if (this->commRank() == 0)
    {
        this->norm_inc_stream_ = std::make_unique<std::ofstream>(
            fname, this->restart_enabled_ ? std::ios::app : std::ios::trunc);
    }
    else
    {
        this->norm_inc_stream_ =
            std::make_unique<::linearSolver::onullstream>();
    }

    if (!this->restart_enabled_)
    {
        auto& fout = *this->norm_inc_stream_;
        auto now = std::chrono::system_clock::now();
        auto in_time_t = std::chrono::system_clock::to_time_t(now);

        fout << COMMENT << "CVFEM solver timestamp: "
             << std::put_time(std::localtime(&in_time_t), "%c\n");
        fout << COMMENT << "Version: " << accel::PROJECT_VERSION << '\n';
        fout << COMMENT << "Git hash: " << accel::GIT_HASH << '\n';
        fout << COMMENT << "Git describe: " << accel::GIT_DESCRIBE << '\n';
        fout << COMMENT << '\n';
        fout << COMMENT << "Utopia norm increment:\n";
        fout << COMMENT << "global_iter" << '\t' << "iter" << '\t'
             << "sim_time[s]" << '\t' << "norm_increment" << '\n';
    }

    // main IO objects
    io_ = utopia::make_unique<IO>(space_structure_);
    io_restart_ = utopia::make_unique<IO>(space_structure_);
    assert(io_);
    assert(io_restart_);

    io_->set_domain_decomposition_method(
        ctrl.solverRef()
            .solverControl_.advancedOptions_.domainDecomposition_.method_);
    io_restart_->set_domain_decomposition_method(
        ctrl.solverRef()
            .solverControl_.advancedOptions_.domainDecomposition_.method_);

    // restart load database
    if (this->restart_enabled_)
    {
        restart_db_ = utopia::make_unique<IO>(space_structure_);
        assert(restart_db_);

        fs::path restart_fluid =
            ctrl.solverRef().restartControl_.inputFilePath_;
        assert(!restart_fluid.empty());
        fs::path dir = restart_fluid.parent_path();
        auto restart_structure = dir / restart_filename;

        controls::hasSTKDatabase(restart_structure);
        restart_db_->set_input_path(restart_structure.c_str());
        restart_db_->set_domain_decomposition_method(
            ctrl.solverRef()
                .solverControl_.advancedOptions_.domainDecomposition_.method_);

        // Remove the `decomposition_method` field to avoid MeshIO::read_minimal
        // overriding it later in the Utopia API.
        YAML::Node mesh_node = fsi_node["mesh"];
        mesh_node.remove("decomposition_method");
    }

#ifdef UTOPIA_ENABLE_SFEM
    utopia::register_sfem_stk_contact();
#endif // UTOPIA_ENABLE_SFEM

#ifdef FLUYA_UTOPIA_YAML_ENABLED

    std::unique_ptr<utopia::YAMLInput> in =
        std::make_unique<utopia::YAMLInput>(fsi_node);
    this->read(*in);

#else

    auto path = fsi_node.template as<std::string>();
    this->cout() << "path for JSON file: " << path << std::endl;

    if (!path.empty())
    {
        auto in = utopia::open_istream(path);
        if (!in)
        {
            utopia::Utopia::Abort("[Error] Invalid path: " + path + "!");
        }

        this->read(*in);
    }

#endif // FLUYA_UTOPIA_YAML_ENABLED

    if (valid_)
    {
        auto mesh = std::make_shared<utopia::stk::Mesh>();
        mesh->wrap(utopia::make_ref(data.mesh_meta_data()),
                   utopia::make_ref(data));

        if (scale_factor_ != 1.0)
        {
            mesh->scale(scale_factor_);
        }
        // space_fluid_.set_print_map(true);
        space_fluid_.init(mesh);
        // space_fluid_.set_n_var(space_structure.n_var());

        // Add fe variables to fluid space for easier handling
        utopia::FEVar var_velocity;
        var_velocity.name = "velocity";
        var_velocity.n_components = space_fluid_.mesh().spatial_dimension();

        assert(BLOCKSIZE >= SPATIAL_DIM);
        space_fluid_.add_variable(var_velocity);
        space_fluid_.initialize(false);

        valid_ = transfer_.init(utopia::make_ref(space_structure_),
                                utopia::make_ref(space_fluid_));

        if (scale_factor_ != 1.0)
        {
            mesh->scale(1.0 / scale_factor_);
        }

        if (valid_)
        {
            valid_ = coupleStructureWithFluid_();
            assert(valid_);
            if (utopia_verbose_)
            {
                space_fluid_.describe(this->cout());
            }
        }
        else
        {
            utopia::err() << "[Error] INVALID transfer" << std::endl;
            utopia::Utopia::Abort("INVALID transfer");
        }
    }
    else
    {
        utopia::err()
            << "[Error] ContextUtopia::init called on undefined structure!\n";
    }

    if (BLOCKSIZE >= SPATIAL_DIM + 1)
    {
        utopia::FEVar pressure_var;
        pressure_var.name = "pressure";
        pressure_var.n_components = 1;
        space_fluid_.add_variable(pressure_var);
    }

    if (!this->restart_enabled_)
    {
        velocity_.data().set(0.0);
    }

    // output data (results)
    auto output_path = ctrl.getResultsDirectory() / structure_filename_;
    if (!output_relative_)
    {
        output_path = structure_filename_;
        fs::create_directories(output_path.parent_path());
    }
    io_->set_output_path(output_path.c_str());
    if (ctrl.solverRef().outputControl_.writeMode_ == stk::io::APPEND_RESULTS)
    {
        controls::hasSTKDatabase(output_path);
        io_->set_output_mode("APPEND");
    }
    else
    {
        io_->set_output_mode("OVERWRITE");
    }

    // restart output data
    auto restart_path = ctrl.getRestartDirectory() / restart_filename;
    io_restart_->set_output_path(restart_path.c_str());
    io_restart_->set_output_mode("RESTART");

    // register output generators
    io_restart_->set_restart_cycle_count(
        ctrl.solverRef().restartControl_.keepNRestartSnapshots_);

    integrator_ = dynamic_cast<utopia::NewmarkIntegrator<FunctionSpace>*>(
        function_.get());
    assert(integrator_);

    integrator_->register_output(*io_);
    sim.registerComplementResult([this](const label /* step */, const scalar t)
    {
        this->integrator_->update_output(*this->io_);
        this->io_->write(-1 /* dummy */, t);
    });

    integrator_->register_output(*io_restart_);
    io_restart_->register_output_field(this->velocity_);
    sim.registerComplementRestart([this](const label /* step */, const scalar t)
    {
        this->integrator_->update_output(*this->io_restart_);
        this->io_restart_->update_output_field(this->velocity_);
        this->io_restart_->write(-1 /* dummy */, t);
    });

    // clear restart database (no longer needed at this point)
    restart_db_ = nullptr;

    UTOPIA_TRACE_REGION_END("ContextUtopia::init");
}

template <size_t N>
void ContextUtopia<N>::read(utopia::Input& in)
{
    in.get("transfer", transfer_);
    in.get("scaleFactor", scale_factor_);
    in.get("verbose", utopia_verbose_);
    valid_ = true;

    in.get("use_strong_velocity_bc", use_strong_velocity_bc_);
    in.get("import_field_data", import_field_data_);
    in.get("cache_hessian_nnz", cache_hessian_nnz_);
    in.get("hessian_use_solid_diagonal", hessian_use_solid_diagonal_);
    in.get("skip_hessian_solid_contribution", skip_hessian_solid_contribution_);
    in.get("inverse_diagonal_scaling", inverse_diagonal_scaling_);
    in.get("use_initial_guess_projection", use_initial_guess_projection_);
    in.get("use_initial_guess_zero", use_initial_guess_zero_);
    in.get("stol", stol_);
    in.get("check_matrix", check_matrix_);
    in.get("reuse_matrix", reuse_matrix_);
    in.get("output",
           [&](utopia::Input& node)
    {
        node.get("relative_to_fluid_output", output_relative_);
        node.get("path", structure_filename_);
    });

    if (reuse_matrix_)
    {
        // In serial petsc performance is bad for reuse_matrix_ == true
        reuse_matrix_ = utopia::mpi_world_size() != 1;
    }

    if (this->restart_enabled_)
    {
        restart_db_->import_all_field_data(true);

        valid_ = restart_db_->open_restart(in);
        restart_db_->load_time_step(stime_->get());
    }
    else
    {
        if (import_field_data_)
        {
            io_->import_all_field_data(true);
            io_->enable_interpolation_mode();

            valid_ = io_->open_input(in);
            io_->load_time_step(stime_->get());
        }
        else
        {
            io_->read(in);
            valid_ = true;
        }
    }

    assert(valid_);
    if (!valid_)
    {
        utopia::Utopia::Abort("MovingObstacle::read: invalid setup");
        return;
    }

    // Velocity
    space_structure_.create_nodal_vector_field(
        space_structure_.mesh().spatial_dimension(), velocity_);
    velocity_.data().set(0);

    space_structure_.create_nodal_vector_field(
        space_structure_.mesh().spatial_dimension(), fluid_velocity_);
    fluid_velocity_.data().set(0);

    {
        auto temp = utopia::ObstacleFEFunctionFactory<FunctionSpace>::make(
                        utopia::make_ref(space_structure_), in)
                        .release();
        function_.reset(dynamic_cast<FEFunction*>(temp));
    }

    if (!function_)
    {
        utopia::Utopia::Abort("Unable to instantiate function correctly!");
    }

    function_->read(in);
    function_->set_time(stime_);
    UVector initial_displacement;
    function_->space()->create_vector(initial_displacement);

    if (this->restart_enabled_)
    {
        function_->setup_IVP(*restart_db_);
        initial_displacement = function_->solution();
        function_->space()->set_overwrite_vector(initial_displacement);
        function_->space()->apply_constraints(initial_displacement);
    }
    else
    {
        if (import_field_data_)
        {
            space_structure_.nodal_field_to_global_vector(initial_displacement);
        }
        else
        {
            initial_displacement.set(0.0);
        }

        function_->space()->set_overwrite_vector(initial_displacement);
        function_->space()->apply_constraints(initial_displacement);
        function_->setup_IVP(initial_displacement);
    }

    function_->create_solution_vector(velocity_.data());

    if (this->restart_enabled_)
    {
        restart_db_->read_nodal(velocity_);
    }
    else
    {
        velocity_.data().set(0.0);

        if (use_initial_guess_projection_)
        {
            function_->initial_guess_for_solver(velocity_.data());
        }
    }

    assert(function_->is_time_dependent());

    if (!function_->is_time_dependent())
    {
        utopia::Utopia::Abort("Utopia: Function must  be time-dependent!");
    }
}

template <size_t N>
void ContextUtopia<N>::filterVelocityBoundaryNodeID(const mesh& mesh,
                                                    velocity& U)
{
    if (!use_strong_velocity_bc_)
    {
        return;
    }

    UTOPIA_TRACE_REGION_BEGIN("ContextUtopia::filterBoundaryNodeID");

    zero_dirichlet_vel_idx_.clear();
    std::vector<bool> visited(mesh.nNodes(), false);

    // enable corrected boundary values in velocity field U
    U.setCorrectedBoundaryNodeValues(true);

    auto filter_local_id = [&](const stk::mesh::BucketVector& buckets)
    {
        const auto& bulk = mesh.bulkDataRef();
        for (const auto bucket : buckets)
        {
            for (const auto& node : *bucket)
            {
                const Index lid = bulk.local_id(node);
                if (!visited[lid])
                {
                    zero_dirichlet_vel_idx_.push_back(lid);
                    visited[lid] = true;
                }
            }
        }
    };

    for (label izone = 0; izone < mesh.nZones(); izone++)
    {
        const zone* zone = mesh.zonePtr(izone);

        for (label ibdr = 0; ibdr < zone->nBoundaries(); ibdr++)
        {
            const boundary* boundary = zone->boundaryPtr(ibdr);
            const boundaryPhysicalType physicalBCType = boundary->type();
            const boundaryConditionType velocityBCType =
                U.boundaryConditionRef(izone, ibdr).type();

            const stk::mesh::Selector bdr_sel =
                mesh.metaDataRef().locally_owned_part() &
                stk::mesh::selectUnion(boundary->parts());

            const stk::mesh::BucketVector& bdr_buckets =
                mesh.bulkDataRef().get_buckets(stk::topology::NODE_RANK,
                                               bdr_sel);

            switch (physicalBCType)
            {
                case boundaryPhysicalType::wall:
                    filter_local_id(bdr_buckets);
                    break;

                case boundaryPhysicalType::inlet:
                    {
                        switch (velocityBCType)
                        {
                            case boundaryConditionType::normalSpeed:
                            case boundaryConditionType::massFlowRate:
                            case boundaryConditionType::specifiedValue:
                                filter_local_id(bdr_buckets);
                                break;

                            default:
                                break;
                        }
                    }
                    break;

                default:
                    break;
            }
        }
    }

    UTOPIA_TRACE_REGION_END("ContextUtopia::filterBoundaryNodeID");
}

} /* namespace fsi */
} /* namespace accel */
