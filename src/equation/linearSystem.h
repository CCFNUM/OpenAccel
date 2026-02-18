// File : linearSystem.h
// Created : Fri Jan 26 2024 09:35:23 (+0100)
// Author : Fabian Wermelinger
// Description: Abstract base for a linear system of equation(s)
// Copyright (c) 2024 CCFNUM, Lucerne University of Applied Sciences and Arts.
// SPDX-License-Identifier: BSD-3-Clause

#ifndef LINEARSYSTEM_H
#define LINEARSYSTEM_H

// code
#include "assembler.h"
#include "git_revision.h"
#include "linearSolver.h"
#include "macros.h"
#include "mesh.h"
#include "messager.h"
#include "residual.h"
#include "simulation.h"
#include "types.h"

// supported linear solvers
#include "HYPRESolver.h"
#include "PETScSolver.h"
#include "TrilinosSolver.h"

namespace accel
{

template <size_t N>
class linearSystem
{
public:
    static constexpr size_t BLOCKSIZE = N;
    using LinearSolver = ::linearSolver::Base<BLOCKSIZE>;
    using Context = typename LinearSolver::Context;

    linearSystem() = delete;
    linearSystem(const linearSystem& c) = delete;
    linearSystem& operator=(const linearSystem& c) = delete;
    linearSystem(linearSystem&& c) = delete;
    linearSystem& operator=(linearSystem&& c) = delete;

    linearSystem(simulation& sim)
        : max_iterations_(std::numeric_limits<int>::max()), min_iterations_(1),
          convergence_tolerance_(1.0e-9), lsolver_(nullptr), sim_(sim),
          system_name_("undefined"), is_converged_(false), write_system_(false),
          normalize_(false), diagonalScale_(false)
    {
        // default equation names
        for (int i = 0; i < LinearSolver::BLOCKSIZE; i++)
        {
            equation_name_.push_back("Eq. " + std::to_string(i + 1));
        }
    }

    virtual ~linearSystem()
    {
        if (lsolver_)
        {
            delete lsolver_;
            lsolver_ = nullptr;
        }
    }

    bool isConverged() const
    {
        return is_converged_;
    }

    bool hasSolver() const
    {
        return !(lsolver_ == nullptr);
    }

    ::linearSolver::GraphLayout setupSolver(const std::string system_name);
    std::shared_ptr<Context> setupSolver(const std::string system_name,
                                         mesh& mesh);
    void setEquationName(const std::vector<std::string>& name);

    simulation& simulationRef()
    {
        return sim_;
    }

    const simulation& simulationRef() const
    {
        return sim_;
    }

    std::shared_ptr<Context> getContext()
    {
        assert(this->lsolver_);
        return lsolver_->getContext();
    }

    std::shared_ptr<const Context> getContext() const
    {
        assert(this->lsolver_);
        return lsolver_->getContext();
    }

    void solve()
    {
        solveSystem_();
    }

protected:
    int max_iterations_;           // outer iterations
    int min_iterations_;           // outer iterations
    double convergence_tolerance_; // final tolerance for control residual
    residualType residual_type_ = residualType::RMS;
    std::vector<std::string> equation_name_; // descriptive names for equations

    LinearSolver* lsolver_;

    virtual void setResidualScales_() = 0;
    virtual void convergenceReport_();
    void solveSystem_();

    // residual I/O
    static constexpr char COMMENT[] = "# ";
    std::shared_ptr<std::ostream> residual_stream_;
    std::string residual_file_name_;
    void initializeHistory_();
    void writeHistory_() const;

private:
    simulation& sim_;
    std::string system_name_;
    bool is_converged_;

    // DEBUG: write system coefficients
    bool write_system_;
    std::string write_system_name_;

    // additional instructions
    bool normalize_;
    bool diagonalScale_;

    std::string canonicalSystemName_() const
    {
        std::string name = this->system_name_;
        std::replace(name.begin(), name.end(), ' ', '_');
        return name;
    }

    void writeLinearSystem_() const;
    void writeVector_(const std::string fname,
                      const typename LinearSolver::Vector& v) const;

    // residual field for debugging
    stk::mesh::Field<scalar>* stkResidualFieldPtr_{nullptr};
    mesh* mesh_{nullptr};

    void writeSTKInitialResiduals_();
};

template <size_t N>
::linearSolver::GraphLayout
linearSystem<N>::setupSolver(const std::string system_name)
{
    const YAML::Node sim_root = sim_.getYAMLSimulationNode();

    // set the name of the system solver is being setup for
    this->system_name_ = system_name;

    if (!sim_root["solver"])
    {
        errorMsg("linearSystem: `solver` node not present in "
                 "`simulation` node (YAML)");
    }
    const YAML::Node& solver_lookup = sim_root["solver"];

    // setup equation system specific config (outer iterations and target
    // tolerance of control residual)
    if (!solver_lookup["solver_control"]["basic_settings"]
                      ["convergence_controls"])
    {
        // NOTE: This error is enforced
        // because changes in yaml file structure may cause the following
        // options to be silently ignored and thus unset. This may cause
        // unintended configuration of simulation runs.
        errorMsg("linearSystem: `convergence_controls` node not present in "
                 "`solver_control` node (YAML)");
    }
    const YAML::Node& convergenceControls =
        solver_lookup["solver_control"]["basic_settings"]
                     ["convergence_controls"];
    if (convergenceControls["max_iterations"])
    {
        max_iterations_ =
            convergenceControls["max_iterations"].template as<int>();
    }
    if (convergenceControls["min_iterations"])
    {
        min_iterations_ =
            convergenceControls["min_iterations"].template as<int>();
    }

    const YAML::Node& convergenceCriteria =
        solver_lookup["solver_control"]["basic_settings"]
                     ["convergence_criteria"];
    if (convergenceCriteria["residual_target"])
    {
        convergence_tolerance_ =
            convergenceCriteria["residual_target"].template as<double>();

        // not used yet BODGEEEEE
        residual_type_ = convertResidualTypeFromString(
            convergenceCriteria["residual_type"].template as<std::string>());
    }

    // setup linear solver (required)
    std::string solver_name;
    std::string solver_type;
    YAML::Node solver;

    // check for linear solver settings in solver control
    std::string equationName = canonicalSystemName_();

    // Truncate the equationName to before the padded dash
    {
        // Find the position of " - "
        size_t pos = equationName.find(" - ");

        if (pos != std::string::npos)
        {
            equationName = equationName.substr(0, pos);
        }
    }

    std::replace(equationName.begin(), equationName.end(), '-', '_');
    ::accel::tolower(equationName);

    const auto& solverControl = sim_.getYAMLSolverControlNode();

    auto setSolver = [&solver_name, &solverControl, &solver_lookup, &solver](
                         const std::string& key)
    {
        if (solverControl["advanced_options"]["linear_solver_settings"][key]
                         ["lookup"])
        {
            solver_name = solverControl["advanced_options"]
                                       ["linear_solver_settings"][key]["lookup"]
                                           .template as<std::string>(); // case
            // sensitive
            if (!solver_lookup[solver_name])
            {
                errorMsg("linearSystem: `" + solver_name +
                         "` not a valid linear solver request for equation '" +
                         key + "'");
            }
            solver.reset(YAML::Clone(solver_lookup[solver_name]));
        }
        else
        {
            solver.reset(
                YAML::Clone(solverControl["advanced_options"]
                                         ["linear_solver_settings"][key]));
        }
    };

    if (solverControl["advanced_options"] &&
        solverControl["advanced_options"]["linear_solver_settings"] &&
        solverControl["advanced_options"]["linear_solver_settings"]
                     [equationName])
    {
        setSolver(equationName);

        // copy equation name to solver name (for print out purposes)
        solver_name = equationName;
        solver_type = equationName;
    }
    else if (solverControl["advanced_options"] &&
             solverControl["advanced_options"]["linear_solver_settings"] &&
             solverControl["advanced_options"]["linear_solver_settings"]
                          ["default"])
    {
        setSolver("default");

        // copy equation name to solver name (for print out purposes)
        solver_name = equationName;
        solver_type = "default";
    }
    else
    {
        errorMsg("linearSystem: no `solver` defined for equation '" +
                 equationName + "'");
    }

    if (solver["write_system"])
    {
        write_system_ = solver["write_system"].template as<bool>();
        write_system_name_ = canonicalSystemName_();
    }

    if (solver["normalize_matrix"])
    {
        normalize_ = solver["normalize_matrix"].template as<bool>();
    }

    if (solver["diagonal_scaling"])
    {
        diagonalScale_ = solver["diagonal_scaling"].template as<bool>();
    }

    std::string family_type;
    if (solver["family"])
    {
        family_type =
            solver["family"].template as<std::string>(); // case insensitive
    }
    else
    {
        errorMsg("linearSystem: `family` node not present in `" + solver_name +
                 "` solver node (YAML)");
    }
    ::accel::tolower(family_type);

    ::linearSolver::GraphLayout layout = ::linearSolver::GraphLayout::UNDEFINED;

    if (family_type == "petsc")
    {
#ifdef HAS_PETSC
        lsolver_ = new ::linearSolver::PETSc<N>(solver);
        layout = ::linearSolver::GraphLayout::ColumnIndexOrder__Global;
#else
        errorMsg("linearSystem: executable does not support PETSc (" +
                 this->system_name_ + ")");
#endif /* HAS_PETSC */
    }
    else if (family_type == "hypre")
    {
#ifdef HAS_HYPRE
        lsolver_ = new ::linearSolver::HYPRE<N>(solver);
        layout = ::linearSolver::GraphLayout::ColumnIndexOrder__Global;
#else
        errorMsg("linearSystem: executable does not support HYPRE (" +
                 this->system_name_ + ")");
#endif /* HAS_HYPRE */
    }
    else if (family_type == "trilinos")
    {
#ifdef HAS_TRILINOS
        lsolver_ = new ::linearSolver::Trilinos<N>(solver);
        layout = ::linearSolver::GraphLayout::ColumnIndexOrder__Global;
#else
        errorMsg("linearSystem: executable does not support Trilinos (" +
                 this->system_name_ + ")");
#endif /* HAS_TRILINOS */
    }
    else
    {
        errorMsg("linearSystem: linear solver family `" + family_type +
                 "` is not supported");
    }

    return layout;
}

template <size_t N>
std::shared_ptr<typename linearSystem<N>::Context>
linearSystem<N>::setupSolver(const std::string system_name, mesh& mesh)
{
    // instantiate solver and set system_name
    std::shared_ptr<typename linearSystem<N>::Context> ctx;
    switch (this->setupSolver(system_name) &
            (::linearSolver::GraphLayout::ColumnIndexOrder__Local |
             ::linearSolver::GraphLayout::ColumnIndexOrder__Global))
    {
        case ::linearSolver::GraphLayout::ColumnIndexOrder__Local:
        case ::linearSolver::GraphLayout::ColumnIndexOrder__Local |
            ::linearSolver::GraphLayout::ColumnIndexOrder__Global:
            ctx = lsolver_->createContext(this->system_name_,
                                          mesh.getLocalOrderGraphPtr());
            break;
        case ::linearSolver::GraphLayout::ColumnIndexOrder__Global:
            ctx = lsolver_->createContext(this->system_name_,
                                          mesh.getGlobalOrderGraphPtr());
            break;
    }

    // create STK residual field
    stkResidualFieldPtr_ = &(mesh.metaDataPtr()->template declare_field<scalar>(
        stk::topology::NODE_RANK, canonicalSystemName_() + "_residual", 1));
    mesh_ = &mesh;
    assert(mesh_);
    for (auto zone : mesh_->zoneVector())
    {
        stk::mesh::put_field_on_mesh(
            *stkResidualFieldPtr_,
            stk::mesh::selectUnion(zone->interiorParts()),
            N,
            nullptr);
    }
    stk::io::create_named_suffix_field_output_type("solverResiduals",
                                                   equation_name_);
    stk::io::set_named_suffix_field_output_type(*stkResidualFieldPtr_,
                                                "solverResiduals");
    return ctx;
}

template <size_t N>
void linearSystem<N>::setEquationName(const std::vector<std::string>& name)
{
    if (name.size() != LinearSolver::BLOCKSIZE)
    {
        errorMsg("linearSystem: equation name vector size does not "
                 "match number of unknowns!");
    }
    equation_name_ = name;
}

template <size_t N>
void linearSystem<N>::solveSystem_()
{
    if (lsolver_)
    {
        if (normalize_)
        {
            lsolver_->getContext()->getCoefficients().normalize();
        }

        if (diagonalScale_)
        {
            lsolver_->getContext()->getCoefficients().diagonalScale();
        }

        if (lsolver_->getContext()->getCallCount() == 0)
        {
            initializeHistory_();
        }

        const auto& outputCtrl = sim_.controlsRef().solverRef().outputControl_;
        const auto& outFreq = outputCtrl.outputFrequency_;
        const bool write_now =
            (sim_.controlsRef().isTransient() &&
             outFreq.option_ == outputFrequencyType::timeInterval)
                ? sim_.writeNowTime(outFreq.timeInterval_)
                : sim_.writeNow(outFreq.timestepInterval_);
        if (write_system_ && write_now)
        {
            this->writeLinearSystem_();
        }
        if (write_now)
        {
            this->writeSTKInitialResiduals_();
        }

        sim_.getProfiler().push("linear_system_solve");
        lsolver_->solve();
        sim_.getProfiler().pop();

        this->convergenceReport_();
    }
    else
    {
        errorMsg("linearSystem: no solver available to solve system!");
    }
}

template <size_t N>
void linearSystem<N>::convergenceReport_()
{
    using Array = typename LinearSolver::Array;

    setResidualScales_();

    assert(lsolver_);
    auto ctx = lsolver_->getContext();
    typename LinearSolver::ControlData cdata = ctx->getControlData();
    const Array residual_scales = ctx->getResidualScales();
    Array& control_residual = ctx->getControlResidual();
    Array residual = {0};
    Array rms_rate = {0};
    std::vector<std::string> convergence_hint;

    for (int i = 0; i < LinearSolver::BLOCKSIZE; i++)
    {
        residual[i] = cdata.scaled_initial_res[i] * residual_scales[i];
        rms_rate[i] = control_residual[i] / residual[i];
        // clip
        cdata.solver_final_res[i] = std::max(cdata.solver_final_res[i], 0.0);

        const double residual_rate =
            cdata.solver_initial_res[i] / cdata.solver_final_res[i];
        const double overall_rate = rms_rate[i];

        std::string hint = "ok";
        if (this->is_converged_)
        {
            hint = "CONV";
        }
        else if (std::abs(1.0 - overall_rate) < 1.0e-6 ||
                 cdata.n_iterations == 0)
        {
            hint = "STAG";
        }
        else if (overall_rate < 1.0 && residual_rate < 1.0)
        {
            hint = "FAIL";
        }
        else if (overall_rate > 1.2 && overall_rate < 2.0 &&
                 residual_rate >= 5.0)
        {
            hint = "good";
        }
        else if (overall_rate >= 2.0 && overall_rate < 10.0 &&
                 residual_rate >= 10.0)
        {
            hint = "GOOD";
        }
        else if (overall_rate >= 10.0 && residual_rate >= 10.0)
        {
            hint = "WOW!";
        }
        convergence_hint.push_back(hint);
    }

    auto& cout = ctx->cout();
    // clang-format off
    cout << '\n';
    cout << "Physics:        " << ctx->getSystemName() << '\n';
    cout << "Solver context: " << ctx->getFamily() << '\n';
    cout << "Iterations:     " << cdata.n_iterations << "/" << ctx->maxIterations() << '\n';
    cout << std::fixed << std::scientific;
    cout << "+----------------+---------+---------+-------------------+---------+------+\n";
    cout << "|    Equation    |  Rate   |   RMS   | LinSol: Start/End |  Drop   | Conv.|\n";
    cout << "+----------------+---------+---------+-------------------+---------+------+\n";
    for (int i = 0; i < LinearSolver::BLOCKSIZE; i++) {
        const size_t strlen = equation_name_[i].size();
        const double rate = (rms_rate[i] > 9.9e99) ? 9.9e99 : rms_rate[i];
        double drop = cdata.solver_final_res[i] / cdata.solver_initial_res[i];
        drop = (drop > 9.9e99) ? 9.9e99 : drop;
        cout <<  "| " << std::setw(14) << std::left << equation_name_[i].substr(0, strlen > 14 ? 14 : strlen)
             << " | " << std::setw(7)  << std::left << std::setprecision(1) << rate
             << " | " << std::setw(7)  << std::left << std::setprecision(1) << residual[i]
             << " | " << std::left << std::setprecision(2) << cdata.solver_initial_res[i] << "/" << std::setprecision(2) << cdata.solver_final_res[i]
             << " | " << std::setw(6) << std::left << std::setprecision(1) << drop
             << " | " << std::setw(4)  << std::left << convergence_hint[i]
             << " |\n";
    }
    cout.unsetf(std::ios::floatfield);
    cout << "+----------------+---------+---------+-------------------+---------+------+";
    cout << std::endl; // flush
    // clang-format on

    // additional reporting from linear solver
    lsolver_->report();

    // update residual and check convergence
    control_residual = residual;
    bool required_iterations = sim_.getIterationCount() > min_iterations_;
    const bool iteration_limit = sim_.getIterationCount() > max_iterations_;
    for (int i = 0; i < LinearSolver::BLOCKSIZE; i++)
    {
        required_iterations =
            required_iterations && (residual[i] < convergence_tolerance_);
    }

    is_converged_ = false;
    if (required_iterations || iteration_limit)
    {
        is_converged_ = true;
    }

    // log residuals
    this->writeHistory_();
}

template <size_t N>
void linearSystem<N>::initializeHistory_()
{
    if (lsolver_)
    {
        const auto ctx = lsolver_->getContext();

        residual_file_name_ = this->canonicalSystemName_();
        residual_file_name_ = sim_.getResidualDirectory() / residual_file_name_;
        const size_t len = residual_file_name_.size();
        int nfiles = 0;
        for (const auto& entry :
             fs::directory_iterator{sim_.getResidualDirectory()})
        {
            if (residual_file_name_.compare(0, len, entry.path(), 0, len) == 0)
            {
                ++nfiles;
            }
        }
        std::ostringstream oss;
        oss << residual_file_name_ << "_" << std::setfill('0') << std::setw(3)
            << ++nfiles << ".out";
        residual_file_name_ = oss.str();

        residual_stream_ = ctx->makeFileStream(residual_file_name_);
        assert(residual_stream_ != nullptr);

        // write header
        auto now = std::chrono::system_clock::now();
        auto in_time_t = std::chrono::system_clock::to_time_t(now);

        auto& fout = *residual_stream_;
        fout << COMMENT << "cN-CVFEM solver timestamp: "
             << std::put_time(std::localtime(&in_time_t), "%c\n");
        fout << COMMENT << "Git revision: " << accel::git_revision << '\n';
        fout << COMMENT << '\n';
        ctx->info(fout, COMMENT);
        fout << COMMENT << '\n';
        fout << COMMENT << "System residual history:\n";
        // clang-format off
        fout << COMMENT
             << "global_iterations" /* 0 */ << '\t'
             << "inner_iterations"  /* 1 */ << '\t'
             << "solver_calls"      /* 2 */ << '\t'
             << "solver_iterations" /* 3 */ << '\t'
             << "sim_time[s]"       /* 4 */;
        constexpr unsigned int residual_start = 5; // residual start column
        // clang-format on
        for (const auto& eq_name : equation_name_)
        {
            fout << '\t' << eq_name << "[control_residual]";
        }
        fout << std::endl;

        // append to region plot list
        unsigned int k = 0;
        for (const auto& eq_name : equation_name_)
        {
            simulation::residualPlotItem item;
            item.data_file = residual_file_name_;
            item.legend_name = eq_name;
            item.xdata_idx = 0;                    // global iteration
            item.ydata_idx = residual_start + k++; // residual component(s)
            sim_.addPlotItem(item);
        }
    }
}

template <size_t N>
void linearSystem<N>::writeHistory_() const
{
    if (lsolver_)
    {
        assert(residual_stream_ != nullptr);

        const auto ctx = lsolver_->getContext();
        const typename LinearSolver::Array control_residual =
            ctx->getControlResidual();
        const typename LinearSolver::ControlData ctrl = ctx->getControlData();

        auto& fout = *residual_stream_;
        fout << sim_.getGlobalIterationCount();
        fout << '\t' << sim_.getIterationCount();
        fout << '\t' << ctx->getCallCount();
        fout << '\t' << ctrl.n_iterations;
        fout << '\t' << std::setprecision(3) << std::scientific
             << sim_.getSimulationTime();
        for (int i = 0; i < LinearSolver::BLOCKSIZE; i++)
        {
            fout << '\t' << std::setprecision(6) << std::scientific
                 << control_residual[i];
        }
        fout << std::endl; // flush
    }
}

template <size_t N>
void linearSystem<N>::writeLinearSystem_() const
{
    assert(lsolver_);
    const auto ctx = lsolver_->getContext();
    assert(ctx);

    std::ostringstream oss;
    oss << write_system_name_ << "_" << std::setfill('0') << std::setw(4)
        << ctx->getCallCount();
    ctx->getAMatrix().writeMatrix(oss.str());
    writeVector_(oss.str() + "_b.bin", ctx->getBVector());
}

template <size_t N>
void linearSystem<N>::writeVector_(const std::string fname,
                                   const typename LinearSolver::Vector& v) const
{
    using DataType = typename LinearSolver::Matrix::DataType;

    assert(lsolver_);
    const auto ctx = lsolver_->getContext();

    const int size = ctx->commSize();
    const int rank = ctx->commRank();
    const MPI_Comm comm = ctx->getCommunicator();

    const int local_n = static_cast<int>(v.size());
    std::vector<int> local_sizes(size, 0);
    MPI_Gather(&local_n, 1, MPI_INT, local_sizes.data(), 1, MPI_INT, 0, comm);

    std::vector<int> displ;
    int offset = 0;
    for (const int k : local_sizes)
    {
        displ.push_back(offset);
        offset += k;
    }

    typename LinearSolver::Vector global_v(offset);
    MPI_Gatherv(v.data(),
                local_n,
                linearSolver::MPIDataType<DataType>::type(),
                global_v.data(),
                local_sizes.data(),
                displ.data(),
                linearSolver::MPIDataType<DataType>::type(),
                0,
                comm);

    if (0 == rank)
    {
        std::ofstream out(fname, std::ios::binary);
        out.write((const char*)global_v.data(),
                  global_v.size() * sizeof(DataType));
        out.close();
    }
}

template <size_t N>
void linearSystem<N>::writeSTKInitialResiduals_()
{
    if (stkResidualFieldPtr_)
    {
        assert(mesh_);
        stk::mesh::BulkData& bulk = mesh_->bulkDataRef();

        // compute residual
        auto ctx = lsolver_->getContext();
        auto& coeff = ctx->getCoefficients();
        linearSolver::residual::compute(coeff);

        const auto& r = ctx->getRVector();
        for (auto zone : mesh_->zoneVector())
        {
            for (const auto bucket : bulk.get_buckets(
                     stk::topology::NODE_RANK,
                     mesh_->metaDataRef().locally_owned_part() &
                         stk::mesh::selectUnion(zone->interiorParts())))
            {
                scalar* stk_r =
                    stk::mesh::field_data(*stkResidualFieldPtr_, *bucket);
                for (unsigned i = 0; i < bucket->size(); i++)
                {
                    const size_t idx = bulk.local_id(bucket->operator[](i));
                    assert(N * idx < r.size());
                    const scalar* src = &r[N * idx];
                    for (int j = 0; j < N; j++)
                    {
                        stk_r[N * i + j] = src[j];
                    }
                }
            }
        }
    }
}

} /* namespace accel */

#endif // LINEARSYSTEM_H
