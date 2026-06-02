# vim: ft=yaml
# This is a 2D case and must be run with a 2D-compiled binary.
mesh:
    file_path: mesh.e
    automatic_decomposition_type: rcb
simulation:
    verbose: 1
    physical_analysis:
        analysis_type:
            option: steady_state
        domains:
        - name: default_domain
          location: [fluid]
          materials: [fluid_1]
          type: fluid
          domain_models:
            reference_pressure: 101325
          fluid_models:
            turbulence:
                option: shear_stress_transport
          boundaries:
          - name: wall
            type: wall
            location: [wall]
          - name: inlet
            type: inlet
            location: [inlet]
            boundary_details:
                mass_and_momentum:
                    option: velocity_components
                    u: 10
                    v: 0
                turbulence:
                    option: k_and_omega
                    k: 0.375
                    omega: 440.15
          - name: outlet
            type: outlet
            location: [outlet]
            boundary_details:
                mass_and_momentum:
                    option: static_pressure
                    relative_pressure: 0
          initialization:
            velocity:
                option: value
                velocity: [0, 0]
            pressure:
                option: value
                pressure: 0
            turbulent_kinetic_energy:
                option: value
                turbulent_kinetic_energy: 0.375
            turbulent_eddy_frequency:
                option: value
                turbulent_eddy_frequency: 440.15
    solver:
        solver_control:
            basic_settings:
                advection_scheme: high_resolution
                turbulence_numerics: upwind
                convergence_controls:
                    min_iterations: 1
                    max_iterations: 2500
                    physical_timescale: 1
                    relaxation_parameters:
                        velocity_relaxation_factor: 0.5
                        turbulence_relaxation_factor: 0.5
                convergence_criteria:
                    residual_type: RMS
                    residual_target: 1e-6
            advanced_options:
                linear_solver_settings:
                    default:
                        family: PETSc
                        min_iterations: 3
                        max_iterations: 20
                        rtol: 1.0e-1
                        atol: 1.0e-12
                        options:
                            ksp_type: fgmres
                            pc_type: bjacobi
                    pressure_correction:
                        family: HYPRE
                        min_iterations: 3
                        max_iterations: 20
                        rtol: 1.0e-3
                        atol: 1.0e-12
                        options:
                            solver:
                                type: GMRES
                            precond:
                                type: BoomerAMG
                                coarsen_type: 10
                                interp_type: 6
                                relax_type: 18
                                strong_threshold: 0.25
                                num_sweeps: 1
                                max_levels: 20
                                aggressive_levels: 1
                                trunc_factor: 0.3
            expert_parameters:
                relax_gradients: false
                incremental_gradient_change: false
                consistent: true
                wall_distance_method: mesh_wave
        output_control:
            file_path: results.e
            output_frequency: 50
            output_fields: [velocity, pressure, turbulent_kinetic_energy, turbulent_eddy_frequency, minimum_distance_to_wall, turbulent_viscosity]
    material_library:
      - name: fluid_1
        thermodynamic_properties:
            equation_of_state:
                option: value
                density: 1.0
        transport_properties:
            dynamic_viscosity:
                option: value
                dynamic_viscosity: 1e-05
