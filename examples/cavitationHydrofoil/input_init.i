#vim : ft = yaml
#This is a 2D case and must be run with a 2D - compiled binary
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
          type: fluid
          location: [fluid]
          materials: [water, water_vapor]
          domain_models:
            reference_pressure: 0
          fluid_models:
            turbulence:
                option: k_epsilon
            multiphase:
                homogeneous: true
                free_surface_model:
                    option: standard
          boundaries:
          - name: hydrofoil
            type: wall
            location: [hydrofoil]
          - name: slip_walls
            type: symmetry
            location: [slip_walls]
          - name: inlet
            type: inlet
            location: [inlet]
            boundary_details:
                mass_and_momentum:
                    option: velocity_components
                    u: 16.91
                    v: 0
                turbulence:
                    option: intensity_and_length_scale
                    fractional_intensity: 0.03
                    eddy_length_scale: 0.0076
            fluid_values:
                water_vapor:
                    volume_fraction:
                        option: value
                        volume_fraction: 0
                water:
                    volume_fraction:
                        option: value
                        volume_fraction: 1
          - name: outlet
            type: outlet
            location: [outlet]
            boundary_details:
                mass_and_momentum:
                    option: static_pressure
                    relative_pressure: 51957
          initialization:
            velocity:
                option: value
                velocity: [16.91, 0]
            pressure:
                option: value
                pressure: 51957
            turbulent_kinetic_energy:
                option: value
                turbulent_kinetic_energy: 0.386
            turbulent_dissipation_rate:
                option: value
                turbulent_dissipation_rate: 5.185
            fluid_specific_initialization:
                water_vapor:
                    volume_fraction:
                        option: value
                        volume_fraction: 0
                water:
                    volume_fraction:
                        option: value
                        volume_fraction: 1
    solver:
        solver_control:
            basic_settings:
                advection_scheme: high_resolution
                turbulence_numerics: upwind
                convergence_controls:
                    min_iterations: 1
                    max_iterations: 300
                    physical_timescale: 1
                    relaxation_parameters:
                        velocity_relaxation_factor: 0.7
                        pressure_relaxation_factor: 0.3
                convergence_criteria:
                    residual_type: RMS
                    residual_target: 1.0e-5
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
        output_control:
            file_path: results_noCav.e
            output_frequency: 10
            restart_file_name: restart_noCav.e
            restart_frequency: 10
            output_fields: [velocity, pressure, density, volume_fraction.water, volume_fraction.water_vapor, turbulent_kinetic_energy, turbulent_dissipation_rate, turbulent_viscosity]
    material_library:
    - name: water
      thermodynamic_properties:
        equation_of_state:
            option: value
            density: 996.2
      transport_properties:
        dynamic_viscosity:
            option: value
            dynamic_viscosity: 8.899e-4
    - name: water_vapor
      thermodynamic_properties:
        equation_of_state:
            option: value
            density: 0.02308
      transport_properties:
        dynamic_viscosity:
            option: value
            dynamic_viscosity: 9.8e-6
