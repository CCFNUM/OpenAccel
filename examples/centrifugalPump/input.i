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
        - name: rotor_domain
          location: [rotor2D-hex]
          materials: [water]
          type: fluid
          domain_models:
            reference_pressure: 0
            domain_motion:
                option: rotating
                origin: [0, 0]
                axis: [0, 0, 1]
                angular_velocity: -209
          fluid_models:
            turbulence:
                option: k_epsilon
          boundaries:
          - name: inlet
            type: inlet
            location: [INLET]
            frame_type: stationary
            boundary_details:
                mass_and_momentum:
                    option: normal_speed
                    normal_speed: 11.4
                turbulence:
                    option: k_and_epsilon
                    k: 0.48735
                    epsilon: 138.342
          - name: blade_rot
            type: wall
            location: [BLADE_ROT]
            frame_type: rotating
          initialization:
            velocity:
                option: value
                velocity: [0, 0]
            pressure:
                option: value
                pressure: 0
            turbulent_kinetic_energy:
                option: value
                turbulent_kinetic_energy: 0.48735
            turbulent_dissipation_rate:
                option: value
                turbulent_dissipation_rate: 138.342
        - name: stator_domain
          location: [stator2D-hex]
          materials: [water]
          type: fluid
          domain_models:
            reference_pressure: 0
          fluid_models:
            turbulence:
                option: k_epsilon
          boundaries:
          - name: outlet
            type: outlet
            location: [OUTLET]
            boundary_details:
                mass_and_momentum:
                    option: average_static_pressure
                    relative_pressure: 0
                    pressure_profile_blend: 0.05
          - name: blade_stat
            type: wall
            location: [BLADE_STAT]
          initialization:
            velocity:
                option: value
                velocity: [0, 0]
            pressure:
                option: value
                pressure: 0
            turbulent_kinetic_energy:
                option: value
                turbulent_kinetic_energy: 0.48735
            turbulent_dissipation_rate:
                option: value
                turbulent_dissipation_rate: 138.342
        interfaces:
        - name: ggi_rotor_stator
          option: general_connection
          search_tolerance: 0.005
          type: fluid_fluid
          side1:
            domain: rotor_domain
            region_list: [GGI_INT]
          side2:
            domain: stator_domain
            region_list: [GGI_EXT]
    solver:
        solver_control:
            basic_settings:
                advection_scheme: high_resolution
                turbulence_numerics: upwind
                convergence_controls:
                    min_iterations: 1
                    max_iterations: 3000
                    physical_timescale: 1e-3
                    relaxation_parameters:
                        relax_mass: 0.3
                        velocity_relaxation_factor: 0.3
                        pressure_relaxation_factor: 0.3
                        turbulence_relaxation_factor: 0.3
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
        output_control:
            file_path: results.e
            output_frequency: 10
            output_fields: [velocity, relative_velocity, pressure, turbulent_kinetic_energy, turbulent_dissipation_rate, turbulent_viscosity, total_pressure, minimum_distance_to_wall]
    material_library:
    - name: water
      thermodynamic_properties:
        equation_of_state:
            option: value
            density: 1000
      transport_properties:
        dynamic_viscosity:
            option: value
            dynamic_viscosity: 0.01
