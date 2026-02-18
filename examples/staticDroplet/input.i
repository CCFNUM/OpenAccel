mesh:
    file_path: mesh.e
simulation:
    verbose: 1
    physical_analysis:
        analysis_type:
            option: transient
            total_time: 0.5
            time_steps:
                option: adaptive
                initial_timestep: 1.0e-4
                timestep_update_frequency: 1
                timestep_adaptation:
                    option: max_courant
                    courant_number: 0.5
                    min_timestep: 1.0e-7
                    max_timestep: 0.1
                    timestep_decrease_factor: 0.8
                    timestep_increase_factor: 1.06
        domains:
          - name: default_domain
            location: [fluid-hex]
            materials: [water, air]
            type: fluid
            domain_models:
                reference_pressure: 0
            fluid_models:
                turbulence:
                    option: laminar
                multiphase:
                    homogeneous: true
                    free_surface_model:
                        option: standard
            fluid_pair_models:
              - pair: [water, air]
                surface_tension:
                    option: continuum_surface_force
                    surface_tension_coefficient: 73.0
            boundaries:
              - name: walls
                type: symmetry
                location: [walls]
              - name: front_and_back
                type: symmetry
                location: [front_and_back]
            initialization:
                velocity:
                    option: value
                    velocity: [0, 0, 0]
                pressure:
                    option: value
                    pressure: 0
                fluid_specific_initialization:
                    water:
                        volume_fraction:
                            option: value
                            input_type: expression
                            volume_fraction: "if ((x-0.5)*(x-0.5) + (y-0.5)*(y-0.5) < 0.04, 1, 0)"
                    air:
                        volume_fraction:
                            option: value
                            input_type: expression
                            volume_fraction: "if ((x-0.5)*(x-0.5) + (y-0.5)*(y-0.5) < 0.04, 0, 1)"
    solver:
        solver_control:
            basic_settings:
                advection_scheme: high_resolution
                transient_scheme: first_order_backward_euler
                reduced_stencil: true
                convergence_controls:
                    min_iterations: 1
                    max_iterations: 10
                convergence_criteria:
                    residual_type: RMS
                    residual_target: 1e-6
                interpolation_scheme:
                    velocity_interpolation_type: linear_linear
            expert_parameters:
                body_force_redistribution: true
            advanced_options:
                pressure_level_information:
                    option: cartesian_coordinates
                    cartesian_coordinates: [0,0,0]
                    relative_pressure_level: 0                
                equation_controls:
                    volume_fraction_smoothing:
                        smooth_volume_fraction: true
                        smoothing_iterations: 3
                        fourier_number: 0.25
                linear_solver_settings:
                    default:
                        family: PETSc
                        max_iterations: 20
                        rtol: 1.0e-1
                        atol: 1.0e-12
                        options:
                            ksp_type: fgmres
                            pc_type: bjacobi                
                    pressure_correction:
                        family: Trilinos
                        max_iterations: 150 # Cap iterations; if it hasn't found it by 150, it's drifting
                        rtol: 1.0e-6
                        atol: 1.0e-10
                        options:
                            belos_solver: gmres
                            preconditioner: riluk
                            preconditioner_parameters:
                                "fact: iluk level-of-fill": 0  # 0 is much faster than 2
                                "fact: absolute threshold": 1.0e-2
                                "fact: relative threshold": 1.0
                            # Belos-specific speed boosts
                            solver_parameters:
                                "Num Blocks": 50            # Size of Krylov subspace
                                "Maximum Restarts": 2       # Prevents infinite loops
                                "Orthogonalization": "ICGS" # Fast orthogonalization
        output_control:
            file_path: results.e
            output_frequency:
                option: time_interval
                time_interval: 0.01
            write_timestep_info: true
            output_fields: [velocity, pressure, volume_fraction.water, curvature.water_air, body_forces, pressure_gradient]
    material_library:
      - name: water
        thermodynamic_properties:
            equation_of_state:
                option: value
                density: 1000
        transport_properties:
            dynamic_viscosity:
                option: value
                dynamic_viscosity: 1e-3
      - name: air
        thermodynamic_properties:
            equation_of_state:
                option: value
                density: 1000
        transport_properties:
            dynamic_viscosity:
                option: value
                dynamic_viscosity: 1e-3
