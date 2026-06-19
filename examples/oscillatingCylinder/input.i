# vim: ft=yaml
# This is a 2D case and must be run with a 2D-compiled binary.
mesh:
    file_path: mesh.e
    automatic_decomposition_type: rcb
simulation:
    verbose: 1
    physical_analysis:
        analysis_type:
            option: transient
            total_time: 49.5072
            time_steps:
                option: constant
                timestep: 0.006876
        domains:
        - name: fluid
          location: [fluid]
          materials: [water]
          type: fluid
          domain_models:
            reference_pressure: 101325
            mesh_deformation:
                option: regions_of_motion_specified
                mesh_motion_model:
                    option: displacement_diffusion
                    mesh_stiffness:
                        option: blended_distance_and_small_volumes
                        distance_exponent: 0.5
                        volume_exponent: 2.0
                displacement_relative_to: initial_mesh
          fluid_models:
            turbulence:
                option: laminar
          boundaries:
          - name: cylinder
            type: wall
            location: [CYLINDER]
            boundary_details:
                mesh_motion:
                    option: periodic_displacement
                    displacement:
                        frequency: 0.20200
                        value: [-0.00795775, 0]
          - name: farfield
            type: wall
            location: [FARFIELD]
          initialization:
            velocity:
                option: value
                velocity: [0, 0]
            pressure:
                option: value
                pressure: 0
    solver:
        solver_control:
            basic_settings:
                advection_scheme: high_resolution
                transient_scheme: second_order_backward_euler
                convergence_controls:
                    min_iterations: 1
                    max_iterations: 20
                    relaxation_parameters:
                        velocity_relaxation_factor: 0.6
                        pressure_relaxation_factor: 0.4
                convergence_criteria:
                    residual_type: RMS
                    residual_target: 1.0e-6
            advanced_options:
                pressure_level_information:
                    option: cartesian_coordinates
                    cartesian_coordinates: [0.55, 0]
                    relative_pressure_level: 0
                linear_solver_settings:
                    default:
                        family: Trilinos
                        min_iterations: 3
                        max_iterations: 50
                        rtol: 1.0e-3
                        atol: 1.0e-12
                        options:
                            belos_solver: gmres
                            preconditioner: riluk
                            preconditioner_parameters:
                                "fact: iluk level-of-fill": 2
                                "fact: drop tolerance": 1.0e-3
                                "fact: absolute threshold": 1.0e-6
                                "fact: relative threshold": 1.0
        output_control:
            file_path: results.e
            output_frequency:
                option: timestep_interval
                timestep_interval: 10
            output_fields: [velocity, pressure, displacement_mesh, velocity_gradient]
            post_process:
            - name: cylinder_forces
              type: force
              options:
                calculate_moment: true
                moment_center: [0, 0]
              location: [CYLINDER]
              frequency: 1
              write_to_file: true
    material_library:
    - name: water
      thermodynamic_properties:
        equation_of_state:
            option: value
            density: 998.2
      transport_properties:
        dynamic_viscosity:
            option: value
            dynamic_viscosity: 1.00818e-3
