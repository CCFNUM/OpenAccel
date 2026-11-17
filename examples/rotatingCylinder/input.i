mesh:
    file_path: mesh.e
    automatic_decomposition_type: rcb

simulation:
    verbose: 1
    physical_analysis:
        analysis_type:
            option: transient
            total_time: 100
            time_steps:
                option: constant
                timestep: 0.005
        domains:
          - name: stator
            location: [stator]
            materials: [fluid_1]
            type: fluid
            domain_models:
                reference_pressure: 0
            fluid_models:
                turbulence:
                    option: laminar
            boundaries:
              - name: inlet
                type: inlet
                location: [inlet]
                boundary_details:
                    mass_and_momentum:
                        option: velocity_components
                        u: 1.0
                        v: 0
              - name: outlet
                type: outlet
                location: [outlet]
                boundary_details:
                    mass_and_momentum:
                        option: static_pressure
                        relative_pressure: 0
              - name: top
                type: symmetry
                location: [top]
              - name: bottom
                type: symmetry
                location: [bottom]
            initialization:
                velocity:
                    option: value
                    velocity: [1, 0]
                pressure:
                    option: value
                    pressure: 0
          - name: rotor
            location: [rotor]
            materials: [fluid_1]
            type: fluid
            domain_models:
                reference_pressure: 0
                domain_motion:
                    option: rotating
                    origin: [0, 0]
                    axis: [0, 0, 1]
                    angular_velocity: 2.0 # rotation speed: alpha = omega*D/(2*U) = omega/2 (here alpha = 1)
            fluid_models:
                turbulence:
                    option: laminar
            boundaries:
              - name: cylinder
                type: wall
                location: [cylinder]
                frame_type: rotating
            initialization:
                velocity:
                    option: value
                    velocity: [1, 0]
                pressure:
                    option: value
                    pressure: 0
        interfaces:
          - name: AMI
            option: general_connection
            search_tolerance: 0.0005
            gauss_lobatto_quadrature: true
            force_nonconformal_treatment: true
            penalty_factor: 1.0
            type: fluid_fluid
            side1:
                domain: stator
                region_list: [AMI_stator]
            side2:
                domain: rotor
                region_list: [AMI_rotor]

    solver:
        solver_control:
            basic_settings:
                transient_scheme: second_order_backward_euler
                advection_scheme: high_resolution
                convergence_controls:
                    min_iterations: 2
                    max_iterations: 40
                    relaxation_parameters:
                        velocity_relaxation_factor: 0.9
                        pressure_relaxation_factor: 0.5
                convergence_criteria:
                    residual_type: RMS
                    residual_target: 1e-05
            advanced_options:
                equation_controls:
                    mesh_motion:
                        freeze_per_timestep: true
                    sub_iterations:
                        pressure_correction: 1
                linear_solver_settings:
                    default:
                        family: trilinos
                        min_iterations: 1
                        max_iterations: 200
                        rtol: 0.0001
                        atol: 1.0e-12
                        options:
                            belos_solver: gmres
                            preconditioner: relaxation
                            preconditioner_parameters:
                                "relaxation: type": Symmetric Gauss-Seidel
                                "relaxation: sweeps": 1
                    pressure_correction:
                        family: HYPRE
                        min_iterations: 1
                        max_iterations: 200
                        rtol: 0.001
                        atol: 1.0e-12
                        options:
                            solver:
                                type: GMRES
                            precond:
                                type: BoomerAMG
            expert_parameters:
                unlimited_advection: true
                linear_flow_predictor: true

        output_control:
            file_path: results.e
            output_frequency:
                option: timestep_interval
                timestep_interval: 200
            output_fields: [velocity, pressure]
            post_process:
            - name: cylinder_forces
              type: force
              options:
                calculate_moment: true
                moment_center: [0, 0]
              location: [cylinder]
              frequency: 1
              write_to_file: true

    material_library:
      - name: fluid_1
        thermodynamic_properties:
            equation_of_state:
                option: value
                density: 1.0
        transport_properties:
            dynamic_viscosity:
                option: value
                dynamic_viscosity: 0.005
