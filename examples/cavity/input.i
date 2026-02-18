mesh:
    file_path: mesh.e
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
                    option: laminar
            boundaries:
              - name: top
                type: wall
                location: [top]
                boundary_details:
                    mass_and_momentum:
                        wall_velocity:
                            option: cartesian_components
                            wall_velocity: [1,0,0]
              - name: sides
                type: wall
                location: [sides]
              - name: front_and_back
                type: symmetry
                location: [frontandback]
            initialization:
                velocity:
                    option: value
                    velocity: [0,0,0]
                pressure:
                    option: value
                    pressure: 0
    solver:
        solver_control:
            basic_settings:
                advection_scheme: upwind
                reduced_stencil: true
                convergence_controls:
                    min_iterations: 1
                    max_iterations: 100
                    physical_timescale: 1
                    relaxation_parameters:
                        velocity_relaxation_factor: 0.8
                        pressure_relaxation_factor: 0.2
                convergence_criteria:
                    residual_type: RMS
                    residual_target: 1e-10
                interpolation_scheme:
                    velocity_interpolation_type: linear_linear
            advanced_options:
                pressure_level_information:
                    option: cartesian_coordinates
                    cartesian_coordinates: [0,0,0]
                    relative_pressure_level: 0
                linear_solver_settings:
                    default:
                        family: Trilinos
                        max_iterations: 200
                        rtol: 1.0e-6
                        atol: 1.0e-12
                        options:
                            belos_solver: fgmres
                            preconditioner: riluk
                    pressure_correction:
                        family: Trilinos
                        max_iterations: 200
                        rtol: 1.0e-6
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
            output_frequency: 10
            output_fields: [velocity, pressure]
    material_library:
      - name: fluid_1
        thermodynamic_properties:
            equation_of_state:
                option: value
                density: 1
        transport_properties:
            dynamic_viscosity:
                option: value
                dynamic_viscosity: 0.01
