mesh:
    file_path: mesh.e
simulation:
    verbose: 1
    physical_analysis:
        analysis_type:
            option: steady_state
        domains:
          - name: default_domain
            location: [unspecified-2-wedge]
            materials: [air]
            type: fluid
            domain_models:
                reference_pressure: 101325
            fluid_models:
                turbulence:
                    option: laminar
            boundaries:
              - name: wall
                type: wall
                location: [wall]
              - name: inlet1
                type: inlet
                location: [inlet1]
                boundary_details:
                    mass_and_momentum:
                        option: velocity_components
                        u: 1
                        v: 0
                        w: 0
              - name: inlet2
                type: inlet
                location: [inlet2]
                boundary_details:
                    mass_and_momentum:
                        option: velocity_components
                        u: 0
                        v: 3
                        w: 0
              - name: outlet
                type: outlet
                location: [outlet]
                boundary_details:
                    mass_and_momentum:
                        option: static_pressure
                        relative_pressure: 0
              - name: front_and_back
                type: symmetry
                location: [front,back]
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
                    max_iterations: 250
                    physical_timescale: 1   
                    relaxation_parameters:
                        velocity_relaxation_factor: 0.8
                        pressure_relaxation_factor: 0.2              
                convergence_criteria:
                    residual_type: RMS
                    residual_target: 1e-6
                interpolation_scheme:
                    velocity_interpolation_type: linear_linear
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
        output_control:
            file_path: results.e
            output_frequency: 10
            output_fields: [velocity, pressure]
    material_library:
      - name: air
        thermodynamic_properties:
            equation_of_state:
                option: value
                density: 1.185
        transport_properties:
            dynamic_viscosity:
                option: value
                dynamic_viscosity: 1.831e-5
