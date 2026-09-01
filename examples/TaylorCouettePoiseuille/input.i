# vim: ft=yaml
# This is a 3D case and must be run with a 3D-compiled binary.
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
          materials: [debug_gas]
          type: fluid
          domain_models:
            reference_pressure: 0
            domain_motion:
                option: rotating
                origin: [0, 0, 0]
                axis: [0, 0, 1]
                angular_velocity: 13.5
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
                    input_type: expression
                    velocity: ["4.5*y - 18.0*y/(x*x + y*y)", "-4.5*x + 18.0*x/(x*x + y*y)", "-47.40*(x*x + y*y) + 205.16*log(sqrt(x*x + y*y)) + 47.40"]
          - name: outlet
            type: outlet
            location: [outlet]
            boundary_details:
                mass_and_momentum:
                    option: static_pressure
                    relative_pressure: 0
          - name: r_wall
            type: wall
            location: [r_wall]
            frame_type: rotating
          - name: s_wall
            type: wall
            location: [s_wall]
            frame_type: rotating
            boundary_details:
                mass_and_momentum:
                    wall_velocity:
                        option: counter_rotating_wall
          initialization:
            velocity:
                option: value
                velocity: [0, 0, 10]
            pressure:
                option: value
                pressure: 0
        interfaces:
        - name: interface_1
          option: rotational_periodicity
          type: fluid_fluid
          search_tolerance: 0.001
          gauss_lobatto_quadrature: false
          rotation_axis: [0, 0, 1]
          axis_location: [0, 0, 0]
          side1:
            domain: default_domain
            region_list: [per1]
          side2:
            domain: default_domain
            region_list: [per2]
    solver:
        solver_control:
            basic_settings:
                advection_scheme: high_resolution
                convergence_controls:
                    min_iterations: 1
                    max_iterations: 250
                    physical_timescale: 1
                    relaxation_parameters:
                        velocity_relaxation_factor: 0.9
                convergence_criteria:
                    residual_type: RMS
                    residual_target: 1e-8
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
                                coarsentype: 10
                                interptype: 6
                                relaxtype: 18
                                strongthreshold: 0.25
                                numsweeps: 1
                                maxlevels: 20
                                aggnumlevels: 1
                                truncfactor: 0.3
            expert_parameters:
                consistent: true
        output_control:
            file_path: results.e
            output_frequency: 10
            output_fields: [velocity, pressure]
            corrected_boundary_values: true
    material_library:
    - name: debug_gas
      thermodynamic_properties:
        equation_of_state:
            option: value
            density: 5.85
      transport_properties:
        dynamic_viscosity:
            option: value
            dynamic_viscosity: 1.3185
