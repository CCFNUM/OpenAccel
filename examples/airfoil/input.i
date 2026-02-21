mesh:
    file_path: mesh.e
simulation:
    verbose: 1
    physical_analysis:
        analysis_type:
            option: steady_state
        domains:
          - name: default_domain
            location: [in-hex,out-hex]
            materials: [air]
            type: fluid
            domain_models:
                reference_pressure: 101325
            fluid_models:
                turbulence:
                    option: laminar
            boundaries:
                - name: airfoil
                  type: wall
                  location: [airfoil]
                - name: inlet
                  type: inlet
                  location: [left]
                  boundary_details:
                      mass_and_momentum:
                          option: velocity_components
                          u: 50
                          v: 0
                          w: 0
                - name: outlet
                  type: outlet
                  location: [right]
                  boundary_details:
                      mass_and_momentum:
                          option: static_pressure
                          relative_pressure: 0
                - name: symmetry
                  type: symmetry
                  location: [front_back_in,front_back_out,top,bot]
            initialization:
                velocity:
                    option: value
                    velocity: [50,0,0]
                pressure:
                    option: value
                    pressure: 0
        interfaces:
          - name: interface1
            option: general_connection
            search_tolerance: 0.01
            gauss_lobatto_quadrature: true
            type: fluid_fluid
            side1:
                domain: default_domain
                region_list: [iface_s1]
            side2:
                domain: default_domain
                region_list: [iface_s2]                    
    solver:
        solver_control:
            basic_settings:
                advection_scheme: high_resolution
                convergence_controls:
                    min_iterations: 1
                    max_iterations: 1000
                    physical_timescale: 0.001
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
                        family: Trilinos
                        min_iterations: 3
                        max_iterations: 20
                        rtol: 1.0e-2
                        atol: 1.0e-12
                        options:
                            belos_solver: gmres
                            preconditioner: ilu  
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
                                coarsen_type: 10      # HMIS (Excellent parallel scaling)
                                interp_type: 6        # Extended+i (Robust for stretched grids)
                                relax_type: 18        # L1-Gauss-Seidel (Stable/Smooth)
                                strong_threshold: 0.25
                                num_sweeps: 1
                                max_levels: 20
                                aggressive_levels: 1  # Reduces memory overhead
                                trunc_factor: 0.3     # Keeps the solver lean                              
        output_control:
            file_path: results.e
            output_frequency: 50
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
                dynamic_viscosity: 1e-24 #inviscid
