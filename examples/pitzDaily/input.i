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
            location: [fluid-hex]
            materials: [air]
            type: fluid
            domain_models:
                reference_pressure: 101325
            fluid_models:
                turbulence:     
                    option: shear-stress-transport
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
                        w: 0
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
                        velocity_relaxation_factor: 0.7
                        pressure_relaxation_factor: 0.3
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
                    pressure_correction:
                        family: HYPRE
                        min_iterations: 3
                        max_iterations: 20
                        rtol: 1.0e-7
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
            output_frequency: 10
            output_fields: [velocity, pressure, turbulent_kinetic_energy, turbulent_eddy_frequency, turbulent_viscosity, total_pressure, minimum_distance_to_wall, wall_scale]
            corrected_boundary_values: true
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
