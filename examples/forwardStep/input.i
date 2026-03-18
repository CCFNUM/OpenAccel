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
            materials: [air_ideal_gas]
            type: fluid
            domain_models:
                reference_pressure: 101325
            fluid_models:
                heat_transfer:     
                    option: total_energy
            boundaries:
              - name: inlet
                type: inlet
                location: [inlet]
                boundary_details:
                    option: supersonic
                    mass_and_momentum:
                        option: velocity_components_and_static_pressure
                        u: 1041.8447718486013
                        v: 0
                        relative_pressure: 0
                    heat_transfer:
                        option: static_temperature
                        static_temperature: 300
              - name: outlet
                type: outlet
                location: [outlet]
                boundary_details:
                    option: supersonic
              - name: obstacle
                type: wall
                location: [obstacle]
              - name: symmetry
                type: symmetry
                location: [bottom,top]
            initialization:
                velocity:
                    option: value
                    velocity: [1041.8447718486013,0]
                pressure:
                    option: value
                    pressure: 0
                temperature:
                    option: value
                    temperature: 300
    solver:
        solver_control:
            basic_settings:
                advection_scheme: upwind
                convergence_controls:
                    min_iterations: 1
                    max_iterations: 1200
                    physical_timescale: 1e-2
                    relaxation_parameters:
                        velocity_relaxation_factor: 0.7
                convergence_criteria:
                    residual_type: RMS
                    residual_target: 1e-8
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
            expert_parameters:  
                consistent: true
        output_control:
            file_path: results.e
            output_frequency: 10
            output_fields: [velocity, pressure, temperature, density, mach_number, total_pressure, total_temperature]
    material_library:
      - name: air_ideal_gas
        thermodynamic_properties:
            equation_of_state:
                option: ideal_gas
                molar_mass: 28.96
            specific_heat_capacity:
                option: value
                specific_heat_capacity: 1004.4                
        transport_properties:
            dynamic_viscosity:
                option: value
                dynamic_viscosity: 1e-16 
            thermal_conductivity:
                option: value
                thermal_conductivity: 2.61e-2                
