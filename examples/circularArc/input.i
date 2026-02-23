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
                    option: thermal_energy
            boundaries:
              - name: inlet
                type: inlet
                location: [inlet]
                boundary_details:
                    mass_and_momentum:
                        option: total_pressure
                        relative_pressure: 18871.405
                    flow_direction:
                        option: cartesian_components
                        x: 1
                        y: 0
                        z: 0
                    heat_transfer:
                        option: total_temperature
                        total_temperature: 315.0095
              - name: outlet
                type: outlet
                location: [outlet]
                boundary_details:
                    mass_and_momentum:
                        option: static_pressure
                        relative_pressure: 0
              - name: slipWall
                type: wall
                location: [slipWall,top]
                boundary_details:
                    mass_and_momentum:
                        option: free_slip_wall                      
              - name: symmetry
                type: symmetry
                location: [symmetry]
            initialization:
                velocity:
                    option: value
                    velocity: [173.64,0,0]
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
                    max_iterations: 1000
                    physical_timescale: 1e-2
                    relaxation_parameters:
                        velocity_relaxation_factor: 0.9
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
