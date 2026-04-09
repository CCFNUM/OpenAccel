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
                    option: shear_stress_transport
                    transitional_turbulence: true
            boundaries:
              - name: inlet
                type: inlet
                location: [inlet]
                boundary_details:
                    mass_and_momentum:
                        option: velocity_components
                        u: 5.4
                        v: 0
                        w: 0
                    turbulence:
                        option: intensity_and_eddy_viscosity_ratio
                        fractional_intensity: 0.033
                        eddy_viscosity_ratio: 12
              - name: outlet
                type: outlet
                location: [outlet]
                boundary_details:
                    mass_and_momentum:
                        option: static_pressure 
                        relative_pressure: 0.0
              - name: plate
                type: wall
                location: [plate]
              - name: slip
                type: wall
                location: [above, top]
                boundary_details:
                    mass_and_momentum:
                        option: free_slip_wall 
              - name: symmetry
                type: symmetry
                location: [defaultfaces]
            initialization:
                velocity:
                    option: value
                    velocity: [5.4,0,0]
                pressure:
                    option: value
                    pressure: 0
                turbulent_kinetic_energy:
                    option: value
                    turbulent_kinetic_energy: 0.047633
                turbulent_eddy_frequency:
                    option: value
                    turbulent_eddy_frequency: 264.63
    solver:
        solver_control:
            basic_settings:
              advection_scheme: high_resolution
              turbulence_numerics: upwind              
              convergence_controls:
                  min_iterations: 1
                  max_iterations: 1000
                  physical_timescale: 1
                  relaxation_parameters:
                      velocity_relaxation_factor: 0.9
                      turbulence_relaxation_factor: 0.9
              convergence_criteria:
                  residual_type: RMS
                  residual_target: 1e-12
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
                                coarsen_type: 10
                                interp_type: 6
                                relax_type: 18
                                strong_threshold: 0.25
                                num_sweeps: 1
                                max_levels: 20
                                aggressive_levels: 1
                                trunc_factor: 0.3 
            expert_parameters:  
                consistent: true
        output_control:
            file_path: results.e
            output_frequency: 50
            output_fields: [velocity, pressure, turbulent_kinetic_energy, turbulent_eddy_frequency, turbulent_viscosity, wall_friction_velocity, minimum_distance_to_wall, density, velocity_gradient, transition_onset_reynolds_number, turbulent_intermittency, wall_shear_stress]
            corrected_boundary_values: true
    material_library:
      - name: air
        thermodynamic_properties:
            equation_of_state:
                option: value
                density: 1.0
        transport_properties:
            dynamic_viscosity:
                option: value
                dynamic_viscosity: 1.5e-5
