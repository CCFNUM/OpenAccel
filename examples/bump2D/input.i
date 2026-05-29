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
            materials: [fluid_1]
            type: fluid
            domain_models:
                reference_pressure: 101325
            fluid_models:
                turbulence:     
                    option: shear-stress-transport
            boundaries:
              - name: bump
                type: wall
                location: [bump]
              - name: inlet
                type: inlet
                location: [inlet]
                boundary_details:
                    mass_and_momentum:
                        option: velocity_components
                        u: 69.44
                        v: 0
                        w: 0
                    turbulence:
                        option: k_and_omega
                        k: 1.08e-3
                        omega: 5220.8
              - name: outlet
                type: outlet
                location: [outlet]
                boundary_details:
                    mass_and_momentum:
                        option: static_pressure
                        relative_pressure: 0
              - name: symmetry
                type: symmetry
                location: [frontandback,symup,symdown,top]
            initialization:
                velocity:
                    option: value
                    velocity: [69.44,0,0]
                pressure:
                    option: value
                    pressure: 0
                turbulent_kinetic_energy:
                    option: value
                    turbulent_kinetic_energy: 1.08e-3
                turbulent_eddy_frequency:
                    option: value
                    turbulent_eddy_frequency: 5220.8    
    solver:
        solver_control:
            basic_settings:
                advection_scheme: high_resolution
                turbulence_numerics: high_resolution
                convergence_controls:
                    min_iterations: 1
                    max_iterations: 100
                    physical_timescale: 0.1
                convergence_criteria:
                    residual_type: RMS
                    residual_target: 1e-12
                interpolation_scheme:
                    pressure_interpolation_type: linear_linear
                    velocity_interpolation_type: linear_linear
                    wall_scale_interpolation_type: linear_linear
                    turbulent_kinetic_energy_interpolation_type: linear_linear
                    turbulent_eddy_frequency_interpolation_type: linear_linear
            advanced_options:  
                linear_solver_settings:
                    default:
                        family: AMGSolver
                        min_iterations: 3
                        max_iterations: 20
                        rtol: 1.0e-1
                        atol: 1.0e-12
                        options:
                            smoother_type: dilu
                            cycle_type: no-cycle
                    coupled_flow:
                        family: AMGSolver
                        min_iterations: 3
                        max_iterations: 20
                        rtol: 1.0e-1
                        atol: 1.0e-12
                        options:
                            smoother_type: dilu
                            cycle_type: r-cycle  
                    wall_scale:
                        family: AMGSolver
                        min_iterations: 3
                        max_iterations: 20
                        rtol: 1.0e-1
                        atol: 1.0e-12
                        options:
                            smoother_type: dilu
                            cycle_type: v-cycle 
        output_control:
            file_path: results.e
            output_frequency: 10
            output_fields: [velocity, pressure, turbulent_kinetic_energy, turbulent_eddy_frequency, wall_scale, minimum_distance_to_wall, turbulent_viscosity]
    material_library:
      - name: fluid_1
        thermodynamic_properties:
            equation_of_state:
                option: value
                density: 1.0
        transport_properties:
            dynamic_viscosity:
                option: value
                dynamic_viscosity: 2.31e-5 
