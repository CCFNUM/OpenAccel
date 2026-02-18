mesh:
    file_path: mesh.e
simulation:
    verbose: 1
    physical_analysis:
        analysis_type:
            option: transient
            total_time: 0.5
            time_steps:
                option: adaptive
                initial_timestep: 0.001
                timestep_update_frequency: 1
                timestep_adaptation:
                    option: max_courant
                    courant_number: 1.0
                    min_timestep: 1.0e-6
                    max_timestep: 1
                    timestep_decrease_factor: 0.8
                    timestep_increase_factor: 1.06
        domains:
          - name: default_domain
            location: [fluid-hex]
            materials: [water,air]
            type: fluid
            domain_models:
                reference_pressure: 101325
                buoyancy_model:
                    option: buoyant
                    gravity: [0,-9.81,0]
                    buoyancy_reference_density: 1
            fluid_models:
                turbulence:     
                    option: laminar
                multiphase:
                    homogeneous: true
                    free_surface_model: 
                        option: standard 
                        flux_corrected_transport: true
                        interface_compression_level: 2
                        n_alpha_corrections: 2
            fluid_pair_models:
              - pair: [water, air]
                surface_tension:
                    option: continuum_surface_force
                    surface_tension_coefficient: 0.07
            boundaries:
              - name: walls
                type: wall
                location: [bot,dam,left,right]
              - name: top
                type: opening
                location: [top]
                boundary_details:
                    mass_and_momentum:
                        option: opening_pressure
                        relative_pressure: 0
                    flow_direction:
                        option: cartesian_components
                        x: 0
                        y: -1
                        z: 0
                fluid_values:
                    water:
                        volume_fraction:
                            option: value
                            volume_fraction: 0
                    air:
                        volume_fraction:
                            option: value
                            volume_fraction: 1
              - name: front_and_back
                type: symmetry
                location: [sym]
            initialization:
                velocity:
                    option: value
                    velocity: [0,0,0]
                pressure:
                    option: value
                    pressure: 0                   
                fluid_specific_initialization:
                    water:
                        volume_fraction:
                            option: value
                            input_type: expression
                            volume_fraction: "if (x<=0.152348 and y<=0.290476, 1, 0)"
                    air:
                        volume_fraction:
                            option: value
                            input_type: expression
                            volume_fraction: "if (x<=0.152348 and y<=0.290476, 0, 1)"
    solver:
        solver_control:
            basic_settings:
                advection_scheme: high_resolution
                transient_scheme: first_order_backward_euler
                reduced_stencil: true
                convergence_controls:
                    min_iterations: 1
                    max_iterations: 5    
                convergence_criteria:
                    residual_type: RMS
                    residual_target: 1e-6
                interpolation_scheme:
                    velocity_interpolation_type: linear_linear             
            advanced_options: 
                equation_controls:
                    volume_fraction_smoothing:
                        smooth_volume_fraction: true
                        smoothing_iterations: 3
                        fourier_number: 0.25
                linear_solver_settings:
                    default:
                        family: PETSc
                        max_iterations: 20
                        rtol: 1.0e-1
                        atol: 1.0e-12
                        options:
                            ksp_type: fgmres
                            pc_type: bjacobi                
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
            output_frequency:
                option: time_interval
                time_interval: 0.05
            write_timestep_info: true
            output_fields: [velocity, pressure, volume_fraction.water]
    material_library:
      - name: water
        thermodynamic_properties:
            equation_of_state:
                option: value
                density: 1000
        transport_properties:
            dynamic_viscosity:
                option: value
                dynamic_viscosity: 1e-3  
      - name: air
        thermodynamic_properties:
            equation_of_state:
                option: value
                density: 1
        transport_properties:
            dynamic_viscosity:
                option: value
                dynamic_viscosity: 1.48e-5                
