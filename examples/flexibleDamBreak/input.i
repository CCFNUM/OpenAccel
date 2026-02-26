mesh:
    file_path: mesh.e
simulation:
    verbose: 1
    physical_analysis:
        analysis_type:
            option: transient
            total_time: 0.3
            time_steps:
                option: adaptive
                initial_timestep: 0.001
                timestep_update_frequency: 1
                timestep_adaptation:
                    option: max_courant
                    courant_number: 0.5
                    min_timestep: 1.0e-6
                    max_timestep: 1
                    timestep_decrease_factor: 0.8
                    timestep_increase_factor: 1.06
        domains:
          - name: fluid
            location: [fluid-hex]
            materials: [water,air]
            type: fluid
            domain_models:
                reference_pressure: 101325
                buoyancy_model:
                    option: buoyant
                    gravity: [0,-9.81,0]
                    buoyancy_reference_density: 1
                mesh_deformation:
                    option: regions_of_motion_specified  
                    mesh_motion_model:
                        option: displacement_diffusion
                        mesh_stiffness: 
                            option: increase_near_small_volumes
                            model_exponent: 2
                    displacement_relative_to: initial_mesh
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
            boundaries:
              - name: walls
                type: wall
                location: [f_bot,f_left,f_right]
              - name: top
                type: opening
                location: [f_top]
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
              - name: f_symmetry
                type: symmetry
                location: [f_front_and_back]
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
                            volume_fraction: "if (x<=0.146 and y<=0.292, 1, 0)"
                    air:
                        volume_fraction:
                            option: value
                            input_type: expression
                            volume_fraction: "if (x<=0.146 and y<=0.292, 0, 1)"
          - name: solid
            location: [solid-hex]
            materials: [solid_1]
            type: solid
            solid_models:
                solid_mechanics:
                    option: linear_elastic
                    formulation: total_lagrangian
                    plane_stress: false
            boundaries:
              - name: bottom
                type: wall
                location: [s_bot]   
                boundary_details:
                    solid_mechanics:
                        option: fixed
              - name: s_symmetry
                type: symmetry
                location: [s_front_and_back]
            initialization:
                displacement:
                    option: value
                    displacement: [0,0,0]
        interfaces:
          - name: interface
            option: general_connection
            type: fluid_solid
            search_tolerance: 0.0005
            side1:
                domain: fluid
                region_list: [f_int]
            side2:
                domain: solid
                region_list: [s_int]
    solver:
        solver_control:
            basic_settings:
                advection_scheme: high_resolution
                transient_scheme: first_order_backward_euler
                convergence_controls:
                    min_iterations: 1
                    max_iterations: 5                
                convergence_criteria:
                    residual_type: RMS
                    residual_target: 1e-8                  
                interpolation_scheme:
                    velocity_interpolation_type: linear_linear             
            advanced_options: 
                equation_controls:
                    volume_fraction_smoothing:
                        smooth_volume_fraction: true
                        smoothing_iterations: 3
                        fourier_number: 0.25
                    sub_iterations:
                        pressure_correction: 4
                    acceleration:
                        solid_displacement:
                            option: aitken
                            initial_omega: 0.1
                            omega_min: 0.01
                            omega_max: 0.5
                    mesh_motion:
                        freeze_per_timestep: false
                        max_smoothing_iters: 25
                interface_transfer:
                    verbose: 1
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
            expert_parameters:
                body_force_redistribution: false
        output_control:
            file_path: results.e
            output_frequency:
                option: timestep_interval
                timestep_interval: 10
            write_timestep_info: true
            output_fields: [velocity, pressure, volume_fraction.water, displacement_mesh, velocity_mesh, density, dynamic_viscosity]
            corrected_boundary_values: true
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
      - name: solid_1
        thermodynamic_properties:
            equation_of_state:
                option: value
                density: 2.7e3
        mechanical_properties:
            young_modulus:
                option: value
                young_modulus: 1e6
            poisson_ratio:
                option: value
                poisson_ratio: 0     
