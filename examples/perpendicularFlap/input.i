mesh:
    file_path: mesh.e
simulation:
    verbose: 1
    physical_analysis:
        analysis_type:
            option: transient
            total_time: 5
            time_steps:
                option: constant
                timestep: 0.01
        domains:
          - name: fluid
            location: [fluid]
            materials: [fluid_1]
            type: fluid
            domain_models:
                reference_pressure: 101325
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
            boundaries:
              - name: inlet
                type: inlet
                location: [fluid_inlet]
                boundary_details:
                    mass_and_momentum:
                        option: velocity_components
                        u: 10
                        v: 0
                        w: 0
              - name: outlet
                type: outlet
                location: [fluid_outlet]
                boundary_details:
                    mass_and_momentum:
                        option: static_pressure
                        relative_pressure: 0
              - name: walls
                type: wall
                location: [fluid_lowerwall,fluid_upperwall]
              - name: fluid_symmetry
                type: symmetry
                location: [fluid_frontandback]
            initialization:
                velocity:
                    option: value
                    velocity: [10,0,0]
                pressure:
                    option: value
                    pressure: 0
          - name: solid
            location: [solid]
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
                location: [solid_bottom]   
                boundary_details:
                    solid_mechanics:
                        option: fixed
              - name: solid_symmetry
                type: symmetry
                location: [solid_frontandback]
            initialization:
                displacement:
                    option: value
                    displacement: [0,0,0]
        interfaces:
          - name: interface
            option: general_connection
            type: fluid_solid
            search_tolerance: 0.005
            side1:
                domain: fluid
                region_list: [fluid_flap]
            side2:
                domain: solid
                region_list: [solid_flap]
    solver:
        solver_control:
            basic_settings:
                transient_scheme: first_order_backward_euler
                advection_scheme: upwind
                convergence_controls:
                    min_iterations: 1
                    max_iterations: 50
                    relaxation_parameters:
                        velocity_relaxation_factor: 0.7
                        pressure_relaxation_factor: 0.3
                convergence_criteria:
                    residual_type: RMS
                    residual_target: 1e-6
            advanced_options:
                equation_controls:
                    sub_iterations:
                        segregated_flow: 10
                        solid_displacement: 30
                    acceleration:
                        solid_displacement:
                            option: aitken
                            initial_omega: 0.5
                            omega_min: 0.1
                            omega_max: 1.0
                    mesh_motion:
                        freeze_per_timestep: false
                        max_smoothing_iters: 3
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
        output_control:
            file_path: results.e
            output_frequency:
                option: timestep_interval
                timestep_interval: 10
            output_fields: [velocity, pressure, density, total_pressure, young_modulus, poisson_ratio, velocity_mesh, displacement_mesh]
            corrected_boundary_values: true
            post_process:
              - name: forces
                type: force
                options:
                    calculate_moment: false
                location: [fluid_flap]
                frequency: 1
                write_to_file: true            
    material_library:
      - name: fluid_1
        thermodynamic_properties:
            equation_of_state:
                option: value
                density: 1
        transport_properties:
            dynamic_viscosity:
                option: value
                dynamic_viscosity: 1
      - name: solid_1
        thermodynamic_properties:
            equation_of_state:
                option: value
                density: 3000
        mechanical_properties:
            young_modulus:
                option: value
                young_modulus: 4e6
            poisson_ratio:
                option: value
                poisson_ratio: 0.3
