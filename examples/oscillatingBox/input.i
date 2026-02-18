mesh:
    file_path: mesh.e
    automatic_decomposition_type: rcb
simulation:
    verbose: 1
    physical_analysis:
        analysis_type:
            option: transient
            total_time: 1
            time_steps:
                option: constant
                timestep: 0.1
        domains:
          - name: fluid
            location: [fluid]
            materials: [air]
            type: fluid
            domain_models:
                reference_pressure: 101325
                mesh_deformation:
                    option: regions_of_motion_specified  
                    mesh_motion_model:
                        option: displacement_diffusion
                        mesh_stiffness: 
                            option: blended_distance_and_small_volumes
                            distance_exponent: 0.5
                            volume_exponent: 2.0
                    displacement_relative_to: initial_mesh              
            fluid_models:
                turbulence:
                    option: laminar
            boundaries:
              - name: oscillating_square
                type: wall
                location: [square]   
                boundary_details:
                    mesh_motion:
                        option: periodic_displacement
                        displacement:
                            frequency: 1
                            value: [0,0.1,0]
              - name: walls
                type: wall
                location: [left,right,bottom,top]   
              - name: symmetry
                type: symmetry
                location: [symmetry]
            initialization:
                velocity:
                    option: value
                    velocity: [0,0,0]
                pressure:
                    option: value
                    pressure: 0
    solver:
        solver_control:
            basic_settings:
                advection_scheme: upwind
                transient_scheme: first_order_backward_euler
                reduced_stencil: true
                convergence_controls:
                    min_iterations: 1
                    max_iterations: 50
                convergence_criteria:
                    residual_type: RMS
                    residual_target: 1e-6
            advanced_options:
                pressure_level_information:
                    option: cartesian_coordinates
                    cartesian_coordinates: [-0.5,-0.3,0]
                    relative_pressure_level: 0                
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
                option: timestep_interval
                timestep_interval: 1
            output_fields: [velocity, pressure, displacement_mesh]                
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
