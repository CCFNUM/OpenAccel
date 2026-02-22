mesh:
    file_path: mesh.e
simulation:
    verbose: 1
    physical_analysis:
        analysis_type:
            option: transient
            time_steps:
                option: constant
                timestep: 1
            total_time: 1000
        domains:
          - name: default_domain
            location: [fluid-hex]
            materials: [fluid_1]
            type: fluid
            domain_models:
                reference_pressure: 0                
                buoyancy_model:
                    option: buoyant
                    gravity: [0,-9.81,0]
                    buoyancy_reference_temperature: 300
            fluid_models:
                turbulence:
                    option: laminar
                heat_transfer:
                    option: thermal_energy
            boundaries:
              - name: ceiling
                type: wall
                location: [ceiling]
                boundary_details:
                    heat_transfer:
                        option: temperature
                        fixed_temperature: 300
              - name: floor
                type: wall
                location: [floor]
                boundary_details:
                    heat_transfer:
                        option: temperature
                        fixed_temperature: 301    
              - name: sideWalls
                type: wall
                location: [sideWalls]                   
              - name: frontandback
                type: symmetry
                location: [frontandback]
            initialization:
                velocity:
                    option: value
                    velocity: [1e-4,0,0]
                pressure:
                    option: value
                    pressure: 0
                temperature:
                    option: value
                    temperature: 300
    solver:
        solver_control:
            basic_settings:
                transient_scheme: first_order_backward_euler
                advection_scheme: upwind
                convergence_controls:
                    min_iterations: 1
                    max_iterations: 1
                convergence_criteria:
                    residual_type: RMS
                    residual_target: 1e-10
            advanced_options:
                pressure_level_information:
                    option: cartesian_coordinates
                    cartesian_coordinates: [0,0,0]
                    relative_pressure_level: 0
                linear_solver_settings:
                    default:
                        family: Trilinos
                        min_iterations: 3
                        max_iterations: 50
                        rtol: 1.0e-2
                        atol: 1.0e-12
                        options:
                            belos_solver: gmres
                            preconditioner: ilu                 
        output_control:
            file_path: results.e
            output_frequency:
                option: timestep_interval
                timestep_interval: 50
            output_fields: [velocity, pressure, temperature]
    material_library:
      - name: fluid_1
        thermodynamic_properties:
            equation_of_state:
                option: value
                density: 1
            specific_heat_capacity:
                option: value
                specific_heat_capacity: 1
        transport_properties:
            dynamic_viscosity:
                option: value
                dynamic_viscosity: 1e-03
            thermal_conductivity:
                option: value
                thermal_conductivity: 1e-03
        buoyancy_properties:
            thermal_expansivity:
                option: value     
                thermal_expansivity: 1e-03         
