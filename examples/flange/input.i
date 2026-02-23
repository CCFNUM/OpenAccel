mesh:
    file_path: mesh.e
simulation:
    verbose: 1
    physical_analysis:
        analysis_type:
            option: transient
            total_time: 3
            time_steps:
                option: constant
                timestep: 0.005
        domains:
          - name: default_domain
            location: [solid-hex,solid-pri]
            materials: [undefined]
            type: solid
            solid_models:
                heat_transfer:
                    option: thermal_energy
            sources:
                energy:
                    option: total_source
                    total_source: 0
            boundaries:
              - name: patch1
                type: wall
                location: [patch1]
              - name: patch2
                type: wall
                location: [patch2]
                boundary_details:
                    heat_transfer:
                        option: temperature
                        fixed_temperature: 273
              - name: patch3
                type: wall
                location: [patch3]
              - name: patch4
                type: wall
                location: [patch4]   
                boundary_details:
                    heat_transfer:
                        option: temperature
                        fixed_temperature: 573                                             
            initialization:
                temperature:
                    option: value
                    temperature: 273
    solver:
        solver_control:
            basic_settings:
                transient_scheme: first_order_backward_euler
                convergence_controls:
                    min_iterations: 1
                    max_iterations: 5               
                convergence_criteria:
                    residual_type: RMS
                    residual_target: 1e-6
                interpolation_scheme:
                    temperature_interpolation_type: linear_linear
            advanced_options:
                linear_solver_settings:
                    default:
                        family: PETSc
                        max_iterations: 20
                        rtol: 1.0e-1
                        atol: 1.0e-12
                        options:
                            ksp_type: fgmres
                            pc_type: bjacobi             
        output_control:
            file_path: results.e
            output_frequency:
                option: timestep_interval
                timestep_interval: 20
            output_fields: [temperature]
    material_library:
      - name: undefined
        thermodynamic_properties:
            equation_of_state:
                option: value
                density: 1
            specific_heat_capacity:
                option: value
                specific_heat_capacity: 1           
        transport_properties:
            thermal_conductivity:
                option: value
                thermal_conductivity: 4e-05
