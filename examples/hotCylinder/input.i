# vim: ft=yaml
# This is a 2D case and must be run with a 2D-compiled binary.
mesh:
    file_path: mesh.e
    automatic_decomposition_type: rcb
simulation:
    verbose: 1
    physical_analysis:
        analysis_type:
            option: steady_state
        domains:
        - name: solid
          location: [fluid-hex]
          materials: [steel]
          type: solid
          solid_models:
            heat_transfer:
                option: thermal_energy
            solid_mechanics:
                option: linear_elastic
                plane_stress: false
                reference_temperature: 100
          boundaries:
          - name: inside
            type: wall
            location: [inside]
            boundary_details:
                heat_transfer:
                    option: temperature
                    fixed_temperature: 200
          - name: outside
            type: wall
            location: [outside]
            boundary_details:
                heat_transfer:
                    option: temperature
                    fixed_temperature: 100
          - name: symmetry
            type: symmetry
            location: [left, bottom]
          initialization:
            displacement:
                option: value
                displacement: [0, 0]
            temperature:
                option: value
                temperature: 150
    solver:
        solver_control:
            basic_settings:
                convergence_controls:
                    min_iterations: 1
                    max_iterations: 2000
                    physical_timescale: 1e4
                convergence_criteria:
                    residual_type: RMS
                    residual_target: 1e-10
            advanced_options:
                linear_solver_settings:
                    default:
                        family: PETSc
                        max_iterations: 20
                        rtol: 1.0e-12
                        atol: 1.0e-20
                        options:
                            ksp_type: preonly
                            pc_type: lu
            expert_parameters:
                bandwidth_reduction: false
        output_control:
            file_path: results.e
            output_frequency: 5
            output_fields: [displacement, stress, strain, temperature, stress_von_mises]
    material_library:
    - name: steel
      thermodynamic_properties:
        equation_of_state:
            option: value
            density: 7854
        specific_heat_capacity:
            option: value
            specific_heat_capacity: 434
      transport_properties:
        thermal_conductivity:
            option: value
            thermal_conductivity: 250
      thermal_properties:
        thermal_expansion_coefficient:
            option: value
            thermal_expansion_coefficient: 1e-5
      mechanical_properties:
        young_modulus:
            option: value
            young_modulus: 2e11
        poisson_ratio:
            option: value
            poisson_ratio: 0.3
