# vim: ft=yaml
# This is a 3D case and must be run with a 3D-compiled binary.
mesh:
    file_path: mesh.e
    automatic_decomposition_type: rcb
simulation:
    verbose: 1
    physical_analysis:
        analysis_type:
            option: transient
            time_steps:
                option: constant
                timestep: 0.05
            total_time: 5
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
                reference_temperature: 300
          boundaries:
          - name: inside
            type: wall
            location: [inside]
            boundary_details:
                heat_transfer:
                    option: temperature
                    input_type: expression
                    fixed_temperature: "300 + 40*t"
                solid_mechanics:
                    option: traction
                    input_type: expression
                    pressure: "0.2e6*t"
                    shear: ["0", "0", "0"]
          - name: outside
            type: wall
            location: [outside]
            boundary_details:
                heat_transfer:
                    option: heat_transfer_coefficient
                    heat_transfer_coefficient: 90
                    outside_temperature: 300
                solid_mechanics:
                    option: traction
                    pressure: 0
                    shear: [0, 0, 0]
          - name: symmetry
            type: symmetry
            location: [symmx, symmy, symmz]
          initialization:
            displacement:
                option: value
                displacement: [0, 0, 0]
            temperature:
                option: value
                temperature: 300
    solver:
        solver_control:
            basic_settings:
                transient_scheme: first_order_backward_euler
                convergence_controls:
                    min_iterations: 1
                    max_iterations: 50
                convergence_criteria:
                    residual_type: RMS
                    residual_target: 1e-7
            advanced_options:
                linear_solver_settings:
                    default:
                        family: PETSc
                        max_iterations: 200
                        rtol: 1.0e-6
                        atol: 1.0e-14
                        options:
                            ksp_type: fgmres
                            pc_type: bjacobi
                    solid_displacement:
                        family: PETSc
                        max_iterations: 200
                        rtol: 1.0e-2
                        atol: 1.0e-14
                        options:
                            ksp_type: gmres
                            pc_type: hypre
                            pc_hypre_type: boomeramg
                            pc_hypre_boomeramg_numfunctions: 3
                            pc_hypre_boomeramg_strong_threshold: 0.7
                            pc_hypre_boomeramg_coarsen_type: HMIS
                            pc_hypre_boomeramg_interp_type: ext+i
                            pc_hypre_boomeramg_agg_nl: 1
                            pc_hypre_boomeramg_p_max: 1
                            pc_hypre_boomeramg_truncfactor: 0.3
            expert_parameters:
                bandwidth_reduction: false
        output_control:
            file_path: results.e
            output_frequency:
                option: timestep_interval
                timestep_interval: 20
            output_fields: [displacement, stress, strain, temperature, stress_von_mises]
    material_library:
    - name: steel
      thermodynamic_properties:
        equation_of_state:
            option: value
            density: 7750
        specific_heat_capacity:
            option: value
            specific_heat_capacity: 486
      transport_properties:
        thermal_conductivity:
            option: value
            thermal_conductivity: 20
      thermal_properties:
        thermal_expansion_coefficient:
            option: value
            thermal_expansion_coefficient: 9.7e-6
      mechanical_properties:
        young_modulus:
            option: value
            young_modulus: 190e9
        poisson_ratio:
            option: value
            poisson_ratio: 0.305
