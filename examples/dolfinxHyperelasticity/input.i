# OpenAccel reproduction of the DOLFINx hyperelasticity tutorial final load.

mesh:
    file_path: mesh.e
    automatic_decomposition_type: rcb
simulation:
    verbose: 1
    physical_analysis:
        analysis_type:
            option: steady_state
        domains:
        - name: beam
          location: [beam]
          materials: [material]
          type: solid
          solid_models:
            solid_mechanics:
                option: neo_hookean
                formulation: total_lagrangian
                plane_stress: false
          boundaries:
          - name: left
            type: wall
            location: [left]
            boundary_details:
                solid_mechanics:
                    option: fixed
          - name: right
            type: wall
            location: [right]
            boundary_details:
                solid_mechanics:
                    option: traction
                    pressure: 0
                    shear: [0, 0, -10]
          - name: traction_free
            type: wall
            location: [bottom, top, front, back]
            boundary_details:
                solid_mechanics:
                    option: traction
                    pressure: 0
                    shear: [0, 0, 0]
          initialization:
            displacement:
                option: value
                displacement: [0, 0, 0]
    solver:
        solver_control:
            basic_settings:
                convergence_controls:
                    min_iterations: 1
                    max_iterations: 200
                    physical_timescale: 1
                    relaxation_parameters:
                        solid_displacement_relaxation_factor: 0.3
                convergence_criteria:
                    residual_type: RMS
                    residual_target: 1e-8
            advanced_options:
                equation_controls:
                    sub_iterations:
                        solid_displacement: 10
                linear_solver_settings:
                    default:
                        family: PETSc
                        min_iterations: 0
                        max_iterations: 10
                        rtol: 1.0e-8
                        atol: 1.0e-12
                        options:
                            ksp_type: preonly
                            pc_type: lu
        output_control:
            file_path: results.e
            output_frequency: 1
            output_fields: [displacement]
            corrected_boundary_values: true
    material_library:
    - name: material
      thermodynamic_properties:
        equation_of_state:
            option: value
            density: 1
      mechanical_properties:
        young_modulus:
            option: value
            young_modulus: 1e4
        poisson_ratio:
            option: value
            poisson_ratio: 0.3
