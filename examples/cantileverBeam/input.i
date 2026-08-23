# vim: ft=yaml
# 3D Neo-Hookean cantilever configured to match the solids4foam tutorial.
mesh:
    file_path: mesh.e
    automatic_decomposition_type: rcb
simulation:
    verbose: 1
    physical_analysis:
        analysis_type:
            # The solids4foam case uses steadyState ddt and d2dt2 schemes.
            option: steady_state
        domains:
        - name: solid
          location: [solid]
          materials: [steel]
          type: solid
          solid_models:
            solid_mechanics:
                option: neo_hookean
                formulation: total_lagrangian
                plane_stress: false
          boundaries:
          - name: fixed
            type: wall
            location: [fixed]
            boundary_details:
                solid_mechanics:
                    option: fixed
          - name: free
            type: wall
            location: [free]
            boundary_details:
                solid_mechanics:
                    option: traction
                    pressure: 0
                    shear: [0, 0]
          - name: top
            type: wall
            location: [top]
            boundary_details:
                solid_mechanics:
                    option: traction
                    pressure: 0
                    shear: [1, 0]
          - name: bottom
            type: wall
            location: [bottom]
            boundary_details:
                solid_mechanics:
                    option: traction
                    pressure: 0
                    shear: [0, 0]
          - name: left
            type: wall
            location: [left]
            boundary_details:
                solid_mechanics:
                    option: traction
                    pressure: 0
                    shear: [0, 0]
          - name: right
            type: wall
            location: [right]
            boundary_details:
                solid_mechanics:
                    option: traction
                    pressure: 0
                    shear: [0, 0]
          initialization:
            displacement:
                option: value
                displacement: [0, 0]
    solver:
        solver_control:
            basic_settings:
                convergence_controls:
                    min_iterations: 1
                    max_iterations: 50
                    physical_timescale: 1
                    relaxation_parameters:
                        solid_displacement_relaxation_factor: 0.3
                convergence_criteria:
                    residual_type: RMS
                    residual_target: 1e-8
            advanced_options:
                equation_controls:
                    sub_iterations:
                        solid_displacement: 1
                linear_solver_settings:
                    default:
                        family: PETSc
                        min_iterations: 0
                        max_iterations: 20
                        rtol: 1.0e-4
                        atol: 1.0e-10
                        options:
                            ksp_type: preonly
                            pc_type: lu
        output_control:
            file_path: results.e
            output_frequency: 1
            output_fields: [displacement,stress,strain]
            corrected_boundary_values: true
    material_library:
    - name: steel
      thermodynamic_properties:
        equation_of_state:
            option: value
            density: 1000
      mechanical_properties:
        young_modulus:
            option: value
            young_modulus: 3e6
        poisson_ratio:
            option: value
            poisson_ratio: 0.3
