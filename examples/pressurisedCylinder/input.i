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
          materials: [cylinder_material]
          type: solid
          solid_models:
            solid_mechanics:
                option: mooney_rivlin
                plane_stress: false
          boundaries:
          - name: internal_pressure
            type: wall
            location: [inner-wall]
            boundary_details:
                solid_mechanics:
                    option: traction
                    pressure: 100e6
                    pressure_type: follower
                    shear: [0, 0]
          - name: external_surface
            type: wall
            location: [outer-wall]
            boundary_details:
                solid_mechanics:
                    option: traction
                    pressure: 0      
                    shear: [0, 0]
          - name: symmetry_planes
            type: symmetry
            location: [symmetry-x, symmetry-y] 
          initialization:
            displacement:
                option: value
                displacement: [0, 0]
    solver:
        solver_control:
            basic_settings:
                convergence_controls:
                    min_iterations: 1
                    max_iterations: 200
                    physical_timescale: 1
                    relaxation_parameters:
                       solid_displacement_relaxation_factor: 0.8
                convergence_criteria:
                    residual_type: RMS
                    residual_target: 1e-16
            advanced_options:
                linear_solver_settings:
                    default:
                        family: PETSc
                        max_iterations: 100
                        rtol: 1.0e-7
                        atol: 1.0e-12
                        options:
                            ksp_type: preonly
                            pc_type: lu
        output_control:
            file_path: results.e
            output_frequency: 10
            output_fields: [displacement, stress, strain]
    material_library:
    - name: cylinder_material
      thermodynamic_properties:
        equation_of_state:
            option: value
            density: 1000 # Dummy density for steady-state
      mechanical_properties:
        c1: 80e6         # 80 MPa
        c2: 20e6         # 20 MPa
        kappa: 4970e6  
        #young_modulus:
          #option: value
          #young_modulus: 1e10 
        #poisson_ratio:
          #option: value
          #poisson_ratio: 0.3

