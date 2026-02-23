mesh:
    file_path: mesh.e
simulation:
    verbose: 1
    physical_analysis:
        analysis_type:
            option: steady_state
        domains:
          - name: solid
            location: [fluid-hex]
            materials: [structural_steel]
            type: solid
            solid_models:
                solid_mechanics:
                    option: linear_elastic
                    plane_stress: true
            boundaries:
              - name: free
                type: wall
                location: [hole,up]
              - name: right
                type: wall
                location: [right]
                boundary_details:
                    solid_mechanics:
                        option: traction
                        pressure: 0  
                        shear: [10000,0,0]
              - name: symmetry
                type: symmetry
                location: [left,down,front_and_back]
            initialization:
                displacement:
                    option: value
                    displacement: [0,0,0]
    solver:
        solver_control:
            basic_settings:
                convergence_controls:
                    min_iterations: 1
                    max_iterations: 500
                    physical_timescale: 1 # dummy
                convergence_criteria:
                    residual_type: RMS
                    residual_target: 1e-8
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
            output_frequency: 20
            output_fields: [displacement, stress, strain]
    material_library:
      - name: structural_steel
        thermodynamic_properties:
            equation_of_state:
                option: value
                density: 7854
        mechanical_properties:
            young_modulus:
                option: value
                young_modulus: 2e11
            poisson_ratio:
                option: value
                poisson_ratio: 0.3
