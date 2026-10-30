mesh:
  file_path: ./box.e
  automatic_decomposition_type: rcb

simulation:
  verbose: 1

  physical_analysis:
    analysis_type:
      option: steady_state

    domains:
      - name: fluid
        location: [background]
        materials: [water]
        type: fluid

        domain_models:
          reference_pressure: 101325

        fluid_models:
          turbulence:
            option: shear_stress_transport

        boundaries:
          - name: inlet
            type: inlet
            location: [inlet]
            boundary_details:
              mass_and_momentum:
                option: velocity_components
                u: 50
                v: 0
                w: 0
              turbulence:
                option: intensity_and_eddy_viscosity_ratio
                fractional_intensity: 0.01
                eddy_viscosity_ratio: 10

          - name: outlet
            type: outlet
            location: [outlet]
            boundary_details:
              mass_and_momentum:
                option: static_pressure
                relative_pressure: 0

          - name: sidewalls
            type: wall
            location: [symmetry]
            boundary_details:
              mass_and_momentum:
                option: free_slip_wall

        initialization:
          velocity:
            option: value
            velocity: [50, 0, 0]

          pressure:
            option: value
            pressure: 0

          turbulent_kinetic_energy:
            option: value
            turbulent_kinetic_energy: 0.24

          turbulent_eddy_frequency:
            option: value
            turbulent_eddy_frequency: 1642.79

    overrides:
      masking:
        - name: sphere
          domains: [fluid]
          forced_layer: 1
          reference_layer: 2
          retain_rings: [1, 2]
          shape:
            option: mesh
            path: sphere.e

  solver:
    solver_control:
      basic_settings:
        advection_scheme: high_resolution
        turbulence_numerics: upwind

        convergence_controls:
          min_iterations: 1
          max_iterations: 1000
          physical_timescale: 0.025
          relaxation_parameters:
            velocity_relaxation_factor: 0.5
            pressure_relaxation_factor: 0.3
            turbulence_relaxation_factor: 0.5

        convergence_criteria:
          residual_type: RMS
          residual_target: 1.0e-6

      advanced_options:
        linear_solver_settings:
          default:
            family: PETSc
            min_iterations: 3
            max_iterations: 50
            rtol: 0.1
            atol: 1.0e-12
            options:
              ksp_type: fgmres
              pc_type: bjacobi

          pressure_correction:
            family: HYPRE
            min_iterations: 3
            max_iterations: 50
            rtol: 0.001
            atol: 1.0e-12
            options:
              solver:
                type: GMRES
              precond:
                type: BoomerAMG
                coarsen_type: 10
                interp_type: 6
                relax_type: 18
                strong_threshold: 0.25
                num_sweeps: 1
                max_levels: 20
                aggressive_levels: 1
                trunc_factor: 0.3

      expert_parameters:
        disable_physics: true
        force_wall_distance_calculation: true
        wall_distance_method: signed_distance_function

    output_control:
      file_path: results.e
      output_frequency: 1
      output_fields: [signed_distance]
      corrected_boundary_values: true

  material_library:
    - name: water
      thermodynamic_properties:
        equation_of_state:
          option: value
          density: 1000

      transport_properties:
        dynamic_viscosity:
          option: value
          dynamic_viscosity: 1.0e-3
