# vim: ft=yaml
# This is a 2D case and must be run with a 2D-compiled binary.
mesh:
  file_path: mesh.e
simulation:
  verbose: 1
  physical_analysis:
    analysis_type:
      option: steady_state
    domains:
    - name: default_domain
      location: [fluid]
      materials: [air]
      type: fluid
      domain_models:
        reference_pressure: 101325
      fluid_models:
        turbulence:
          option: laminar
      boundaries:
      - name: wall
        type: wall
        location: [wall]
      - name: inlet1
        type: inlet
        location: [inlet1]
        boundary_details:
          mass_and_momentum:
            option: velocity_components
            u: 1
            v: 0
      - name: inlet2
        type: inlet
        location: [inlet2]
        boundary_details:
          mass_and_momentum:
            option: velocity_components
            u: 0
            v: 3
      - name: outlet
        type: outlet
        location: [outlet]
        boundary_details:
          mass_and_momentum:
            option: static_pressure
            relative_pressure: 0
      initialization:
        velocity:
          option: value
          velocity: [0, 0]
        pressure:
          option: value
          pressure: 0
  solver:
    solver_control:
      basic_settings:
        advection_scheme: high_resolution
        convergence_controls:
          min_iterations: 1
          max_iterations: 100
          physical_timescale: 10
          relaxation_parameters:
              velocity_relaxation_factor: 0.7
              pressure_relaxation_factor: 0.3
        convergence_criteria:
          residual_type: RMS
          residual_target: 1e-6
      advanced_options:
        linear_solver_settings:
          default:
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
      expert_parameters:
        relax_gradients: false                
    output_control:
      file_path: results.e
      output_frequency: 10
      output_fields: [velocity, pressure]
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
