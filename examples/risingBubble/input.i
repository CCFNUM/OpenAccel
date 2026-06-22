# vim: ft=yaml
# Hysing et al. (2009) 2D rising-bubble benchmark, Case 1.
# Domain [0,1]x[0,2]; bubble r=0.25 at (0.5,0.5); rho1/rho2=1000/100,
# mu1/mu2=10/1, g=0.98, sigma=24.5. No-slip top/bottom, free-slip left/right.
# Reference (t=3): center-of-mass y~1.081, max rise vel~0.2417, min circularity~0.9013.
mesh:
    file_path: mesh.e
    automatic_decomposition_type: rcb
simulation:
    verbose: 1
    physical_analysis:
        analysis_type:
            option: transient
            total_time: 5.0
            time_steps:
                option: adaptive
                initial_timestep: 1.0e-3
                timestep_update_frequency: 1
                timestep_adaptation:
                    option: max_courant
                    courant_number: 0.25
                    min_timestep: 1.0e-6
                    max_timestep: 0.02
                    timestep_decrease_factor: 0.8
                    timestep_increase_factor: 1.05
        domains:
        - name: default_domain
          location: [fluid]
          # bubble (fluid2) listed first => secondary phase (gets smoothed alpha,
          # needed for curvature); ambient (fluid1) last => primary phase.
          materials: [fluid2, fluid1]
          type: fluid
          domain_models:
            reference_pressure: 0
            buoyancy_model:
                option: buoyant
                gravity: [0, -0.98]
                buoyancy_reference_density: 1000
          fluid_models:
            turbulence:
                option: laminar
            multiphase:
                homogeneous: true
                free_surface_model:
                    option: standard
                    flux_corrected_transport: true
                    interface_compression_level: 2
                    n_alpha_corrections: 2
          fluid_pair_models:
          - pair: [fluid2, fluid1]
            surface_tension:
                option: continuum_surface_force
                surface_tension_coefficient: 24.5
          boundaries:
          - name: bottom
            type: wall
            location: [bottom]
          - name: top
            type: wall
            location: [top]
          - name: left
            type: symmetry
            location: [left]
          - name: right
            type: symmetry
            location: [right]
          initialization:
            velocity:
                option: value
                velocity: [0, 0]
            pressure:
                option: value
                pressure: 0
            fluid_specific_initialization:
                fluid2:
                    volume_fraction:
                        option: value
                        input_type: expression
                        volume_fraction: "if ((x-0.5)*(x-0.5) + (y-0.5)*(y-0.5) < 0.0625, 1, 0)"
                fluid1:
                    volume_fraction:
                        option: value
                        input_type: expression
                        volume_fraction: "if ((x-0.5)*(x-0.5) + (y-0.5)*(y-0.5) < 0.0625, 0, 1)"
    solver:
        solver_control:
            basic_settings:
                advection_scheme: high_resolution
                transient_scheme: first_order_backward_euler
                convergence_controls:
                    min_iterations: 1
                    max_iterations: 10
                convergence_criteria:
                    residual_type: RMS
                    residual_target: 1e-6
            advanced_options:
                pressure_level_information:
                    option: cartesian_coordinates
                    cartesian_coordinates: [0, 0]
                    relative_pressure_level: 0
                equation_controls:
                    volume_fraction_smoothing:
                        smooth_volume_fraction: true
                        smoothing_iterations: 3
                        fourier_number: 0.25
                        curvature_smoothing_iterations: 6
                linear_solver_settings:
                    default:
                        family: PETSc
                        min_iterations: 3
                        max_iterations: 20
                        rtol: 1.0e-1
                        atol: 1.0e-12
                        options:
                            ksp_type: fgmres
                            pc_type: bjacobi
                    pressure_correction:
                        family: HYPRE
                        min_iterations: 3
                        max_iterations: 20
                        rtol: 1.0e-3
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
                body_force_redistribution: true
                relax_gradients: true
                incremental_gradient_change: false
        output_control:
            file_path: results.e
            output_frequency:
                option: time_interval
                time_interval: 0.05
            write_timestep_info: true
            output_fields: [velocity, pressure, volume_fraction.fluid2, curvature.fluid2_fluid1, body_forces]
    material_library:
    - name: fluid1
      thermodynamic_properties:
        equation_of_state:
            option: value
            density: 1000
      transport_properties:
        dynamic_viscosity:
            option: value
            dynamic_viscosity: 10
    - name: fluid2
      thermodynamic_properties:
        equation_of_state:
            option: value
            density: 100
      transport_properties:
        dynamic_viscosity:
            option: value
            dynamic_viscosity: 1
