# vim: ft=yaml
# Euler-Euler batch settling column.
#
# A dilute suspension of heavy particles (alpha = 0.05, rho = 2500, d = 50 um)
# in water settles under gravity in a closed 4 mm x 50 mm column. Only drag
# couples the phases -- no lift, no virtual mass, no mass transfer.
#
# This is the first case that exercises gravity, buoyancy, the alpha*grad(p)
# split, Schiller-Naumann drag and genuine volume-fraction transport, none of
# which slipRelaxation touches.
#
# Steady-state analysis. The column is closed, so the net volumetric flux
# vanishes:  alpha_d U_d + alpha_c U_c = 0. Writing the two momentum balances,
# dividing each by its own alpha and subtracting eliminates grad(p):
#
#     (rho_d - rho_c) g = K s (1/alpha_d + 1/alpha_c)
#
# with s = U_d - U_c the slip. For Schiller-Naumann,
# K = 0.75 C_D rho_c alpha_c alpha_d |s| / d, so every alpha cancels and
#
#     0.75 C_D(Re) rho_c s^2 / d = (rho_d - rho_c) g ,   Re = rho_c |s| d / mu
#
# i.e. the slip is exactly the single-particle terminal velocity u_t, and
#
#     U_d = -(1 - alpha_d) u_t     (settling, downward)
#     U_c = +alpha_d u_t           (return flow, upward)
#
# For these numbers u_t = 1.982957e-3 m/s at Re = 0.0991, which is 2.97 %
# below the Stokes value -- so the case genuinely tests C_D(Re), not just
# Stokes drag. The clear-water interface descends at |U_d| = 1.883809e-3 m/s.
#
# NOTE ON RICHARDSON-ZAKI: this model yields U_d/u_t = (1 - alpha_d)^1, i.e. a
# hindered-settling exponent of exactly 1. The empirical Richardson-Zaki
# exponent is n ~ 4.65 in the Stokes regime. The difference is a missing swarm
# / hindrance correction to the drag, not a solver defect: plain
# Schiller-Naumann is a single-particle law. Do not treat disagreement with
# Richardson-Zaki as a bug -- see README.md.
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
                option: constant
                timestep: 2.0e-3
        domains:
        - name: fluid
          location: [fluid]
          # particles listed first => secondary phase (transported alpha);
          # water last => primary phase, from the closure sum(alpha_k) = 1.
          materials: [particles, water]
          type: fluid
          domain_models:
            reference_pressure: 0
            buoyancy_model:
                option: buoyant
                gravity: [0, 0]
                buoyancy_reference_density: 1000
          fluid_models:
            turbulence:
                option: laminar
            multiphase:
                model: euler_euler
                homogeneous: false
                primary_phase: water
          fluid_pair_models:
          - pair: [particles, water]
            drag:
                option: schiller_naumann
                dispersed_phase: particles
                diameter: 5.0e-5
          boundaries:
          # Closed column: impermeable top and bottom, symmetry on the sides
          # so the problem stays one-dimensional.
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
                particles:
                    velocity:
                        option: value
                        velocity: [0, 0]
                    volume_fraction:
                        option: value
                        volume_fraction: 0.05
                water:
                    velocity:
                        option: value
                        velocity: [0, 0]
                    volume_fraction:
                        option: value
                        volume_fraction: 0.95
    solver:
        solver_control:
            basic_settings:
                advection_scheme: high_resolution
                transient_scheme: first_order_backward_euler
                convergence_controls:
                    # SIMPLE under-relaxation is required; without it the outer
                    # loop over-corrects and the pressure diverges.
                    relaxation_parameters:
                        velocity_relaxation_factor: 0.7
                        pressure_relaxation_factor: 0.3
                    min_iterations: 2
                    max_iterations: 30
                convergence_criteria:
                    residual_type: RMS
                    residual_target: 1e-9
            advanced_options:
                # Closed domain: no boundary sets the pressure level.
                pressure_level_information:
                    option: cartesian_coordinates
                    cartesian_coordinates: [0.002, 0.025]
                    relative_pressure_level: 0
                linear_solver_settings:
                    default:
                        family: PETSc
                        min_iterations: 3
                        max_iterations: 50
                        rtol: 1.0e-4
                        atol: 1.0e-14
                        options:
                            ksp_type: fgmres
                            pc_type: bjacobi
                    # Poisson-like system on a 12.5:1 aspect-ratio column:
                    # fgmres+bjacobi stalls here (100/100 iters, residual drop
                    # ~0.96). AMG is the appropriate preconditioner.
                    pressure_correction:
                        family: PETSc
                        min_iterations: 3
                        max_iterations: 200
                        rtol: 1.0e-10
                        atol: 1.0e-14
                        options:
                            ksp_type: fgmres
                            pc_type: bjacobi
            expert_parameters:
                body_force_redistribution: true
        output_control:
            file_path: results.e
            output_frequency:
                option: time_interval
                time_interval: 0.1
            write_timestep_info: true
            output_fields: [velocity.particles, velocity.water, pressure,
                            volume_fraction.particles, volume_fraction.water]
    material_library:
    - name: particles
      thermodynamic_properties:
        equation_of_state:
            option: value
            density: 2500
      transport_properties:
        dynamic_viscosity:
            option: value
            dynamic_viscosity: 1.0e-3
    - name: water
      thermodynamic_properties:
        equation_of_state:
            option: value
            density: 1000
      transport_properties:
        dynamic_viscosity:
            option: value
            dynamic_viscosity: 1.0e-3
