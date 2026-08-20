# vim: ft=yaml
# Euler-Euler slip-relaxation verification case.
#
# Two incompressible phases fill a small periodic box with uniform, constant
# volume fractions and no gravity, no pressure gradient and no mass transfer.
# The phases start at different uniform velocities and are coupled only by a
# constant interphase drag coefficient K. Every field stays spatially uniform,
# so the discretization contributes no error and the case is a direct test of
# the phasic momentum assembly and the drag exchange term.
#
# With m_k = alpha_k * rho_k the momentum equations reduce to
#     m_1 dU_1/dt = -K (U_1 - U_2)
#     m_2 dU_2/dt = +K (U_1 - U_2)
# so the mixture momentum m_1 U_1 + m_2 U_2 is conserved exactly and the slip
# velocity decays as u_r(t) = u_r(0) exp(-t/tau) with
#     tau = 1 / ( K * (1/m_1 + 1/m_2) ).
#
# This case:
#     alpha_1 = 0.4, rho_1 = 1000  ->  m_1 = 400
#     alpha_2 = 0.6, rho_2 =  100  ->  m_2 =  60
#     U_1(0) = (1, 0), U_2(0) = (0, 0),  K = 500
#     tau      = 1200 / 11500      = 0.10434783 s
#     U_inf    = 400 / 460         = 0.86956522 m/s
#     U_1(t)   = U_inf + (m_2/(m_1+m_2)) exp(-t/tau)
#     U_2(t)   = U_inf - (m_1/(m_1+m_2)) exp(-t/tau)
#
# Total time 0.6 s is about 5.75 tau, i.e. the slip decays by ~99.7 %.
# Verify with `python3 verify.py` after running.
mesh:
    file_path: mesh.e
    automatic_decomposition_type: rcb
simulation:
    verbose: 1
    physical_analysis:
        analysis_type:
            option: transient
            total_time: 0.6
            time_steps:
                option: constant
                timestep: 2.5e-4
        domains:
        - name: fluid
          location: [fluid]
          # fluid2 listed first => secondary phase, gets the transported alpha
          # equation; fluid1 last => primary phase, obtained from the closure
          # sum(alpha_k) = 1.
          materials: [fluid2, fluid1]
          type: fluid
          domain_models:
            reference_pressure: 0
            # no buoyancy_model block => no gravity, as the case requires
          fluid_models:
            turbulence:
                option: laminar
            multiphase:
                model: euler_euler
                homogeneous: false
                primary_phase: fluid1
          fluid_pair_models:
          - pair: [fluid2, fluid1]
            drag:
                option: constant
                coefficient: 500.0
          boundaries:
          # per1/per2 carry the translational-periodicity interface below and
          # are therefore not declared as boundaries. The y-normal faces are
          # symmetry planes: the flow is along x with U_y = 0 everywhere, so
          # they are exactly satisfied and inject nothing.
          - name: bottom
            type: symmetry
            location: [bottom]
          - name: top
            type: symmetry
            location: [top]
          initialization:
            velocity:
                option: value
                velocity: [0, 0]
            pressure:
                option: value
                pressure: 0
            fluid_specific_initialization:
                fluid2:
                    velocity:
                        option: value
                        velocity: [0, 0]
                    volume_fraction:
                        option: value
                        volume_fraction: 0.6
                fluid1:
                    velocity:
                        option: value
                        velocity: [1, 0]
                    volume_fraction:
                        option: value
                        volume_fraction: 0.4
        interfaces:
        - name: periodicX
          option: translational_periodicity
          search_tolerance: 1.0e-6
          gauss_lobatto_quadrature: true
          type: fluid_fluid
          side1:
            domain: fluid
            region_list: [per1]
          side2:
            domain: fluid
            region_list: [per2]
    solver:
        solver_control:
            basic_settings:
                advection_scheme: high_resolution
                transient_scheme: first_order_backward_euler
                convergence_controls:
                    # SIMPLE under-relaxation. Without pressure
                    # under-relaxation the outer loop over-corrects and the
                    # pressure field diverges even though momentum, drag and
                    # the mass fluxes are all correct.
                    relaxation_parameters:
                        velocity_relaxation_factor: 0.7
                        pressure_relaxation_factor: 0.3
                    min_iterations: 2
                    max_iterations: 30
                convergence_criteria:
                    residual_type: RMS
                    residual_target: 1e-10
            advanced_options:
                # The domain is periodic in x and symmetric in y, so no
                # boundary sets the pressure level: pin it at one point.
                pressure_level_information:
                    option: cartesian_coordinates
                    # interior node: keep the reference off the periodic
                    # per1/per2 faces
                    cartesian_coordinates: [0.5, 0.5]
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
                    pressure_correction:
                        family: PETSc
                        min_iterations: 3
                        max_iterations: 50
                        rtol: 1.0e-6
                        atol: 1.0e-14
                        options:
                            ksp_type: fgmres
                            pc_type: bjacobi
        output_control:
            file_path: results.e
            output_frequency:
                option: time_interval
                # ~20 samples per tau, and fine enough to resolve the very
                # early transient where the analytic comparison is sharpest
                time_interval: 0.005
            write_timestep_info: true
            output_fields: [velocity.fluid1, velocity.fluid2, pressure,
                            volume_fraction.fluid1, volume_fraction.fluid2,
                            # per-phase drag source: sum over phases must be 0
                            interphase_momentum_source.fluid1,
                            interphase_momentum_source.fluid2,
                            # phasic mass-divergence residuals: must be ~0
                            divergence.fluid1, divergence.fluid2]
    material_library:
    - name: fluid1
      thermodynamic_properties:
        equation_of_state:
            option: value
            density: 1000
      transport_properties:
        dynamic_viscosity:
            option: value
            dynamic_viscosity: 1.0e-3
    - name: fluid2
      thermodynamic_properties:
        equation_of_state:
            option: value
            density: 100
      transport_properties:
        dynamic_viscosity:
            option: value
            dynamic_viscosity: 1.0e-4
