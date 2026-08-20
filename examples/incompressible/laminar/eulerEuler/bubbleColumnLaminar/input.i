# vim: ft=yaml
# Euler-Euler bubble column -- OpenAccel equivalent of the OpenFOAM tutorial
#   $FOAM_TUTORIALS/multiphaseEuler/bubbleColumnLaminar
#
# Matched to the OpenFOAM case:
#   geometry   0.15 m x 1.0 m          (blockMeshDict vertices)
#   mesh       12 x 80 QUAD4           (included refined comparison mesh)
#   phases     air dispersed in water, d = 3 mm
#   drag       Schiller-Naumann, no lift / turbulent dispersion
#   turbulence laminar for both phases (momentumTransport.* simulationType)
#   gravity    (0, -9.81)              (constant/g)
#   initial    smoothed water/air interface centred at y = 0.701
#              U.air = (0, 0.02), U.water = (0, 0)
#   inlet      bottom: alpha.air = 0.5, U.air = (0, 0.02), U.water = (0, 0)
#   outlet     top: opening total pressure, pressure-driven velocity
#   walls      left/right: no-slip
#
# This setup has been run to 5 s on 3x20, 6x40, and 12x80 meshes. First-order
# upwind is intentional: the current high-resolution CVFEM limiter develops
# isolated phasic-velocity spikes under refinement, while bounded upwind gives
# mesh-stable bulk profiles close to the matched OpenFOAM calculation.
mesh:
    file_path: mesh.e
    automatic_decomposition_type: rcb
simulation:
    verbose: 1
    physical_analysis:
        analysis_type:
            option: transient
            total_time: 5
            time_steps:
                option: constant
                # Segregated Euler-Euler with buoyancy-driven phase slip is
                # transient-stability limited, not just accuracy limited: the
                # air phase sees alpha*(rho_air - rho_ref)*g against a
                # transient coefficient of only alpha*rho_air, so its
                # acceleration is O(1e3 m/s^2) until drag balances it. 2.5e-3
                # diverges; 1e-3 is stable and runs 0.1 s in a few seconds.
                timestep: 0.001
        domains:
        - name: fluid
          location: [fluid]
          # air listed first => secondary phase (transported alpha);
          # water last => primary phase, from the closure sum(alpha_k) = 1.
          materials: [air, water]
          type: fluid
          domain_models:
            reference_pressure: 0
            buoyancy_model:
                option: buoyant
                gravity: [0, -9.81]
                buoyancy_reference_density: 1000
          fluid_models:
            turbulence:
                option: laminar
            multiphase:
                model: euler_euler
                homogeneous: false
                primary_phase: water
                residual_volume_fraction: 1.0e-6
          fluid_pair_models:
          - pair: [air, water]
            drag:
                option: schiller_naumann
                dispersed_phase: air
                diameter: 3.0e-3
                blending:
                    inverted_diameter: 1.0e-4
                    # Match OpenFOAM's `blending segregated` thresholds. Their
                    # overlap activates the Marschall interface-drag branch.
                    min_partly_continuous: 0.5
                    min_fully_continuous: 0.7
                    segregated_m: 0.5
                    segregated_n: 8.0

          boundaries:
          - name: wall_l
            type: wall
            location: [wall_l]
          - name: wall_r
            type: wall
            location: [wall_r]
          - name: inlet
            type: inlet
            location: [inlet]
            boundary_details:
                mass_and_momentum:
                    option: velocity_components
                    u: 0
                    v: 0.02
            fluid_values:
                air:
                    velocity:
                        option: value
                        velocity: [0, 0.02]
                    volume_fraction:
                        option: value
                        volume_fraction: 0.5
                water:
                    velocity:
                        option: value
                        velocity: [0, 0]
                    volume_fraction:
                        option: value
                        volume_fraction: 0.5
          - name: outlet
            type: opening
            location: [outlet]
            boundary_details:
                mass_and_momentum:
                    option: opening_pressure
                    relative_pressure: 0
                    # Matches OpenFOAM's prghTotalPressure fields U.air,
                    # phi.air, and rho.air during reverse flow.
                    reference_phase: air
                flow_direction:
                    option: cartesian_components
                    x: 0
                    y: -1
            fluid_values:
                air:
                    volume_fraction:
                        option: value
                        volume_fraction: 1
                water:
                    volume_fraction:
                        option: value
                        volume_fraction: 0
          initialization:
            velocity:
                option: value
                velocity: [0, 0]
            pressure:
                option: value
                pressure: 0
            fluid_specific_initialization:
                air:
                    velocity:
                        option: value
                        velocity: [0, 0.02]
                    volume_fraction:
                        option: value
                        input_type: expression
                        # Smooth over about two cells and avoid pure-phase
                        # singularities while diagnosing solver stability.
                        volume_fraction: "0.5 + 0.49 * tanh((y - 0.701) / 0.02)"
                water:
                    velocity:
                        option: value
                        velocity: [0, 0]
                    volume_fraction:
                        option: value
                        input_type: expression
                        volume_fraction: "0.5 - 0.49 * tanh((y - 0.701) / 0.02)"
    solver:
        solver_control:
            basic_settings:
                advection_scheme: upwind
                transient_scheme: first_order_backward_euler
                convergence_controls:
                    relaxation_parameters:
                        # Damp the localized phasic-momentum mode at the
                        # moving alpha interface. This acts on each momentum
                        # correction, rather than clipping the final velocity.
                        velocity_relaxation_factor: 0.5
                        pressure_relaxation_factor: 0.2
                    min_iterations: 2
                    # Match OpenFOAM's nOuterCorrectors = 3. Continuing a
                    # transient corrector to 30 iterations can pass the
                    # residual minimum and excite a nonphysical phase mode.
                    max_iterations: 3
                convergence_criteria:
                    residual_type: RMS
                    residual_target: 1e-8
            advanced_options:
                pressure_level_information:
                    option: cartesian_coordinates
                    cartesian_coordinates: [0.075, 0.5]
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
                        family: HYPRE
                        min_iterations: 3
                        max_iterations: 100
                        rtol: 1.0e-6
                        atol: 1.0e-14
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
        output_control:
            file_path: results.e
            output_frequency:
                option: time_interval
                # 501 frames over 5 s so ParaView can animate the transient
                time_interval: 0.01
            write_timestep_info: true
            output_fields: [velocity.air, velocity.water, pressure,
                            volume_fraction.air, volume_fraction.water]
    material_library:
    - name: air
      thermodynamic_properties:
        equation_of_state:
            option: value
            density: 1.2
      transport_properties:
        dynamic_viscosity:
            option: value
            dynamic_viscosity: 1.84e-5
    - name: water
      thermodynamic_properties:
        equation_of_state:
            option: value
            density: 1000
      transport_properties:
        dynamic_viscosity:
            option: value
            dynamic_viscosity: 3.645e-4
