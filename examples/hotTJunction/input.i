# vim: ft=yaml
# This is a 2D case and must be run with a 2D-compiled binary.
mesh:
    file_path: mesh.e
    automatic_decomposition_type: rcb
simulation:
    verbose: 1
    physical_analysis:
        analysis_type:
            option: transient
            total_time: 0.2
            time_steps:
                option: constant
                timestep: 0.001
        domains:
        - name: fluid
          location: [fluid-hex]
          materials: [air]
          type: fluid
          domain_models:
            reference_pressure: 0
            buoyancy_model:
                option: buoyant
                gravity: [0, -9.81]
                buoyancy_reference_temperature: 303
            mesh_deformation:
                option: regions_of_motion_specified
                mesh_motion_model:
                    option: displacement_diffusion
                    mesh_stiffness:
                        option: increase_near_small_volumes
                        model_exponent: 2
                displacement_relative_to: initial_mesh
          fluid_models:
            turbulence:
                option: laminar
            heat_transfer:
                option: thermal_energy
                include_viscous_work: false
          boundaries:
          - name: inlet1
            type: inlet
            location: [inlet1]
            boundary_details:
                mass_and_momentum:
                    option: velocity_components
                    u: 0.1
                    v: 0
                heat_transfer:
                    option: static_temperature
                    static_temperature: 280
          - name: inlet2
            type: inlet
            location: [inlet2]
            boundary_details:
                mass_and_momentum:
                    option: velocity_components
                    input_type: expression
                    velocity: ["0", "-10*min(t, 0.2 - t)"]
                heat_transfer:
                    option: static_temperature
                    static_temperature: 320
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
            temperature:
                option: value
                temperature: 280
        - name: solid
          location: [solid-hex]
          materials: [wall]
          type: solid
          solid_models:
            heat_transfer:
                option: thermal_energy
            solid_mechanics:
                option: linear_elastic
                plane_stress: false
                reference_temperature: 280
          boundaries:
          - name: fixed
            type: wall
            location: [fixed]
            boundary_details:
                solid_mechanics:
                    option: fixed
                heat_transfer:
                    option: temperature
                    fixed_temperature: 280
          - name: free
            type: wall
            location: [free]
            boundary_details:
                solid_mechanics:
                    option: traction
                    pressure: 0
                    shear: [0, 0]
          initialization:
            displacement:
                option: value
                displacement: [0, 0]
            temperature:
                option: value
                temperature: 280
        interfaces:
        - name: interface
          option: general_connection
          type: fluid_solid
          search_tolerance: 0.000005
          side1:
            domain: fluid
            region_list: [fluid_interface]
          side2:
            domain: solid
            region_list: [solid_interface]
    solver:
        solver_control:
            basic_settings:
                advection_scheme: high_resolution
                transient_scheme: first_order_backward_euler
                convergence_controls:
                    min_iterations: 1
                    max_iterations: 50
                    relaxation_parameters:
                        velocity_relaxation_factor: 0.7
                        pressure_relaxation_factor: 0.3
                convergence_criteria:
                    residual_type: RMS
                    residual_target: 1e-9
                    physics_convergence:
                        enabled: true
                        criteria: [fsi_interface_residual]
                        targets:
                            fsi_interface_residual: 1e-5
            advanced_options:
                equation_controls:
                    sub_iterations:
                        segregated_flow: 20
                        solid_displacement: 1
                    acceleration:
                        solid_displacement:
                            option: aitken
                            initial_omega: 0.5
                            omega_min: 0.1
                            omega_max: 1
                    mesh_motion:
                        freeze_per_timestep: false
                        max_smoothing_iters: 6
                linear_solver_settings:
                    default:
                        family: Trilinos
                        min_iterations: 3
                        max_iterations: 200
                        rtol: 1.0e-6
                        atol: 1.0e-14
                        options:
                            belos_solver: gmres
                            preconditioner: ilu
                    pressure_correction:
                        family: PETSc
                        max_iterations: 200
                        rtol: 1.0e-6
                        atol: 1.0e-14
                        options:
                            ksp_type: gmres
                            pc_type: hypre
                            pc_hypre_type: boomeramg
                    solid_displacement:
                        family: PETSc
                        max_iterations: 1
                        rtol: 1.0e-12
                        atol: 1.0e-20
                        options:
                            ksp_type: preonly
                            pc_type: lu
            expert_parameters:
                bandwidth_reduction: false
                solid_assembler_type: sfem
        output_control:
            file_path: results.e
            output_frequency:
                option: timestep_interval
                timestep_interval: 10
            output_fields: [velocity, pressure, temperature, displacement, stress_von_mises, stress, displacement_mesh]
    material_library:
    - name: air
      thermodynamic_properties:
        equation_of_state:
            option: value
            density: 1
        specific_heat_capacity:
            option: value
            specific_heat_capacity: 1010
      transport_properties:
        dynamic_viscosity:
            option: value
            dynamic_viscosity: 2.1e-5
        thermal_conductivity:
            option: value
            thermal_conductivity: 0.03
      buoyancy_properties:
        thermal_expansivity:
            option: value
            thermal_expansivity: 2.85e-3
    - name: wall
      thermodynamic_properties:
        equation_of_state:
            option: value
            density: 1
        specific_heat_capacity:
            option: value
            specific_heat_capacity: 1010
      transport_properties:
        thermal_conductivity:
            option: value
            thermal_conductivity: 0.04
      thermal_properties:
        thermal_expansion_coefficient:
            option: value
            thermal_expansion_coefficient: 1e-5
      mechanical_properties:
        young_modulus:
            option: value
            young_modulus: 5e3
        poisson_ratio:
            option: value
            poisson_ratio: 0.3
