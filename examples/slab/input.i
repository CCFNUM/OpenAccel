mesh:
    file_path: mesh.e
simulation:
    verbose: 1
    physical_analysis:
        analysis_type:
            option: steady_state
        domains:
          - name: fluid
            location: [fluid]
            materials: [air]
            type: fluid
            domain_models:
                reference_pressure: 101325
            fluid_models:
                turbulence:
                    option: laminar
                heat_transfer:
                    option: thermal_energy
                    include_viscous_work: false
            boundaries:
              - name: inlet
                type: inlet
                location: [inlet]
                boundary_details:
                    mass_and_momentum:
                        option: velocity_components
                        u: 0.02083
                        v: 0
                        w: 0
                    heat_transfer:
                        option: static_temperature
                        static_temperature: 300
              - name: outlet
                type: outlet
                location: [outlet]
                boundary_details:
                    mass_and_momentum:
                        option: static_pressure
                        relative_pressure: 0
              - name: wall
                type: wall
                location: [wall]
              - name: slip
                type: wall
                location: [slip, top]
                boundary_details:
                    mass_and_momentum:
                        option: free_slip_wall
              - name: symmetry
                type: symmetry
                location: [fluidfrontandback]
            initialization:
                velocity:
                    option: value
                    velocity: [0.02083,0,0]
                pressure:
                    option: value
                    pressure: 0
                temperature:
                    option: value
                    temperature: 300
          - name: solid
            location: [solid]
            materials: [copper]
            type: solid
            solid_models:
                heat_transfer:
                    option: thermal_energy
            boundaries:
              - name: outer
                type: wall
                location: [solidbottom]   
                boundary_details:
                    heat_transfer:
                        option: temperature
                        fixed_temperature: 400
              - name: leftandright
                type: wall
                location: [solidsides]
              - name: symmetry
                type: symmetry
                location: [solidfrontandback]
            initialization:
                temperature:
                    option: value
                    temperature: 400
        interfaces:
          - name: interface
            option: general_connection
            type: fluid_solid
            search_tolerance: 0.0001
            gauss_lobatto_quadrature: false
            side1:
                domain: fluid
                region_list: [intf]
            side2:
                domain: solid
                region_list: [ints]
    solver:
        solver_control:
            basic_settings:
                advection_scheme: high_resolution
                convergence_controls:
                    min_iterations: 1
                    max_iterations: 1000
                    physical_timescale: 100
                    relaxation_parameters:
                        velocity_relaxation_factor: 0.8
                        pressure_relaxation_factor: 0.2
                convergence_criteria:
                    residual_type: RMS
                    residual_target: 1e-10
            advanced_options:
                linear_solver_settings:
                    default:
                        family: Trilinos
                        max_iterations: 200
                        rtol: 1.0e-6
                        atol: 1.0e-12
                        options:
                            belos_solver: gmres
                            preconditioner: ilu
        output_control:
            file_path: results.e
            output_frequency: 10
            output_fields: [velocity, pressure, temperature]
    material_library:
      - name: air
        thermodynamic_properties:
            equation_of_state:
                option: value
                density: 1.2
            specific_heat_capacity:
                option: value
                specific_heat_capacity: 1006.43
        transport_properties:
            dynamic_viscosity:
                option: value
                dynamic_viscosity: 1e-5  
            thermal_conductivity:
                option: value
                thermal_conductivity: 0.01006
      - name: copper
        thermodynamic_properties:
            equation_of_state:
                option: value
                density: 8800.0
            specific_heat_capacity:
                option: value
                specific_heat_capacity: 420.0
        transport_properties:
            thermal_conductivity:
                option: value
                thermal_conductivity: 0.20129
