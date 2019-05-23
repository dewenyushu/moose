T_room = 300
T_melt = 800
T_ambient = 300

[Mesh]
  [./mesh]
    type = FileMeshGenerator
    # file = ./input_mesh/solid_cylinder_17_substrate_100x100_fine_v0.e
    file = ./input_mesh/solid_cylinder_17_substrate_100x100.e
  [../]
  [./add_set1]
    type = SubdomainBoundingBoxGenerator
    input = mesh
    block_id = 3
    bottom_left = '-50 -50 0'
    top_right = '50 50 5'
  [../]
  [./add_set2]
    type = SubdomainBoundingBoxGenerator
    input = add_set1
    block_id = 1
    bottom_left = '-50 -50 5'
    top_right = '50 50 25'
  [../]
  [./add_set3]
    type = SubdomainBoundingBoxGenerator
    input = add_set2
    block_id = 2
    bottom_left = '0 0 5'
    top_right = '2 2 5'
  [../]
  [./moving_boundary]
    type = SideSetsAroundSubdomainGenerator
    input = add_set3
    block = 2
    new_boundary = 'moving_boundary'
  []
  [./middle]
    type = SideSetsAroundSubdomainGenerator
    input = moving_boundary
    block = 3
    new_boundary = 'middle'
    normal = '0 0 1'
  []
[]

[Problem]
  kernel_coverage_check = false
  # material_coverage_check = false
[]

[Variables]
  [./temp]
    initial_condition = ${T_room}
  [../]
[]

[AuxVariables]
  [./disp_x]
    block = '2 3'
  [../]
  [./disp_y]
    block = '2 3'
  [../]
  [./disp_z]
    block = '2 3'
  [../]
  [./processor_id]
    order = CONSTANT
    family = MONOMIAL
  [../]
[]

[Kernels]
  [./time]
    type = ADHeatConductionTimeDerivative
    variable = temp
  [../]
  [./heat_conduct]
    type = ADHeatConduction
    variable = temp
    use_displaced_mesh = true
    thermal_conductivity = thermal_conductivity
  [../]
  [./heatsource]
    type = ADMatHeatSource
    material_property = volumetric_heat
    variable = temp
    scalar = 1
    use_displaced_mesh = true
  [../]
[]

[AuxKernels]
  [./processor_id_aux]
    type = ProcessorIDAux
    variable = processor_id
    execute_on = timestep_begin
  [../]
[]

[BCs]
  [./bottom_temp]
    type = ADDirichletBC
    variable = temp
    boundary = 1
    value = ${T_room}
  [../]
  [./convective_substrate]
    # type = ConvectiveFluxFunction # Convective flux, e.g. q'' = h*(Tw - Tf)
    type = ADConvectiveHeatFluxBC
    variable = temp
    boundary = '2'
    # coefficient = 2e-5
    heat_transfer_coefficient = 2e-4
    T_infinity = ${T_ambient}
  [../]
  [./convective_middle]
    # type = ConvectiveFluxFunction # Convective flux, e.g. q'' = h*(Tw - Tf)
    type = ADConvectiveHeatFluxBC
    variable = temp
    boundary = 'middle'
    # coefficient = 2e-5
    heat_transfer_coefficient = 2e-4
    T_infinity = ${T_ambient}
  [../]
[]

[MultiApps]
  [thermo_mech]
    type = TransientMultiApp
    positions = '0.0 0.0 0.0'
    input_files = coarse_sub_app.i

    catch_up = true
    max_catch_up_steps = 2
    max_failures = 2
    keep_solution_during_restore = true
    execute_on = 'timestep_end'
  []
[]

[Transfers]
  [to_mech]
    type = MultiAppCopyTransfer
    direction = to_multiapp
    execute_on = 'timestep_end'
    multi_app = thermo_mech
    source_variable = temp
    variable = temp_aux
  []
  [to_disp_x]
    type = MultiAppCopyTransfer
    direction = from_multiapp
    execute_on = 'timestep_end'
    multi_app = thermo_mech
    source_variable = disp_x
    variable = disp_x
  []
  [to_disp_y]
    type = MultiAppCopyTransfer
    direction = from_multiapp
    execute_on = 'timestep_end'
    multi_app = thermo_mech
    source_variable = disp_y
    variable = disp_y
  []
  [to_disp_z]
    type = MultiAppCopyTransfer
    direction = from_multiapp
    execute_on = 'timestep_end'
    multi_app = thermo_mech
    source_variable = disp_z
    variable = disp_z
  []
[]

[Functions]
  [./heat_source_x]
    type = PiecewiseLinear
    data_file = ./input_mat_params/path_eq_t_x_coarse.csv
    format = columns
    scale_factor = 1.0
  [../]
  [./heat_source_y]
    type = PiecewiseLinear
    data_file = ./input_mat_params/path_eq_t_y_coarse.csv
    format = columns
    scale_factor = 1.0
  [../]
  [./heat_source_z]
    type = PiecewiseLinear
    data_file = ./input_mat_params/path_eq_t_z_coarse.csv
    format = columns
    scale_factor = 1.0
  [../]
  [./specific_heat]
    type = PiecewiseLinear
    data_file = ./input_mat_params/Specific_Heat.csv
    format = columns
    scale_factor = 1.0
  [../]
  [./thermal_conductivity]
    type = PiecewiseLinear
    data_file = ./input_mat_params/Thermal_Conductivity.csv
    format = columns
    scale_factor = 0.05e-3
  [../]
[]

[Materials]
  [./heat]
    type = ADHeatConductionMaterial
    specific_heat_temperature_function = specific_heat
    thermal_conductivity_temperature_function = thermal_conductivity
    temp = temp
  [../]
  [./volumetric_heat]
    type = FunctionPathEllipsoidHeatSource
    a = 4
    b = 4
    c = 4
    power = 2000
    efficienty = 0.36
    factor = 1.0
    function_x= heat_source_x
    function_y= heat_source_y
    function_z= heat_source_z
  [../]
  [./density]
    type = ADDensity
    density = 7.609e-6
  [../]
[]

[UserObjects]
  [./activated_elem_uo]
    type = ActivateElementTemp
    execute_on = timestep_begin
    function_x= heat_source_x
    function_y= heat_source_y
    function_z= heat_source_z
    active_subdomain_id = 2
    inactive_subdomain_id = 3
    variable_activation = true
    coupled_var = temp
    activate_value = ${T_melt}
    expand_boundary_name = 'moving_boundary'
  [../]
[]

[Preconditioning]
  [./smp]
    type = SMP
    full = true
  [../]
[]

[Executioner]
  type = Transient
  solve_type = 'NEWTON'

  # automatic_scaling = true

  petsc_options_iname = '-ksp_type -pc_type -pc_factor_mat_solver_package'
  petsc_options_value = 'preonly lu       superlu_dist'

  # petsc_options_iname = '-pc_type -ksp_type -pc_factor_shift_type -pc_factor_shift_amount'
  # petsc_options_value = 'lu  preonly NONZERO 1e-10'

  line_search = 'none'

  l_max_its = 100
  nl_max_its = 20
  nl_rel_tol = 1e-6
  nl_abs_tol = 1e-6

  start_time = 0.0
  end_time = 1.0 # 49
  dt = 0.1
  dtmin = 1e-2

  auto_advance = true # cut time-step when subapp fails
[]

[Outputs]
  file_base = './coarse_output/AM_combined_master'
  [./exodus]
    type = Exodus
    execute_on = 'INITIAL TIMESTEP_END'
  [../]
[]
