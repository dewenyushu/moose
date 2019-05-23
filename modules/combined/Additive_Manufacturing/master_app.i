[Mesh]
  [./mesh]
    type = FileMeshGenerator
    file = ./input_mesh/solid_cylinder_17_substrate_100x100_fine_v0.e
  [../]
[]

[Variables]
  [./temp]
    initial_condition = 293
  [../]
[]

[AuxVariables]
  [./disp_x]
    initial_condition = 0
  [../]
  [./disp_y]
    initial_condition = 0
  [../]
  [./disp_z]
    initial_condition = 0
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

[BCs]
  [./bottom_temp]
    type = DirichletBC
    variable = temp
    boundary = 1
    value = 293
  [../]
[]

[MultiApps]
  [thermo_mech]
    type = TransientMultiApp
    positions = '0.0 0.0 0.0'
    # input_files = sub_app.i
    input_files = linear_sub_app.i
    # sub_cycling = true
    # steady_state_tol = 1e-6
    # detect_steady_state = true
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
 # [to_temp]
 #   type = MultiAppCopyTransfer
 #   direction = from_multiapp
 #   execute_on = 'timestep_end'
 #   multi_app = thermo_mech
 #   source_variable = temp
 #   variable = temp
 # []
[]

[Functions]
  [./heat_source_x]
    type = PiecewiseLinear
    data_file = ./input_mat_params/path_eq_t_x.csv
    format = columns
    scale_factor = 1.0
  [../]
  [./heat_source_y]
    type = PiecewiseLinear
    data_file = ./input_mat_params/path_eq_t_y.csv
    format = columns
    scale_factor = 1.0
  [../]
  [./heat_source_z]
    type = PiecewiseLinear
    data_file = ./input_mat_params/path_eq_t_z.csv
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
    a = 4.0
    b = 4.0
    c = 4.0
    power = 4500
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

[Preconditioning]
  [./smp]
    type = SMP
    full = true
  [../]
[]

# [Adaptivity]
#   [./Indicators]
#     [./jump]
#       type = ValueJumpIndicator
#       variable = temp
#       outputs = exodus
#     [../]
#   [../]
#   [./Markers]
#     [./error]
#       type = ErrorFractionMarker
#       indicator = jump
#       coarsen = 0.1
#       refine = 0.95
#     [../]
#   [../]
#   marker = error
#   max_h_level = 2
#   cycles_per_step = 2
# []

[Executioner]
  type = Transient
  solve_type = 'NEWTON'

  petsc_options_iname = '-ksp_type -pc_type -pc_factor_mat_solver_package'
  petsc_options_value = 'preonly lu       superlu_dist'

  l_max_its = 100
  nl_rel_tol = 1e-6
  nl_abs_tol = 1e-6

  start_time = 0.0
  end_time = 0.5 # 49
  dt = 0.05
  dtmin = 1e-4

  auto_advance = false # cut time-step when subapp fails
[]

[Outputs]
  exodus = true
[]
