v = 10
r = 0.75

[Mesh]
  [./gen]
    type = GeneratedMeshGenerator
    dim = 3
    xmin =0
    xmax =10
    ymin =0
    ymax =2
    zmin =0
    zmax =1
    nx=100
    ny=20
    nz=10
  [../]
[]


[Variables]
  [./temp]
  [../]
[]

[ICs]
  [./temp_substrate]
    type = ConstantIC
    variable = temp
    value = 300
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
    thermal_conductivity = thermal_conductivity
  [../]
  [./heatsource]
    type = ADMatHeatSource
    material_property = volumetric_heat
    variable = temp
    scalar = 1
  [../]
[]

[Functions]
  [./heat_source_x]
    type = ParsedFunction
    value ='v*t'
    vars = 'v'
    vals = ${v}
  [../]
  [./heat_source_y]
    type = ParsedFunction
    value = 1
  [../]
  [./heat_source_z]
    type = ParsedFunction
    value = 1
  [../]
  [./specific_heat]
    type = PiecewiseLinear
    data_file = ./../AM_Brick/AM_Brick_parameters/Specific_Heat.csv
    format = columns
    scale_factor = 1.0
  [../]
  [./thermal_conductivity]
    type = PiecewiseLinear
    data_file = ./../AM_Brick/AM_Brick_parameters/Thermal_Conductivity.csv
    format = columns
    scale_factor = 0.05e-3
  [../]
[]

[BCs]
  [./temp_bottom_fix]
    type = ADDirichletBC
    variable = temp
    boundary = back
    value = 300
  [../]

  # [./convective_substrate]
  #   type = ADConvectiveHeatFluxBC
  #   variable = temp
  #   boundary = 'bottom top left right'
  #   heat_transfer_coefficient = 0.016
  #   T_infinity = 300
  # [../]
[]

[Materials]
  [./volumetric_heat]
    type = FunctionPathEllipsoidHeatSource
    r = ${r}
    power = 350
    efficiency = 0.36
    factor = 1
    function_x= heat_source_x
    function_y= heat_source_y
    function_z= heat_source_z
  [../]
  [./density]
    type = ADDensity
    density = 7.609e-6
  [../]
  [./heat]
    type = ADHeatConductionMaterial
    specific_heat_temperature_function = specific_heat
    thermal_conductivity_temperature_function = thermal_conductivity
    temp = temp
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

  automatic_scaling = true

  #Preconditioned JFNK (default)
  solve_type = 'NEWTON'

  # # petsc_options = '-snes_ksp'
  # petsc_options_iname = '-pc_type -ksp_type'
  # petsc_options_value = 'lu  preonly'

  petsc_options_iname = '-ksp_type -pc_type -pc_factor_mat_solver_package'
  petsc_options_value = 'preonly lu       superlu_dist'

  line_search = 'none'

  l_max_its = 100
  nl_max_its = 20
  nl_rel_tol = 1e-5
  nl_abs_tol = 1e-6
  # l_tol = 1e-5

  start_time = 0.0
  end_time = 1.0
  dt = 0.02
  dtmin = 1e-4
[]

[Outputs]
  file_base = './output_r${r}/heat_source_out'
  [./exodus]
    type = Exodus
    interval = 1
  [../]
  # [checkpoint]
  #   type = Checkpoint
  #   num_files = 2
  # []
[]
