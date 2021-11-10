# sample #40

speed = 10e-3
power = 300e-3
r = 200e-3
dt = 1

[Mesh]
  [gen]
    type = GeneratedMeshGenerator
    dim = 3
    xmax = 1
    ymax = 5
    zmax = 0.5
    nx = 10
    ny = 50
    nz = 5
  []
[]

[Variables]
  [temp]
  []
[]

[ICs]
  [temp_substrate]
    type = ConstantIC
    variable = temp
    value = 300
  []
[]

[Kernels]
  [time]
    type = ADHeatConductionTimeDerivative
    variable = temp
  []
  [heat_conduct]
    type = ADHeatConduction
    variable = temp
    thermal_conductivity = thermal_conductivity
  []
  [heatsource]
    type = ADMatHeatSource
    material_property = volumetric_heat
    variable = temp
    scalar = 1
  []
[]

[Functions]
  [heat_source_x]
    type = ConstantFunction
    value = '0.5'
  []
  [heat_source_y]
    type = ParsedFunction
    value = '0.5+${speed}*t'
  []
  [heat_source_z]
    type = ConstantFunction
    value = '0.5'
  []
  [specific_heat]
    type = PiecewiseLinear
    x='197.79 298.46 600.31 1401.01 1552.59 1701.44 '
    y='426.69 479.77  549.54 676.94 695.14 726.99'
    format = columns
    scale_factor = 1.0
  []
  [thermal_conductivity]
    type = PiecewiseLinear
    x = '198 298 500 801 1001 1400 1601'
    y = '247.72 285.64 358.44 446.41 491.91 554.09 569.26'
    format = columns
    scale_factor = 0.05e-6
  []
[]

[BCs]
  [temp_bottom_fix]
    type = ADDirichletBC
    variable = temp
    boundary = back
    value = 300
  []

  [convective_substrate]
    type = ADConvectiveHeatFluxBC
    variable = temp
    boundary = 'bottom top left right front'
    heat_transfer_coefficient = 16e-6 # W/m^2/K ->
    T_infinity = 300
  []
[]

[Materials]
  [volumetric_heat]
    type = FunctionPathEllipsoidHeatSource
    r = ${r}
    power = ${power}
    efficiency = 0.36
    factor = 0.5
    function_x = heat_source_x
    function_y = heat_source_y
    function_z = heat_source_z
    heat_source_type = 'line'
    threshold_length = 0.1 #mm
    number_time_integration = 10
  []
  [density]
    type = ADDensity
    density = 7609e-9 # kg/m^3 -> 1e-9 kg/mm^3
  []
  [heat]
    type = ADHeatConductionMaterial
    specific_heat_temperature_function = specific_heat
    thermal_conductivity_temperature_function = thermal_conductivity
    temp = temp
  []
[]

[Preconditioning]
  [smp]
    type = SMP
    full = true
  []
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

  start_time = 0.0
  end_time = '${fparse 5/speed}'
  dt = ${dt}
  num_steps = 40
  dtmin = 1e-6
[]

[Outputs]
  file_base = 'output/no_activation_amr/out'
  [exodus]
    type = Exodus
    # execute_on = 'INITIAL TIMESTEP_END'
    interval = 1
  []
[]

[Adaptivity]
  steps = 1
  marker = marker
  initial_marker = marker
  max_h_level = 1
  [Indicators/indicator]
    type = GradientJumpIndicator
    variable = temp
    scale_by_flux_faces = true
  []
  [Markers/marker]
    type = ErrorFractionMarker
    indicator = indicator
    coarsen = 0.2
    refine = 0.5
  []
[]
