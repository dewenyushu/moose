v = 100e-3 # 100 mm/s = 100e-3 mm/ms
r = 50e-3 # 50 um = 50e-3 mm
power = 45e-3 # 45W = 45kg*m^2/s^3 = 45e-3 kg*mm^2/ms^3
scan_type = mixed

[Mesh]
  [gen]
    type = GeneratedMeshGenerator
    dim = 3
    xmin = 0
    xmax = 5
    ymin = 0
    ymax = 1
    zmin = 0
    zmax = 0.5
    nx = 100
    ny = 20
    nz = 10
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
    type = ParsedFunction
    value = 'v*t'
    vars = 'v'
    vals = ${v}
  []
  [heat_source_y]
    type = ParsedFunction
    value = 0.5
  []
  [heat_source_z]
    type = ParsedFunction
    value = 0.5
  []
  [specific_heat]
    type = PiecewiseLinear
    data_file = Specific_Heat.csv # kg*m^2/s^2 -> kg*mm^2/ms^2
    format = columns
    scale_factor = 1.0
  []
  [thermal_conductivity]
    type = PiecewiseLinear
    data_file = Thermal_Conductivity.csv #W/mk -> kg*m^2/s^3/m/K ->kg*m/s^3/K -> 1e-6 kg*mm/ms^3/K
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

  # [./convective_substrate]
  #   type = ADConvectiveHeatFluxBC
  #   variable = temp
  #   boundary = 'bottom top left right front'
  #   heat_transfer_coefficient = 16e-6
  #   T_infinity = 300
  # [../]
[]

[Materials]
  [volumetric_heat]
    type = FunctionPathEllipsoidHeatSource
    r = ${r}
    power = ${power}
    efficiency = 0.7
    factor = 1
    function_x = heat_source_x
    function_y = heat_source_y
    function_z = heat_source_z
    heat_source_type = ${scan_type}
    threshold_length = 0.01
  []
  [density]
    type = ADDensity
    density = 7.609e-6
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
  # l_tol = 1e-5

  start_time = 0.0
  end_time = 50 # ms
  dt = 10 # ms
  dtmin = 1e-4
[]

# [Adaptivity]
#   steps = 1
#   marker = marker
#   initial_marker = marker
#   max_h_level = 1
#   [Indicators/indicator]
#     type = GradientJumpIndicator
#     variable = temp
#     scale_by_flux_faces = true
#   []
#   [Markers/marker]
#     type = ErrorFractionMarker
#     indicator = indicator
#     coarsen = 0.5
#     refine = 0.7
#   []
# []

[Outputs]
  file_base = './output_r${r}/heat_source_out_${power}W_${scan_type}'
  [exodus]
    type = Exodus
    interval = 1
  []
  # [checkpoint]
  #   type = Checkpoint
  #   num_files = 2
  # []
[]
