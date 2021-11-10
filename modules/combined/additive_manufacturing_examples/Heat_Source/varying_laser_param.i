# sample #40

v = 10.58e-3 # 10 mm/s = 10e-3 mm/ms
r = 100e-3 # 200 um = 100e-3 mm
power = 300e-3 # 300W = kg*m^2/s^3 = 300e-3 kg*mm^2/ms^3
scan_type = mixed
dt = 50

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
    data_file = Specific_Heat.csv # J/kg/K -> kg*m^2/s^2/kg/K -> mm^2/ms^2/K
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
    heat_source_type = ${scan_type}
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
  # l_tol = 1e-5

  start_time = 0.0
  end_time = 500 # ms
  dt = ${dt} # ms
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

[VectorPostprocessors]
  [line]
    type = LineValueSampler
    variable = 'temp'
    start_point = '0 0.5 0.5'
    end_point = '5 0.5 0.5'
    num_points = 100
    sort_by = x
    outputs = vpp
    execute_on = 'INITIAL TIMESTEP_END'
  []
[]

[Outputs]
  [exodus]
    type = Exodus
    interval = 1
    file_base = './output_r${r}/heat_source_out_${power}W_${scan_type}_dt${dt}'
  []
  [vpp]
    type = CSV
    delimiter = ','
    sync_times = '300'
    sync_only = true
    time_data = true
    file_base = './output_r${r}/heat_source_out_${power}W_${scan_type}_dt${dt}_vpp/center_temperature'
  []
[]
