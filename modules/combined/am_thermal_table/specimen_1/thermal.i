T_room = 303
T_ambient = 303
T_melt = 1563

#speed = 25e-3 # 25 mm/s = 25e-3 mm/ms
r = 2 # 2 mm, TBD
dt = 6 #'${fparse 0.3*r/speed}' # ms
factor = 1.6

refine = 0

[Mesh]
  [mesh]
    type = FileMeshGenerator
    file = SCAN_TRACKS2.e
  []

  [add_set3]
    type = GeneratedMeshGenerator
    dim = 3
    xmax = 0.001
    ymax = 0.001
    zmin = -0.001
    subdomain_ids = 2
  []
  [moving_boundary]
    type = SideSetsAroundSubdomainGenerator
    input = add_set3
    block = 2
    new_boundary = 'moving_boundary'
  []

  [cmbn]
    type = CombinerGenerator
    inputs = 'mesh moving_boundary'
  []

  uniform_refine = ${refine}
[]

[Problem]
  material_coverage_check = false
[]

[Variables]
  [temp]
    block = '1 2 3 4'
  []
[]

[ICs]
  [temp_substrate]
    type = ConstantIC
    variable = temp
    value = ${T_room}
    block = '1 3 4'
  []
  [temp_product]
    type = FunctionIC
    variable = temp
    function = temp_ic
    block = '2'
  []
[]

[AuxVariables]
  [processor_id]
    order = CONSTANT
    family = MONOMIAL
  []
[]

[Kernels]
  [time]
    type = ADHeatConductionTimeDerivative
    variable = temp
  []
  [heat_conduc]
    type = ADHeatConduction
    variable = temp
    # use_displaced_mesh = true
    thermal_conductivity = thermal_conductivity
  []
  [heatsource]
    type = ADMatHeatSource
    material_property = volumetric_heat
    variable = temp
    scalar = 1
    # use_displaced_mesh = true
  []
[]

[AuxKernels]
  [processor_id_aux]
    type = ProcessorIDAux
    variable = processor_id
    execute_on = timestep_begin
  []
[]

[BCs]
  [bottom_temp]
    type = ADDirichletBC
    variable = temp
    boundary = 'fix_temp_table'
    value = ${T_room}
  []
  [convective]
    # type = ConvectiveFluxFunction # Convective flux, e.g. q'' = h*(Tw - Tf)
    type = ADConvectiveHeatFluxBC
    variable = temp
    boundary = 'convective_table moving_boundary'
    heat_transfer_coefficient = 1e-5 # W/m^2/K ->
    T_infinity = ${T_ambient}
  []
[]

[Functions]
  [heat_source_x]
    type = PiecewiseLinear
    data_file = 'SCAN_TRACKS2_X_coord.csv'
    format = columns
  []
  [heat_source_y]
    type = PiecewiseLinear
    data_file = 'SCAN_TRACKS2_Y_coord.csv'
    format = columns
  []
  [heat_source_z]
    type = PiecewiseLinear
    data_file = 'SCAN_TRACKS2_Z_coord.csv'
    format = columns
  []
  [effective_power]
    type = PiecewiseConstant
    data_file = 'SCAN_TRACKS2_Power.csv'
    direction = RIGHT_INCLUSIVE
    format = columns
    scale_factor = 1e-3 #3000W = kg*m^2/s^3 = 300e-3 kg*mm^2/ms^3
  []
  [specific_heat_alloy]
    type = PiecewiseLinear
    data_file = 'Specific_heat.txt'
    format = columns
    scale_factor = 1.0
  []
  [thermal_conductivity_alloy]
    type = PiecewiseLinear
    data_file = 'Thermal_conductivity.txt'
    format = columns
    scale_factor = 1.0e-6
  []
  [density_alloy]
    type = PiecewiseLinear
    data_file = 'Density.txt'
    format = columns
    scale_factor = 1.0e-9
  []
  [specific_heat_substrate]
    type = ConstantFunction
    value = 450.0
  []
  [thermal_conductivity_substrate]
    type = ConstantFunction
    value = 16e-6
  []
  [temp_ic]
    type = ParsedFunction
    expression = 'if(t<=0, temp_room, temp_melt)'
    symbol_names = 'temp_room temp_melt'
    symbol_values = '${T_room} ${T_melt}'
  []
[]

[Materials]
  [heat_alloy]
    type = ADHeatConductionMaterial
    specific_heat_temperature_function = specific_heat_alloy
    thermal_conductivity_temperature_function = thermal_conductivity_alloy
    temp = temp
    block = '1 2'
  []
  [heat_substrate]
    type = ADHeatConductionMaterial
    specific_heat_temperature_function = specific_heat_substrate
    thermal_conductivity_temperature_function = thermal_conductivity_substrate
    temp = temp
    block = '3 4'
  []
  [volumetric_heat_alloy] # TODO: need to separate?
    type = FunctionPathGaussianHeatSource
    r = ${r}
    power = effective_power
    efficiency = 1.0
    factor = ${factor}
    function_x = heat_source_x
    function_y = heat_source_y
    function_z = heat_source_z
    heat_source_type = 'line'
    threshold_length = 2.0 #mm
    number_time_integration = 10
    block = '1 2 3 4'
  []
  [density_alloy]
    type = ADCoupledValueFunctionMaterial
    function = density_alloy
    v = temp
    prop_name = "density"
    block = '1 2'
  []
  [density_substrate]
    type = ADDensity
    density = 7894e-9 # kg/m^3 -> 1e-9 kg/mm^3
    block = '3 4'
  []
[]

[UserObjects]
  [activated_elem_uo_beam]
    type = CoupledVarThresholdElementSubdomainModifier
    execute_on = 'TIMESTEP_BEGIN'
    coupled_var = temp
    block = 1
    subdomain_id = 2
    criterion_type = ABOVE
    threshold = ${T_melt}
    moving_boundary_name = 'moving_boundary'
  []
[]

[Adaptivity]
  marker = marker
  initial_marker = marker
  max_h_level = 4
  [Indicators]
    [indicator]
      type = GradientJumpIndicator
      variable = temp
    []
  []
  [Markers]
    [efm]
      type = ErrorFractionMarker
      indicator = indicator
      coarsen = 0.1
      refine = 0.5
      outputs = none
    []
    [marker]
      type = BoundaryPreservedMarker
      preserved_boundary = moving_boundary
      marker = efm
    []
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
  solve_type = 'PJFNK'

  # automatic_scaling = true

  petsc_options_iname = '-ksp_type -pc_type -pc_factor_mat_solver_package -pc_factor_shift_type -pc_factor_shift_amount'
  petsc_options_value = 'preonly lu       superlu_dist NONZERO 1e-10'

  line_search = 'none'

  l_max_its = 100
  nl_max_its = 15
  nl_rel_tol = 1e-8
  nl_abs_tol = 1e-10

  start_time = 0.0
  end_time = 60910 #125600
  dt = 40 # ms
  dtmin = 1e-6

  auto_advance = true # cut time-step when subapp fails

  error_on_dtmin = false
[]

[Outputs]
  file_base = 'output/Line_thermal_r_${r}_dt_${dt}'
  csv = true
  [exodus]
    type = Exodus
    file_base = 'output/Exodus_r_${r}_dt_${dt}/Thermal'
    # execute_on = 'INITIAL TIMESTEP_END'
    interval = 20
  []
[]

[Postprocessors]
  [x_coord]
    type = FunctionValuePostprocessor
    function = heat_source_x
    outputs = 'console csv'
  []
  [y_coord]
    type = FunctionValuePostprocessor
    function = heat_source_y
    outputs = 'console csv'
  []
  [z_coord]
    type = FunctionValuePostprocessor
    function = heat_source_z
    outputs = 'console csv'
  []
  [power]
    type = FunctionValuePostprocessor
    function = effective_power
    outputs = 'console csv'
  []
[]

[VectorPostprocessors]
  [temperature]
    type = PointValueSampler
    variable = 'temp'
    points = '24.0 160.0 12.0
              24.0 140.0 12.0
              24.0 120.0 12.0
              24.0 100.0 12.0
              24.0 80.0 12.0
              24.0 60.0 12.0
              24.0 40.0 12.0
              24.0 160.0 14.0
              24.0 140.0 14.0
              24.0 120.0 14.0
              24.0 100.0 14.0
              24.0 80.0 14.0
              24.0 60.0 14.0
              24.0 40.0 14.0
              24.0 160.0 16.0
              24.0 140.0 16.0
              24.0 120.0 16.0
              24.0 100.0 16.0
              24.0 80.0 16.0
              24.0 60.0 16.0
              24.0 40.0 16.0
              24.0 160.0 18.0
              24.0 140.0 18.0
              24.0 120.0 18.0
              24.0 100.0 18.0
              24.0 80.0 18.0
              24.0 60.0 18.0
              24.0 40.0 18.0
              24.0 160.0 20.0
              24.0 140.0 20.0
              24.0 120.0 20.0
              24.0 100.0 20.0
              24.0 80.0 20.0
              24.0 60.0 20.0
              24.0 40.0 20.0
              24.0 160.0 22.0
              24.0 140.0 22.0
              24.0 120.0 22.0
              24.0 100.0 22.0
              24.0 80.0 22.0
              24.0 60.0 22.0
              24.0 40.0 22.0
              '
    sort_by = id
  []
[]
