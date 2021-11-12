T_room = 300
T_ambient = 300
T_melt = 1700

elevate_z = 0 # 200e-3 mm

# speed = 2 # mm/s
speed = 10.58e-3 # 10 mm/s = 10e-3 mm/ms
power = 300e-3 # 300W = kg*m^2/s^3 = 300e-3 kg*mm^2/ms^3
r = 300e-3 # 400 um = 400e-3 mm
dt = 1 # ms

[Mesh]
  [mesh]
    type = GeneratedMeshGenerator
    dim = 3
    xmin =-1
    xmax = 1
    ymin = -2.5
    ymax = 2.5
    zmin = 0
    zmax = 2
    nx = 40
    ny = 100
    nz = 40
  []
  [add_set1]
    type = SubdomainBoundingBoxGenerator
    input = mesh
    block_id = 3
    bottom_left = '-50 -50 0'
    top_right = '50 50 1'
  []
  [add_set2]
    type = SubdomainBoundingBoxGenerator
    input = add_set1
    block_id = 1
    bottom_left = '-50 -50 1'
    top_right = '50 50 2'
  []
  [add_set3]
    type = SubdomainBoundingBoxGenerator
    input = add_set2
    block_id = 2
    bottom_left = '-0.1 -2 1'
    top_right = '0 -1.9 1.1'
  []
  [add_set4]
    type = SubdomainBoundingBoxGenerator
    input = add_set3
    block_id = 4
    bottom_left = '-0.1 -2 0.9'
    top_right = '0 -1.9 1'
  []
  [moving_boundary]
    type = SideSetsAroundSubdomainGenerator
    input = add_set4
    block = 2
    new_boundary = 'moving_boundary'
  []
  [middle]
    type = SideSetsAroundSubdomainGenerator
    input = moving_boundary
    block = 3
    new_boundary = 'middle'
    normal = '0 0 1'
  []
[]

[Problem]
  kernel_coverage_check = false
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
    type = ConstantIC
    variable = temp
    value = ${T_melt}
    block = '2'
  []
[]

[AuxVariables]
  [processor_id]
    order = CONSTANT
    family = MONOMIAL
  []
  [power_aux]
    order = CONSTANT
    family = MONOMIAL
  []
  [speed_aux]
    order = CONSTANT
    family = MONOMIAL
  []
  [x_coord]
    order = FIRST
    family = LAGRANGE
  []
  [y_coord]
    order = FIRST
    family = LAGRANGE
  []
  [z_coord]
    order = FIRST
    family = LAGRANGE
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
    use_displaced_mesh = true
    thermal_conductivity = thermal_conductivity
  []
  [heatsource]
    type = ADMatHeatSource
    material_property = volumetric_heat
    variable = temp
    scalar = 1
    use_displaced_mesh = true
  []
[]

[AuxKernels]
  [processor_id_aux]
    type = ProcessorIDAux
    variable = processor_id
    execute_on = timestep_begin
  []
  [power]
    type = ConstantAux
    variable = power_aux
    value = ${power}
    # execute_on = 'INITIAL TIMESTEP_BEGIN'
  []
  [speed]
    type = ConstantAux
    variable = speed_aux
    value = ${speed}
    # execute_on = 'INITIAL TIMESTEP_BEGIN'
  []
  [x_coord]
    type = FunctionAux
    function = x_coord
    variable = x_coord
  []
  [y_coord]
    type = FunctionAux
    function = y_coord
    variable = y_coord
  []
  [z_coord]
    type = FunctionAux
    function = z_coord
    variable = z_coord
  []
[]

[BCs]
  [bottom_temp]
    type = ADDirichletBC
    variable = temp
    boundary = 1
    value = ${T_room}
  []
  [convective]
    # type = ConvectiveFluxFunction # Convective flux, e.g. q'' = h*(Tw - Tf)
    type = ADConvectiveHeatFluxBC
    variable = temp
    boundary = '2'
    # coefficient = 2e-5
    heat_transfer_coefficient = 2e-5 # W/m^2/K ->
    T_infinity = ${T_ambient}
  []
  [convective_air]
    # type = ConvectiveFluxFunction # Convective flux, e.g. q'' = h*(Tw - Tf)
    type = ADConvectiveHeatFluxBC
    variable = temp
    boundary = '4'
    # coefficient = 2e-5
    heat_transfer_coefficient = 2e-5 # W/m^2/K ->
    T_infinity = ${T_ambient}
  []
  # [./convective_middle]
  #   # type = ConvectiveFluxFunction # Convective flux, e.g. q'' = h*(Tw - Tf)
  #   type = ADConvectiveHeatFluxBC
  #   variable = temp
  #   boundary = 'middle'
  #   # coefficient = 2e-5
  #   heat_transfer_coefficient = 2e-5
  #   T_infinity = ${T_ambient}
  # [../]
[]

[Functions]
  [heat_source_x]
    type = ConstantFunction
    value = 0
  []
  [heat_source_y]
    type = ParsedFunction
    value = '-2 + ${speed}*t '
  []
  [heat_source_z]
    type = ParsedFunction
    value = '1+${elevate_z}'
  []
  [specific_heat_metal]
    type = PiecewiseLinear
    x = '-1e7 197.79  298.46  600.31 1401.01 1552.59 1701.44 1e7'
    y = '426.69 426.69 479.77 549.54 676.94  695.14 726.99 726.99'
    scale_factor = 1.0
  []
  [thermal_conductivity_metal]
    type = PiecewiseLinear
    x = '-1e7 198.84  298.10  398.75 500.76 601.40 700.64 801.27 901.89 1001.12 1098.98 1200.96 '
        '1301.56 1400.78 1501.37 1601.96 1e7'
    y = '247.72 247.72 285.64 323.55 358.44 390.29 417.59 446.41 469.16 491.91 510.11 528.31 540.44 '
        '554.09 561.67 569.26 569.26'
    format = columns
    scale_factor = 0.05e-6
  []
  [specific_heat_air]
    # type = PiecewiseLinear # make sure we do not have big jumps in the air-metal interface
    # x = '-1e-7 0 300 1701.44 1e7'
    # y = '1.008e3 1.008e3 1.008e3 726.99  726.99'
    # scale_factor = 1.0
    type = PiecewiseLinear
    data_file = AirSpecificHeat.csv
    format = columns
    scale_factor = 1.0

  []
  [thermal_conductivity_air]
    # type = PiecewiseLinear # make sure we do not have big jumps in the air-metal interface
    # x = '-1e7 0 300 1601.96 1e7'
    # y = '0.025e-6 0.025e-6 0.025e-6 28.463e-6 28.463e-6'
    # scale_factor = 1.0
    type = PiecewiseLinear
    data_file = AirConductivity.csv
    format = columns
    scale_factor = 1e-6
  []
  # for monitoring the deposited material geometry
  [scan_length_y]
    type = ParsedFunction
    value = '${speed}*t '
  []
  [x_coord]
    type = ParsedFunction
    value = 'x'
  []
  [y_coord]
    type = ParsedFunction
    value = 'y'
  []
  [z_coord]
    type = ParsedFunction
    value = 'z'
  []
[]

[Materials]
  [heat_metal]
    type = ADHeatConductionMaterial
    specific_heat_temperature_function = specific_heat_metal
    thermal_conductivity_temperature_function = thermal_conductivity_metal
    temp = temp
    block = '2 3 4'
  []
  [heat_air]
    type = ADHeatConductionMaterial
    specific_heat_temperature_function = specific_heat_air #1.008e3 #1.008KJ/Kg*K -> same for J
    thermal_conductivity_temperature_function = thermal_conductivity_air # 0.025e-6
    temp = temp
    block = '1'
  []
  [volumetric_heat_air]
    type = FunctionPathEllipsoidHeatSource
    r = ${r}
    power = ${power}
    efficiency = 0.36
    factor = 0.8e-3
    function_x = heat_source_x
    function_y = heat_source_y
    function_z = heat_source_z
    heat_source_type = 'line'
    threshold_length = 0.1 #mm
    number_time_integration = 10
    block = '1'
  []
  [volumetric_heat_metal]
    type = FunctionPathEllipsoidHeatSource
    r = ${r}
    power = ${power}
    efficiency = 0.36
    factor = 0.9
    function_x = heat_source_x
    function_y = heat_source_y
    function_z = heat_source_z
    heat_source_type = 'line'
    threshold_length = 0.1 #mm
    number_time_integration = 10
    block = '2 3 4'
  []
  [density_metal]
    type = ADDensity
    density = 7609e-9 # kg/m^3 -> 1e-9 kg/mm^3
    block = '2 3 4'
  []
  [density_air]
    type = ADDensity
    density = 1.1644e-9 # 1.1644 kg/m^3 -> 1.1644e-9 kg/mm^3
    block = '1'
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
  [activated_elem_uo_melt]
    type = CoupledVarThresholdElementSubdomainModifier
    execute_on = 'TIMESTEP_BEGIN'
    coupled_var = temp
    block = 3
    subdomain_id = 4
    criterion_type = ABOVE
    threshold = ${T_melt}
    # moving_boundary_name = 'moving_boundary'
  []
  # [activated_elem_uo]
  #   type = ActivateElementsCoupled
  #   coupled_var = temp
  #   activate_type = above
  #   active_subdomain_id = '2'
  #   expand_boundary_name= 'moving_boundary'
  #   activate_value= ${T_melt}
  #   execute_on = 'TIMESTEP_END'
  # []
[]

# [Adaptivity]
#   steps = 1
#   marker = marker
#   initial_marker = marker
#   max_h_level = 2
#   [Indicators/indicator]
#     type = GradientJumpIndicator
#     variable = temp
#   []
#   [Markers/marker]
#     type = ErrorFractionMarker
#     indicator = indicator
#     coarsen = 0
#     refine = 0.5
#   []
# []

[Preconditioning]
  [smp]
    type = SMP
    full = true
  []
[]

[Executioner]
  type = Transient
  solve_type = 'PJFNK'

  automatic_scaling = true

  petsc_options_iname = '-ksp_type -pc_type -pc_factor_mat_solver_package'
  petsc_options_value = 'preonly lu       superlu_dist'

  # petsc_options_iname = '-pc_type -ksp_type -pc_factor_shift_type -pc_factor_shift_amount'
  # petsc_options_value = 'lu  preonly NONZERO 1e-10'

  line_search = 'none'

  l_max_its = 100
  nl_max_its = 40
  nl_rel_tol = 1e-5
  nl_abs_tol = 1e-6

  start_time = 0.0
  end_time = '${fparse 4/speed}'
  dt = ${dt} # ms
  dtmin = 1e-6

  auto_advance = true # cut time-step when subapp fails

  error_on_dtmin = false
[]

[Outputs]
  file_base = 'output_laser_factor/Line_master_speed_${speed}_power_${power}_r_${r}'
  csv = true
  [exodus]
    type = Exodus
    file_base = 'output_laser_factor/Exodus_speed_${speed}_power_${power}_r_${r}/Line_thermal_melt'
    # execute_on = 'INITIAL TIMESTEP_END'
    interval = 1
  []
[]

[Postprocessors]
  [bead_max_temperature]
    type = ElementExtremeValue
    variable = temp
    value_type = max
    block = '2'
    outputs = 'csv'
  []
  [bead_min_temperature]
    type = ElementExtremeValue
    variable = temp
    value_type = min
    block = '2'
    outputs = 'csv'
  []
  [bead_volume]
    type = VolumePostprocessor
    block = '2'
    use_displaced_mesh = true
    outputs = 'csv console'
  []
  [melt_volume]
    type = VolumePostprocessor
    block = '4'
    use_displaced_mesh = true
    outputs = 'csv console'
  []
  [pp_power]
    type = ElementAverageValue
    variable = power_aux
    outputs = 'csv'
  []
  [pp_speed]
    type = ElementAverageValue
    variable = speed_aux
    outputs = 'csv'
  []
  [bead_x_coord_max]
    type = NodalExtremeValue
    variable = x_coord
    value_type = max
    block = '2'
    outputs = 'csv console'
  []
  [bead_x_coord_min]
    type = NodalExtremeValue
    variable = x_coord
    value_type = min
    block = '2'
    outputs = 'csv'
  []
  [bead_y_coord_max]
    type = NodalExtremeValue
    variable = y_coord
    value_type = max
    block = '2'
    outputs = 'csv console'
  []
  [bead_y_coord_min]
    type = NodalExtremeValue
    variable = y_coord
    value_type = min
    block = '2'
    outputs = 'csv console'
  []
  [bead_z_coord_max]
    type = NodalExtremeValue
    variable = z_coord
    value_type = max
    block = '2'
    outputs = 'csv console'
  []
  [bead_z_coord_min]
    type = NodalExtremeValue
    variable = z_coord
    value_type = min
    block = '2'
    outputs = 'csv'
  []
  [melt_x_coord_max]
    type = NodalExtremeValue
    variable = x_coord
    value_type = max
    block = '4'
    outputs = 'csv console'
  []
  [melt_x_coord_min]
    type = NodalExtremeValue
    variable = x_coord
    value_type = min
    block = '4'
    outputs = 'csv'
  []
  [melt_y_coord_max]
    type = NodalExtremeValue
    variable = y_coord
    value_type = max
    block = '4'
    outputs = 'csv console'
  []
  [melt_y_coord_min]
    type = NodalExtremeValue
    variable = y_coord
    value_type = min
    block = '4'
    outputs = 'csv console'
  []
  [melt_z_coord_max]
    type = NodalExtremeValue
    variable = z_coord
    value_type = max
    block = '4'
    outputs = 'csv'
  []
  [melt_z_coord_min]
    type = NodalExtremeValue
    variable = z_coord
    value_type = min
    block = '4'
    outputs = 'csv console'
  []
[]
