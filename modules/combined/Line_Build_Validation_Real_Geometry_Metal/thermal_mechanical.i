T_room = 300
T_ambient = 300
T_melt = 1700

elevate_z = 0 # 200e-3 mm

# speed = 2 # mm/s
speed = 10.58e-3 # 10 mm/s = 10e-3 mm/ms
power = 300e-3 # 300W = kg*m^2/s^3 = 300e-3 kg*mm^2/ms^3
r = 300e-3 # 400 um = 400e-3 mm
dt = 2 #'${fparse 0.3*r/speed}' # ms

refine = 0

[GlobalParams]
  displacements = 'disp_x disp_y disp_z'
  volumetric_locking_correction = true
[]

[Problem]
  kernel_coverage_check = false
  material_coverage_check = false
[]

[Mesh]
  [mesh]
    type = GeneratedMeshGenerator
    dim = 3
    xmin = -1
    xmax = 1
    ymin = -2.5
    ymax = 2.5
    zmin = 0
    zmax = 2
    nx = 20
    ny = 50
    nz = 20
  []
  [add_set1]
    type = SubdomainBoundingBoxGenerator
    input = mesh
    block_id = 3
    bottom_left = '-50 -50 0'
    top_right = '50 50 1.0'
  []
  [add_set2]
    type = SubdomainBoundingBoxGenerator
    input = add_set1
    block_id = 1
    bottom_left = '-50 -50 1.0'
    top_right = '50 50 2'
  []
  [add_set3]
    type = SubdomainBoundingBoxGenerator
    input = add_set2
    block_id = 2
    bottom_left = '-0.1 -2 1.0'
    top_right = '0 -1.9 1.1'
  []
  [moving_boundary]
    type = SideSetsAroundSubdomainGenerator
    input = add_set3
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
  displacements = 'disp_x disp_y disp_z'

  uniform_refine = ${refine}
[]

[Variables]
  [temp]
    block = '1 2 3'
  []
  [disp_x]
    block = '2 3'
  []
  [disp_y]
    block = '2 3'
  []
  [disp_z]
    block = '2 3'
  []
[]

[ICs]
  [temp_substrate]
    type = ConstantIC
    variable = temp
    value = ${T_room}
    block = '1 3'
  []
  [temp_product]
    type = FunctionIC
    variable = temp
    function = temp_ic
    block = '2'
  []
[]

[AuxVariables]
  [von_mises]
    order = CONSTANT
    family = MONOMIAL
    block = '2 3'
  []
  [plastic_strain_eff]
    order = CONSTANT
    family = MONOMIAL
    block = '2 3'
  []
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

[Modules/TensorMechanics/Master]
  strain = FINITE
  incremental = true
  add_variables = true
  generate_output = 'stress_xx stress_yy stress_zz stress_xy stress_yz stress_xz strain_yy strain_xx '
                    'strain_zz strain_xy strain_xz strain_yz'
  use_automatic_differentiation = true
  [product]
    block = '2'
    eigenstrain_names = 'thermal_eigenstrain_product'
    use_automatic_differentiation = true
  []
  [substrate]
    block = '3'
    eigenstrain_names = 'thermal_eigenstrain_substrate'
    use_automatic_differentiation = true
  []
[]

[Kernels]
  [time]
    type = ADHeatConductionTimeDerivative
    variable = temp
    use_displaced_mesh = false
  []
  [heat_conduc]
    type = ADHeatConduction
    variable = temp
    use_displaced_mesh = false
    thermal_conductivity = thermal_conductivity
  []
  [heatsource]
    type = ADMatHeatSource
    material_property = volumetric_heat
    variable = temp
    scalar = 1
    use_displaced_mesh = false
  []
[]

[AuxKernels]
  [processor_id_aux]
    type = ProcessorIDAux
    variable = processor_id
    execute_on = timestep_begin
  []
  [von_mises_kernel]
    type = ADRankTwoScalarAux
    variable = von_mises
    rank_two_tensor = stress
    execute_on = timestep_end
    scalar_type = VonMisesStress
    block = '2 3'
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
    boundary = 'back'
    value = ${T_room}
  []
  [convective]
    # type = ConvectiveFluxFunction # Convective flux, e.g. q'' = h*(Tw - Tf)
    type = ADConvectiveHeatFluxBC
    variable = temp
    boundary = 'bottom front left right top'
    # coefficient = 2e-5
    heat_transfer_coefficient = 2e-5 # W/m^2/K ->
    T_infinity = ${T_ambient}
  []
  [ux_bottom_fix]
    type = ADDirichletBC
    variable = disp_x
    boundary = 'back'
    value = 0.0
  []
  [uy_bottom_fix]
    type = ADDirichletBC
    variable = disp_y
    boundary = 'back'
    value = 0.0
  []
  [uz_bottom_fix]
    type = ADDirichletBC
    variable = disp_z
    boundary = 'back'
    value = 0.0
  []
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
  [temp_ic]
    type = ParsedFunction
    value = 'if(t<=0, temp_room, temp_melt)'
    vars = 'temp_room temp_melt'
    vals = '${T_room} ${T_melt}'
  []
[]

[Materials]
  [E]
    type = ADPiecewiseLinearInterpolationMaterial
    x = '0 294.994  1671.48  1721.77 1e7'
    y = '201.232e3 201.232e3 80.0821e3 6.16016e3 6.16016e3' #MPa # 10^9 Pa = 10^9 kg/m/s^2 = kg/mm/ms^2
    # y = '6.16016e3 6.16016e3 6.16016e3 6.16016e3 6.16016e3'
    property = youngs_modulus
    variable = temp
    extrapolation = false
    block = '2 3'
  []
  [nu]
    type = ADPiecewiseLinearInterpolationMaterial
    x = '0 294.994 1669.62 1721.77 1e7'
    y = '0.246407 0.246407   0.36961  0.36961 0.36961' #''0.513347 0.513347'
    property = poissons_ratio
    variable = temp
    extrapolation = false
    block = '2 3'
  []
  [elasticity_tensor]
    type = ADComputeVariableIsotropicElasticityTensor
    youngs_modulus = youngs_modulus
    poissons_ratio = poissons_ratio
    block = '2 3'
  []
  [thermal_expansion_strain_product]
    type = ADComputeThermalExpansionEigenstrain
    stress_free_temperature = ${T_melt}
    # thermal_expansion_coeff = 1.72e-5
    thermal_expansion_coeff = 6.72e-6 #1.72e-5 /K
    temperature = temp
    eigenstrain_name = thermal_eigenstrain_product
    block = '2'
  []
  [thermal_expansion_strain_substrate]
    type = ADComputeThermalExpansionEigenstrain
    stress_free_temperature = ${T_room}
    # thermal_expansion_coeff = 1.72e-5
    thermal_expansion_coeff = 6.72e-6 #1.72e-5 /K
    temperature = temp
    eigenstrain_name = thermal_eigenstrain_substrate
    block = '3'
  []

  [stress]
    type = ADComputeFiniteStrainElasticStress
    block = '2 3'
  []

  [heat_metal]
    type = ADHeatConductionMaterial
    specific_heat_temperature_function = specific_heat_metal
    thermal_conductivity_temperature_function = thermal_conductivity_metal
    temp = temp
    block = '1 2 3'
  []
  [volumetric_heat_metal]
    type = FunctionPathEllipsoidHeatSource
    r = ${r}
    power = ${power}
    efficiency = 0.36
    factor = 1.6
    function_x = heat_source_x
    function_y = heat_source_y
    function_z = heat_source_z
    heat_source_type = 'line'
    threshold_length = 0.1 #mm
    number_time_integration = 10
    block = '1 2 3'
  []
  # [density_metal]
  #   type = ADDensity
  #   density = 7609e-9 # kg/m^3 -> 1e-9 kg/mm^3
  #   block = '1 2 3'
  # []
  [density_metal]
    type = ADGenericConstantMaterial
    prop_names = 'density'
    prop_values = 7609e-9 # kg/m^3 -> 1e-9 kg/mm^3
    block = '1 2 3'
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
  max_h_level = 1
  [Indicators/indicator]
    type = GradientJumpIndicator
    variable = temp
  []
  [Markers/marker]
    type = ErrorFractionMarker
    indicator = indicator
    coarsen = 0.1
    refine = 0.5
    check_subdomain_consistent_for_coarsen = true
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

  automatic_scaling = true

  petsc_options_iname = '-ksp_type -pc_type -pc_factor_mat_solver_package -pc_factor_shift_type '
                        '-pc_factor_shift_amount'
  petsc_options_value = 'preonly lu       superlu_dist NONZERO 1e-10'

  line_search = 'none'

  l_max_its = 100
  nl_max_its = 15
  nl_rel_tol = 1e-8
  nl_abs_tol = 1e-10

  start_time = 0.0
  end_time = '${fparse 3/speed}'
  dt = ${dt} # ms
  dtmin = 1e-6

  auto_advance = true # cut time-step when subapp fails

  error_on_dtmin = false
[]

[Outputs]
  file_base = 'output_thermal_mech/thermal_mech'
  csv = true
  [exodus]
    type = Exodus
    file_base = 'output_thermal_mech/thermal_mech/thermal_mech'
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
    # use_displaced_mesh = true
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
  [max_von_mises_stress]
    type = ElementExtremeValue
    variable = von_mises
    value_type = max
    block = '2'
  []
  [min_von_mises_stress]
    type = ElementExtremeValue
    variable = von_mises
    value_type = min
    block = '2'
  []
[]
