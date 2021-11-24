T_room = 300
T_ambient = 300
T_melt = 1000

speed = 10e-3
power = 300e-3
r = 200e-3
dt = 1

[Mesh]
  [mesh]
    type = GeneratedMeshGenerator
    dim = 3
    xmax = 1
    ymax = 1
    zmax = 1
    nx = 10
    ny = 10
    nz = 10
  []
  [add_set1]
    type = SubdomainBoundingBoxGenerator
    input = mesh
    block_id = 3
    bottom_left = '0 0 0'
    top_right = '1 5 0.5'
  []
  [add_set2]
    type = SubdomainBoundingBoxGenerator
    input = add_set1
    block_id = 1
    bottom_left = '0 0 0.5'
    top_right = '1 5 1'
  []
  [add_set3]
    type = SubdomainBoundingBoxGenerator
    input = add_set2
    block_id = 2
    bottom_left = '0.5 0.0 0.5'
    top_right = '0.6 0.1 0.6'
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
  # displacements = 'disp_x disp_y disp_z'
[]

[Problem]
  kernel_coverage_check = false
  # material_coverage_check = false
[]

[Variables]
  [temp]
    block = '1 2 3'
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
    type = ConstantIC
    variable = temp
    value = ${T_melt}
    block = '2'
  []
[]

[Kernels]
  [time]
    type = ADHeatConductionTimeDerivative
    variable = temp
  []
  [heat_conduct_metal]
    type = ADHeatConduction
    variable = temp
    use_displaced_mesh = true
    thermal_conductivity = thermal_conductivity
    block = '2 3'
  []
  [heat_conduct_air]
    type = ADHeatConduction
    variable = temp
    use_displaced_mesh = true
    thermal_conductivity = thermal_conductivity # 0.025 W/(m·K) 1e-6 for mm and ms
    block = '1'
  []
  [heatsource]
    type = ADMatHeatSource
    material_property = volumetric_heat
    variable = temp
    scalar = 1
    use_displaced_mesh = true
  []
[]

[BCs]
  [bottom_temp]
    type = ADDirichletBC
    variable = temp
    boundary = back
    value = ${T_room}
  []
  [convective]
    type = ADConvectiveHeatFluxBC
    variable = temp
    boundary = 'left right top bottom front'
    heat_transfer_coefficient = 2e-5
    T_infinity = ${T_ambient}
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
  [specific_heat_metal]
    type = PiecewiseLinear
    x='197.79 298.46 600.31 1401.01 1552.59 1701.44 '
    y='426.69 479.77  549.54 676.94 695.14 726.99'
    format = columns
    scale_factor = 1.0
  []
  [thermal_conductivity_metal]
    type = PiecewiseLinear
    x = '198 298 500 801 1001 1400 1601'
    y = '247.72 285.64 358.44 446.41 491.91 554.09 569.26'
    format = columns
    scale_factor = 0.05e-6
  []
[]

[Materials]
  [heat_metal]
    type = ADHeatConductionMaterial
    specific_heat_temperature_function = specific_heat_metal
    thermal_conductivity_temperature_function = thermal_conductivity_metal
    temp = temp
    block = '1 2 3'
  []
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
  [density_metal]
    type = ADDensity
    density = 7609e-9 # kg/mm^3
    block = '1 2 3'
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

  petsc_options_iname = '-ksp_type -pc_type -pc_factor_mat_solver_package'
  petsc_options_value = 'preonly lu       superlu_dist'

  line_search = 'none'

  l_max_its = 100
  nl_max_its = 20
  nl_rel_tol = 1e-5
  nl_abs_tol = 1e-6

  start_time = 0.0
  end_time = '${fparse 4/speed}'
  dt = ${dt}
  num_steps = 40
  dtmin = 1e-6
[]

[UserObjects]
  [activated_elem_uo]
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

[Outputs]
  file_base = 'output/thermal_activation_amr/out'
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
