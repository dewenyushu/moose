T_melt = 600
T_room = 300

[Problem]
  kernel_coverage_check = false
[]


[Mesh]
  [./mesh]
    type = GeneratedMeshGenerator
    dim = 3
    xmin =0
    xmax =10
    ymin =0
    ymax =10
    zmin =0
    zmax =0.5
    nx=20
    ny=20
    nz=1
  [../]
  [./add_set1]
    type = SubdomainBoundingBoxGenerator
    input = mesh
    block_id = 1
    bottom_left = '0 0 0'
    top_right = '5 10 0.5'
  [../]
  [./add_set2]
    type = SubdomainBoundingBoxGenerator
    input = add_set1
    block_id = 2
    bottom_left = '5 0 0'
    top_right = '10 10 0.5'
  [../]

[]

[Variables]
  [./temp]
    initial_condition = ${T_room}
    block = 1
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

[Functions]
  [./heat_source_x]
    type = ParsedFunction
    value= '5.25'
  [../]
  [./heat_source_y]
    type = ParsedFunction
    value= '2.5*t'
  [../]
  [./heat_source_z]
    type = ParsedFunction
    value= '0.25'
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

[BCs]
  [./temp_bottom_fix]
    type = ADDirichletBC
    variable = temp
    boundary = 0
    value = ${T_room}
  [../]
[]

[Materials]
  [./thermal_expansion_strain]
    type = ADComputeThermalExpansionEigenstrain
    stress_free_temperature = 293
    thermal_expansion_coeff = 0.5e-5
    temperature = temp
    eigenstrain_name = thermal_eigenstrain
    # activated_elem_aux = activated_elem
    melt_temperature = ${T_melt}
    block = 1
  [../]
  [./volumetric_heat]
    type = FunctionPathEllipsoidHeatSource
    a = 2
    b = 2
    c = 2
    power = 1000
    efficienty = 1.0
    factor = 2
    function_x= heat_source_x
    function_y= heat_source_y
    function_z= heat_source_z
  [../]
  [./density]
    type = ADDensity
    density = 7.609e-6
    block = 1
  [../]
  [./heat]
    type = ADHeatConductionMaterial
    specific_heat_temperature_function = specific_heat
    thermal_conductivity_temperature_function = thermal_conductivity
    temp = temp
    block = 1
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

  # petsc_options = '-snes_ksp'
  petsc_options_iname = '-pc_type -ksp_type'
  petsc_options_value = 'lu  preonly'

  line_search = 'none'

  l_max_its = 10
  nl_max_its = 20
  nl_rel_tol = 1e-4
  nl_abs_tol = 1e-12
  # l_tol = 1e-5

  start_time = 0.0
  end_time = 2.0
  dt = 1e-1
  dtmin = 1e-4
[]

# [AuxVariables]
#   [./activated_elem]
#     order = CONSTANT
#     family = MONOMIAL
#   [../]
# []

# [AuxKernels]
#   [./activated_elem]
#     type = ActivatedElementsMarker
#     melt_temperature = 480
#     temp_aux = temp
#     variable = activated_elem
#     execute_on = timestep_begin
#   [../]
# []

[UserObjects]
  [./activated_elem_uo]
    type = ActivateElementTemp
    execute_on = timestep_begin
    function_x= heat_source_x
    function_y= heat_source_y
    function_z= heat_source_z
    active_subdomain_id = 1
  [../]
[]

# [Postprocessors]
#   [./temperature]
#     type = ElementAverageValue
#     variable = temp
#   [../]
# []

[Outputs]
  file_base = 'temp_out'
  [./exodus]
    type = Exodus
  [../]
[]
