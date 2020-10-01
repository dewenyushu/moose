T_room = 300

[Problem]
  kernel_coverage_check = false
[]

[Mesh]
  [./gen]
    type = GeneratedMeshGenerator
    dim = 3
    xmin =0
    xmax =2
    ymin =0
    ymax =2
    zmin =0
    zmax =2
    nx=2
    ny=2
    nz=2
  [../]
  [./add_set1]
    type = SubdomainBoundingBoxGenerator
    input = gen
    block_id = 1
    bottom_left = '0 0 0'
    top_right = '2 2 1'
  [../]
  [./add_set2]
    type = SubdomainBoundingBoxGenerator
    input = add_set1
    block_id = 2
    bottom_left = '0 0 1'
    top_right = '2 2 2'
  [../]
  [./rename_bnd]
    type = RenameBoundaryGenerator
    input = add_set2
    # bottom = 0, front = 1, right = 2, back = 3, left = 4, top = 5
    old_boundary_id='0 1 2 3 4 5'
    new_boundary_name='bottom front right back left top'
  [../]
  [./middle]
    type = SideSetsAroundSubdomainGenerator
    input = rename_bnd
    block = 1
    new_boundary = 'moving_boundary'
    normal = '0 0 1'
  []
[]

[Variables]
  [./temp]
    initial_condition = ${T_room}
    block = '1'
  [../]
[]

[Kernels]
  [./time]
    type = ADHeatConductionTimeDerivative
    variable = temp
    block = '1'
  [../]
  [./heat_conduct]
    type = ADHeatConduction
    variable = temp
    use_displaced_mesh = true
    thermal_conductivity = thermal_conductivity
    block = '1'
  [../]
  [./heatsource]
    type = ADMatHeatSource
    material_property = volumetric_heat
    variable = temp
    scalar = 1
    use_displaced_mesh = true
    block = '1'
  [../]

[]

[Functions]
  [./heat_source_x]
    type = ParsedFunction
    value= '0.5'
  [../]
  [./heat_source_y]
    type = ParsedFunction
    value= '7.5*t'
  [../]
  [./heat_source_z]
    type = ParsedFunction
    value= '1.5'
  [../]
[]

[BCs]
  [./temp_bottom_fix]
    type = ADDirichletBC
    variable = temp
    boundary = 'bottom'
    value = ${T_room}
  [../]

  [./temp_moving_fix]
    type = ADDirichletBC
    variable = temp
    boundary = 'moving_boundary'
    value = 500
  [../]
  # [./convective]
  #   type = ConvectiveFluxFunction # Convective flux, e.g. q'' = h*(Tw - Tf)
  #   variable = temp
  #   boundary = 'moving_boundary'
  #   coefficient = 2e-5
  #   T_infinity = ${T_room}
  # [../]
[]

[Materials]
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
    density = 4.43e-6
  [../]
  [./heat]
    type = ADHeatConductionMaterial
    specific_heat = 603
    thermal_conductivity = 10e-2
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
  end_time = 0.2
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
    # activate_tol=1e-2
    active_subdomain_id = 1
    expand_boundary_name = 'moving_boundary'
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
