T_room = 1000

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
    nx=4
    ny=4
    nz=4
  [../]
  [./add_set1]
    type = SubdomainBoundingBoxGenerator
    input = gen
    block_id = 1
    bottom_left = '0 0 0'
    top_right = '2 2 0.5'
  [../]
  [./add_set2]
    type = SubdomainBoundingBoxGenerator
    input = add_set1
    block_id = 2
    bottom_left = '0 0 0.5'
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
  [../]
[]

[ICs]
  [./temp_substrate]
    type = ConstantIC
    variable = temp
    value = 293
    block = '2'
  [../]
  [./temp_product]
    type = ConstantIC
    variable = temp
    value = ${T_room}
    block = '1'
  [../]
[]

[AuxVariables]
  [./prop1]
    order = CONSTANT
    family = MONOMIAL
    block = '1'
  [../]
  [./ad_prop1]
    order = CONSTANT
    family = MONOMIAL
    block = '1'
  [../]
  [./density]
    order = CONSTANT
    family = MONOMIAL
    block = '1'
  [../]
  [./specific_heat]
    order = CONSTANT
    family = MONOMIAL
    block = '1'
  [../]
  [./thermal_conductivity]
    order = CONSTANT
    family = MONOMIAL
    block = '1'
  [../]
[]

[Kernels]
  # [./time]
  #   type = ADHeatConductionTimeDerivative
  #   variable = temp
  #   block = '1'
  # [../]
  # [./heat_conduct]
  #   type = ADHeatConduction
  #   variable = temp
  #   use_displaced_mesh = true
  #   thermal_conductivity = thermal_conductivity
  #   block = '1'
  # [../]
  # # [./heat_source]
  # #   type = HeatSource
  # #   function = 10
  # #   variable = temp
  # #   block = '1'
  # # [../]
  # [./heatsource]
  #   type = ADMatHeatSource
  #   material_property = volumetric_heat
  #   variable = temp
  #   scalar = 1
  #   use_displaced_mesh = true
  #   block = '1'
  # [../]
[]

[AuxKernels]
  [./prop1_output]
    type = MaterialRealAux
    variable = prop1
    property = test_mat_property
    block = '1'
  [../]
  [./ad_prop1_output]
    type = ADMaterialRealAux
    variable = ad_prop1
    property = ad_test_mat_property
    block = '1'
  [../]
  [./density_output]
    type = ADMaterialRealAux
    variable = density
    property = density
    block = '1'
  [../]
  [./sp_heat_output]
    type = ADMaterialRealAux
    variable = specific_heat
    property = specific_heat
    block = '1'
  [../]
  [./thermal_conductivity_output]
    type = ADMaterialRealAux
    variable = thermal_conductivity
    property = thermal_conductivity
    block = '1'
  [../]
[]

[Functions]
  # [./heat_source_x]
  #   type = PiecewiseLinear
  #   x ='0 0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9 1.0 1.1 1.2'
  #   y ='0.75 0.75 1.25 1.25 0.75 0.75 1.25 1.25 0.75 0.75 1.25 1.25 0.75'
  # [../]
  # [./heat_source_y]
  #   type = PiecewiseLinear
  #   x ='0 0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9 1.0 1.1 1.2'
  #   y ='0.75 0.75 0.75 1.25 1.25 0.75 0.75 1.25 1.25 0.75 0.75 1.25 1.25'
  # [../]
  # [./heat_source_z]
  #   type = PiecewiseLinear
  #   x ='0 0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9 1.0 1.1 1.2'
  #   y ='0.75 0.75 0.75 0.75 0.75 1.25 1.25 1.25 1.25 1.75 1.75 1.75 1.75'
  # [../]
  [./heat_source_x]
    type = PiecewiseLinear
    x ='0 0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9 1.0 1.1 1.2'
    y ='1 1 1.5 1.5 1 1 1.5 1.5 1 1 1.5 1.5 1'
  [../]
  [./heat_source_y]
    type = PiecewiseLinear
    x ='0 0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9 1.0 1.1 1.2'
    y ='1 1 1 1.5 1.5 1 1 1.5 1.5 1 1 1.5 1.5'
  [../]
  [./heat_source_z]
    type = PiecewiseLinear
    x ='0 0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9 1.0 1.1 1.2'
    y ='1 1 1 1 1 1.5 1.5 1.5 1.5 2 2 2 2'
  [../]
[]

[BCs]
  # [./temp_bottom_fix]
  #   type = ADDirichletBC
  #   variable = temp
  #   boundary = 'bottom'
  #   value = ${T_room}
  # [../]

  # [./temp_moving_fix]
  #   type = ADDirichletBC
  #   variable = temp
  #   boundary = 'moving_boundary'
  #   value = 500
  # [../]

  # [./convective]
  #   type = ADConvectiveHeatFluxBC # Convective flux, e.g. q'' = h*(Tw - Tf)
  #   variable = temp
  #   boundary = 'moving_boundary'
  #   heat_transfer_coefficient = 2e-5
  #   T_infinity = ${T_room}
  # [../]
[]

[Materials]
  # [./volumetric_heat]
  #   type = FunctionPathEllipsoidHeatSource
  #   a = 2
  #   b = 2
  #   c = 2
  #   power = 1000
  #   efficienty = 1.0
  #   factor = 2
  #   function_x= heat_source_x
  #   function_y= heat_source_y
  #   function_z= heat_source_z
  # [../]
  [./volumetric_heat]
    type = ADGenericFunctionMaterial
    prop_names = 'volumetric_heat'
    prop_values = 0.1
  [../]
  [./density]
    type = ADDensity
    density = 4.43e-6
    block = '1'
  [../]
  [./heat]
    type = ADHeatConductionMaterial
    specific_heat = 603
    thermal_conductivity = 10e-2
    block = '1'
  [../]
  [./stateful]
    type = StatefulTest
    prop_names = test_mat_property
    prop_values = 1.0
    block = '1'
  [../]
  [./ad_stateful]
    type = ADStatefulTest
    prop_names = ad_test_mat_property
    prop_values = 1.0
    block = '1'
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
  petsc_options_iname = '-pc_type -ksp_type'
  petsc_options_value = 'lu  preonly'

  line_search = 'none'

  l_max_its = 10
  nl_max_its = 20
  nl_rel_tol = 1e-4
  nl_abs_tol = 1e-12
  # l_tol = 1e-5

  start_time = 0.0
  end_time =0.4
  dt = 1e-1
  dtmin = 1e-4
[]

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

[Outputs]
  file_base = 'output/out_'
  [./exodus]
    type = Exodus
  [../]
  [checkpoint]
    type = Checkpoint
    num_files = 2
  []
[]
