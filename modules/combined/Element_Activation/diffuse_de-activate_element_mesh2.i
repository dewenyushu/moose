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
    zmax =1
    nx=2
    ny=2
    nz=1
  [../]
  [./subdomain_id]
    input = gen
    type = ElementSubdomainIDGenerator
    subdomain_ids = '1 2
                     1 2'
  [../]
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
    value= '1.5'
  [../]
  [./heat_source_y]
    type = ParsedFunction
    value= 'if (t<=0.1 , 0.5, -1.0)'
  [../]
  [./heat_source_z]
    type = ParsedFunction
    value= '0.25'
  [../]
  [./heat_source_y_deactivate]
    type = ParsedFunction
    value= 'if (t>0.1 , 0.5, -1.0)'
  [../]
[]

[BCs]
  [./temp_bottom_fix]
    type = ADDirichletBC
    variable = temp
    boundary = 'back'
    value = ${T_room}
  [../]

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

[UserObjects]
  [./activated_elem_uo]
    type = ActivateElementTemp
    execute_on = timestep_begin
    function_x= heat_source_x
    function_y= heat_source_y
    function_z= heat_source_z
    # activate_tol=1e-2
    active_subdomain_id = 1
  [../]
  [./de-ctivated_elem_uo]
    type = ActivateElementTemp
    execute_on = timestep_begin
    function_x= heat_source_x
    function_y= heat_source_y_deactivate
    function_z= heat_source_z
    # activate_tol=1e-2
    active_subdomain_id = 2
  [../]
[]

[AuxVariables]
  [./processor_id]
    order = CONSTANT
    family = MONOMIAL
  [../]
[]

[AuxKernels]
  [./processor_id_aux]
    type = ProcessorIDAux
    variable = processor_id
    execute_on = timestep_begin
  [../]
[]

[Outputs]
  file_base = 'temp_out'
  [./exodus]
    type = Exodus
  [../]
[]
