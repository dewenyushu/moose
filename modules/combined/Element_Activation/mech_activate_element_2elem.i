T_melt = 600
T_room = 300

# [GlobalParams]
#   displacements = 'disp_x disp_y disp_z'
# []

[Problem]
  kernel_coverage_check = false
  # material_coverage_check = false
[]


[Mesh]
  [./mesh]
    type = GeneratedMeshGenerator
    dim = 3
    xmin =0
    xmax =2
    ymin =0
    ymax =1
    zmin =0
    zmax =1
    nx=2
    ny=1
    nz=1
  [../]
  [./add_set1]
    type = SubdomainBoundingBoxGenerator
    input = mesh
    block_id = 1
    bottom_left = '0 0 0'
    top_right = '1 1 1'
  [../]
  [./add_set2]
    type = SubdomainBoundingBoxGenerator
    input = add_set1
    block_id = 2
    bottom_left = '1 0 0'
    top_right = '2 1 1'
  [../]
  displacements='disp_x disp_y disp_z'
[]

[Variables]
  [./disp_x]
    block = 1
  [../]
  [./disp_y]
    block = 1
  [../]
  [./disp_z]
    block = 1
  [../]
  [./temp]
    initial_condition = ${T_room}
    block = 1
  [../]
[]

[Kernels]
  [./time]
    type = ADHeatConductionTimeDerivative
    variable = temp
    block =1
  [../]
  [./heat_conduct]
    type = ADHeatConduction
    variable = temp
    use_displaced_mesh = true
    thermal_conductivity = thermal_conductivity
    block =1
  [../]
  [./heatsource]
    type = ADMatHeatSource
    material_property = volumetric_heat
    variable = temp
    scalar = 1
    use_displaced_mesh = true
    block =1
  [../]
  [./disp_x]
    type = ADStressDivergenceTensors
    component = 0
    variable = disp_x
    use_displaced_mesh = true
    block = 1
    displacements = 'disp_x disp_y disp_z'
  [../]
  [./disp_y]
    type = ADStressDivergenceTensors
    component = 1
    variable = disp_y
    use_displaced_mesh = true
    block = 1
    displacements = 'disp_x disp_y disp_z'
  [../]
  [./disp_z]
    type = ADStressDivergenceTensors
    component = 2
    variable = disp_z
    use_displaced_mesh = true
    block = 1
    displacements = 'disp_x disp_y disp_z'
  [../]
[]


[Functions]
  [./heat_source_x]
    type = ParsedFunction
    value= '0.5 + 7.0*t'
  [../]
  [./heat_source_y]
    type = ParsedFunction
    value= '0.5'
  [../]
  [./heat_source_z]
    type = ParsedFunction
    value= '0.5'
  [../]
[]

[BCs]
  [./ux_bottom_fix]
    type = ADDirichletBC
    variable = disp_x
    boundary = 'back'
    value = 0.0
  [../]
  [./uy_bottom_fix]
    type = ADDirichletBC
    variable = disp_y
    boundary = 'back'
    value = 0.0
  [../]
  [./uz_bottom_fix]
    type = ADDirichletBC
    variable = disp_z
    boundary = 'back'
    value = 0.0
  [../]
  [./temp_bottom_fix]
    type = ADDirichletBC
    variable = temp
    boundary = 0
    value = ${T_room}
  [../]
[]

[Materials]
  [./elasticity_tensor]
    type = ADComputeVariableIsotropicElasticityTensor
    youngs_modulus = 2e5
    poissons_ratio = 0.3
    block = 1
  [../]
  [./strain]
    # type = ADComputeFiniteStrain
    type = ADComputeSmallStrain
    eigenstrain_names = thermal_eigenstrain
    block = 1
    displacements = 'disp_x disp_y disp_z'
  [../]
  [./stress]
    type = ADComputeLinearElasticStress
    block = 1
  [../]
  [./thermal_expansion_strain]
    type = ADComputeThermalExpansionEigenstrain
    stress_free_temperature = 293
    thermal_expansion_coeff = 1e-6
    temperature = temp
    eigenstrain_name = thermal_eigenstrain
    melt_temperature = ${T_melt}
    block = 1
  [../]
  [./volumetric_heat]
    type = FunctionPathEllipsoidHeatSource
    a = 2
    b = 2
    c = 2
    power = 10
    efficienty = 1.0
    factor = 2
    function_x= heat_source_x
    function_y= heat_source_y
    function_z= heat_source_z
  [../]
  [./density]
    type = ADDensity
    density = 7.609e-6
    displacements='disp_x disp_y disp_z'
    block = 1
  [../]
  [./heat]
    type = ADHeatConductionMaterial
    specific_heat = 600
    thermal_conductivity = 1e-2
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
