T_room = 300
T_melt = 1000

[GlobalParams]
  displacements = 'disp_x disp_y disp_z'
[]

[Problem]
  kernel_coverage_check = false
[]

[Mesh]
  [./gen]
    type = GeneratedMeshGenerator
    dim = 3
    xmin =0
    xmax =4
    ymin =0
    ymax =4
    zmin =0
    zmax =1
    nx=4
    ny=4
    nz=1
  [../]
  [./add_set1]
    type = SubdomainBoundingBoxGenerator
    input = gen
    block_id = 1
    bottom_left = '0 0 0'
    top_right = '2 4 1'
  [../]
  [./add_set2]
    type = SubdomainBoundingBoxGenerator
    input = add_set1
    block_id = 2
    bottom_left = '2 0 0'
    top_right = '4 4 1'
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
    normal = '1 0 0'
  []
[]

[Variables]
  [./disp_x]
    block = '1'
  [../]
  [./disp_y]
    block = '1'
  [../]
  [./disp_z]
    block = '1'
  [../]
  [./temp]
    initial_condition = ${T_room}
    block = '1'
  [../]
[]

[Modules/TensorMechanics/Master]
    strain = FINITE
    incremental = true
    add_variables = true
    generate_output = 'stress_xx stress_yy stress_zz stress_xy stress_yz stress_xz strain_yy strain_xx strain_zz strain_xy strain_xz strain_yz'
    use_automatic_differentiation = true
  [product]
    block = '1'
    eigenstrain_names = 'thermal_eigenstrain_product'
    use_automatic_differentiation = true
  []
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

[AuxVariables]
[]

[AuxKernels]
[]

[Functions]
  [./heat_source_x]
    type = ParsedFunction
    value= '2.25'
  [../]
  [./heat_source_y]
    type = ParsedFunction
    value= '5.0*t'
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
    boundary = 'bottom'
    value = 0.0
  [../]
  [./uy_bottom_fix]
    type = ADDirichletBC
    variable = disp_y
    boundary = 'bottom'
    value = 0.0
  [../]
  [./uz_bottom_fix]
    type = ADDirichletBC
    variable = disp_z
    boundary = 'bottom'
    value = 0.0
  [../]
  [./temp_bottom_fix]
    type = ADDirichletBC
    variable = temp
    boundary = 'bottom'
    value = ${T_room}
  [../]

  # [./temp_moving_fix]
  #   type = ADDirichletBC
  #   variable = temp
  #   boundary = 'moving_boundary'
  #   value = ${T_room}
  # [../]


  # [./convective]
  #   # type = ConvectiveFluxFunction # Convective flux, e.g. q'' = h*(Tw - Tf)
  #   type = ADConvectiveHeatFluxBC
  #   variable = temp
  #   boundary = 'moving_boundary'
  #   # coefficient = 2e-5
  #   heat_transfer_coefficient = 2e-5
  #   T_infinity = 300
  # [../]
[]

[Materials]
  [./elasticity_tensor]
    type = ADComputeVariableIsotropicElasticityTensor
    youngs_modulus = 300
    poissons_ratio = 0.3
    block = '1'
  [../]

  [./stress]
    type = ADComputeLinearElasticStress
    block = '1'
  [../]

  [./thermal_expansion_strain_product]
    type = ADComputeThermalExpansionEigenstrain
    stress_free_temperature = ${T_melt}
    thermal_expansion_coeff = 0.5e-5
    temperature = temp
    eigenstrain_name = thermal_eigenstrain_product
    block = '1'
  [../]

  [./volumetric_heat]
    type = FunctionPathEllipsoidHeatSource
    a = 2
    b = 2
    c = 2
    power = 100
    efficienty = 1.0
    factor = 2
    function_x= heat_source_x
    function_y= heat_source_y
    function_z= heat_source_z
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
  [../]
[]

[UserObjects]
  [./activated_elem_uo]
    type = ActivateElementTemp
    execute_on = timestep_begin
    function_x= heat_source_x
    function_y= heat_source_y
    function_z= heat_source_z
    # activate_distance = 1.3
    active_subdomain_id = 1
    expand_boundary_name = 'moving_boundary'
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

  petsc_options_iname = '-ksp_type -pc_type -pc_factor_mat_solver_package'
  petsc_options_value = 'preonly lu       superlu_dist'

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



[Outputs]
  file_base = 'mech_output/out'
  [./exodus]
    type = Exodus
  [../]
  [checkpoint]
    type = Checkpoint
    num_files = 2
  []
[]
