T_melt = 1000
T_room = 300

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
    xmax =2
    ymin =0
    ymax =2
    zmin =0
    zmax =2
    nx=4
    ny=4
    nz=4
  [../]
  [./rename_bnd]
    type = RenameBoundaryGenerator
    input = gen
    # bottom = 0, front = 1, right = 2, back = 3, left = 4, top = 5
    old_boundary_id='0 1 2 3 4 5'
    new_boundary_name='bottom front right back left top'
  [../]
  [./add_set1]
    type = SubdomainBoundingBoxGenerator
    input = rename_bnd
    block_id = 3
    bottom_left = '0 0 0'
    top_right = '2 2 0.5'
  [../]
  [./add_set2]
    type = SubdomainBoundingBoxGenerator
    input = add_set1
    block_id = 1
    bottom_left = '0 0 0.5'
    top_right = '2 2 2'
  [../]
  [./add_set3]
    type = SubdomainBoundingBoxGenerator
    input = add_set2
    block_id = 2
    bottom_left = '0 0 0.5'
    top_right = '1.0 1.0 1.0'
  [../]
  [./product_top]
    type = SideSetsAroundSubdomainGenerator
    input = add_set3
    block = 2
    new_boundary = 'product_top'
    normal = '0 0 1'
  [../]
  [./moving_boundary]
    type = SideSetsAroundSubdomainGenerator
    input = product_top
    block = 2
    new_boundary = 'moving_boundary'
  []
  [./middle]
    type = SideSetsAroundSubdomainGenerator
    input = moving_boundary
    block = 3
    new_boundary = 'moving_boundary'
    normal = '0 0 1'
  []
  displacements='disp_x disp_y disp_z'
[]

[Variables]
  [./disp_x]
    block = '2 3'
  [../]
  [./disp_y]
    block = '2 3'
  [../]
  [./disp_z]
    block = '2 3'
  [../]
  [./temp]
    block = '2 3'
  [../]
[]

[ICs]
  [./temp_substrate]
    type = ConstantIC
    variable = temp
    value = ${T_room}
    block = '3'
  [../]
  [./temp_product]
    type = ConstantIC
    variable = temp
    value = ${T_melt}
    block = '2'
  [../]
[]

[Modules/TensorMechanics/Master]
  # [./all]
    strain = FINITE
    incremental = true
    add_variables = true
    generate_output = 'stress_xx stress_yy stress_zz stress_xy stress_yz stress_xz strain_yy strain_xx strain_zz strain_xy strain_xz strain_yz'
    use_automatic_differentiation = true
    # eigenstrain_names = 'thermal_eigenstrain'
    # block = '2 3'
  # [../]
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
  [./time]
    type = ADHeatConductionTimeDerivative
    variable = temp
    block = '2 3'
  [../]
  [./heat_conduct]
    type = ADHeatConduction
    variable = temp
    use_displaced_mesh = true
    thermal_conductivity = thermal_conductivity
    block = '2 3'
  [../]
  [./heatsource]
    type = ADMatHeatSource
    material_property = volumetric_heat
    variable = temp
    scalar = 1
    use_displaced_mesh = true
    block = '2 3'
  [../]
[]

[AuxVariables]
  # [./processor_id]
  #   order = CONSTANT
  #   family = MONOMIAL
  # [../]
  # [./von_mises]
  #   order = CONSTANT
  #   family = MONOMIAL
  #   block = '1'
  # [../]
  # [./plastic_strain_eff]
  #   order = CONSTANT
  #   family = MONOMIAL
  #   block = '1'
  # [../]
[]

[AuxKernels]
  # [./processor_id_aux]
  #   type = ProcessorIDAux
  #   variable = processor_id
  #   execute_on = timestep_begin
  # [../]
  # [./von_mises_kernel]
  #   type = ADRankTwoScalarAux
  #   variable = von_mises
  #   rank_two_tensor = stress
  #   execute_on = timestep_end
  #   scalar_type = VonMisesStress
  #   block = '1'
  # [../]
  # [./eff_plastic_strain_kernel]
  #   type = ADRankTwoScalarAux
  #   variable = plastic_strain_eff
  #   rank_two_tensor = plastic_strain
  #   execute_on = timestep_end
  #   scalar_type = EffectiveStrain
  #   block = '1'
  # [../]
[]

[Functions]
  [./heat_source_x]
    type = PiecewiseConstant
    x ='0   0.5   0.6  0.7  0.8  0.9 1.0'
    y ='0.5 1.25  1.25 1.75 1.75 1.25 1.25'
  [../]
  [./heat_source_y]
    type = PiecewiseConstant
    x ='0    0.5   0.6  0.7  0.8  0.9  1.0'
    y ='0.5  0.25  0.25 0.25 0.25 0.75 0.75'
  [../]
  [./heat_source_z]
    type = PiecewiseConstant
    x ='0 0.5 1.0'
    y ='0.5 0.75 0.75'
  [../]
  [./specific_heat]
    type = PiecewiseLinear
    data_file = ./input_params/Specific_Heat.csv
    format = columns
    scale_factor = 1.0
  [../]
  [./thermal_conductivity]
    type = PiecewiseLinear
    data_file = ./input_params/Thermal_Conductivity.csv
    format = columns
    scale_factor = 0.05e-3
  [../]

  [./pull_top]
    type = ParsedFunction
    value = 'if(t>0.5, 0.2, 0.4*t)'
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

  [./uz_top_pull]
    type = ADFunctionDirichletBC
    variable = disp_z
    boundary = 'product_top'
    function = pull_top
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


  [./convective]
    # type = ConvectiveFluxFunction # Convective flux, e.g. q'' = h*(Tw - Tf)
    type = ADConvectiveHeatFluxBC
    variable = temp
    boundary = 'moving_boundary'
    # coefficient = 2e-5
    heat_transfer_coefficient = 2e-5
    T_infinity = 300
  [../]
[]

[Materials]
  [./E]
    type = ADPiecewiseLinearInterpolationMaterial
    x = '294.994  1671.48  1721.77'
    y = '201.232e3 80.0821e3 6.16016e3' # [N/mm^3], i.e., MPa
    property = youngs_modulus
    variable = temp
    extrapolation = false
    block = '2 3'
  [../]
  [./nu]
    type = ADPiecewiseLinearInterpolationMaterial
    x = '294.994 1669.62 1721.77'
    y = '0.246407   0.36961  0.513347'
    property = poissons_ratio
    variable = temp
    extrapolation = false
    block = '2 3'
  [../]
  [./elasticity_tensor]
    type = ADComputeVariableIsotropicElasticityTensor
    youngs_modulus = youngs_modulus
    poissons_ratio = poissons_ratio
    block = '2 3'
  [../]

  # [./stress]
  #   type = ADComputeLinearElasticStress
  #   block = '1'
  # [../]

  [./radial_return_stress]
    type = ADComputeMultipleInelasticStress
    inelastic_models = 'rate_temp_plas'
    block = '2 3'
  [../]
  [./rate_temp_plas]
    type = ADRateTempDependentStressUpdate
    temperature = temp
    Y0 = 5.264e03 # [MPa]
    Rd1 = 8.565e-4 # [MPa]
    hxi = 1.670e-06 # [mm/(s*MPa)]
    Ex = '294.994  1671.48  1721.77'
    Ey = '201.232e3 80.0821e3 6.16016e3' # [N/mm^3], i.e., MPa
    nux = '294.994 1669.62 1721.77'
    nuy = '0.246407   0.36961  0.513347'
    # n2 = ${T_melt}
    absolute_tolerance = 1e-8
    block = '2 3'

    internal_solve_full_iteration_history = true
    internal_solve_output_on = "on_error"
  [../]

  [./thermal_expansion_strain_product]
    type = ADComputeThermalExpansionEigenstrain
    stress_free_temperature = ${T_melt}
    thermal_expansion_coeff = 0.5e-5
    temperature = temp
    eigenstrain_name = thermal_eigenstrain_product
    block = '2'
  [../]

  [./thermal_expansion_strain_substrate]
    type = ADComputeThermalExpansionEigenstrain
    stress_free_temperature = ${T_room}
    thermal_expansion_coeff = 0.5e-5
    temperature = temp
    eigenstrain_name = thermal_eigenstrain_substrate
    block = '3'
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
    density = 7.609e-6
    block = '2 3'
  [../]
  [./heat]
    type = ADHeatConductionMaterial
    specific_heat_temperature_function = specific_heat
    thermal_conductivity_temperature_function = thermal_conductivity
    temp = temp
    block = '2 3'
  [../]
[]

[UserObjects]
  [./activated_elem_uo]
    type = ActivateElementTemp
    execute_on = timestep_begin
    function_x= heat_source_x
    function_y= heat_source_y
    function_z= heat_source_z
    active_subdomain_id = 2
    inactive_subdomain_id = 3
    variable_activation = false
    # activate_distance = 1.4
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
  end_time = 1.0
  dt = 0.1
  dtmin = 1e-4
[]



[Outputs]
  file_base = 'output_4by4/out'
  [./exodus]
    type = Exodus
  [../]
  # [checkpoint]
  #   type = Checkpoint
  #   num_files = 2
  # []
[]
