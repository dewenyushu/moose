T_room = 300
T_melt = 1000
T_ambient = 300

[GlobalParams]
  displacements = 'disp_x disp_y disp_z'
[]

[Problem]
  kernel_coverage_check = false
[]

[Mesh]
  [./mesh]
    type = FileMeshGenerator
    file = ./input_mesh/solid_cylinder_17_substrate_100x100.e
  [../]
  [./add_set1]
    type = SubdomainBoundingBoxGenerator
    input = mesh
    block_id = 3
    bottom_left = '-50 -50 0'
    top_right = '50 50 5'
  [../]
  [./add_set2]
    type = SubdomainBoundingBoxGenerator
    input = add_set1
    block_id = 1
    bottom_left = '-50 -50 5'
    top_right = '50 50 25'
  [../]
  [./add_set3]
    type = SubdomainBoundingBoxGenerator
    input = add_set2
    block_id = 2
    bottom_left = '-2 -2 5'
    top_right = '2 2 6'
  [../]
  [./moving_boundary]
    type = SideSetsAroundSubdomainGenerator
    input = add_set3
    block = 2
    new_boundary = 'moving_boundary'
  []
  [./middle]
    type = SideSetsAroundSubdomainGenerator
    input = moving_boundary
    block = 3
    new_boundary = 'middle'
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

# [AuxVariables]
# []

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

# [AuxKernels]
#   [./ad_prop1_output]
#     type = ADMaterialRealAux
#     variable = ad_prop1
#     property = ad_test_mat_property
#     block = '2 3'
#   [../]
# []

[Functions]
  [./heat_source_x]
    type = PiecewiseLinear
    data_file = ./input_mat_params/path_eq_t_x.csv
    format = columns
    scale_factor = 1.0
  [../]
  [./heat_source_y]
    type = PiecewiseLinear
    data_file = ./input_mat_params/path_eq_t_y.csv
    format = columns
    scale_factor = 1.0
  [../]
  [./heat_source_z]
    type = PiecewiseLinear
    data_file = ./input_mat_params/path_eq_t_z.csv
    format = columns
    scale_factor = 1.0
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
  [./ux_bottom_fix]
    type = ADDirichletBC
    variable = disp_x
    boundary = 1
    value = 0.0
  [../]
  [./uy_bottom_fix]
    type = ADDirichletBC
    variable = disp_y
    boundary = 1
    value = 0.0
  [../]
  [./uz_bottom_fix]
    type = ADDirichletBC
    variable = disp_z
    boundary = 1
    value = 0.0
  [../]
  [./temp_bottom_fix]
    type = ADDirichletBC
    variable = temp
    boundary = 1
    value = ${T_room}
  [../]


  [./convective]
    # type = ConvectiveFluxFunction # Convective flux, e.g. q'' = h*(Tw - Tf)
    type = ADConvectiveHeatFluxBC
    variable = temp
    boundary = 'moving_boundary'
    # coefficient = 2e-5
    heat_transfer_coefficient = 2e-4
    T_infinity = ${T_ambient}
  [../]
  # [./convective_substrate]
  #   # type = ConvectiveFluxFunction # Convective flux, e.g. q'' = h*(Tw - Tf)
  #   type = ADConvectiveHeatFluxBC
  #   variable = temp
  #   boundary = 2
  #   # coefficient = 2e-5
  #   heat_transfer_coefficient = 2e-4
  #   T_infinity = ${T_ambient}
  # [../]
  # [./convective_middle]
  #   # type = ConvectiveFluxFunction # Convective flux, e.g. q'' = h*(Tw - Tf)
  #   type = ADConvectiveHeatFluxBC
  #   variable = temp
  #   boundary = 'middle'
  #   # coefficient = 2e-5
  #   heat_transfer_coefficient = 2e-4
  #   T_infinity = ${T_ambient}
  # [../]
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
    a = 4.0
    b = 4.0
    c = 4.0
    power = 450
    efficienty = 0.36
    factor = 1.0
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
  [./ad_stateful]
    type = ADStatefulTest
    prop_names = ad_test_mat_property
    prop_values = 1.0
    block = '2 3'
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
  # petsc_options_iname = '-pc_type -ksp_type'
  # petsc_options_value = 'lu  preonly'

  petsc_options_iname = '-ksp_type -pc_type -pc_factor_mat_solver_package'
  petsc_options_value = 'preonly lu       superlu_dist'

  line_search = 'none'

  l_max_its = 100
  nl_max_its = 20
  nl_rel_tol = 1e-6
  nl_abs_tol = 1e-6
  # l_tol = 1e-5

  start_time = 0.0
  end_time = 0.2
  dt = 0.05
  dtmin = 1e-2
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
    expand_boundary_name = 'moving_boundary'
    activate_distance = 2.0
  [../]
[]

[Outputs]
  file_base = './output/EA_thermal_out'
  [./exodus]
    type = Exodus
  [../]
  [checkpoint]
    type = Checkpoint
    num_files = 2
  []
[]
