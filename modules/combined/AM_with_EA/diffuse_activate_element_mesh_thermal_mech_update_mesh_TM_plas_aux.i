T_melt = 1700
T_room = 300

[GlobalParams]
  displacements = 'disp_x disp_y disp_z'
[]

[Problem]
  kernel_coverage_check = false
  # material_coverage_check = false
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
  displacements='disp_x disp_y disp_z'
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
  [./all]
    strain = FINITE
    incremental = true
    add_variables = true
    generate_output = 'stress_xx stress_yy stress_zz stress_xy stress_yz stress_xz strain_yy strain_xx strain_zz strain_xy strain_xz strain_yz'
    use_automatic_differentiation = true
    eigenstrain_names = thermal_eigenstrain
    block = 1
  [../]
[]

[Kernels]
  [./time]
    type = ADHeatConductionTimeDerivative
    variable = temp
    block = 1
  [../]
  [./heat_conduct]
    type = ADHeatConduction
    variable = temp
    use_displaced_mesh = true
    thermal_conductivity = thermal_conductivity
    block=1
  [../]
  [./heatsource]
    type = ADMatHeatSource
    material_property = volumetric_heat
    variable = temp
    scalar = 1
    use_displaced_mesh = true
    block=1
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
    boundary = 'back'
    value = ${T_room}
  [../]
[]

[Materials]
  [./E]
    type = ADPiecewiseLinearInterpolationMaterial
    x = '294.994  1671.48  1721.77'
    y = '201.232e3 80.0821e3 6.16016e3' # [N/mm^3], i.e., MPa
    property = youngs_modulus
    variable = temp
    block = 1
    extrapolation = true
  [../]
  [./nu]
    type = ADPiecewiseLinearInterpolationMaterial
    x = '294.994 1669.62 1721.77'
    y = '0.246407   0.36961  0.513347'
    property = poissons_ratio
    variable = temp
    block = 1
    extrapolation = true
  [../]
  [./elasticity_tensor]
    type = ADComputeVariableIsotropicElasticityTensor
    youngs_modulus = youngs_modulus
    poissons_ratio = poissons_ratio
    block = 1
  [../]
  # [./stress]
  #   type = ADComputeLinearElasticStress
  #   block = 1
  # [../]

  # [./radial_return_stress]
  #   type = ADComputeMultipleInelasticStress
  #   inelastic_models = 'power_law_creep'
  #   block = 1
  # [../]
  # [./power_law_creep]
  #   type = ADPowerLawCreepStressUpdate
  #   coefficient = 1.0e-15
  #   n_exponent = 4
  #   activation_energy = 3.0e5
  #   temperature = temp
  #   block = 1
  # [../]

  [./radial_return_stress]
    type = ADComputeMultipleInelasticStress
    inelastic_models = 'rate_temp_plas'
    block = 1
  [../]
  [./rate_temp_plas]
    type = ADRateTempDependentStressUpdate
    temperature = temp
    theta_melt = 2000 # set to a large value to avoid fluid part
    Y0 = 5.264e03 # [MPa]
    Rd1 = 8.565e-4 # [MPa]
    hxi = 1.670e-06 # [mm/(s*MPa)]
    Ex = '294.994  1671.48  1721.77'
    Ey = '201.232e3 80.0821e3 6.16016e3' # [N/mm^3], i.e., MPa
    nux = '294.994 1669.62 1721.77'
    nuy = '0.246407   0.36961  0.513347'
    block = 1
  [../]

  [./thermal_expansion_strain]
    type = ADComputeThermalExpansionEigenstrain
    stress_free_temperature = ${T_melt}
    thermal_expansion_coeff = 1e-6
    temperature = temp
    eigenstrain_name = thermal_eigenstrain
    # activated_elem_aux = activated_elem
    # melt_temperature = ${T_melt}
    block = 1
  [../]
  [./volumetric_heat]
    type = FunctionPathEllipsoidHeatSource
    a = 2
    b = 2
    c = 2
    power = 400
    efficienty = 0.36
    factor = 1
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
    specific_heat_temperature_function = specific_heat
    thermal_conductivity_temperature_function = thermal_conductivity
    temp = temp
    block = 1
  [../]
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

[AuxVariables]
  [./von_mises]
    order = CONSTANT
    family = MONOMIAL
    block = 1
  [../]
  [./plastic_strain_eff]
    order = CONSTANT
    family = MONOMIAL
    block = 1
  [../]
  # [./stress_free_temp]
  #   [./InitialCondition]
  #     type = FunctionIC
  #     function = init_temp
  #   [../]
  # [../]
[]

[AuxKernels]
  [./von_mises_kernel]
    type = ADRankTwoScalarAux
    variable = von_mises
    rank_two_tensor = stress
    execute_on = timestep_end
    scalar_type = VonMisesStress
    block = 1
  [../]
  [./eff_plastic_strain_kernel]
    type = ADRankTwoScalarAux
    variable = plastic_strain_eff
    rank_two_tensor = plastic_strain
    execute_on = timestep_end
    scalar_type = EffectiveStrain
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
  # petsc_options_iname = '-pc_type -ksp_type -pc_factor_shift_type -pc_factor_shift_amount'
  # petsc_options_value = 'lu  preonly NONZERO 1e-10'

  petsc_options_iname = '-pc_type -ksp_type'
  petsc_options_value = 'lu  preonly'

  line_search = 'none'

  l_max_its = 10
  nl_max_its = 20
  nl_rel_tol = 1e-6
  nl_abs_tol = 1e-12
  # l_tol = 1e-5

  start_time = 0.0
  end_time = 0.2
  dt = 1e-1
  dtmin = 1e-4
[]


# [Postprocessors]
#   [./temperature]
#     type = ElementAverageValue
#     variable = temp
#   [../]
# []

[Outputs]
  file_base = 'temp_out_TM'
  [./exodus]
    type = Exodus
  [../]
[]
