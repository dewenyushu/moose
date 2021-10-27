T_room = 500
T_melt = 1500
T_ambient = 500

speed = 1 # mm/s
power = 1 # w

[GlobalParams]
  displacements = 'disp_x disp_y disp_z'
[]

[Problem]
  kernel_coverage_check = false
[]

[Mesh]
  [mesh]
    type = FileMeshGenerator
    file = box_10x20x4.e
  []
  [add_set1]
    type = SubdomainBoundingBoxGenerator
    input = mesh
    block_id = 3
    bottom_left = '-50 -50 0'
    top_right = '50 50 4'
  []
  [add_set2]
    type = SubdomainBoundingBoxGenerator
    input = add_set1
    block_id = 1
    bottom_left = '-50 -50 4'
    top_right = '50 50 24'
  []
  [add_set3]
    type = SubdomainBoundingBoxGenerator
    input = add_set2
    block_id = 2
    bottom_left = '-0.5 -8 4'
    # top_right = '-1.8 -1.8 4.2'
    # top_right = '-1.5 -1.5 4.5'
    top_right = '0 -7.5 4.5'
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
  displacements = 'disp_x disp_y disp_z'
[]

[Variables]
  [disp_x]
    block = '2 3'
  []
  [disp_y]
    block = '2 3'
  []
  [disp_z]
    block = '2 3'
  []
  [temp]
    block = '2 3'
  []
[]

[AuxVariables]
  [temp_aux]
    order = FIRST
    family = LAGRANGE
  []
  [von_mises]
    order = CONSTANT
    family = MONOMIAL
    block = '2 3'
  []
  [plastic_strain_eff]
    order = CONSTANT
    family = MONOMIAL
    block = '2 3'
  []
  [power_aux]
    order = CONSTANT
    family = MONOMIAL
  []
  [speed_aux]
    order = CONSTANT
    family = MONOMIAL
  []
[]

[ICs]
  [temp_substrate]
    type = ConstantIC
    variable = temp
    value = ${T_room}
    block = '3'
  []
  [temp_product]
    type = ConstantIC
    variable = temp
    value = ${T_room}
    block = '2'
  []
[]

[Modules/TensorMechanics/Master]
  # [./all]
  strain = FINITE
  incremental = true
  add_variables = true
  generate_output = 'stress_xx stress_yy stress_zz stress_xy stress_yz stress_xz strain_yy strain_xx '
                    'strain_zz strain_xy strain_xz strain_yz'
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
  [time]
    type = ADHeatConductionTimeDerivative
    variable = temp
    block = '2 3'
  []
  [heat_conduct]
    type = ADHeatConduction
    variable = temp
    use_displaced_mesh = true
    thermal_conductivity = thermal_conductivity
    block = '2 3'
  []
  [heatsource]
    type = ADMatHeatSource
    material_property = volumetric_heat
    variable = temp
    scalar = 1
    use_displaced_mesh = true
    block = '2 3'
  []
[]

[AuxKernels]
  [von_mises_kernel]
    type = ADRankTwoScalarAux
    variable = von_mises
    rank_two_tensor = stress
    execute_on = timestep_end
    scalar_type = VonMisesStress
    block = '2 3'
  []
  [power]
    type = ConstantAux
    variable = power_aux
    value = ${power}
    # execute_on = 'INITIAL TIMESTEP_BEGIN'
  []
  [speed]
    type = ConstantAux
    variable = speed_aux
    value = ${speed}
    # execute_on = 'INITIAL TIMESTEP_BEGIN'
  []
  # [./eff_plastic_strain_kernel]
  #   type = ADRankTwoScalarAux
  #   variable = plastic_strain_eff
  #   rank_two_tensor = plastic_strain
  #   execute_on = timestep_end
  #   scalar_type = EffectiveStrain
  #   block = '2 3'
  # [../]
[]

[Functions]
  [heat_source_x]
    type = ConstantFunction
    value = 0
  []
  [heat_source_y]
    type = ParsedFunction
    value = '-8 + ${speed}*t '
  []
  [scan_length_y]
    type = ParsedFunction
    value = '${speed}*t '
  []
  [heat_source_z]
    type = ConstantFunction
    value = 4
  []
  [specific_heat]
    type = PiecewiseLinear
    data_file = Specific_Heat.csv
    format = columns
    scale_factor = 1.0
  []
  [thermal_conductivity]
    type = PiecewiseLinear
    data_file = Thermal_Conductivity.csv
    format = columns
    scale_factor = 0.05e-3
  []
[]

[BCs]
  [ux_bottom_fix]
    type = ADDirichletBC
    variable = disp_x
    boundary = 1
    value = 0.0
  []
  [uy_bottom_fix]
    type = ADDirichletBC
    variable = disp_y
    boundary = 1
    value = 0.0
  []
  [uz_bottom_fix]
    type = ADDirichletBC
    variable = disp_z
    boundary = 1
    value = 0.0
  []
  [temp_bottom_fix]
    type = ADDirichletBC
    variable = temp
    boundary = 1
    value = ${T_room}
  []

  [convective]
    # type = ConvectiveFluxFunction # Convective flux, e.g. q'' = h*(Tw - Tf)
    type = ADConvectiveHeatFluxBC
    variable = temp
    boundary = 'moving_boundary'
    # coefficient = 2e-5
    heat_transfer_coefficient = 0.016
    T_infinity = ${T_ambient}
  []
  [convective_substrate]
    # type = ConvectiveFluxFunction # Convective flux, e.g. q'' = h*(Tw - Tf)
    type = ADConvectiveHeatFluxBC
    variable = temp
    boundary = 2
    # coefficient = 2e-5
    heat_transfer_coefficient = 0.016
    T_infinity = ${T_ambient}
  []
  # [./convective_middle]
  #   # type = ConvectiveFluxFunction # Convective flux, e.g. q'' = h*(Tw - Tf)
  #   type = ADConvectiveHeatFluxBC
  #   variable = temp
  #   boundary = 'middle'
  #   # coefficient = 2e-5
  #   heat_transfer_coefficient = 0.016
  #   T_infinity = ${T_ambient}
  # [../]
[]

[Materials]
  [E]
    type = ADPiecewiseLinearInterpolationMaterial
    x = '294.994  1671.48  1721.77'
    y = '201.232e3 80.0821e3 6.16016e3' # [N/mm^3], i.e., MPa
    property = youngs_modulus
    variable = temp
    extrapolation = false
    block = '2 3'
  []
  [nu]
    type = ADPiecewiseLinearInterpolationMaterial
    x = '294.994 1669.62 1721.77'
    y = '0.246407   0.36961  0.513347'
    property = poissons_ratio
    variable = temp
    extrapolation = false
    block = '2 3'
  []
  [elasticity_tensor]
    type = ADComputeVariableIsotropicElasticityTensor
    youngs_modulus = youngs_modulus
    poissons_ratio = poissons_ratio
    block = '2 3'
  []

  [stress]
    type = ADComputeFiniteStrainElasticStress
    block = '2 3'
  []

  # [./radial_return_stress]
  #   type = ADComputeMultipleInelasticStress
  #   inelastic_models = 'power_law_hardening'
  #   block = '2 3'
  # [../]
  #
  # [power_law_hardening]
  #   type = ADIsotropicPowerLawHardeningStressUpdate
  #   strength_coefficient = 847 #K
  #   strain_hardening_exponent = 0.06 #n
  #   relative_tolerance = 1e-6
  #   absolute_tolerance = 1e-6
  #   temperature = temp
  #   block = '2 3'
  # []

  # [./rate_temp_plas]
  #   type = ADRateTempDependentStressUpdate
  #   temperature = temp
  #   Y0 = 5.264e03 # [MPa]
  #   Rd1 = 8.565e-4 # [MPa]
  #   hxi = 1.670e-06 # [mm/(s*MPa)]
  #   Ex = '294.994  1671.48  1721.77'
  #   Ey = '201.232e3 80.0821e3 6.16016e3' # [N/mm^3], i.e., MPa
  #   nux = '294.994 1669.62 1721.77'
  #   nuy = '0.246407   0.36961  0.513347'
  #   # n2 = ${T_melt}
  #   absolute_tolerance = 1e-8
  #   block = '2 3'
  #   use_substep = true
  #   max_inelastic_increment = 0.02
  # [../]

  [thermal_expansion_strain_product]
    type = ADComputeThermalExpansionEigenstrain
    stress_free_temperature = ${T_room}
    # thermal_expansion_coeff = 1.72e-5
    thermal_expansion_coeff = 1.72e-6
    temperature = temp
    eigenstrain_name = thermal_eigenstrain_product
    block = '2'
  []
  [thermal_expansion_strain_substrate]
    type = ADComputeThermalExpansionEigenstrain
    stress_free_temperature = ${T_room}
    # thermal_expansion_coeff = 1.72e-5
    thermal_expansion_coeff = 1.72e-6
    temperature = temp
    eigenstrain_name = thermal_eigenstrain_substrate
    block = '3'
  []

  [volumetric_heat]
    type = FunctionPathEllipsoidHeatSource
    r = 2.0
    power = ${power}
    efficiency = 0.8
    factor = 1
    function_x = heat_source_x
    function_y = heat_source_y
    function_z = heat_source_z
  []
  [density]
    type = ADDensity
    density = 7.609e-6
    block = '2 3'
  []
  [heat]
    type = ADHeatConductionMaterial
    specific_heat_temperature_function = specific_heat
    thermal_conductivity_temperature_function = thermal_conductivity
    temp = temp
    block = '2 3'
  []
[]

[UserObjects]
  [activated_elem_uo]
    type = CoupledVarThresholdElementSubdomainModifier
    execute_on = 'TIMESTEP_BEGIN'
    coupled_var = temp_aux
    block = 1
    subdomain_id = 2
    criterion_type = ABOVE
    threshold = ${T_melt}
    moving_boundary_name = 'moving_boundary'
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
  #Preconditioned JFNK (default)
  solve_type = 'PJFNK'

  # automatic_scaling = true

  petsc_options_iname = '-ksp_type -pc_type -pc_factor_mat_solver_package'
  petsc_options_value = 'preonly lu       superlu_dist'

  # petsc_options_iname = '-pc_type -ksp_type -pc_factor_shift_type -pc_factor_shift_amount'
  # petsc_options_value = 'lu  preonly NONZERO 1e-10'

  line_search = 'none'

  l_max_its = 100
  nl_max_its = 20
  nl_rel_tol = 1e-5
  nl_abs_tol = 1e-6

  start_time = 0.0
  end_time = '${fparse 12/speed}'
  dt = 0.1
  dtmin = 1e-4

  error_on_dtmin = false
[]

[Outputs]
  file_base = 'output/Line_sub_speed_${speed}_power_${power}'
  csv = true
  # [exodus]
  #   type = Exodus
  #   file_base = './output/Exodus_speed_${speed}_power_${power}/Line_sub'
  #   # execute_on = 'INITIAL TIMESTEP_END'
  #   interval = 1
  # []
[]

[Postprocessors]
  [max_temperature]
    type = ElementExtremeValue
    variable = temp
    value_type = max
    block = '2'
  []
  [min_temperature]
    type = ElementExtremeValue
    variable = temp
    value_type = min
    block = '2'
  []
  [max_von_mises_stress]
    type = ElementExtremeValue
    variable = von_mises
    value_type = max
    block = '2'
  []
  [min_von_mises_stress]
    type = ElementExtremeValue
    variable = von_mises
    value_type = min
    block = '2'
  []
  [deposition_volume]
    type = VolumePostprocessor
    block = '2'
    use_displaced_mesh = true
  []
  [deposition_length]
    type = FunctionValuePostprocessor
    function = scan_length_y
  []
  [pp_power]
    type = ElementAverageValue
    variable = power_aux
  []
  [pp_speed]
    type = ElementAverageValue
    variable = speed_aux
  []
[]
