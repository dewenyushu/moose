# Interpolated material property data from paper
Ex = '294.994  1671.48  1721.77'
Ey = '201.232e9 80.0821e9 6.16016e9'
nux = '294.994 1669.62 1721.77'
nuy = '0.246407   0.36961  0.513347'

T_init = 1000
T_melt = 2000

T_room = 300

[GlobalParams]
  displacements = 'disp_x disp_y disp_z'
[]

[Mesh]
  [./mesh]
    type = FileMeshGenerator
    # file = cylinder_product.e
    file = ./input_mesh/solid_cylinder_product.e
  [../]
  [./add_set1]
    type = SubdomainBoundingBoxGenerator
    input = mesh
    block_id = 2
    bottom_left = '-25 -25 0'
    top_right = '25 25 5'
  [../]
  [./add_set2]
    type = SubdomainBoundingBoxGenerator
    input = add_set1
    block_id = 3
    bottom_left = '25 25 5'
    top_right = '25 25 20'
  [../]
  [./add_bnd]
    type = SideSetsFromAllElementFaces
    input = add_set2
    block = 1
    new_boundary = 'moving_boundary'
  [../]
[]

[Variables]
  [./disp_x]
    initial_condition = 0
  [../]
  [./disp_y]
    initial_condition = 0
  [../]
  [./disp_z]
    initial_condition = 0
  [../]
  [./temp]
    initial_condition = ${T_room}
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
  # [./heat_source_x]
  #   type = ParsedFunction
  #   value= '8.0*cos(2.0*pi/3.0*t)'
  # [../]
  # [./heat_source_y]
  #   type = ParsedFunction
  #   value= '8.0*sin(2.0*pi/3.0*t)'
  # [../]
  [./heat_source_x]
    type = PiecewiseLinear
    data_file = ./input_mat_params/path_t_x.csv
    format = columns
    scale_factor = 1.0
  [../]
  [./heat_source_y]
    type = PiecewiseLinear
    data_file = ./input_mat_params/path_t_y.csv
    format = columns
    scale_factor = 1.0
  [../]
  [./heat_source_z]
    type = PiecewiseLinear
    data_file = ./input_mat_params/path_t_z.csv
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
    type = ConvectiveFluxFunction # Convective flux, e.g. q'' = h*(Tw - Tf)
    variable = temp
    boundary = 'moving_boundary'
    coefficient = 1e-1
    T_infinity = ${T_room}
    marker_uo = activated_elem_uo
  [../]
  [./convective_substrate]
    type = ConvectiveFluxFunction # Convective flux, e.g. q'' = h*(Tw - Tf)
    variable = temp
    boundary = 2
    coefficient = 1e-1
    T_infinity = ${T_room}
  [../]
[]

[Materials]
  [./E]
    type = ADPiecewiseLinearInterpolationMaterial
    x = ${Ex}
    y = ${Ey}
    property = youngs_modulus
    variable = temp
    extrapolation = true
  [../]
  [./nu]
    type = ADPiecewiseLinearInterpolationMaterial
    x = ${nux}
    y = ${nuy}
    property = poissons_ratio
    variable = temp
    extrapolation = true
  [../]
  [./elasticity_tensor]
    type = ADComputeVariableIsotropicElasticityTensor
    youngs_modulus = youngs_modulus
    poissons_ratio = poissons_ratio
  [../]
  [./radial_return_stress]
    type = ADComputeMultipleInelasticStress
    inelastic_models = 'rate_temp_plas'
  [../]
  [./rate_temp_plas]
    type = ADRateTempDependentStressUpdate
    temperature = temp
    theta_melt= ${T_melt}
    Ex = ${Ex}
    Ey = ${Ey}
    nux = ${nux}
    nuy = ${nuy}
  [../]
  [./thermal_expansion_strain]
    type = ADComputeThermalExpansionEigenstrain
    stress_free_temperature = ${T_room}
    thermal_expansion_coeff = 5e-6
    temperature = temp
    eigenstrain_name = thermal_eigenstrain
  [../]
  [./volumetric_heat]
    type = FunctionPathEllipsoidHeatSource
    a = 5
    b = 5
    c = 3
    power = 1000
    efficienty = 0.5
    factor = 2
    function_x= heat_source_x
    function_y= heat_source_y
    function_z= heat_source_z
  [../]
  [./density]
    type = ADDensity
    density = 7.609e-6
  [../]
  [./heat]
    type = ADHeatConductionMaterial
    specific_heat_temperature_function = specific_heat
    thermal_conductivity_temperature_function = thermal_conductivity
    temp = temp
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
  end_time = 0.4
  dt = 2e-1
  dtmin = 1e-4
[]

[AuxVariables]
  [./von_mises]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./plastic_strain_eff]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./strain_eff]
    order = CONSTANT
    family = MONOMIAL
  [../]
  # [./stress_free_temp]
  #   [./InitialCondition]
  #     type = FunctionIC
  #     function = init_temp
  #   [../]
  # [../]
  [./activated_elem]
    order = CONSTANT
    family = MONOMIAL
  [../]
[]

[AuxKernels]
  [./von_mises_kernel]
    type = ADRankTwoScalarAux
    variable = von_mises
    rank_two_tensor = stress
    execute_on = timestep_end
    scalar_type = VonMisesStress
  [../]
  [./eff_plastic_strain_kernel]
    type = ADRankTwoScalarAux
    variable = plastic_strain_eff
    rank_two_tensor = plastic_strain
    execute_on = timestep_end
    scalar_type = EffectiveStrain
  [../]
  [./activated_elem]
    type = ActivatedElementsMarker
    melt_temperature = 480
    temp_aux = temp
    variable = activated_elem
    execute_on = timestep_begin
  [../]
[]

[UserObjects]
  # [./activated_elem_uo]
  #   type = ActivatedElementsMarkerUO
  #   melt_temperature = 480
  #   temp_aux = temp
  #   execute_on = timestep_begin
  # [../]
  [./activated_elem_uo]
    type = ActivatedElementsMarkerUO
    melt_temperature = 480
    temp_aux = temp
    execute_on = timestep_begin
  [../]
[]

[Postprocessors]
  # [./strain_yy]
  #   type = ElementAverageValue
  #   variable = strain_yy
  # [../]
  # [./strain_xy]
  #   type = ElementAverageValue
  #   variable = strain_xy
  # [../]
  # [./stress_yy]
  #   type = ElementAverageValue
  #   variable = stress_yy
  # [../]
  # [./stress_xy]
  #   type = ElementAverageValue
  #   variable = stress_xy
  # [../]
  [./temperature]
    type = ElementAverageValue
    variable = temp
  [../]
[]

[Outputs]
  file_base = 'stress_free_T${T_init}_out'
  [./exodus]
    type = Exodus
  [../]
  [./csv]
    type = CSV
    # file_base = 'tension_ramp_T${T_init}_out'
    execute_on = 'INITIAL TIMESTEP_END'
  [../]
[]
