# Interpolation data from paper
Ex = '294.994  1671.48  1721.77'
Ey = '201.232e9 80.0821e9 6.16016e9'
nux = '294.994 1669.62 1721.77'
nuy = '0.246407   0.36961  0.513347'

T_init = 1000
T_melt = 2000

T_room = 800

refine=0

[GlobalParams]
  displacements = 'disp_x disp_y disp_z'
[]

[Mesh]
  type = GeneratedMesh
  dim = 3
  xmin =-5
  xmax =5
  ymin =0
  ymax =0.5
  zmin =-5
  zmax = 5
  nx = 20
  ny = 1
  nz = 20
  uniform_refine = ${refine}
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
    [./InitialCondition]
      type = FunctionIC
      function = init_temp
    [../]
    # initial_condition = ${T_room}
  [../]
[]

[Modules/TensorMechanics/Master]
  [./all]
    strain = FINITE
    incremental = true
    add_variables = true
    generate_output = 'stress_xx stress_yy stress_zz stress_xy stress_yz stress_xz elastic_strain_yy plastic_strain_yy strain_yy strain_xx strain_zz strain_xy strain_xz strain_yz'
    use_automatic_differentiation = true
    eigenstrain_names = thermal_eigenstrain
    use_displaced_mesh = true
  [../]
[]

[Kernels]
  # [./heat]
  #   type = ADHeatConduction
  #   variable = temp
  #   thermal_conductivity = thermal_conductivity
  # [../]
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
  # [./heatsource]
  #   type = ADMatHeatSource
  #   material_property = volumetric_heat
  #   variable = temp
  #   scalar = 1
  #   use_displaced_mesh = true
  # [../]
[]

[Functions]
  # [./init_temp]
  #   type = ParsedFunction
  #   value = 'if(x>-0.6 & x<0.6 & z>-0.6 & z<0.6, ${T_init}, ${T_room})'
  # [../]
  [./init_temp]
    type = ParsedFunction
    value = 'if( y>0.4 & x*x + z*z< 2.0, ${T_init}*exp(-(x*x+z*z)/5.68), ${T_room})'
  [../]
[]

[BCs]
  [./u_y_bottom_fix]
    type = ADDirichletBC
    variable = disp_y
    boundary = bottom
    value = 0.0
  [../]
  [./u_x_bottom_fix]
    type = ADDirichletBC
    variable = disp_x
    boundary = bottom
    value = 0.0
  [../]
  [./u_z_bottom_fix]
    type = ADDirichletBC
    variable = disp_z
    boundary = bottom
    value = 0.0
  [../]
  [./u_xy_back_fix]
    type = ADDirichletBC
    variable = disp_z
    boundary = back
    value = 0.0
  [../]
  [./u_xy_front_fix]
    type = ADDirichletBC
    variable = disp_z
    boundary = front
    value = 0.0
  [../]
  [./u_yz_left_fix]
    type = ADDirichletBC
    variable = disp_x
    boundary = left
    value = 0.0
  [../]
  [./u_yz_right_fix]
    type = ADDirichletBC
    variable = disp_x
    boundary = left
    value = 0.0
  [../]



  [./temp_bottom_fix]
    type = ADDirichletBC
    variable = temp
    boundary = bottom
    value = ${T_room}
  [../]
  # [./temp_side_fix]
  #   type = ADDirichletBC
  #   variable = temp
  #   boundary = 'back front right left'
  #   value = ${T_room}
  # [../]
  # [./temp_xy_front_fix]
  #   type = ADDirichletBC
  #   variable = temp
  #   boundary = front
  #   value = ${T_room}
  # [../]
  # [./temp_top_ramp]
  #   type = FunctionDirichletBC
  #   variable = temp
  #   boundary = top
  #   function = temp_ramp
  # [../]
[]

[Materials]
  [./E]
    type = ADPiecewiseLinearInterpolationMaterial
    x = ${Ex}
    y = ${Ey}
    property = youngs_modulus
    variable = temp
  [../]
  [./nu]
    type = ADPiecewiseLinearInterpolationMaterial
    x = ${nux}
    y = ${nuy}
    property = poissons_ratio
    variable = temp
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
    stress_free_temperature = stress_free_temp
    thermal_expansion_coeff = 1e-4
    temperature = temp
    eigenstrain_name = thermal_eigenstrain
  [../]
  [./volumetric_heat]
    type = ADGenericConstantMaterial
    prop_names = 'volumetric_heat'
    prop_values = '10'
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
  end_time = 0.1
  dt = 1e-2
  dtmin = 1e-4
[]

[AuxVariables]
  [./von_mises]
    #Dependent variable used to visualize the Von Mises stress
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./stress_free_temp]
    [./InitialCondition]
      type = FunctionIC
      function = init_temp
    [../]
  [../]
[]

[AuxKernels]
  [./von_mises_kernel]
    #Calculates the von mises stress and assigns it to von_mises
    type = ADRankTwoScalarAux
    variable = von_mises
    rank_two_tensor = stress
    execute_on = timestep_end
    scalar_type = VonMisesStress
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
  [./stress_yy]
    type = ElementAverageValue
    variable = stress_yy
  [../]
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
