# Interpolated material property data from paper
Ex = '294.994  1671.48  1721.77'
Ey = '201.232e9 80.0821e9 6.16016e9'
nux = '294.994 1669.62 1721.77'
nuy = '0.246407   0.36961  0.513347'

T_init = 1000
T_melt = 2000

T_room = 300

# refine=0

[GlobalParams]
  displacements = 'disp_x disp_y disp_z'
[]

[Mesh]
  [./block1]
    type = GeneratedMeshGenerator
    dim = 3
    xmin =0
    xmax =10
    ymin =0
    ymax =10
    zmin =0
    zmax =0.5
    nx = 20
    ny = 20
    nz = 1
  [../]
  [./block1_id]
    type = SubdomainIDGenerator
    input = block1
    subdomain_id = 1
  [../]
  [./block1_rename_bnd]
    type = RenameBoundaryGenerator
    input = block1_id
    # back = 0, bottom = 1, right = 2, top = 3, left = 4, front = 5
    old_boundary_id='0 1 2 3 4 5'
    new_boundary_id='10 11 12 13 14 15'
  [../]
  # [./product]
  #   type = SideSetsFromAllElementFaces
  #   input = block1_rename_bnd
  #   block = 1
  #   new_boundary = 'moving_boundary'
  # [../]

  [./block2]
    type = GeneratedMeshGenerator
    dim = 3
    xmin =0
    xmax =10
    ymin =0
    ymax =10
    zmin =-0.5
    zmax =0
    nx = 20
    ny = 20
    nz = 1
  [../]
  [./block2_rename_bnd]
    type = RenameBoundaryGenerator
    input = block2
    # back = 0, bottom = 1, right = 2, top = 3, left = 4, front = 5
    old_boundary_id='0 1 2 3 4 5'
    new_boundary_id='20 21 22 23 24 25'
  [../]
  [./substrate]
    type = SubdomainIDGenerator
    input = block2_rename_bnd
    subdomain_id = 2
  [../]

  [./combined]
    type = MeshCollectionGenerator
    # inputs = 'product substrate'
    inputs = 'block1_rename_bnd substrate'
  [../]
  # uniform_refine = ${refine}
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
    use_displaced_mesh = true
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
  #   value= '5.0*t'
  # [../]
  # [./heat_source_y]
  #   type = ParsedFunction
  #   value= '0.25'
  # [../]
  # [./heat_source_z]
  #   type = ParsedFunction
  #   value= '0.25'
  # [../]
  [./heat_source_x]
    type = PiecewiseLinear
    x='0.0 2.0 2.2 4.2 4.4 6.4 6.6 8.6 8.8 10.8 11.0 13.0 13.2 15.2 15.4 17.4 17.6 19.6 19.8 21.8 22.0 24.0'
    y='0.0 10.0 10.0 0.0 0.0 10.0 10.0 0.0 0.0 10.0 10.0 0.0 0.0 10.0 10.0 0.0 0.0 10.0 10.0 0.0 0.0 10.0'
  [../]
  [./heat_source_y]
    type = PiecewiseLinear
    x='0.0 2.0 2.2 4.2 4.4 6.4 6.6 8.6 8.8 10.8 11.0 13.0 13.2 15.2 15.4 17.4 17.6 19.6 19.8 21.8 22.0 24.0'
    y='0.0 0.0 1.0 1.0 2.0 2.0 3.0 3.0 4.0 4.0 5.0 5.0 6.0 6.0 7.0 7.0 8.0 8.0 9.0 9.0 10.0 10.0'
  [../]
  [./heat_source_z]
    type = ConstantFunction
    value = 0.25  # half the height
  [../]
[]

[BCs]
  [./ux_bottom_fix]
    type = ADDirichletBC
    variable = disp_x
    boundary = 20
    value = 0.0
  [../]
  [./uy_bottom_fix]
    type = ADDirichletBC
    variable = disp_y
    boundary = 20
    value = 0.0
  [../]
  [./uz_bottom_fix]
    type = ADDirichletBC
    variable = disp_z
    boundary = 20
    value = 0.0
  [../]
  # [./uy_back_fix]
  #   type = ADDirichletBC
  #   variable = disp_y
  #   boundary = 23
  #   value = 0.0
  # [../]
  # [./uy_front_fix]
  #   type = ADDirichletBC
  #   variable = disp_y
  #   boundary = 21
  #   value = 0.0
  # [../]
  # [./ux_left_fix]
  #   type = ADDirichletBC
  #   variable = disp_x
  #   boundary = 24
  #   value = 0.0
  # [../]
  # [./ux_right_fix]
  #   type = ADDirichletBC
  #   variable = disp_x
  #   boundary = 22
  #   value = 0.0
  # [../]

  # [./temp_bottom_fix]
  #   type = ADDirichletBC
  #   variable = temp
  #   boundary = 20
  #   value = ${T_room}
  # [../]

  # [./convective]
  #   type = ConvectiveFluxFunction # Convective flux, e.g. q'' = h*(Tw - Tf)
  #   variable = temp
  #   boundary = 'moving_boundary'
  #   coefficient = 1e-1
  #   T_infinity = ${T_room}
  #   marker_uo = activated_elem_uo
  # [../]
  [./convective_substrate]
    type = ConvectiveFluxFunction # Convective flux, e.g. q'' = h*(Tw - Tf)
    variable = temp
    boundary = '20 21 22 23 24'
    coefficient = 1e-1
    T_infinity = 300
  [../]
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
    stress_free_temperature = ${T_room}
    thermal_expansion_coeff = 1e-5
    temperature = temp
    eigenstrain_name = thermal_eigenstrain
  [../]
  [./volumetric_heat]
    type = FunctionPathEllipsoidHeatSource
    a = 2
    b = 2
    c = 2
    power = 1000
    efficienty = 1.0
    factor = 2
    velocity = 5
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
  end_time = 2
  dt = 1e-1
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
