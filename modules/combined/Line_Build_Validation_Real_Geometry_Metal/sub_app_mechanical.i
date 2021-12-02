T_room = 300
# T_ambient = 300
T_melt = 1700

# speed = 2 # mm/s
speed = 10.58e-3 # 10 mm/s = 10e-3 mm/ms
power = 300e-3 # 300W = kg*m^2/s^3 = 300e-3 kg*mm^2/ms^3
r = 300e-3 # 400 um = 400e-3 mm
dt = 2 #'${fparse 0.3*r/speed}' # ms

refine = 0

[GlobalParams]
  displacements = 'disp_x disp_y disp_z'
  volumetric_locking_correction = true
[]

[Problem]
  kernel_coverage_check = false
  material_coverage_check = false
[]

[Mesh]
  [mesh]
    type = GeneratedMeshGenerator
    dim = 3
    xmin = -1
    xmax = 1
    ymin = -2.5
    ymax = 2.5
    zmin = 0
    zmax = 2
    nx = 20
    ny = 50
    nz = 20
  []
  [add_set1]
    type = SubdomainBoundingBoxGenerator
    input = mesh
    block_id = 3
    bottom_left = '-50 -50 0'
    top_right = '50 50 1.0'
  []
  [add_set2]
    type = SubdomainBoundingBoxGenerator
    input = add_set1
    block_id = 1
    bottom_left = '-50 -50 1.0'
    top_right = '50 50 2'
  []
  [add_set3]
    type = SubdomainBoundingBoxGenerator
    input = add_set2
    block_id = 2
    bottom_left = '-0.1 -2 1.0'
    top_right = '0 -1.9 1.1'
  []
  # [add_set4]
  #   type = SubdomainBoundingBoxGenerator
  #   input = add_set3
  #   block_id = 4
  #   bottom_left = '-0.1 -2 0.9'
  #   top_right = '0 -1.9 1.0'
  # []
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

  uniform_refine = ${refine}
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
[]

[AuxVariables]
  [temp_aux]
    order = FIRST
    family = LAGRANGE
    block = '1 2 3'
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
  [x_coord]
    order = FIRST
    family = LAGRANGE
  []
  [y_coord]
    order = FIRST
    family = LAGRANGE
  []
  [z_coord]
    order = FIRST
    family = LAGRANGE
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
  # [x_coord]
  #   type = FunctionAux
  #   function = x_coord
  #   variable = x_coord
  # []
  # [y_coord]
  #   type = FunctionAux
  #   function = y_coord
  #   variable = y_coord
  # []
  # [z_coord]
  #   type = FunctionAux
  #   function = z_coord
  #   variable = z_coord
  # []
[]

[Functions]
  # [heat_source_x]
  #   type = ConstantFunction
  #   value = 0
  # []
  # [heat_source_y]
  #   type = ParsedFunction
  #   value = '-2 + ${speed}*t '
  # []
  # [heat_source_z]
  #   type = ConstantFunction
  #   value = 0.5
  # []
  # [scan_length_y]
  #   type = ParsedFunction
  #   value = '${speed}*t '
  # []
  # for monitoring the deposited material geometry
  # [x_coord]
  #   type = ParsedFunction
  #   value = 'x'
  # []
  # [y_coord]
  #   type = ParsedFunction
  #   value = 'y'
  # []
  # [z_coord]
  #   type = ParsedFunction
  #   value = 'z'
  # []
[]

[BCs]
  [ux_bottom_fix]
    type = ADDirichletBC
    variable = disp_x
    boundary = 'back'
    value = 0.0
  []
  [uy_bottom_fix]
    type = ADDirichletBC
    variable = disp_y
    boundary = 'back'
    value = 0.0
  []
  [uz_bottom_fix]
    type = ADDirichletBC
    variable = disp_z
    boundary = 'back'
    value = 0.0
  []
[]

[Materials]
  [E]
    type = ADPiecewiseLinearInterpolationMaterial
    x = '0 294.994  1671.48  1721.77 1e7'
    y = '201.232e3 201.232e3 80.0821e3 6.16016e3 6.16016e3' #MPa # 10^9 Pa = 10^9 kg/m/s^2 = kg/mm/ms^2
    # y = '6.16016e3 6.16016e3 6.16016e3 6.16016e3 6.16016e3'
    property = youngs_modulus
    variable = temp_aux
    extrapolation = false
    block = '2 3'
  []
  [nu]
    type = ADPiecewiseLinearInterpolationMaterial
    x = '0 294.994 1669.62 1721.77 1e7'
    y = '0.246407 0.246407   0.36961  0.36961 0.36961' #''0.513347 0.513347'
    property = poissons_ratio
    variable = temp_aux
    extrapolation = false
    block = '2 3'
  []
  [elasticity_tensor]
    type = ADComputeVariableIsotropicElasticityTensor
    youngs_modulus = youngs_modulus
    poissons_ratio = poissons_ratio
    block = '2 3'
  []

  [thermal_expansion_strain_product]
    type = ADComputeThermalExpansionEigenstrain
    stress_free_temperature = ${T_melt}
    # thermal_expansion_coeff = 1.72e-5
    thermal_expansion_coeff = 6.72e-6 #1.72e-5 /K
    temperature = temp_aux
    eigenstrain_name = thermal_eigenstrain_product
    block = '2'
  []
  [thermal_expansion_strain_substrate]
    type = ADComputeThermalExpansionEigenstrain
    stress_free_temperature = ${T_room}
    # thermal_expansion_coeff = 1.72e-5
    thermal_expansion_coeff = 6.72e-6 #1.72e-5 /K
    temperature = temp_aux
    eigenstrain_name = thermal_eigenstrain_substrate
    block = '3'
  []

  [stress]
    type = ADComputeFiniteStrainElasticStress
    block = '2 3'
  []

  # [./radial_return_stress]
  #   type = ADComputeMultipleInelasticStress
  #   inelastic_models = 'rate_temp_plas'
  #   block = '2 3 4'
  # [../]

  # [power_law_hardening]
  #   type = ADIsotropicPowerLawHardeningStressUpdate
  #   strength_coefficient = 847 #K
  #   strain_hardening_exponent = 0.06 #n
  #   relative_tolerance = 1e-6
  #   absolute_tolerance = 1e-6
  #   temperature = temp_aux
  #   block = '2 3'
  # []

  # [./rate_temp_plas]
  #   type = ADRateTempDependentStressUpdate
  #   temperature = temp_aux
  #   Y0 = 5.264 # 5.264e03 [MPa]
  #   Rd1 = 8.565e-7 #8.565e-4 [MPa]
  #   hxi = 1.670e6
  #   Ex = '294.994  1671.48  1721.77'
  #   Ey = '201.232 80.0821 6.16016' # GPa
  #   f1 = 9.178e-05 #[1/ms]
  #   nux = '294.994 1669.62 1721.77'
  #   nuy = '0.246407   0.36961  0.513347'
  #   # n2 = ${T_melt}
  #   # absolute_tolerance = 1e-8
  #   # relative_tolerance = 1e-6
  #   block = '2 3 4'
  #   # use_substep = true
  #   # max_inelastic_increment = 0.02
  # [../]
[]

[UserObjects]
  [activated_elem_uo_beam]
    type = CoupledVarThresholdElementSubdomainModifier
    execute_on = 'TIMESTEP_BEGIN'
    coupled_var = temp_aux
    block = 1
    subdomain_id = 2
    criterion_type = ABOVE
    threshold = ${T_melt}
    moving_boundary_name = 'moving_boundary'
    apply_initial_conditions = false
  []
  # [activated_elem_uo_melt]
  #   type = CoupledVarThresholdElementSubdomainModifier
  #   execute_on = 'TIMESTEP_BEGIN'
  #   coupled_var = temp_aux
  #   block = 3
  #   subdomain_id = 4
  #   criterion_type = ABOVE
  #   threshold = ${T_melt}
  #   # moving_boundary_name = 'moving_boundary'
  #   apply_initial_conditions = false
  # []
  # [activated_elem_uo]
  #   type = ActivateElementsCoupled
  #   coupled_var = temp_aux
  #   activate_type = above
  #   active_subdomain_id = '2'
  #   expand_boundary_name= 'moving_boundary'
  #   activate_value= ${T_melt}
  #   execute_on = 'TIMESTEP_BEGIN'
  # []
[]

# [Adaptivity]
#   marker = marker
#   initial_marker = marker
#   max_h_level = 1
#   [Indicators/indicator]
#     type = GradientJumpIndicator
#     variable = temp_aux
#   []
#   [Markers/marker]
#     type = ErrorFractionMarker
#     indicator = indicator
#     coarsen = 0.0
#     refine = 0.5
#     check_subdomain_consistent_for_coarsen = true
#   []
# []

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

  automatic_scaling = true

  petsc_options_iname = '-ksp_type -pc_type -pc_factor_mat_solver_package -pc_factor_shift_type '
                        '-pc_factor_shift_amount'
  petsc_options_value = 'preonly lu       superlu_dist NONZERO 1e-10'

  # petsc_options_iname = '-pc_type -ksp_type -pc_factor_shift_type -pc_factor_shift_amount'
  # petsc_options_value = 'lu  preonly NONZERO 1e-10'

  line_search = 'none'

  l_max_its = 100
  nl_max_its = 15
  nl_rel_tol = 1e-8
  nl_abs_tol = 1e-10

  start_time = 0.0
  end_time = '${fparse 3/speed}'
  dt = ${dt} # ms
  dtmin = 1e-6

  error_on_dtmin = false
[]

[Outputs]
  file_base = 'output_multiapp/Line_sub_speed_${speed}_power_${power}_r_${r}_dt_${dt}'
  csv = true
  [exodus]
    type = Exodus
    file_base = 'output_multiapp/Exodus_speed_${speed}_power_${power}_r_${r}_dt_${dt}/Sub_adapt'
    # execute_on = 'INITIAL TIMESTEP_END'
    interval = 1
  []
[]

[Postprocessors]
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
  [pp_power]
    type = ElementAverageValue
    variable = power_aux
  []
  [pp_speed]
    type = ElementAverageValue
    variable = speed_aux
  []
[]
