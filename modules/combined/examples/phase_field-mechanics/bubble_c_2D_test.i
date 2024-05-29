[Mesh]
  [gmg]
    type = DistributedRectilinearMeshGenerator
    dim = 2
    nx = 300
    ny = 300
    xmax = 3.0
    ymax = 3.0
    partition = square
  []
[]

[GlobalParams]
  # op_num = 10
  # var_name_base = gr
  displacements = 'disp_x disp_y'
  # int_width = 0.03
[]

[Variables]
  [./d]
  [../]
[]


[AuxVariables]
  [./stress_yy]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./force_x]
    order = FIRST
    family = LAGRANGE
  [../]
  [./force_y]
    order = FIRST
    family = LAGRANGE
  [../]
  [./c]
    order = FIRST
    family = LAGRANGE
  [../]
  [./C1111]
    order = CONSTANT
    family = MONOMIAL
    block = 0
  [../]
  [./strain_yy]
    family = MONOMIAL
    order = CONSTANT
  [../]
  [./bounds_dummy]
    order = FIRST
    family = LAGRANGE
  [../]
  [disp_x]
  []
  [disp_y]
  []
[]

[ICs]
  [c]
    type = SmoothCircleIC
    variable = c
    x1 = 1.5
    y1 = 1.5
    radius = 0.5
    invalue = 1.0
    outvalue = 0.0
    int_width = 0.03
  []
[]

[Functions]
  [./tfunc]
    type = ParsedFunction
    value = '0*t'
  [../]
  [./pressure]
    type = PiecewiseLinear
    x = '0 150 200'
    y = '80 200 200'
  [../]
[]

[Bounds]
  [./d_upper_bound]
    type = ConstantBoundsAux
    variable = bounds_dummy
    bounded_variable = d
    bound_type = upper
    bound_value = 1.0
  [../]
  [./d_lower_bound]
    type = VariableOldValueBoundsAux
    variable = bounds_dummy
    bounded_variable = d
    bound_type = lower
  [../]
[]

[Kernels]
  [./pfbulk]
    type = AllenCahn
    variable = d
    mob_name = L
    f_name = F
  [../]
  [./dcdt]
    type = TimeDerivative
    variable = d
  [../]
  [./acint]
    type = ACInterface
    variable = d
    mob_name = L
    kappa_name = kappa_op
  [../]
[]

[AuxKernels]
  [./C1111]
    type = RankFourAux
    variable = C1111
    rank_four_tensor = matrix_elasticity_tensor
    index_l = 0
    index_j = 0
    index_k = 0
    index_i = 0
    execute_on = timestep_end
  [../]
[]

[Materials]
  [./pfbulkmat]
    type = GenericConstantMaterial
    prop_names = 'l visco'
    prop_values = '0.01 1e-3'
  [../]
  [pressure]
    type = GenericFunctionMaterial
    block = 0
    prop_names = fracture_pressure
    prop_values = pressure
    # factor = 1e-6
  []
  [./gc]
    type = GenericConstantMaterial
    prop_names = gc_prop
    prop_values = '0.0012'
  [../]

  [./define_mobility]
    type = ParsedMaterial
    material_property_names = 'gc_prop visco'
    f_name = L
    function = '1.0/(gc_prop * visco)'
  [../]
  [./define_kappa]
    type = ParsedMaterial
    material_property_names = 'gc_prop l'
    f_name = kappa_op
    function = 'gc_prop * l / 3.14159 * 2'
  [../]

  [./stress_void]
    type = ComputeLinearElasticStress
    block = 0
    base_name = void
  [../]
  [./strain_void]
    type = ComputeSmallStrain
    block = 0
    base_name = void
  [../]

  [./strain]
    type = ComputeSmallStrain
    block = 0
    base_name = matrix
  [../]

  [./damage_stress]
    type = ComputeLinearElasticPFFractureStress
    c = d
    E_name = 'elastic_energy'
    D_name = 'degradation'
    F_name = 'local_fracture_energy'
    I_name = 'indicator_function'
   decomposition_type = strain_vol_dev
   #decomposition_type = none
    use_snes_vi_solver = true
    base_name = matrix
  [../]
  [./indicator_function]
    type = DerivativeParsedMaterial
    f_name = indicator_function
    args = 'd'
    function = 'd*d'
    derivative_order = 2
  [../]
  [./degradation]
    type = DerivativeParsedMaterial
    f_name = degradation
    args = 'd'
    function = '((1.0-d)^2+eta)/((1.0-d)^2+d*(1-0.5*d)*(4/3.14159/l*E*gc_prop/sigma^2))'
    material_property_names = 'gc_prop l'
    constant_names       = 'E sigma eta'
    constant_expressions = '385000 130 1e-4'
    derivative_order = 2
  [../]
  [./fracture_energy]
    type = DerivativeParsedMaterial
    f_name = local_fracture_energy
    args = 'd'
    material_property_names = 'gc_prop l'
    function = 'gc_prop/l/3.14159*(2*d-d^2)'
    derivative_order = 2
  [../]
#  [./fracture_driving_energy]
#    type = DerivativeSumMaterial
#    args = d
#    sum_materials = 'elastic_energy local_fracture_energy'
#    derivative_order = 2
#    f_name = F
#  [../]
  [./fracture_driving_energy]
    type = DerivativeParsedMaterial
    args = 'c d'
    material_property_names = 'elastic_energy(d) local_fracture_energy(d)'
    function = '(1-c)*elastic_energy + local_fracture_energy'
    derivative_order = 2
    f_name = F
  [../]

  [./const_stress]
    type = ComputeExtraStressConstant
    block = 0
    base_name = void
    extra_stress_tensor = '-1 -1 -1 0 0 0'
    prefactor = fracture_pressure
  [../]
  [./global_stress]
    type = TwoPhaseStressMaterial
    base_A = matrix
    base_B = void
  [../]
  [./switching]
    type = SwitchingFunctionMaterial
    eta = c
  [../]
  [./elasticity_tensor]
    type = ComputeElasticityTensor
    C_ijkl = '385000 0.23'
    fill_method = symmetric_isotropic_E_nu
    base_name = matrix
  [../]
  [./elasticity_tensor_void]
    type = ComputeElasticityTensor
    C_ijkl = '3.85 0.23'
    fill_method = symmetric_isotropic_E_nu
    base_name = void
  [../]
[]

[Postprocessors]
  [./ave_stress_top]
    type = SideAverageValue
    variable = stress_yy
    boundary = top
  [../]
  [./disp_y_top]
    type = SideAverageValue
    variable = disp_y
    boundary = top
  [../]
  [./react_y_top]
    type = NodalSum
    variable = force_y
    boundary = top
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
  solve_type = NEWTON
  petsc_options_iname = '-pc_type -ksp_type -snes_type -pc_factor_shift_type -pc_factor_shift_amount '
  petsc_options_value = 'lu preonly  vinewtonrsls NONZERO 1e-10'

#  solve_type = PJFNK
#  petsc_options_iname = '-pc_type -sub_pc_type -snes_type'
#  petsc_options_value = 'asm lu vinewtonrsls'
  nl_rel_tol = 1e-6  ##nonlinear relative tolerance
  nl_abs_tol = 1e-6
  l_max_its = 10   ##max linear iterations Previous:200
  nl_max_its = 20  ##max nonlinear iterations Previous:50
  start_time=0
  line_search = 'none'
  end_time = 2000
  dtmax = 1
  dtmin = 1e-14
  automatic_scaling = true
#  [./TimeStepper]
#    type = IterationAdaptiveDT
#    dt = 1
#    optimal_iterations = 10
#    iteration_window = 0
#    growth_factor = 1.2
#    cutback_factor = 0.5
#  [../]
[]

[Outputs]
  print_linear_converged_reason = false
  print_nonlinear_converged_reason = false
  print_linear_residuals = false
  nemesis = true
  csv = true
#gnuplot = true
[]
