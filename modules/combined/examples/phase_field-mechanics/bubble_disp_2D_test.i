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

[Problem]
  type = ReferenceResidualProblem
      reference_vector = 'ref'
        extra_tag_vectors = 'ref'
[]

[GlobalParams]
  # op_num = 10
  # var_name_base = gr
  # int_width = 0.03
  displacements = 'disp_x disp_y'
[]

[MultiApps]
  [damage]
    type = TransientMultiApp
    input_files = 'bubble_c_2D_test.i'
  []
[]

[Transfers]
  [to_disp_x]
    type = MultiAppCopyTransfer
    multi_app = 'damage'
    direction = to_multiapp
    source_variable = 'disp_x'
    variable = 'disp_x'
  []
  [to_disp_y]
    type = MultiAppCopyTransfer
    multi_app = 'damage'
    direction = to_multiapp
    source_variable = 'disp_y'
    variable = 'disp_y'
  []
  [from_d]
    type = MultiAppCopyTransfer
    multi_app = 'damage'
    direction = from_multiapp
    source_variable = 'd'
    variable = 'd'
  []
[]

[Variables]
  [./disp_x]
  [../]
  [./disp_y]
  [../]
[]

[Modules]
  [./TensorMechanics]
    [./Master]
      displacements = 'disp_x disp_y'
      [./mech]
        add_variables = true
        strain = SMALL
        additional_generate_output = 'stress_yy stress_xy stress_xx strain_xx strain_xy strain_yy strain_zz hydrostatic_stress mid_principal_stress min_principal_stress max_principal_stress'
        decomposition_method = EigenSolution
        save_in = 'force_x force_y'
        displacements = 'disp_x disp_y'
        extra_vector_tags = 'ref'
      [../]
    [../]
  [../]
[]

[AuxVariables]
  [d]
  []
  [./stress_yy]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./force_x]
  [../]
  [./force_y]
  [../]
  [./c]
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
  # [./pressure]
  #   type = ParsedFunction
  #   value = '100'
  # [../]

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

[BCs]
  [./yfix]
    type = DirichletBC
    variable = disp_y
    boundary = bottom
    value = 0
  [../]
  [./xfix]
    type = DirichletBC
    variable = disp_x
    boundary = left
    value = 0
  [../]
  [./Pressure]
    [./coolantPressure]
      boundary = 'top right'
      factor = 0
      function = 1
    [../]
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
  [./ave_stress_right]
    type = SideAverageValue
    variable = stress_xx
    boundary = right
  [../]
  [./disp_x_right]
    type = SideAverageValue
    variable = disp_x
    boundary = right
  [../]
  [./react_x_top]
    type = NodalSum
    variable = force_x
    boundary = right
  [../]
[]

[Preconditioning]
  active = 'smp'
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
#  petsc_options_iname = '-ksp_type -pc_type -sub_pc_type -snes_max_it -sub_pc_factor_shift_type -pc_asm_overlap -snes_type'
#  petsc_options_value = 'gmres asm lu 100 NONZERO 2 vinewtonrsls'
  nl_rel_tol = 1e-6  ##nonlinear relative tolerance
  nl_abs_tol = 1e-6
  l_max_its = 10   ##max linear iterations Previous:200
  nl_max_its = 20  ##max nonlinear iterations Previous:50
  start_time = 0
  line_search = 'none'
  end_time = 200
  num_steps = 1500
  dt = 5
  dtmin = 1e-15
  automatic_scaling = true
#  [./TimeStepper]
#    type = IterationAdaptiveDT
#    dt = 1
#    optimal_iterations = 10
#    iteration_window = 0
#    growth_factor = 1.2
#    cutback_factor = 0.5
#  [../]

  picard_max_its = 20
  picard_rel_tol = 1e-6
  picard_abs_tol = 1e-6
  accept_on_max_picard_iteration = true
[]

[Outputs]
  # [exodus]
  #   type = Exodus
  #   interval = 25
  #   execute_on = 'initial timestep_end'
  # []
  csv = true
  nemesis = true
#gnuplot = true
[]
