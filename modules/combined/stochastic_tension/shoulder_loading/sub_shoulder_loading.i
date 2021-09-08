[Mesh]
  [gmsh]
    type = GmshTensionSampleGenerator
    dim = 2
    right_circle_x = 5.75
    right_circle_y = 0
    right_circle_r = 0.75
    left_circle_x = -5.75
    left_circle_y = 0
    left_circle_r = 0.75
    fillet_radius = 0.3
    spacing = 0.15
    scale = 1000.0 # scale to micron
  []
[]

[GlobalParams]
  displacements = 'disp_x disp_y'
  out_of_plane_strain = strain_zz
[]

[Variables]
  [strain_zz]
  []
[]

[AuxVariables]
  [temp]
    initial_condition = 297
  []
  [nl_strain_zz]
    order = CONSTANT
    family = MONOMIAL
  []
[]

[Modules/TensorMechanics/Master]
  [all]
    planar_formulation = WEAK_PLANE_STRESS
    strain = FINITE
    incremental = true
    add_variables = true
    generate_output = 'stress_xx stress_yy stress_xy strain_xx strain_yy strain_xy elastic_strain_yy '
                      'vonmises_stress'
    use_automatic_differentiation = true
  []
[]

[AuxKernels]
  [strain_zz]
    type = ADRankTwoAux
    rank_two_tensor = total_strain
    variable = nl_strain_zz
    index_i = 2
    index_j = 2
  []
[]

[Materials]
  [elasticity_tensor]
    type = ADComputeIsotropicElasticityTensor
    youngs_modulus = 200e3 # MPa
    poissons_ratio = 0.3
    constant_on = SUBDOMAIN
  []
  [radial_return_stress]
    type = ADComputeMultipleInelasticStress
    inelastic_models = 'power_law_creep'
  []
  [power_law_creep]
    type = ADPowerLawCreepStressUpdate
    coefficient = 1.0e-15
    n_exponent = 4
    activation_energy = 3.0e2
    temperature = temp
    relative_tolerance = 1e-6
    absolute_tolerance = 1e-6
  []
[]

[BCs]
  [left_x]
    type = FunctionDirichletBC
    variable = disp_x
    boundary = '16 17 21 22'
    function = -50*t
  []
  # [./left_y]
  #   type = DirichletBC
  #   variable = disp_y
  #   boundary = '18 20'
  #   value = 0
  # [../]
  [right_x]
    type = FunctionDirichletBC
    variable = disp_x
    boundary = '5 6 10 11'
    function = 50*t
  []
  # [./right_y]
  #   type = DirichletBC
  #   variable = disp_y
  #   boundary = '7 9'
  #   value = 0
  # [../]
[]

[Preconditioning]
  [smp]
    type = SMP
    full = true
  []
[]

[Executioner]
  type = Transient

  solve_type = 'PJFNK'

  petsc_options_iname = '-pc_type -pc_factor_mat_solver_package -pc_factor_shift_type '
                        '-pc_factor_shift_amount'
  petsc_options_value = 'lu    superlu_dist NONZERO 1e-13'
  line_search = 'none'

  [TimeStepper]
    type = IterationAdaptiveDT
    dt = 0.1
    optimal_iterations = 8
    iteration_window = 2
  []
  dtmin = 1e-4
  end_time = 8
[]

[Outputs]
  file_base = ./output
  execute_on = 'initial timestep_end'
  exodus = false
  csv = false
[]

[Postprocessors]
  [stress_xx_center]
    type = PointValue
    variable = stress_xx
    point = '0 0 0'
    use_displaced_mesh = false
  []
  [strain_xx_center]
    type = PointValue
    variable = strain_xx
    point = '0 0 0'
    use_displaced_mesh = false
  []
  [stress_xx_top]
    type = PointValue
    variable = stress_xx
    point = '0 0.6 0'
    use_displaced_mesh = false
  []
  [stress_xx_left]
    type = PointValue
    variable = stress_xx
    point = '-0.6 0 0'
    use_displaced_mesh = false
  []
  [stress_xx_right]
    type = PointValue
    variable = stress_xx
    point = '0.6 0 0'
    use_displaced_mesh = false
  []
  [strain_xx_top]
    type = PointValue
    variable = strain_xx
    point = '0 0.6 0'
    use_displaced_mesh = false
  []
  [stress_xx_bot]
    type = PointValue
    variable = stress_xx
    point = '0 -0.6 0'
    use_displaced_mesh = false
  []
  [strain_xx_bot]
    type = PointValue
    variable = strain_xx
    point = '0 -0.6 0'
    use_displaced_mesh = false
  []
  [strain_xx_left]
    type = PointValue
    variable = strain_xx
    point = '-0.6 0 0'
    use_displaced_mesh = false
  []
  [strain_xx_right]
    type = PointValue
    variable = strain_xx
    point = '0.6 0 0'
    use_displaced_mesh = false
  []
[]

# [Reporters]
#   [acc]
#     type = AccumulateReporter
#     reporters = 'stress_xx/value strain_xx/value'
#   []
# []
