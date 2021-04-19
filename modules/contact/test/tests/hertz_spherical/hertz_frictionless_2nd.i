[GlobalParams]
  # volumetric_locking_correction = true
  displacements = 'disp_x disp_y'
  order = SECOND
  family = LAGRANGE
[]

[Mesh]
  second_order = true
  [input_file]
    type = FileMeshGenerator
    file = fine_spherical_2d.e
  []
  [secondary]
    type = LowerDBlockFromSidesetGenerator
    new_block_id = 10001
    new_block_name = 'secondary_lower'
    sidesets = '3' # pellet_outer_surface
    input = input_file
  []
  [primary]
    type = LowerDBlockFromSidesetGenerator
    new_block_id = 10000
    sidesets = '2' # clad_inside_right
    new_block_name = 'primary_lower'
    input = secondary
  []
[]

[Problem]
  type = ReferenceResidualProblem
  extra_tag_vectors = 'ref'
  reference_vector = 'ref'
[]

[Variables]
  [./disp_x]
  [../]
  [./disp_y]
  [../]
  [frictionless_normal_lm]
    order = SECOND
    family = LAGRANGE
    block = 'secondary_lower'
    use_dual = false
  []
[]

[AuxVariables]
  [./stress_xx]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./stress_yy]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./stress_xy]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./saved_x]
  [../]
  [./saved_y]
  [../]
  [./diag_saved_x]
  [../]
  [./diag_saved_y]
  [../]
  # [./inc_slip_x]
  # [../]
  # [./inc_slip_y]
  # [../]
  # [./accum_slip_x]
  # [../]
  # [./accum_slip_y]
  # [../]
  # [./tang_force_x]
  # [../]
  # [./tang_force_y]
  # [../]
  [./von_mises]
    #Dependent variable used to visualize the Von Mises stress
    order = CONSTANT
    family = MONOMIAL
  [../]
[]

[Functions]
  [./disp_ramp_vert]
    type = PiecewiseLinear
    x = '0. 1. 3.5'
    y = '0. -0.2 -0.2'
  [../]
  [./disp_ramp_horz]
    type = PiecewiseLinear
    x = '0. 1. 3.5'
    y = '0. 0.0 0.0014'
  [../]
[]

[Kernels]
  [./TensorMechanics]
    use_displaced_mesh = true
    # save_in = 'saved_x saved_y'
    extra_vector_tags = 'ref'
    block = '1 2 3 4 5 6 7'
    strain = FINITE
  [../]
[]

[AuxKernels]
  [./stress_xx]
    type = RankTwoAux
    rank_two_tensor = stress
    variable = stress_xx
    index_i = 0
    index_j = 0
    execute_on = timestep_end
    block = '1 2 3 4 5 6 7'
  [../]
  [./stress_yy]
    type = RankTwoAux
    rank_two_tensor = stress
    variable = stress_yy
    index_i = 1
    index_j = 1
    execute_on = timestep_end
    block = '1 2 3 4 5 6 7'
  [../]
  [./stress_xy]
    type = RankTwoAux
    rank_two_tensor = stress
    variable = stress_xy
    index_i = 0
    index_j = 1
    execute_on = timestep_end
    block = '1 2 3 4 5 6 7'
  [../]
  # [./inc_slip_x]
  #   type = PenetrationAux
  #   variable = inc_slip_x
  #   execute_on = timestep_end
  #   boundary = 3
  #   paired_boundary = 2
  # [../]
  # [./inc_slip_y]
  #   type = PenetrationAux
  #   variable = inc_slip_y
  #   execute_on = timestep_end
  #   boundary = 3
  #   paired_boundary = 2
  # [../]
  # [./accum_slip_x]
  #   type = PenetrationAux
  #   variable = accum_slip_x
  #   execute_on = timestep_end
  #   boundary = 3
  #   paired_boundary = 2
  # [../]
  # [./accum_slip_y]
  #   type = PenetrationAux
  #   variable = accum_slip_y
  #   execute_on = timestep_end
  #   boundary = 3
  #   paired_boundary = 2
  # [../]
  # [./tang_force_x]
  #   type = PenetrationAux
  #   variable = tang_force_x
  #   quantity = tangential_force_x
  #   boundary = 3
  #   paired_boundary = 2
  # [../]
  # [./tang_force_y]
  #   type = PenetrationAux
  #   variable = tang_force_y
  #   quantity = tangential_force_y
  #   boundary = 3
  #   paired_boundary = 2
  # [../]
  # [./penetration]
  #   type = PenetrationAux
  #   variable = penetration
  #   boundary = 3
  #   paired_boundary = 2
  # [../]
  [./von_mises_kernel]
    #Calculates the von mises stress and assigns it to von_mises
    type = RankTwoScalarAux
    variable = von_mises
    rank_two_tensor = stress
    execute_on = timestep_end
    scalar_type = VonMisesStress
    block = '1 2 3 4 5 6 7'
  [../]
[]

[Postprocessors]
  [./bot_react_x]
    type = NodalSum
    variable = saved_x
    boundary = 1
  [../]
  [./bot_react_y]
    type = NodalSum
    variable = saved_y
    boundary = 1
  [../]
  [./top_react_x]
    type = NodalSum
    variable = saved_x
    boundary = 4
  [../]
  [./top_react_y]
    type = NodalSum
    variable = saved_y
    boundary = 4
  [../]
  # [./disp_x226]
  #   type = NodalVariableValue
  #   nodeid = 225
  #   variable = disp_x
  # [../]
  # [./disp_y226]
  #   type = NodalVariableValue
  #   nodeid = 225
  #   variable = disp_y
  # [../]
  [./_dt]
    type = TimestepSize
  [../]
  [./num_lin_it]
    type = NumLinearIterations
  [../]
  [./num_nonlin_it]
    type = NumNonlinearIterations
  [../]
  [./total_lin_it]
    type = CumulativeValuePostprocessor
    postprocessor = num_lin_it
  [../]
  [./total_nonlin_it]
    type = CumulativeValuePostprocessor
    postprocessor = num_nonlin_it
  [../]
[]

[BCs]
  [./side_x]
    type = DirichletBC
    variable = disp_y
    boundary = '1 2'
    value = 0.0
  [../]
  [./bot_y]
    type = DirichletBC
    variable = disp_x
    boundary = '1 2'
    value = 0.0
  [../]
  [./top_y_disp]
    type = FunctionDirichletBC
    variable = disp_y
    boundary = 4
    function = disp_ramp_vert
  [../]
  [./top_x_disp]
    type = DirichletBC
    variable = disp_x
    boundary = 4
    value = 0
  [../]
[]

[Materials]
  [./stuff1_elas_tens]
    type = ComputeIsotropicElasticityTensor
    block = '1'
    youngs_modulus = 1e10
    poissons_ratio = 0.0
  [../]
  [./stuff1_strain]
    type = ComputeFiniteStrain
    block = '1'
  [../]
  [./stuff1_stress]
    type = ComputeFiniteStrainElasticStress
    block = '1'
  [../]
  [./stuff2_elas_tens]
    type = ComputeIsotropicElasticityTensor
    block = '2 3 4 5 6 7'
    youngs_modulus = 1e6
    poissons_ratio = 0.3
  [../]
  [./stuff2_strain]
    type = ComputeFiniteStrain
    block = '2 3 4 5 6 7'
  [../]
  [./stuff2_stress]
    type = ComputeFiniteStrainElasticStress
    block = '2 3 4 5 6 7'
  [../]
[]

[Executioner]
  type = Transient
  #Preconditioned JFNK (default)
  solve_type = 'PJFNK'

  petsc_options_iname = '-pc_type -pc_factor_mat_solver_package'
  petsc_options_value = 'lu     superlu_dist'

  line_search = 'none'

  nl_abs_tol = 1e-6
  nl_rel_tol = 1e-5
  l_max_its = 100
  nl_max_its = 200

  start_time = 0.0
  end_time = 3.5
  l_tol = 1e-3
  dt = 0.025
  dtmin = 1e-4
[]

[Preconditioning]
  [./SMP]
    type = SMP
    full = true
  [../]
[]

# [VectorPostprocessors]
#   [./x_disp]
#     type = NodalValueSampler
#     variable = disp_x
#     boundary = '3 4'
#     sort_by = id
#   [../]
#   [./y_disp]
#     type = NodalValueSampler
#     variable = disp_y
#     boundary = '3 4'
#     sort_by = id
#   [../]
#   [./cont_press]
#     type = NodalValueSampler
#     variable = frictionless_normal_lm
#     boundary = '3'
#     sort_by = id
#   [../]
# []

[Outputs]
  print_linear_residuals = true
  perf_graph = true
  exodus = true
  csv = true

  [./console]
    type = Console
    max_rows = 5
  [../]
  # [./chkfile]
  #   type = CSV
  #   show = 'x_disp y_disp cont_press'
  #   start_time = 0.9
  #   execute_vector_postprocessors_on = timestep_end
  # [../]
  # [./chkfile2]
  #   type = CSV
  #   show = 'bot_react_x bot_react_y top_react_x top_react_y'
  #   start_time = 0.9
  #   execute_vector_postprocessors_on = timestep_end
  # [../]
  # [./outfile]
  #   type = CSV
  #   delimiter = ' '
  #   execute_vector_postprocessors_on = none
  # [../]
[]

[Constraints]
  # All constraints below for mechanical contact (Mortar)
  [weighted_gap_lm]
    type = ComputeWeightedGapLMMechanicalContact
    primary_boundary = 2
    secondary_boundary = 3
    primary_subdomain = 10000
    secondary_subdomain = 10001
    variable = frictionless_normal_lm
    disp_x = disp_x
    disp_y = disp_y
    use_displaced_mesh = true
  []
  [ncp_lm]
    type = ApplyPenetrationConstraintLMMechanicalContact
    primary = 2
    secondary = 3
    variable = frictionless_normal_lm
    primary_variable = disp_x
    c = 1.0e8
    tangential_tolerance = 0.01
  []
  [x]
    type = NormalMortarMechanicalContact
    primary_boundary = '2'
    secondary_boundary = '3'
    primary_subdomain = '10000'
    secondary_subdomain = '10001'
    variable = frictionless_normal_lm
    secondary_variable = disp_x
    component = x
    use_displaced_mesh = true
    compute_lm_residuals = false
  []
  [y]
    type = NormalMortarMechanicalContact
    primary_boundary = '2'
    secondary_boundary = '3'
    primary_subdomain = '10000'
    secondary_subdomain = '10001'
    variable = frictionless_normal_lm
    secondary_variable = disp_y
    component = y
    use_displaced_mesh = true
    compute_lm_residuals = false
  []
[]


# [VectorPostprocessors]
#   [contact_post]
#     type = NodalValueSampler
#     variable = frictionless_normal_lm
#     boundary = '3'
#     sort_by = x
#   []
#   [disp_x]
#     type = NodalValueSampler
#     variable = disp_x
#     boundary = '3'
#     sort_by = x
#   []
#   [disp_y]
#     type = NodalValueSampler
#     variable = disp_y
#     boundary = '3'
#     sort_by = x
#   []
# []


# [Contact]
#   [./interface]
#     primary = 2
#     secondary = 3
#     normalize_penalty = true
#     tangential_tolerance = 1e-3
#     penalty = 1e+10
#   [../]
# []
