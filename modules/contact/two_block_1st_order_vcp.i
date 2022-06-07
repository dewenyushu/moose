[GlobalParams]
  displacements = 'disp_x disp_y'
  volumetric_locking_correction = true
[]

theta = 60
velocity = 0.15

[Mesh]
  [left_block]
    type = GeneratedMeshGenerator
    dim = 2
    xmin = -0.35
    xmax = -0.05
    ymin = -1
    ymax = 0
    nx = 1
    ny = 3
    elem_type = QUAD4
  []
  [left_block_sidesets]
    type = RenameBoundaryGenerator
    input = left_block
    old_boundary = '0 1 2 3'
    new_boundary = '10 11 12 13'
  []
  [left_block_sideset_names]
    type = RenameBoundaryGenerator
    input = left_block_sidesets
    old_boundary = '10 11 12 13'
    new_boundary = 'l_bottom l_right l_top l_left'
  []
  [left_block_id]
    type = SubdomainIDGenerator
    input = left_block_sideset_names
    subdomain_id = 1
  []
  [right_block]
    type = GeneratedMeshGenerator
    dim = 2
    xmin = 0
    xmax = 0.3
    ymin = -1
    ymax = 0
    nx = 1
    ny = 2
    elem_type = QUAD4
  []
  [right_block_sidesets]
    type = RenameBoundaryGenerator
    input = right_block
    old_boundary = '0 1 2 3'
    new_boundary = '20 21 22 23'
  []
  [right_block_sideset_names]
    type = RenameBoundaryGenerator
    input = right_block_sidesets
    old_boundary = '20 21 22 23'
    new_boundary = 'r_bottom r_right r_top r_left'
  []
  [right_block_id]
    type = SubdomainIDGenerator
    input = right_block_sideset_names
    subdomain_id = 2
  []

  [combined_mesh]
    type = MeshCollectionGenerator
    inputs = 'left_block_id right_block_id'
  []

  [rotate_mesh]
    type = TransformGenerator
    input = combined_mesh
    transform = ROTATE
    vector_value = '0 0 ${theta}'
  []
[]

[Modules/TensorMechanics/Master]
  [all]
    strain = FINITE
    incremental = true
    add_variables = true
    block = '1 2'
  []
[]

[Functions]
  [horizontal_movement]
    type = ParsedFunction
    value = '${velocity} * t * cos(${theta}/180*pi)'
  []
  [vertical_movement]
    type = ParsedFunction
    value = '${velocity} * t * sin(${theta}/180*pi)'
  []
[]

[BCs]
  [push_left_x]
    type = FunctionDirichletBC
    variable = disp_x
    boundary = 13
    function = horizontal_movement
  []
  [fix_right_x]
    type = DirichletBC
    variable = disp_x
    boundary = 21
    value = 0.0
  []
  [fix_right_y]
    type = DirichletBC
    variable = disp_y
    boundary = 21
    value = 0.0
  []
  [push_left_y]
    type = FunctionDirichletBC
    variable = disp_y
    boundary = 13
    function = vertical_movement
  []
[]

[Materials]
  [elasticity_tensor_left]
    type = ComputeIsotropicElasticityTensor
    block = 1
    youngs_modulus = 1.0e6
    poissons_ratio = 0.3
  []
  [stress_left]
    type = ComputeFiniteStrainElasticStress
    block = 1
  []

  [elasticity_tensor_right]
    type = ComputeIsotropicElasticityTensor
    block = 2
    youngs_modulus = 1.0e6
    poissons_ratio = 0.3
  []
  [stress_right]
    type = ComputeFiniteStrainElasticStress
    block = 2
  []
[]

[Contact]
  [leftright]
    mesh = combined_mesh
    secondary = '11'
    primary = '23'
    formulation = mortar
    model = frictionless
    correct_edge_dropping = true
    use_dual = true
  []
[]

[Preconditioning]
  [vcp]
    type = VCP
    full = true
    lm_variable = 'leftright_normal_lm'
    primary_variable = 'disp_x disp_y'
    preconditioner = 'LU'
    is_lm_coupling_diagonal = true
    adaptive_condensation = true
  []
[]

[Executioner]
  type = Transient
  solve_type = 'NEWTON'

  # petsc_options = '-snes_converged_reason -ksp_converged_reason -snes_view'

  petsc_options_iname = '  -mat_mffd_err -pc_factor_shift_type -pc_factor_shift_amount'
  petsc_options_value = '      1e-5          NONZERO               1e-10'

  line_search = none

  dt = 0.1
  dtmin = 0.1
  end_time = 0.4

  l_max_its = 100

  nl_max_its = 20
  nl_rel_tol = 1e-6
  snesmf_reuse_base = false
[]

[Outputs]
  exodus = true
  file_base = './output/1st_order_${theta}_degree_3-2_out'
  [comp]
    type = CSV
    show = 'tot_lin_it tot_nonlin_it'
    execute_on = 'FINAL'
  []
[]

[Postprocessors]
  [contact]
    type = ContactDOFSetSize
    variable = leftright_normal_lm
    subdomain = leftright_secondary_subdomain
  []
  [normal_lm]
    type = ElementAverageValue
    variable = leftright_normal_lm
    block = leftright_secondary_subdomain
  []
  [avg_disp_x]
    type = ElementAverageValue
    variable = disp_x
    block = '1 2'
  []
  [avg_disp_y]
    type = ElementAverageValue
    variable = disp_y
    block = '1 2'
  []
  [max_disp_x]
    type = ElementExtremeValue
    variable = disp_x
    block = '1 2'
  []
  [max_disp_y]
    type = ElementExtremeValue
    variable = disp_y
    block = '1 2'
  []
  [min_disp_x]
    type = ElementExtremeValue
    variable = disp_x
    block = '1 2'
    value_type = min
  []
  [min_disp_y]
    type = ElementExtremeValue
    variable = disp_y
    block = '1 2'
    value_type = min
  []
  [num_lin_it]
    type = NumLinearIterations
  []
  [num_nonlin_it]
    type = NumNonlinearIterations
  []
  [tot_lin_it]
    type = CumulativeValuePostprocessor
    postprocessor = num_lin_it
  []
  [tot_nonlin_it]
    type = CumulativeValuePostprocessor
    postprocessor = num_nonlin_it
  []
[]
