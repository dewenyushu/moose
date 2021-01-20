[Mesh]
  [./original_file_mesh]
    type = FileMeshGenerator
    file = non_conform_2blocks.e
    # file = dmp_2blocks_coarse.e
  [../]
  [slave]
    input = original_file_mesh
    type = LowerDBlockFromSidesetGenerator
    sidesets = '10'
    new_block_id = '100'
    new_block_name = 'slave_lower'
  []
  [master]
    input = slave
    type = LowerDBlockFromSidesetGenerator
    sidesets = '20'
    new_block_id = '200'
    new_block_name = 'master_lower'
  []
[]

[Functions]
  [./exact_sln]
    type = ParsedFunction
    value = sin(2*pi*x)*sin(2*pi*y)
  [../]
  [./ffn]
    type = ParsedFunction
    value = 8*pi*pi*sin(2*pi*x)*sin(2*pi*y)
  [../]
[]

[Variables]
  [./u]
    order = FIRST
    family = LAGRANGE
    block = '1 2'
  [../]

  [./lm]
    order = FIRST
    family = LAGRANGE
    block = slave_lower
    use_dual = true
  [../]
[]

[Kernels]
  [./diff]
    type = Diffusion
    variable = u
  [../]
  [./ffn]
    type = BodyForce
    variable = u
    function = ffn
  [../]
[]

[Constraints]
  [./ced]
    type = EqualValueConstraint
    variable = lm
    secondary_variable = u
    primary_boundary = 20
    primary_subdomain = 200
    secondary_boundary = 10
    secondary_subdomain = 100
  [../]
[]

[BCs]
  [./all]
    type = DirichletBC
    variable = u
    boundary = '30 40'
    value = 0.0
  [../]
  [./neumann]
    type = FunctionGradientNeumannBC
    exact_solution = exact_sln
    variable = u
    boundary = '50 60'
  [../]
[]

[Postprocessors]
  [./l2_error]
    type = ElementL2Error
    variable = u
    function = exact_sln
    block = '1 2'
    execute_on = 'initial timestep_end'
  [../]
[]

[Preconditioning]
  [./pdmp]
    type = PDMP
    full = true
    variable = 'lm'
    coupled_variable = 'u'
    secondary_subdomain = 1
    secondary_boundary = 10
    primary_subdomain = 2
    primary_boundary = 20
    preconditioner = 'AMG'
  [../]
[]


[Executioner]
  type = Steady
  solve_type = 'NEWTON'

  petsc_options = '-snes_converged_reason -ksp_converged_reason -snes_view'
  # petsc_options_iname = '-pc_factor_shift_type -pc_factor_shift_amount -pc_factor_mat_solver_type'
  # petsc_options_value = 'nonzero 1e-8 superlu_dist'

  # petsc_options = '-snes_converged_reason -ksp_converged_reason -ksp_view_mat -snes_view' # -pc_svd_monitor'
  # petsc_options_iname = '-pc_factor_shift_type -pc_factor_shift_amount -pc_factor_mat_solver_type'
  # petsc_options_value = 'nonzero 1e-8 superlu_dist'

  # petsc_options_iname = '-pc_type -pc_factor_shift_type -pc_factor_shift_amount'
  # petsc_options_value = 'lu nonzero 1e-8'

  # petsc_options_iname = ' -pc_factor_shift_type -pc_factor_shift_amount'
  # petsc_options_value = ' nonzero 1e-8'

  l_max_its = 10
  nl_max_its = 1

  nl_rel_tol = 1e-6
  l_tol = 1e-12
[]

[Outputs]
  csv = true

  # [dofmap]
  #   type = DOFMap
  #   execute_on = 'initial'
  # []
[]
