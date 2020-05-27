refine = 0

# [GlobalParams]
#   displacements = 'disp_x disp_y'
#   volumetric_locking_correction = true
# []

[Mesh]
  second_order = true
  [./original_file_mesh]
    type = FileMeshGenerator
    file = ddm_mortar_non_conform_2blocks.e
    # file = ddm_mortar_2blocks_fine_slave.e
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
  uniform_refine =  ${refine}
[]

[Functions]
  [./exact_sln]
    type = ParsedFunction
    # value = x*x+y*y
    value = sin(2*pi*x)*sin(2*pi*y)
  [../]
  [./exact_sln_lambda]
    type = ParsedFunction
    value = -2*pi*sin(2*pi*y)
  [../]
  [./ffn]
    type = ParsedFunction
    # value = -4
    value = 8*pi*pi*sin(2*pi*x)*sin(2*pi*y)
  [../]
[]

[Variables]
  [./u]
    order = SECOND
    family = LAGRANGE
    block = '1 2'
  [../]

  [./lm]
    order = SECOND
    family = LAGRANGE
    block = slave_lower
    use_dual = false
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
    slave_variable = u
    master_boundary = 20
    master_subdomain = 200
    slave_boundary = 10
    slave_subdomain = 100
  [../]
[]

[BCs]
  [./all]
    type = FunctionDirichletBC
    variable = u
    boundary = '30 40 50 60'
    function = exact_sln
  [../]
[]

[Postprocessors]
  [./l2_u]
    type = ElementL2Error
    variable = u
    function = exact_sln
    block = '1 2'
    execute_on = 'initial timestep_end'
  [../]
  [./l2_lm]
    type = ElementL2Error
    variable = lm
    function = exact_sln_lambda
    block = '100'
    execute_on = 'initial timestep_end'
  [../]
  [h]
     type = AverageElementSize
     execute_on = 'initial timestep_end'
     block = '1 2'
  []
[]

[Preconditioning]
  [./fmp]
    type = SMP
    full = true
  [../]
[]

[Executioner]
  type = Steady
  solve_type = 'NEWTON'

  petsc_options = '-snes_converged_reason -ksp_converged_reason -snes_ksp_ew'
  petsc_options_iname = '-pc_type -pc_factor_shift_type -pc_factor_shift_amount'
  petsc_options_value = 'lu NONZERO   1e-15'

  # petsc_options_iname = '-pc_type -ksp_view_mat'
 	# petsc_options_value = 'svd ascii:matrix_ref${refine}.m:ascii_matlab'

  nl_rel_tol = 1e-12
  l_tol = 1e-12
[]

[Outputs]
  [./exodus]
    file_base = DDM_LU_exodus_ref${refine}
    type = Exodus
  [../]
  # [dof_map]
  #   file_base = DDM_LU_dofmap_ref${refine}
  #   type = DOFMap
  #   execute_on = 'initial'
  # []
  # [./csv]
  #   file_base = DDM_LU_exodus_ref${refine}
  #   type = CSV
  # [../]
[]
