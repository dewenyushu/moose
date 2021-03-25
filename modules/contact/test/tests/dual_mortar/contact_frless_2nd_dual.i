offset = 0.0
vy = 0.1

max_lx = 0.08

[GlobalParams]
  displacements = 'disp_x disp_y'
  volumetric_locking_correction = true
[]

[Mesh]
  [./left_block]
    type = GeneratedMeshGenerator
    dim = 2
    xmin = -1
    xmax = 0
    ymin = -1
    ymax = 0
    nx = 2
    ny = 2
    elem_type = QUAD9
  [../]
  [./left_block_sidesets]
    type = RenameBoundaryGenerator
    input = left_block
    old_boundary = '0 1 2 3'
    new_boundary = '10 11 12 13'
  [../]
  [./left_block_id]
    type = SubdomainIDGenerator
    input = left_block_sidesets
    subdomain_id = 1
  [../]
  [./right_block]
    type = GeneratedMeshGenerator
    dim = 2
    xmin = 0
    xmax = 1
    ymin = -1
    ymax = 1
    nx = 2
    ny = 4
    elem_type = QUAD9
  [../]
  [right_block_sidesets]
    type = RenameBoundaryGenerator
    input = right_block
    old_boundary = '0 1 2 3'
    new_boundary = '20 21 22 23'
  []
  [./right_block_id]
    type = SubdomainIDGenerator
    input = right_block_sidesets
    subdomain_id = 2
  [../]

  [./combined_mesh]
    type = MeshCollectionGenerator
    inputs = 'left_block_id right_block_id'
  [../]
[]

[Variables]
  [./disp_x]
    block = '1 2'
    family = LAGRANGE
    order = SECOND
  [../]
  [./disp_y]
    block = '1 2'
    family = LAGRANGE
    order = SECOND
  [../]
[]

[Modules/TensorMechanics/Master]
  [./all]
    strain = FINITE
    incremental = true
    add_variables = true
    use_automatic_differentiation = true
    block = '1 2'
  [../]
[]

[Functions]
  [./horizontal_movement]
    type = ParsedFunction
    value = '${max_lx}-${max_lx}*exp(-t)'
  [../]
  [./vertical_movement]
    type = ParsedFunction
    value = '${vy}*t+${offset}'
  [../]
[]

[BCs]
  [./push_left_x]
    type = FunctionDirichletBC
    variable = disp_x
    boundary = 13
    function = horizontal_movement
  [../]
  [./fix_right_x]
    type = DirichletBC
    variable = disp_x
    boundary = 21
    value = 0.0
  [../]
  [./fix_right_y]
    type = DirichletBC
    variable = disp_y
    boundary = 21
    value = 0.0
  [../]
  [./push_left_y]
    type = FunctionDirichletBC
    variable = disp_y
    boundary = 13
    function = vertical_movement
  [../]
[]

[Materials]
  [./elasticity_tensor_left]
    type = ADComputeIsotropicElasticityTensor
    block = 1
    youngs_modulus = 1.0e3
    poissons_ratio = 0.3
  [../]
  [./stress_left]
    type = ADComputeFiniteStrainElasticStress
    block = 1
  [../]

  [./elasticity_tensor_right]
    type = ADComputeIsotropicElasticityTensor
    block = 2
    youngs_modulus = 1.0e3
    poissons_ratio = 0.3
  [../]
  [./stress_right]
    type = ADComputeFiniteStrainElasticStress
    block = 2
  [../]
[]

[Contact]
  [leftright]
    mesh = combined_mesh
    secondary = '11'
    primary = '23'

    formulation = mortar
    model = frictionless

    use_dual = true
  [../]
[]

[ICs]
  [./disp_y]
    block = 1
    variable = disp_y
    value = ${offset}
    type = ConstantIC
  [../]
  [./disp_x]
    block = 1
    variable = disp_x
    value = ${offset}
    type = ConstantIC
  [../]
[]

[Preconditioning]
  [./SMP]
    type = SMP
    full = true
  [../]
[]

[Executioner]
  type = Transient
  solve_type = 'NEWTON'

  petsc_options = '-snes_converged_reason -ksp_converged_reason -snes_view'

  petsc_options_iname = '-pc_type -pc_factor_shift_type -pc_factor_shift_amount'
  petsc_options_value = 'lu NONZERO   1e-15'

  dt = 0.1
  dtmin = 0.1
  end_time = 0.1

  l_max_its = 20

  nl_max_its = 10
  nl_rel_tol = 1e-6
  snesmf_reuse_base = false
[]

[Outputs]
  file_base = ./contact_frless_dual_out
  [./comp]
    type = CSV
  [../]
  [./exo]
    type = Exodus
  [../]
[]


[Postprocessors]
  [./total_nl]
    type = CumulativeValuePostprocessor
    postprocessor = nl
  [../]
  [./total_lin]
    type = CumulativeValuePostprocessor
    postprocessor = lin
  [../]
  [nl]
    type = NumNonlinearIterations
  []
  [lin]
    type = NumLinearIterations
  []
  [./contact]
    type = ContactDOFSetSize
    variable = leftright_normal_lm
    subdomain = leftright_secondary_subdomain
  []
  [./normal_lm]
    type = ElementAverageValue
    variable = leftright_normal_lm
    block = leftright_secondary_subdomain
  [../]
  [./avg_disp_x]
    type = ElementAverageValue
    variable = disp_x
    block = '1 2'
  [../]
  [./avg_disp_y]
    type = ElementAverageValue
    variable = disp_y
    block = '1 2'
  [../]
  [./max_disp_x]
    type = ElementExtremeValue
    variable = disp_x
    block = '1 2'
  [../]
  [./max_disp_y]
    type = ElementExtremeValue
    variable = disp_y
    block = '1 2'
  [../]
  [./min_disp_x]
    type = ElementExtremeValue
    variable = disp_x
    block = '1 2'
    value_type = min
  [../]
  [./min_disp_y]
    type = ElementExtremeValue
    variable = disp_y
    block = '1 2'
    value_type = min
  [../]
[]
