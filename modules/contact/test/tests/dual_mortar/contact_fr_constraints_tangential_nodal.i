offset = 0.0
vy = 0.1

max_lx = 0.1

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
    nx = 4
    ny = 4
    elem_type = QUAD4
    boundary_id_offset = 10
    boundary_name_prefix = left
  [../]
  [./left_block_id]
    type = SubdomainIDGenerator
    input = left_block
    subdomain_id = 1
  [../]
  [./right_block]
    type = GeneratedMeshGenerator
    dim = 2
    xmin = 0
    xmax = 1
    ymin = -1
    ymax = 1
    nx = 6
    ny = 12
    elem_type = QUAD4
    boundary_id_offset = 20
    boundary_name_prefix = right
  [../]
  [./right_block_id]
    type = SubdomainIDGenerator
    input = right_block
    subdomain_id = 2
  [../]

  [./combined_mesh]
    type = MeshCollectionGenerator
    inputs = 'left_block_id right_block_id'
  [../]

  [./secondary]
    input = combined_mesh
    type = LowerDBlockFromSidesetGenerator
    sidesets = 'left_right'
    new_block_id = '3'
    new_block_name = 'secondary_lower'
  [../]
  [./primary]
    input = secondary
    type = LowerDBlockFromSidesetGenerator
    sidesets = 'right_left'
    new_block_id = '4'
    new_block_name = 'primary_lower'
  [../]
[]

[Variables]
  [./disp_x]
    block = '1 2'
    family = LAGRANGE
    order = FIRST
  [../]
  [./disp_y]
    block = '1 2'
    family = LAGRANGE
    order = FIRST
  [../]
  [./normal_lm]
    block = 'secondary_lower'
    # family = MONOMIAL
    # order = CONSTANT
  [../]
  [./tangential_lm]
    block = 'secondary_lower'
    # family = MONOMIAL
    # order = CONSTANT
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

[Constraints]
  # All constraints below for mechanical contact (Mortar)
  [normal_lm]
    type = NormalNodalLMMechanicalContact
    primary = right_left
    secondary = left_right
    variable = normal_lm
    primary_variable = disp_x
    disp_y = disp_y
    use_displaced_mesh = true
    ncp_function_type = min
  []
  [normal_x]
    type = NormalMortarMechanicalContact
    primary_boundary = '23'
    secondary_boundary = '11'
    primary_subdomain = 4
    secondary_subdomain = 3
    variable = normal_lm
    secondary_variable = disp_x
    component = x
    use_displaced_mesh = true
    compute_lm_residuals = false
  []
  [normal_y]
    type = NormalMortarMechanicalContact
    primary_boundary = 23
    secondary_boundary = 11
    primary_subdomain = 4
    secondary_subdomain = 3
    variable = normal_lm
    secondary_variable = disp_y
    component = y
    use_displaced_mesh = true
    compute_lm_residuals = false
  []
  [tangential_lm] #nodal
    type = TangentialNodalLMMechanicalContact
    primary = right_left
    secondary = left_right
    variable = tangential_lm
    primary_variable = disp_x
    disp_y = disp_y
    contact_pressure = normal_lm
    ncp_function_type = fb
    mu = 0.3
  []
  [tangential_x]
    type = TangentialMortarMechanicalContact
    primary_boundary = 23
    secondary_boundary = 11
    primary_subdomain = 4
    secondary_subdomain = 3
    variable = tangential_lm
    secondary_variable = disp_x
    component = x
    use_displaced_mesh = true
    compute_lm_residuals = false
  []
  [tangential_y]
    type = TangentialMortarMechanicalContact
    primary_boundary = 23
    secondary_boundary = 11
    primary_subdomain = 4
    secondary_subdomain = 3
    variable = tangential_lm
    secondary_variable = disp_y
    component = y
    use_displaced_mesh = true
    compute_lm_residuals = false
  []
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

  petsc_options = '-snes_converged_reason -ksp_converged_reason -snes_view -snes_test_jacobian'

  petsc_options_iname = '-pc_type -pc_factor_shift_type -pc_factor_shift_amount'
  petsc_options_value = 'lu NONZERO   1e-15'

  dt = 0.1
  dtmin = 1e-4
  end_time = 0.5

  l_max_its = 20

  nl_max_its = 8
  nl_rel_tol = 1e-6
  snesmf_reuse_base = false
[]

[Outputs]
  file_base = ./contact_fr_constraints_out
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
    variable = normal_lm
    subdomain = '3'
  []
  [./normal_lm]
    type = ElementAverageValue
    variable = normal_lm
    block = '3'
  [../]
  [./tangential_lm]
    type = ElementAverageValue
    variable = tangential_lm
    block = '3'
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
