# UserObject J2 test, with hardening, but with rate=0
# apply uniform compression in x direction to give
# trial stress_xx = -5, so sqrt(3*J2) = 5
# with zero Poisson's ratio, this should return to
# stress_xx = -3, stress_yy = -1 = stress_zz,
# for strength = 2
# (note that stress_xx - stress_yy = stress_xx - stress_zz = -2, so sqrt(3*j2) = 2,
#  and that the mean stress remains = -5/3)

[Mesh]
  type = GeneratedMesh
  dim = 2
  nx = 20
  ny = 20
  nz = 0
  xmax = 1000
  ymax = 1000
  zmax = 0
  # elem_type = QUAD4
  # uniform_refine = 2
[]


[Variables]
  [./disp_x]
  [../]
  [./disp_y]
  [../]
  # [./disp_z]
  # [../]
[]

[Kernels]
  [./TensorMechanics]
    displacements = 'disp_x disp_y'
  [../]
[]

[BCs]
  [./top_displacement]
    type = DirichletBC
    variable = disp_y
    boundary = 'top'
    value = 50.0
  [../]
  [./x_anchor]
    type = DirichletBC
    variable = disp_x
    boundary = 'left right'
    value = 0.0
  [../]
  [./y_anchor]
    type = DirichletBC
    variable = disp_y
    boundary = bottom
    value = 0.0
  [../]
[]

[AuxVariables]
  [./stress_xx]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./stress_xy]
    order = CONSTANT
    family = MONOMIAL
  [../]
  # [./stress_xz]
  #   order = CONSTANT
  #   family = MONOMIAL
  # [../]
  [./stress_yy]
    order = CONSTANT
    family = MONOMIAL
  [../]
  # [./stress_yz]
  #   order = CONSTANT
  #   family = MONOMIAL
  # [../]
  # [./stress_zz]
  #   order = CONSTANT
  #   family = MONOMIAL
  # [../]
  [./f]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./iter]
    order = CONSTANT
    family = MONOMIAL
  [../]
[]

[AuxKernels]
  [./stress_xx]
    type = RankTwoAux
    rank_two_tensor = stress
    variable = stress_xx
    index_i = 0
    index_j = 0
  [../]
  [./stress_xy]
    type = RankTwoAux
    rank_two_tensor = stress
    variable = stress_xy
    index_i = 0
    index_j = 1
  [../]
  # [./stress_xz]
  #   type = RankTwoAux
  #   rank_two_tensor = stress
  #   variable = stress_xz
  #   index_i = 0
  #   index_j = 2
  # [../]
  [./stress_yy]
    type = RankTwoAux
    rank_two_tensor = stress
    variable = stress_yy
    index_i = 1
    index_j = 1
  [../]
  # [./stress_yz]
  #   type = RankTwoAux
  #   rank_two_tensor = stress
  #   variable = stress_yz
  #   index_i = 1
  #   index_j = 2
  # [../]
  # [./stress_zz]
  #   type = RankTwoAux
  #   rank_two_tensor = stress
  #   variable = stress_zz
  #   index_i = 2
  #   index_j = 2
  # [../]
  [./f]
    type = MaterialStdVectorAux
    index = 0
    property = plastic_yield_function
    variable = f
  [../]
  [./iter]
    type = MaterialRealAux
    property = plastic_NR_iterations
    variable = iter
  [../]
[]

[Postprocessors]
  [./s_xx]
    type = PointValue
    point = '0 0 0'
    variable = stress_xx
  [../]
  [./s_xy]
    type = PointValue
    point = '0 0 0'
    variable = stress_xy
  [../]
  # [./s_xz]
  #   type = PointValue
  #   point = '0 0 0'
  #   variable = stress_xz
  # [../]
  [./s_yy]
    type = PointValue
    point = '0 0 0'
    variable = stress_yy
  [../]
  # [./s_yz]
  #   type = PointValue
  #   point = '0 0 0'
  #   variable = stress_yz
  # [../]
  # [./s_zz]
  #   type = PointValue
  #   point = '0 0 0'
  #   variable = stress_zz
  # [../]
  [./f]
    type = PointValue
    point = '0 0 0'
    variable = f
  [../]
  [./iter]
    type = PointValue
    point = '0 0 0'
    variable = iter
  [../]
[]

[UserObjects]
  [./str]
    type = TensorMechanicsHardeningConstant
    value = 2
  [../]
  [./j2]
    type = TensorMechanicsPlasticJ2
    yield_strength = str
    yield_function_tolerance = 1E-3
    internal_constraint_tolerance = 1E-9
  [../]
[]

[Materials]
  [./elasticity_tensor]
    type = ComputeElasticityTensor
    block = 0
    fill_method = symmetric_isotropic
    C_ijkl = '0 1E6'
  [../]
  [./strain]
    type = ComputeFiniteStrain
    block = 0
    displacements = 'disp_x disp_y'
  [../]
  [./mc]
    type = ComputeMultiPlasticityStress
    block = 0
    ep_plastic_tolerance = 1E-9
    plastic_models = j2
    debug_fspb = crash
  [../]
[]


[Executioner]
  type = Transient
  scheme = bdf2
  solve_type = PJFNK
  # line_search = none
  petsc_options_iname = '-pc_type -pc_hypre_type -ksp_gmres_restart -pc_hypre_boomeramg_strong_threshold'
  petsc_options_value = 'hypre boomeramg 31 0.7'
  l_tol = 1.0e-3
  l_max_its = 30
  nl_max_its = 25
  nl_rel_tol = 1.0e-6
  nl_abs_tol = 1e-8
  start_time = 0.0
  num_steps = 50
  automatic_scaling = true
  # [./TimeStepper]
  #   type = IterationAdaptiveDT
  #   dt = 1.5
  #   growth_factor = 1.2
  #   cutback_factor = 0.8
  #   optimal_iterations = 8
  # [../]
  # [./Adaptivity]
  #   initial_adaptivity = 2
  #   refine_fraction = 0.8
  #   coarsen_fraction = 0.05
  #   max_h_level = 3
  # [../]
[]


[Outputs]
  # file_base = hard1
  exodus = true
  [./csv]
    type = CSV
    [../]
[]
