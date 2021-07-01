[Mesh]
  [file_mesh]
    type = FileMeshGenerator
    file = Dream3D_149_Crystals.e
  []
  [sidesets]
    type = SideSetsFromNormalsGenerator
    input = file_mesh
    normals = '-1 0 0
               1 0 0
               0 -1 0
               0 1 0
              0 0 -1
               0 0 1'
    fixed_normal = true
    new_boundary = 'left right bottom top back front'
  []
[]

[Variables]
  [disp_x]
  []
  [disp_y]
  []
  [disp_z]
  []
[]

[GlobalParams]
  volumetric_locking_correction = true
[]

[AuxVariables]

  [pk2_zz]
    order = CONSTANT
    family = MONOMIAL
  []

  [vonmises]
    order = CONSTANT
    family = MONOMIAL
  []

  [stress_yy]
    order = CONSTANT
    family = MONOMIAL
  []

  [e_yy]
    order = CONSTANT
    family = MONOMIAL
  []

  [fp_00]
    order = CONSTANT
    family = MONOMIAL
  []
  [fp_01]
    order = CONSTANT
    family = MONOMIAL
  []
  [fp_02]
    order = CONSTANT
    family = MONOMIAL
  []
  [fp_10]
    order = CONSTANT
    family = MONOMIAL
  []
  [fp_11]
    order = CONSTANT
    family = MONOMIAL
  []
  [fp_12]
    order = CONSTANT
    family = MONOMIAL
  []
  [fp_20]
    order = CONSTANT
    family = MONOMIAL
  []
  [fp_21]
    order = CONSTANT
    family = MONOMIAL
  []
  [fp_22]
    order = CONSTANT
    family = MONOMIAL
  []

  [euler1]
    order = CONSTANT
    family = MONOMIAL
  []
  [euler2]
    order = CONSTANT
    family = MONOMIAL
  []
  [euler3]
    order = CONSTANT
    family = MONOMIAL
  []
  [stress_zz]
    order = CONSTANT
    family = MONOMIAL
  []
  [e_zz]
    order = CONSTANT
    family = MONOMIAL
  []

  [rot00]
    order = CONSTANT
    family = MONOMIAL
  []
  [rot01]
    order = CONSTANT
    family = MONOMIAL
  []
  [rot02]
    order = CONSTANT
    family = MONOMIAL
  []
  [rot10]
    order = CONSTANT
    family = MONOMIAL
  []
  [rot11]
    order = CONSTANT
    family = MONOMIAL
  []
  [rot12]
    order = CONSTANT
    family = MONOMIAL
  []
  [rot20]
    order = CONSTANT
    family = MONOMIAL
  []
  [rot21]
    order = CONSTANT
    family = MONOMIAL
  []
  [rot22]
    order = CONSTANT
    family = MONOMIAL
  []

  [crysrot00]
    order = CONSTANT
    family = MONOMIAL
  []
  [crysrot01]
    order = CONSTANT
    family = MONOMIAL
  []
  [crysrot02]
    order = CONSTANT
    family = MONOMIAL
  []
  [crysrot10]
    order = CONSTANT
    family = MONOMIAL
  []
  [crysrot11]
    order = CONSTANT
    family = MONOMIAL
  []
  [crysrot12]
    order = CONSTANT
    family = MONOMIAL
  []
  [crysrot20]
    order = CONSTANT
    family = MONOMIAL
  []
  [crysrot21]
    order = CONSTANT
    family = MONOMIAL
  []
  [crysrot22]
    order = CONSTANT
    family = MONOMIAL
  []
[]

[Functions]
  [tdisp]
    type = ParsedFunction
    value = 0.001*t
  []
[]

[UserObjects]
  [prop_read]
    type = ElementPropertyReadFile
    prop_file_name = 'HexagonalSingleEquiaxed_EulerAngle.txt'
    # Enter file data as prop#1, prop#2, .., prop#nprop
    nprop = 3
    execute_on = timestep_end
    read_type = block
    nblock = 150
  []
[]

[AuxKernels]

  [vonmises]
    type = RankTwoScalarAux
    rank_two_tensor = stress
    variable = vonmises
    scalar_type = VonMisesStress
    execute_on = timestep_end
  []

  [stress_yy]
    type = RankTwoAux
    variable = stress_yy
    rank_two_tensor = stress
    index_j = 1
    index_i = 1
    execute_on = timestep_end
  []
  [e_yy]
    type = RankTwoAux
    variable = e_yy
    rank_two_tensor = total_lagrangian_strain
    index_j = 1
    index_i = 1
    execute_on = timestep_end
  []
  [stress_zz]
    type = RankTwoAux
    variable = stress_zz
    rank_two_tensor = stress
    index_j = 2
    index_i = 2
    execute_on = timestep_end
  []
  [e_zz]
    type = RankTwoAux
    variable = e_zz
    rank_two_tensor = total_lagrangian_strain
    index_j = 2
    index_i = 2
    execute_on = timestep_end
  []

  [fp_00]
    type = RankTwoAux
    variable = fp_00
    rank_two_tensor = plastic_deformation_gradient
    index_j = 0
    index_i = 0
    execute_on = timestep_end
  []

  [fp_01]
    type = RankTwoAux
    variable = fp_01
    rank_two_tensor = plastic_deformation_gradient
    index_j = 1
    index_i = 0
    execute_on = timestep_end
  []

  [fp_02]
    type = RankTwoAux
    variable = fp_02
    rank_two_tensor = plastic_deformation_gradient
    index_j = 2
    index_i = 0
    execute_on = timestep_end
  []

  [fp_10]
    type = RankTwoAux
    variable = fp_10
    rank_two_tensor = plastic_deformation_gradient
    index_j = 0
    index_i = 1
    execute_on = timestep_end
  []

  [fp_11]
    type = RankTwoAux
    variable = fp_11
    rank_two_tensor = plastic_deformation_gradient
    index_j = 1
    index_i = 1
    execute_on = timestep_end
  []

  [fp_12]
    type = RankTwoAux
    variable = fp_12
    rank_two_tensor = plastic_deformation_gradient
    index_j = 2
    index_i = 1
    execute_on = timestep_end
  []

  [fp_20]
    type = RankTwoAux
    variable = fp_20
    rank_two_tensor = plastic_deformation_gradient
    index_j = 0
    index_i = 2
    execute_on = timestep_end
  []

  [fp_21]
    type = RankTwoAux
    variable = fp_21
    rank_two_tensor = plastic_deformation_gradient
    index_j = 1
    index_i = 2
    execute_on = timestep_end
  []

  [fp_22]
    type = RankTwoAux
    variable = fp_22
    rank_two_tensor = plastic_deformation_gradient
    index_j = 2
    index_i = 2
    execute_on = timestep_end
  []

  [euler1]
    type = MaterialRealVectorValueAux
    variable = euler1
    property = updated_Euler_angle
    component = 0
    execute_on = timestep_end
  []
  [euler2]
    type = MaterialRealVectorValueAux
    variable = euler2
    property = updated_Euler_angle
    component = 1
    execute_on = timestep_end
  []
  [euler3]
    type = MaterialRealVectorValueAux
    variable = euler3
    property = updated_Euler_angle
    component = 2
    execute_on = timestep_end
  []
  [pk2_zz]
    type = RankTwoAux
    variable = pk2_zz
    rank_two_tensor = second_piola_kirchhoff_stress
    index_j = 2
    index_i = 2
    execute_on = timestep_end
  []

  [rot00]
    type = RankTwoAux
    variable = rot00
    rank_two_tensor = update_rotation
    index_j = 0
    index_i = 0
    execute_on = timestep_end
  []
  [rot01]
    type = RankTwoAux
    variable = rot01
    rank_two_tensor = update_rotation
    index_j = 1
    index_i = 0
    execute_on = timestep_end
  []
  [rot02]
    type = RankTwoAux
    variable = rot02
    rank_two_tensor = update_rotation
    index_j = 2
    index_i = 0
    execute_on = timestep_end
  []
  [rot10]
    type = RankTwoAux
    variable = rot10
    rank_two_tensor = update_rotation
    index_j = 0
    index_i = 1
    execute_on = timestep_end
  []
  [rot11]
    type = RankTwoAux
    variable = rot11
    rank_two_tensor = update_rotation
    index_j = 1
    index_i = 1
    execute_on = timestep_end
  []
  [rot12]
    type = RankTwoAux
    variable = rot12
    rank_two_tensor = update_rotation
    index_j = 2
    index_i = 1
    execute_on = timestep_end
  []
  [rot20]
    type = RankTwoAux
    variable = rot20
    rank_two_tensor = update_rotation
    index_j = 0
    index_i = 2
    execute_on = timestep_end
  []
  [rot21]
    type = RankTwoAux
    variable = rot21
    rank_two_tensor = update_rotation
    index_j = 1
    index_i = 2
    execute_on = timestep_end
  []
  [rot22]
    type = RankTwoAux
    variable = rot22
    rank_two_tensor = update_rotation
    index_j = 2
    index_i = 2
    execute_on = timestep_end
  []

  [crysrot00]
    type = RankTwoAux
    variable = crysrot00
    rank_two_tensor = crysrot
    index_j = 0
    index_i = 0
    execute_on = timestep_end
  []
  [crysrot01]
    type = RankTwoAux
    variable = crysrot01
    rank_two_tensor = crysrot
    index_j = 1
    index_i = 0
    execute_on = timestep_end
  []
  [crysrot02]
    type = RankTwoAux
    variable = crysrot02
    rank_two_tensor = crysrot
    index_j = 2
    index_i = 0
    execute_on = timestep_end
  []
  [crysrot10]
    type = RankTwoAux
    variable = crysrot10
    rank_two_tensor = crysrot
    index_j = 0
    index_i = 1
    execute_on = timestep_end
  []
  [crysrot11]
    type = RankTwoAux
    variable = crysrot11
    rank_two_tensor = crysrot
    index_j = 1
    index_i = 1
    execute_on = timestep_end
  []
  [crysrot12]
    type = RankTwoAux
    variable = crysrot12
    rank_two_tensor = crysrot
    index_j = 2
    index_i = 1
    execute_on = timestep_end
  []
  [crysrot20]
    type = RankTwoAux
    variable = crysrot20
    rank_two_tensor = crysrot
    index_j = 0
    index_i = 2
    execute_on = timestep_end
  []
  [crysrot21]
    type = RankTwoAux
    variable = crysrot21
    rank_two_tensor = crysrot
    index_j = 1
    index_i = 2
    execute_on = timestep_end
  []
  [crysrot22]
    type = RankTwoAux
    variable = crysrot22
    rank_two_tensor = crysrot
    index_j = 2
    index_i = 2
    execute_on = timestep_end
  []
[]

[BCs]

  # [Periodic]
  #   [all]
  #     variable = 'disp_x'
  #     auto_direction = 'z'
  #   []
  # []

  [fix_x]
    type = DirichletBC
    variable = disp_x
    boundary = 'left'
    value = 0
  []

  [fix_y]
    type = DirichletBC
    variable = disp_y
    boundary = 'bottom'
    value = 0
  []

  [fix_z]
    type = DirichletBC
    variable = disp_z
    boundary = 'back'
    value = 0
  []

  [tdisp]
    type = FunctionDirichletBC
    variable = disp_z
    boundary = 'front'
    function = tdisp
  []
[]

[Materials]
  [strain]
    type = ComputeFiniteStrain
    displacements = 'disp_x disp_y disp_z'
  []
  [elasticity_tensor]
    type = ComputeElasticityTensorCP
    C_ijkl = '1.0675e5 0.6041e5 0.6041e5 1.0675e5 0.6041e5 1.0675e5 0.2834e5 0.2834e5 0.2834e5'
    fill_method = symmetric9
    read_prop_user_object = prop_read
  []
  [stress]
    type = ComputeMultipleCrystalPlasticityStress
    crystal_plasticity_models = 'trial_xtalpl'
    tan_mod_type = none
    maximum_substep_iteration = 5
  []
  [trial_xtalpl]
    type = CrystalPlasticityKalidindiUpdate
    stol = 1e-2
    number_slip_systems = 12
    slip_sys_file_name = slip_sys_fcc.txt
    ao = 0.001
    xm = 0.05
    gss_initial = 31
    r = 1.4
    h = 75
    t_sat = 63
    gss_a = 2.25
  []
  [euler_angle_mag]
    type = UpdatedEulerAngle
  []
[]

[Postprocessors]

  # [stress_zz]
  #   type = ElementAverageValue
  #   variable = stress_zz
  # []
  # [e_zz]
  #   type = ElementAverageValue
  #   variable = e_zz
  # []
  # [vonMises]
  #   type = ElementAverageValue
  #   variable = vonmises
  # []
  # [stress_yy]
  #   type = ElementAverageValue
  #   variable = stress_yy
  # []
  # [e_yy]
  #   type = ElementAverageValue
  #   variable = e_yy
  # []
  #
  # [fp_22]
  #   type = ElementAverageValue
  #   variable = fp_22
  # []
  #
  # [pk2_zz]
  #   type = ElementAverageValue
  #   variable = pk2_zz
  # []
  #
  # [rot22]
  #   type = ElementAverageValue
  #   variable = rot22
  # []
  #
  # [crysrot22]
  #   type = ElementAverageValue
  #   variable = crysrot22
  # []

  [euler1]
    type = ElementAverageValue
    variable = euler1
    execute_on = 'TIMESTEP_END'
  []
  [euler2]
    type = ElementAverageValue
    variable = euler2
    execute_on = 'TIMESTEP_END'
  []
  [euler3]
    type = ElementAverageValue
    variable = euler3
    execute_on = 'TIMESTEP_END'
  []
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

  petsc_options_iname = '-pc_type -pc_factor_mat_solver_package -pc_factor_shift_type -pc_factor_shift_amount'
  petsc_options_value = 'lu    superlu_dist   NONZERO   1e-10'
  line_search = 'none'

  nl_abs_tol = 1e-6
  nl_rel_tol = 1e-6

  l_max_its = 50
  l_tol = 1e-3

  dt = 10
  dtmin = 1e-4
  end_time = 5000
[]

[Outputs]
  file_base = 'cp_mat_based_out'
  exodus = true
  csv = true
  interval = 5
  checkpoint = true
[]

[Kernels]
  [TensorMechanics]
    displacements = 'disp_x disp_y disp_z'
    use_displaced_mesh = true
  []
[]
