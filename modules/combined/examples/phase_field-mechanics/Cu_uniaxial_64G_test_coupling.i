[GlobalParams]
  displacements = 'disp_x disp_y disp_z'
[]

[Mesh]
  displacements = 'disp_x disp_y disp_z'
  [fmg]
    type = GeneratedMeshGenerator
    dim = 3
    nx = 40
    ny = 40
    nz = 1
    xmax = 1000
    ymax = 1000
    zmax = 1.0
  []
  [subdomain]
    input = fmg
    type = SubdomainPerElementGenerator
    element_ids = '0 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26 27 28 29 30 31 32 33 34 35 36 37 38 39'
    subdomain_ids = '0 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26 27 28 29 30 31 32 33 34 35 36 37 38 39'
  []
  [bot_corner]
    type = ExtraNodesetGenerator
    new_boundary = bot_corner
    input = subdomain
    coord = '0 0.0 0.0'
  []
[]

[Modules]

  [TensorMechanics]

    [Master]
      [all]
        strain = FINITE
        add_variables = true
        volumetric_locking_correction = true
        generate_output = 'stress_yy strain_yy'
        save_in = 'resid_x resid_y resid_z'
      []
    []
  []
[]

[AuxVariables]
  [resid_x]
  []
  [resid_y]
  []
  [resid_z]
  []
  [euler_angle_1]
    order = CONSTANT
    family = MONOMIAL
  []
  [euler_angle_2]
    order = CONSTANT
    family = MONOMIAL
  []
  [euler_angle_3]
    order = CONSTANT
    family = MONOMIAL
  []
  # [updated_euler_angle_1]
  #   order = CONSTANT
  #   family = MONOMIAL
  # []
  # [updated_euler_angle_2]
  #   order = CONSTANT
  #   family = MONOMIAL
  # []
  # [updated_euler_angle_3]
  #   order = CONSTANT
  #   family = MONOMIAL
  # []

  [updated_quaternion_x]
    order = CONSTANT
    family = MONOMIAL
  []
  [updated_quaternion_y]
    order = CONSTANT
    family = MONOMIAL
  []
  [updated_quaternion_z]
    order = CONSTANT
    family = MONOMIAL
  []
  [updated_quaternion_w]
    order = CONSTANT
    family = MONOMIAL
  []


  [unique_grains]
    order = CONSTANT
    family = MONOMIAL
  []
[]

[AuxKernels]
  # [updated_euler_angle_1]
  #   type = MaterialRealVectorValueAux
  #   variable = updated_euler_angle_1
  #   property = updated_Euler_angle
  #   component = 0
  #   execute_on = 'TIMESTEP_END'
  # []
  # [updated_euler_angle_2]
  #   type = MaterialRealVectorValueAux
  #   variable = updated_euler_angle_2
  #   property = updated_Euler_angle
  #   component = 1
  #   execute_on = 'TIMESTEP_END'
  # []
  # [updated_euler_angle_3]
  #   type = MaterialRealVectorValueAux
  #   variable = updated_euler_angle_3
  #   property = updated_Euler_angle
  #   component = 2
  #   execute_on = 'TIMESTEP_END'
  # []

  [updated_quaternion_x]
    type = MaterialStdVectorAux
    variable = updated_quaternion_x
    property = updated_quaternion
    index = 0
    execute_on = 'TIMESTEP_END'
  []
  [updated_quaternion_y]
    type = MaterialStdVectorAux
    variable = updated_quaternion_y
    property = updated_quaternion
    index = 1
    execute_on = 'TIMESTEP_END'
  []
  [updated_quaternion_z]
    type = MaterialStdVectorAux
    variable = updated_quaternion_z
    property = updated_quaternion
    index = 2
    execute_on = 'TIMESTEP_END'
  []
  [updated_quaternion_w]
    type = MaterialStdVectorAux
    variable = updated_quaternion_w
    property = updated_quaternion
    index = 3
    execute_on = 'TIMESTEP_END'
  []
[]

[UserObjects]
  [assign_block_id]
    type = VariableValueElementSubdomainModifier
    coupled_var = 'unique_grains'
    execute_on = 'TIMESTEP_BEGIN'
    execution_order_group = -1
  []

  # [avg_ea1]
  #   type = BlockAverage
  #   variable = euler_angle_1
  #   execute_on = 'TIMESTEP_END'
  #   execution_order_group = -1
  # []
  # [avg_ea2]
  #   type = BlockAverage
  #   variable = euler_angle_2
  #   execute_on = 'TIMESTEP_END'
  #   execution_order_group = -1
  # []
  # [avg_ea3]
  #   type = BlockAverage
  #   variable = euler_angle_3
  #   execute_on = 'TIMESTEP_END'
  #   execution_order_group = -1
  # []

  [avg_qx]
    type = BlockAverage
    variable = updated_quaternion_x
    execute_on = 'TIMESTEP_END'
    execution_order_group = -1
  []
  [avg_qy]
    type = BlockAverage
    variable = updated_quaternion_y
    execute_on = 'TIMESTEP_END'
    execution_order_group = -1
  []
  [avg_qz]
    type = BlockAverage
    variable = updated_quaternion_z
    execute_on = 'TIMESTEP_END'
    execution_order_group = -1
  []
  [avg_qw]
    type = BlockAverage
    variable = updated_quaternion_w
    execute_on = 'TIMESTEP_END'
    execution_order_group = -1
  []
[]

[BCs]
  [x_roller]
    type = DirichletBC
    preset = true
    variable = disp_y
    boundary = bottom
    value = 0.0
  []
  [y_roller]
    type = DirichletBC
    preset = true
    variable = disp_x
    boundary = left
    value = 0.0
  []
  [z_roller]
    type = DirichletBC
    preset = true
    variable = disp_z
    boundary = back
    value = 0.0
  []

  [y_pull_function]
    type = FunctionDirichletBC
    variable = disp_z
    boundary = front
    function = '1.0e-4*t'
  []
[]

[Materials]
  [elasticity_tensor_copper]
    type = ComputeElasticityTensorCP
    C_ijkl = '1.684e5 1.214e5 1.214e5 1.684e5 1.214e5 1.684e5 0.754e5 0.754e5 0.754e5'
    fill_method = symmetric9
    euler_angle_variables = 'euler_angle_1 euler_angle_2 euler_angle_3'
  []

  [stress]
    type = ComputeMultipleCrystalPlasticityStress
    crystal_plasticity_models = 'trial_xtalpl_copper'
    tan_mod_type = exact
  []
  [trial_xtalpl_copper]
    type = CrystalPlasticityKalidindiUpdate
    number_slip_systems = 12
    slip_sys_file_name = 'input_slip_sys.txt'
  []
  [updated_euler_angle]
    type = ComputeUpdatedEulerAngle
    radian_to_degree = false
  []
[]

[Postprocessors]
  [stress_zz]
    type = ElementAverageValue
    variable = stress_yy
  []
  [strain_zz]
    type = ElementAverageValue
    variable = strain_yy
  []
[]

[VectorPostprocessors]
  # [updated_grain_ea]
  #   type = BlockAverageFromUserObjects
  #   block_average_uos = "avg_ea1 avg_ea2 avg_ea3"
  #   execute_on = 'TIMESTEP_END'
  # []
  [updated_grain_ea]
    type = BlockOrientationFromQuaternionUserObjects
    quaternion_average_uos = "avg_qx avg_qy avg_qz avg_qw"
    execute_on = 'TIMESTEP_END'
  []
[]

[Preconditioning]
  [SMP]
    type = SMP
    full = true
  []
[]

[Executioner]
  type = Transient

  solve_type = 'PJFNK'
  # petsc_options_iname = '-pc_type -pc_asm_overlap -sub_pc_type -ksp_type -ksp_gmres_restart'
  # petsc_options_value = ' asm      2              lu            gmres     200'

  petsc_options_iname = '-pc_type -pc_factor_mat_solver_package'
  petsc_options_value = 'lu superlu_dist'

  l_tol = 1e-8
  nl_abs_tol = 1e-7
  nl_rel_tol = 1e-6
  nl_max_its = 20
  nl_forced_its = 1
  l_max_its = 100
  num_steps = 5
  dt = 0.1

  # [TimeStepper]
  #   type = IterationAdaptiveDT
  #   dt = 0.1
  #   iteration_window = 2
  #   optimal_iterations = 10
  #   growth_factor = 1.2
  #   cutback_factor = 0.8
  # []
[]

[Outputs]
  execute_on = 'initial timestep_end'
  file_base = 64Cu_test
  csv = true
  print_linear_residuals = true
  perf_graph = true
  exodus = true
[]
