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
    [bot_corner]
      type = ExtraNodesetGenerator
      new_boundary = bot_corner
      input = fmg
      coord = '0 0.0 0.0'
    []
  []

[Modules/TensorMechanics/Master]
    [all]
      strain = FINITE
      add_variables = true
      volumetric_locking_correction = true
      generate_output = 'stress_yy strain_yy'
      save_in = 'resid_x resid_y resid_z'
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
  [unique_grains]
    order = CONSTANT
    family = MONOMIAL
  []
[]

[UserObjects]
     [assign_block_id]
      type = VariableValueElementSubdomainModifier
      coupled_var = 'unique_grains'
      execute_on = 'INITIAL TIMESTEP_BEGIN'
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

[Preconditioning]
  [SMP]
    type = SMP
    full=true
  []
[]


[Executioner]
  type = Transient

  solve_type = 'PJFNK'
  petsc_options_iname = '-pc_type -pc_asm_overlap -sub_pc_type -ksp_type -ksp_gmres_restart'
  petsc_options_value = ' asm      2              lu            gmres     200'

  l_tol = 1e-8
  nl_abs_tol = 1e-7
  nl_rel_tol = 1e-6
  nl_max_its = 20
  nl_forced_its = 1
  l_max_its = 100
  start_time = 0.0
  end_time = 100

  [TimeStepper]
    type = IterationAdaptiveDT
    dt=1e-6
    iteration_window = 2
    optimal_iterations = 10
    growth_factor = 1.2
    cutback_factor = 0.3
  []
[]

[Outputs]
  execute_on = 'initial timestep_end'
  file_base = 64Cu_test
  csv = true
  print_linear_residuals = true
  perf_graph = true
  exodus = true
[]
