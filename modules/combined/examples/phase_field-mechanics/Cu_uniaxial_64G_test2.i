L = 0.4 #mm
# area = ${fparse L*L}

[GlobalParams]
    displacements = 'disp_x disp_y disp_z'
[]

[Mesh]
    # construct_side_list_from_node_list = false
    displacements = 'disp_x disp_y disp_z'
    # construct_side_list_from_node_list = false
    [fmg]
      type = GeneratedMeshGenerator
      dim = 3
      nx = 5
      ny = 5
      nz = 2
      xmax = 100
      ymax = 100
      zmax = 100.0
      # file = 64G_512E_M.e
    []
    [bot_corner]
      type = ExtraNodesetGenerator
      new_boundary = bot_corner
      input = fmg
      coord = '0 0.0 0.0'
    []
    # [add_side_sets]
    #   type = SideSetsFromNormalsGenerator
    #   normals = '1  0  0
    #              0  1  0
    #              0  0  1
    #             -1  0  0
    #              0 -1  0
    #              0  0 -1'
    #   fixed_normal = false
    #   new_boundary = 'xp_face yp_face zp_face xn_face yn_face zn_face'
    #   input=bot_corner
    # []
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
[]

[UserObjects]
    [prop_read]
      type = PropertyReadFile
      prop_file_name = 'EulerAngles.txt'
      #NB: euler angles must be supplied in degrees and internally converted to radians.
      #Random grain orientation.
      read_type = 'element'
      nprop = 3
      # nblock=64
      use_zero_based_block_indexing = 'false'
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
       variable = disp_y
       boundary = top
       function = '1.0e-3*${L}*t'
     []
   []

[Materials]
  [elasticity_tensor_copper]
    type = ComputeElasticityTensorCP
    C_ijkl = '1.684e5 1.214e5 1.214e5 1.684e5 1.214e5 1.684e5 0.754e5 0.754e5 0.754e5'
    fill_method = symmetric9
    read_prop_user_object = prop_read
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
  # [react_z]
  #   type = NodalSum
  #   variable = resid_z
  #   boundary = 'back'
  # []
  # [avg_dispz]
  #   type = AverageNodalVariableValue
  #   variable = disp_z
  #   boundary = 'front'
  # []

  # [avgstress]
  #   type = ParsedPostprocessor
  #   function = '-react_z /${area} '
  #   pp_names = 'react_z'
  # []

  # [avgstrain]
  #   type = ParsedPostprocessor
  #   function = 'avg_dispz /${L} '
  #   pp_names = 'avg_dispz'
  # []

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
  end_time = 5

  [TimeStepper]
    type = IterationAdaptiveDT
    dt=1e-6
    # maxdt = 0.5
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
