
[GlobalParams]
  order = SECOND
  family = LAGRANGE
  displacements = 'disp_x disp_y disp_z'
[]

[Mesh]
  file = two_wires_mesh.e
  # file = twowirehexfiner.e
  patch_update_strategy = iteration
  patch_size = 30
[]

# [Variables]
#   [./disp_x]
#   [../]
#   [./disp_y]
#   [../]
#   [./disp_z]
#   [../]
#
# []

[Preconditioning]
  [./SMP]
    type = SMP
    full = true
  [../]
[]



[AuxVariables]
  [./saved_x]
  [../]
  [./saved_y]
  [../]
  [./saved_z]
  [../]
[]

[Functions]
  [./rampConstantAngle]
    type = PiecewiseLinear
    x = '0. 360.'
    y = '0. 360.'
    scale_factor = 0.5
  [../]
  [./-rampConstantAngle]
    type = PiecewiseLinear
    x = '0. 360.'
    y = '0. 360.'
    scale_factor = -1
  [../]
  [./decreaseZ]
    type = ParsedFunction
    value = 0 #if(t>25,-0.01*(t-25),0)
  [../]
[]

[Modules]
  [./TensorMechanics]
    [./Master]
      [./all]
        strain = FINITE
        add_variables = true
      [../]
      save_in = 'saved_x saved_y saved_z'
      use_displaced_mesh = true
    [../]
  [../]
[]



[BCs]


  [./MainBottom_x]
    type = DirichletBC
    variable = disp_x
    boundary = 2
    value = 0
  [../]
  [./MainBottom_y]
    type = DirichletBC
    variable = disp_y
    boundary = 2
    value = 0
  [../]
  [./MainBottom_z]
    type = DirichletBC
    variable = disp_z
    boundary = 2
    value = 0
  [../]
  [./Maintop_x]
    type = DisplacementAboutAxis
    boundary = 3
    function = rampConstantAngle
    angle_units = degrees
    axis_origin = '0. 0. 0.'
    axis_direction = '0. 0. 1.0'
    component = 0
    variable = disp_x
  [../]
  [./Maintop_y]
    type = DisplacementAboutAxis
    boundary = 3
    function = rampConstantAngle
    angle_units = degrees
    axis_origin = '0. 0. 0.'
    axis_direction = '0. 0. 1.0'
    component = 1
    variable = disp_y
  [../]
  [./Maintop_z]
    type = DirichletBC
    variable = disp_z
    boundary = 3
    value = 0
  [../]




  [./Bottom_x]
    type = DirichletBC
    variable = disp_x
    boundary = 5
    value = 0
  [../]
  [./Bottom_y]
    type = DirichletBC
    variable = disp_y
    boundary = 5
    value = 0
  [../]
  [./Bottom_z]
    type = DirichletBC
    variable = disp_z
    boundary = 5
    value = 0
  [../]
  [./top_x]
    type = DisplacementAboutAxis
    boundary = 6
    function = rampConstantAngle
    angle_units = degrees
    axis_origin = '0. 0. 0.'
    axis_direction = '0. 0. 1.0'
    component = 0
    variable = disp_x
  [../]
  [./top_y]
    type = DisplacementAboutAxis
    boundary = 6
    function = rampConstantAngle
    angle_units = degrees
    axis_origin = '0. 0. 0.'
    axis_direction = '0. 0. 1.0'
    component = 1
    variable = disp_y
  [../]
  [./top_z]
    type = DirichletBC
    variable = disp_z
    boundary = 6
    value = 0
  [../]

[] # BCs


[Materials]
  [./elasticity_tensor]   #Silicone_Rubber
    type = ComputeIsotropicElasticityTensor
    block = '1 2'
    youngs_modulus = 128E9
    poissons_ratio = 0.36
  [../]
  [./elastic_stress]
    type = ComputeFiniteStrainElasticStress
    block = '1 2'
  [../]
[]

[Contact]
  [./contact]
    formulation = kinematic
    master = 1
    slave = 4
    model = coulomb
    order = SECOND
    penalty = 1E7
    tangential_tolerance = .01
    normal_smoothing_distance = .01
    normalize_penalty = true
  [../]

[]

[Executioner]

  type = Transient
  solve_type = 'PJFNK'

  petsc_options = '-snes_ksp_ew'
  petsc_options_iname = '-ksp_gmres_restart -pc_type -pc_factor_mat_solver_package'
  petsc_options_value = '200 lu superlu_dist'

  line_search = 'none'

  l_max_its = 7000
  nl_max_its = 100
  nl_abs_tol = 0.9E-8
  nl_rel_tol = 0.9E-7
  l_tol = 0.9E-3

  start_time = 0.0
  dt = 0.5
  dtmin = 0.05

  end_time = 360
[]

[Outputs]
  file_base = ./torque_reaction_wire_coarse/torque_reaction_wires_out
  [./exodus]
        type = Exodus
  [../]
  [./console]
        type = Console
        max_rows = 5
  [../]
  [./pgraph]
        type = PerfGraphOutput
        execute_on = FINAL
        level = 1
  [../]
[]

[Postprocessors]
  [./nl_its]
    type = NumNonlinearIterations
    execute_on = 'initial timestep_end'
  [../]
  [./lin_its]
    type = NumLinearIterations
    execute_on = 'initial timestep_end'
  [../]
  [./cumul_nl]
    type = CumulativeValuePostprocessor
    postprocessor = nl_its
  [../]
  [./cumul_lin]
    type = CumulativeValuePostprocessor
    postprocessor = lin_its
  [../]
[]
