radius= 0.5
height=5.0
#
# nz=2
# ring_num=2
# sector=4

refine=0

Job=1

end_time=1

[GlobalParams]
  order = FIRST
  family = LAGRANGE
  displacements = 'disp_x disp_y disp_z'
  volumetric_locking_correction = true
[]

[Mesh]
  patch_update_strategy = iteration
  patch_size = 10

  [file]
    type = FileMeshGenerator
    file = two_wires_mesh_4blocks.e
  []
  [slave]
    type = LowerDBlockFromSidesetGenerator
    sidesets = '40'
    new_block_id = 200
    new_block_name = 'slave_lower'
    input = file
  []
  [master]
    type = LowerDBlockFromSidesetGenerator
    sidesets = '10'
    new_block_id = 100
    new_block_name = 'master_lower'
    input = slave
  []

  uniform_refine =${refine}
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
    x = '0. ${end_time}.'
    y = '0. ${end_time}.'
    scale_factor = 0.5
  [../]
  [./-rampConstantAngle]
    type = PiecewiseLinear
    x = '0. ${end_time}.'
    y = '0. ${end_time}.'
    scale_factor = -0.5
  [../]
  [./decreaseZ]
    type = ParsedFunction
    value = 0 #if(t>25,-0.01*(t-25),0)
  [../]
[]

[Modules/TensorMechanics/Master]
  [./action]
    strain = FINITE
    add_variables = true
    block = '1 2 100 200'
    use_automatic_differentiation = false
  [../]
[]



[BCs]


  [./MainBottom_x]
    type = DisplacementAboutAxis
    boundary = 20
    function = -rampConstantAngle
    angle_units = degrees
    axis_origin = '0. 0. 0.'
    axis_direction = '0. 0. 1.0'
    component = 0
    variable = disp_x
  [../]
  [./MainBottom_y]
    type = DisplacementAboutAxis
    boundary = 20
    function = -rampConstantAngle
    angle_units = degrees
    axis_origin = '0. 0. 0.'
    axis_direction = '0. 0. 1.0'
    component = 1
    variable = disp_y
  [../]
  [./MainBottom_z]
    type = DirichletBC
    variable = disp_z
    boundary = 20
    value = 0
  [../]
  [./Maintop_x]
    type = DisplacementAboutAxis
    boundary = 30
    function = rampConstantAngle
    angle_units = degrees
    axis_origin = '0. 0. 0.'
    axis_direction = '0. 0. 1.0'
    component = 0
    variable = disp_x
  [../]
  [./Maintop_y]
    type = DisplacementAboutAxis
    boundary = 30
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
    boundary = 30
    value = 0
  [../]




  [./Bottom_x]
    type = DisplacementAboutAxis
    boundary = 50
    function = -rampConstantAngle
    angle_units = degrees
    axis_origin = '0. 0. 0.'
    axis_direction = '0. 0. 1.0'
    component = 0
    variable = disp_x
  [../]
  [./Bottom_y]
    type = DisplacementAboutAxis
    boundary = 50
    function = -rampConstantAngle
    angle_units = degrees
    axis_origin = '0. 0. 0.'
    axis_direction = '0. 0. 1.0'
    component = 1
    variable = disp_y
  [../]
  [./Bottom_z]
    type = DirichletBC
    variable = disp_z
    boundary = 50
    value = 0
  [../]
  [./top_x]
    type = DisplacementAboutAxis
    boundary = 60
    function = rampConstantAngle
    angle_units = degrees
    axis_origin = '0. 0. 0.'
    axis_direction = '0. 0. 1.0'
    component = 0
    variable = disp_x
  [../]
  [./top_y]
    type = DisplacementAboutAxis
    boundary = 60
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
    boundary = 60
    value = 0
  [../]

[] # BCs


[Materials]
  [./elasticity_tensor]   #Silicone_Rubber
    type = ComputeIsotropicElasticityTensor
    # type = ADComputeIsotropicElasticityTensor
    block = '1 2 100 200'
    youngs_modulus = 128E9
    poissons_ratio = 0.36
  [../]
  [./elastic_stress]
    type = ComputeFiniteStrainElasticStress
    # type = ADComputeFiniteStrainElasticStress
    block = '1 2 100 200'
  [../]
[]

[Contact]
  [./contact]
    mesh = block_sidesets
    formulation = kinematic
    primary = '10'
    secondary = '40'
    model = coulomb
    penalty = 1E7
    tangential_tolerance = .01
    normal_smoothing_distance = .01
    normalize_penalty = true
  [../]
[]

[Preconditioning]
   [./FSP]
     type = FSP
     # It is the starting point of splitting
     topsplit = 'contact_interior' # 'contact_interior' should match the following block name
     [./contact_interior]
       splitting          = 'contact interior'
       splitting_type     = multiplicative
     [../]
     [./interior]
       type = ContactSplit
       vars = 'disp_x disp_y'
       uncontact_primary =  '10'
       uncontact_secondary =  '40'
       # contact_displaced = '20'
       blocks              = '1 2'
       include_all_contact_nodes = 1

       # petsc_options_iname = '-ksp_type -pc_type -pc_hypre_type '
       # petsc_options_value = '  preonly hypre  boomeramg'
       petsc_options_iname = '-ksp_type -pc_sub_type -pc_factor_shift_type  -pc_factor_shift_amount'
       petsc_options_value = '  preonly lu NONZERO 1e-15'
      [../]
      [./contact]
       type = ContactSplit
       vars = 'disp_x disp_y'
       contact_primary =  '10'
       contact_secondary =  '40'
       # contact_displaced = '20'
       include_all_contact_nodes = 1
       blocks = '200'


       petsc_options_iname = '-ksp_type -pc_sub_type -pc_factor_shift_type  -pc_factor_shift_amount'
       petsc_options_value = '  preonly lu NONZERO 1e-15'
     [../]
   [../]
 []



[Executioner]

  type = Transient
  solve_type = 'PJFNK'

  petsc_options=' -snes_monitor -ksp_monitor'

  line_search = 'none'

  l_max_its = 200
  nl_max_its = 20
  nl_abs_tol = 0.9E-8
  nl_rel_tol = 0.9E-7
  l_tol = 0.9E-3

  start_time = 0.0
  dt = 0.25
  dtmin = 0.01

  end_time = ${end_time}

  [./Predictor]
    type = SimplePredictor
    scale = 0.5
  [../]
[]

# [Preconditioning]
#   [./SMP]
#     type = SMP
#     full = true
#   [../]
# []
#
# [Executioner]
#
#   type = Transient
#   solve_type = 'PJFNK'
#
#   petsc_options = '-snes_ksp_ew -log_view'
#   petsc_options_iname = '-ksp_gmres_restart -pc_type -pc_factor_mat_solver_package'
#   petsc_options_value = '200 lu superlu_dist'
#
#   line_search = 'none'
#
#   l_max_its = 7000
#   nl_max_its = 100
#   nl_abs_tol = 0.9E-8
#   nl_rel_tol = 0.9E-7
#   l_tol = 0.9E-3
#
#   start_time = 0.0
#   dt = 0.25
#   dtmin = 0.01
#
#   end_time = ${end_time}
#
#   [./Predictor]
#     type = SimplePredictor
#     scale = 0.5
#   [../]
# []

[Outputs]
  file_base = ./stroing_scaling_output/torque_2wires_height${height}_radii${radius}_${end_time}deg_job${Job}_out
  interval = 4
  # [./exodus]
  #       type = Exodus
  # [../]
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
