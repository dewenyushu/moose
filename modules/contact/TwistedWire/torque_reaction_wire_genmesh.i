radius= 1.0
height=60.0

nz=600
ring_num=10
sector=8

refine=0

Job=1

end_time= 540

order = FIRST

[GlobalParams]
  order = ${order}
  family = LAGRANGE
  displacements = 'disp_x disp_y disp_z'
  volumetric_locking_correction = true
[]

[Mesh]
  patch_update_strategy = iteration
  patch_size = 10
  partitioner = centroid
  centroid_partitioner_direction = z

  [./c1]
    type = ConcentricCircleMeshGenerator
    num_sectors = ${sector}
    radii = '${radius}'
    rings = '${ring_num}'
    has_outer_square = false
    preserve_volumes = off
    smoothing_max_it = 3
  []
  [./c1_rename]
    type = RenameBoundaryGenerator
    input = c1
    old_boundary_name = 'outer'
    new_boundary_name = '10'
  []
  [./cylinder1]
    input = c1_rename
    type = MeshExtruderGenerator
    num_layers = ${nz}
    extrusion_vector = '0 0 ${height}'
    bottom_sideset = '20'
    top_sideset = '30'
  [../]
  [./cylinder1_translate]
    type = TransformGenerator
    transform = TRANSLATE
    vector_value = '-${radius} 0 0'
    input = cylinder1
  [../]
  [./cylinder1_id]
    type = SubdomainIDGenerator
    input = cylinder1_translate
    subdomain_id = 1
  [../]

  [./c2]
    type = ConcentricCircleMeshGenerator
    num_sectors = ${sector}
    radii = '${radius}'
    rings = '${ring_num}'
    has_outer_square = false
    preserve_volumes = off
    smoothing_max_it = 3
  []
  [./c2_rename]
    type = RenameBoundaryGenerator
    input = c2
    old_boundary_name = 'outer'
    new_boundary_name = '40'
  []
  [./cylinder2]
    input = c2_rename
    type = MeshExtruderGenerator
    num_layers = ${nz}
    extrusion_vector = '0 0 ${height}'
    bottom_sideset = '50'
    top_sideset = '60'
  [../]
  [./cylinder2_translate]
    type = TransformGenerator
    transform = TRANSLATE
    vector_value = '${radius} 0 0'
    input = cylinder2
  [../]
  [./cylinder2_id]
    type = SubdomainIDGenerator
    input = cylinder2_translate
    subdomain_id = 2
  [../]

  [./combined]
    type = MeshCollectionGenerator
    inputs = 'cylinder1_id cylinder2_id'
  [../]
  [./block_sidesets]
    type = SideSetsFromPointsGenerator
    input = combined
    points = '-${radius} 0 0
              -${radius} 0 ${height}
              ${radius} 0 0
              ${radius} 0 ${height}
              ${radius} ${radius} 5
              -${radius} ${radius} 5'
    new_boundary = '20 30 50 60 40 10'
  [../]
  uniform_refine =${refine}
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
    block = '1 2'
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
    block = '1 2'
    youngs_modulus = 128E9
    poissons_ratio = 0.36
  [../]
  [./elastic_stress]
    type = ComputeFiniteStrainElasticStress
    # type = ADComputeFiniteStrainElasticStress
    block = '1 2'
  [../]
[]

[Contact]
  [./contact]
    mesh = block_sidesets
    formulation = kinematic
    master = 10
    slave = 40
    model = coulomb
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
  dt = 2.0
  dtmin = 0.1

  end_time = ${end_time}

  [./Predictor]
    type = SimplePredictor
    scale = 0.5
  [../]
[]

[Outputs]
  file_base = ./output/torque_2wires_height${height}_radii${radius}_${end_time}deg_job${Job}_out
  interval = 4
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
