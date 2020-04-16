refine=0

[GlobalParams]
  order = FIRST
  family = LAGRANGE
  displacements = 'disp_x disp_y disp_z'
  volumetric_locking_correction = true
[]

[Mesh]
  [./c1]
    type = ConcentricCircleMeshGenerator
    num_sectors = 4
    radii = '0.5'
    rings = '3'
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
    num_layers = 20
    extrusion_vector = '0 0 10'
    bottom_sideset = '20'
    top_sideset = '30'
  [../]
  [./cylinder1_translate]
    type = TransformGenerator
    transform = TRANSLATE
    vector_value = '0 0 0'
    input = cylinder1
  [../]
  [./cylinder1_id]
    type = SubdomainIDGenerator
    input = cylinder1_translate
    subdomain_id = 1
  [../]

  [./block_sidesets]
    type = SideSetsFromPointsGenerator
    input = cylinder1_id
    points = '0 0 0
              0 0 10
              -0.5 0 5'
    new_boundary = '20 30 10'
  [../]
  uniform_refine =${refine}
[]

[AuxVariables]
  [./saved_x]
  [../]
  [./saved_y]
  [../]
  [./saved_z]
  [../]
  [./stress_xx]      # stress aux variables are defined for output; this is a way to get integration point variables to the output file
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./vonmises]
    order = CONSTANT
    family = MONOMIAL
  [../]
[]

[Functions]
  [./rampConstantAngle]
    type = PiecewiseLinear
    x = '0. 360.'
    y = '0. 360.'
    # scale_factor = 0.1
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

[AuxKernels]
  [./stress_xx]
    type = RankTwoAux
    rank_two_tensor = stress
    variable = stress_xx
    index_i = 0
    index_j = 0
    execute_on = 'timestep_end'
  [../]
  [./vonmises]
    type = RankTwoScalarAux
    rank_two_tensor = stress
    variable = vonmises
    scalar_type = VonMisesStress
    execute_on = timestep_end
  [../]
[]


[BCs]


  [./MainBottom_x]
    type = DirichletBC
    variable = disp_x
    boundary = 20
    value = 0
  [../]
  [./MainBottom_y]
    type = DirichletBC
    variable = disp_y
    boundary = 20
    value = 0
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
    axis_origin = '0. 0. 10.'
    axis_direction = '0. 0. 1.0'
    component = 0
    variable = disp_x
  [../]
  [./Maintop_y]
    type = DisplacementAboutAxis
    boundary = 30
    function = rampConstantAngle
    angle_units = degrees
    axis_origin = '0. 0. 10.'
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

[] # BCs


[Materials]
  [./elasticity_tensor]   #Silicone_Rubber
    type = ComputeIsotropicElasticityTensor
    block = 1
    youngs_modulus = 128E9
    poissons_ratio = 0.36
  [../]
  [./elastic_stress]
    type = ComputeFiniteStrainElasticStress
    block = 1
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

  end_time = 40
[]

[Outputs]
  file_base = ./output/torque_1wire_out
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
