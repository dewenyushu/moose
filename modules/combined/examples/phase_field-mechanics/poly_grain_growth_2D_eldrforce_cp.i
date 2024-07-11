[Mesh]
  type = GeneratedMesh
  dim = 2
  nx = 20
  ny = 20
  nz = 0
  xmax = 1000
  ymax = 1000
  zmax = 0
  elem_type = QUAD4
  uniform_refine = 2
[]

[GlobalParams]
  op_num = 8
  var_name_base = gr
  grain_num = 36
  displacements = 'disp_x disp_y'
[]

[Variables]
  [PolycrystalVariables]
  []
[]

[UserObjects]
  [euler_angle_file]
    type = EulerAngleUpdateFromReporter
    file_name = grn_36_rand_2D.tex
    execute_on = 'initial timestep_begin'
    execution_order_group = -1 # execute before grain tracker

    euler_angle_0_name = updated_ea/ea0
    euler_angle_1_name = updated_ea/ea1
    euler_angle_2_name = updated_ea/ea2
    grain_id_name = updated_ea/grain_id
  []
  # [euler_angle_file]
  #   type = EulerAngleFileReader
  #   file_name = grn_36_rand_2D.tex
  #   execute_on = 'initial timestep_begin'
  #   execution_order_group = -1 # execute before grain tracker
  # []
  [voronoi]
    type = PolycrystalVoronoi
    # coloring_algorithm = bt
  []
  # [elasticity_tensor_copper]
  #   type = ComputeElasticityTensorCP
  #   C_ijkl = '1.684e5 1.214e5 1.214e5 1.684e5 1.214e5 1.684e5 0.754e5 0.754e5 0.754e5'
  #   fill_method = symmetric9
  #   euler_angle_variables = 'euler_angle_1 euler_angle_2 euler_angle_3'
  # []

  [grain_tracker]
    type = GrainTrackerElasticity
    threshold = 0.2
    compute_var_to_feature_map = true
    execute_on = 'initial timestep_begin'
    flood_entity_type = ELEMENTAL

    C_ijkl = '1.27e5 0.708e5 0.708e5 1.27e5 0.708e5 1.27e5 0.7355e5 0.7355e5 0.7355e5'
    fill_method = symmetric9
    euler_angle_provider = euler_angle_file
  []
  # [./str]
  #   type = TensorMechanicsHardeningConstant
  #   value = 1
  # [../]
  # [./j2]
  #   type = TensorMechanicsPlasticJ2
  #   yield_strength = str
  #   yield_function_tolerance = 1E-3
  #   internal_constraint_tolerance = 1E-9
  # [../]
[]

[MultiApps]
  [sub_app]
    type = TransientMultiApp
    # positions = '0 0 0'
    input_files = 'Cu_uniaxial_64G_test_coupling.i'
    # app_type = coupling_xolotlApp
    # execute_on = TIMESTEP_END
    # library_path = 'lib'
  []
[]

[Transfers]
  [from_disp_x]
    type = MultiAppGeometricInterpolationTransfer
    direction = from_multiapp
    multi_app = sub_app
    source_variable = 'disp_x'
    variable = 'disp_x'
  []
  [from_disp_y]
    type = MultiAppGeometricInterpolationTransfer
    direction = from_multiapp
    multi_app = sub_app
    source_variable = 'disp_y'
    variable = 'disp_y'
  []
  # [./from_sub_disloc]
  #   type = MultiAppInterpolationTransfer
  #   direction = from_multiapp
  #   multi_app = sub_app
  #   source_variable = 'disloc'
  #   variable = 'disloc'
  # [../]
  [tosub_euler0]
    type = MultiAppGeometricInterpolationTransfer
    direction = to_multiapp
    multi_app = sub_app
    source_variable = euler_angle0
    variable = euler_angle_1
    num_points = 1 # Set the value to be equal to that of the point nearest to the element (no interpolation)
  []
  [tosub_euler1]
    type = MultiAppGeometricInterpolationTransfer
    direction = to_multiapp
    multi_app = sub_app
    source_variable = euler_angle1
    variable = euler_angle_2
    num_points = 1
  []
  [tosub_euler2]
    type = MultiAppGeometricInterpolationTransfer
    direction = to_multiapp
    multi_app = sub_app
    source_variable = euler_angle2
    variable = euler_angle_3
    num_points = 1
  []
  [tosub_grain]
    type = MultiAppGeometricInterpolationTransfer
    direction = to_multiapp
    multi_app = sub_app
    source_variable = 'unique_grains'
    variable = 'unique_grains'
    num_points = 1
  []
[]

[ICs]
  [PolycrystalICs]
    [PolycrystalColoringIC]
      polycrystal_ic_uo = voronoi
    []
  []
[]

[AuxVariables]
  [bnds]
    order = FIRST
    family = LAGRANGE
  []
  [disp_x]
    order = FIRST
    family = LAGRANGE
  []
  [disp_y]
    order = FIRST
    family = LAGRANGE
  []
  [elastic_strain11]
    order = CONSTANT
    family = MONOMIAL
  []
  [elastic_strain22]
    order = CONSTANT
    family = MONOMIAL
  []
  [elastic_strain12]
    order = CONSTANT
    family = MONOMIAL
  []
  [unique_grains]
    order = CONSTANT
    family = MONOMIAL
  []
  [var_indices]
    order = CONSTANT
    family = MONOMIAL
  []
  # [./vonmises_stress]
  #   order = CONSTANT
  #   family = MONOMIAL
  # [../]
  [C1111]
    order = CONSTANT
    family = MONOMIAL
  []
  [euler_angle0]
    order = CONSTANT
    family = MONOMIAL
  []
  [euler_angle1]
    order = CONSTANT
    family = MONOMIAL
  []
  [euler_angle2]
    order = CONSTANT
    family = MONOMIAL
  []
  # [./f]
  #   order = CONSTANT
  #   family = MONOMIAL
  # [../]
  # [./iter]
  #   order = CONSTANT
  #   family = MONOMIAL
  # [../]
[]

[Kernels]
  [PolycrystalKernel]
  []
  [PolycrystalElasticDrivingForce]
  []
  # [./TensorMechanics]
  #   use_displaced_mesh = true
  #   displacements = 'disp_x disp_y'
  # [../]
[]

# [Modules/TensorMechanics/Master]
#   [all]
#     strain = FINITE
#     add_variables = true
#     volumetric_locking_correction = true
#     # generate_output = 'stress_zz strain_zz'
#     # save_in = 'resid_x resid_y'
#   []
# []

[AuxKernels]
  [BndsCalc]
    type = BndsCalcAux
    variable = bnds
    execute_on = timestep_end
  []
  [elastic_strain11]
    type = RankTwoAux
    variable = elastic_strain11
    rank_two_tensor = elastic_strain
    index_i = 0
    index_j = 0
    execute_on = timestep_end
  []
  [elastic_strain22]
    type = RankTwoAux
    variable = elastic_strain22
    rank_two_tensor = elastic_strain
    index_i = 1
    index_j = 1
    execute_on = timestep_end
  []
  [elastic_strain12]
    type = RankTwoAux
    variable = elastic_strain12
    rank_two_tensor = elastic_strain
    index_i = 0
    index_j = 1
    execute_on = timestep_end
  []
  [unique_grains]
    type = FeatureFloodCountAux
    variable = unique_grains
    execute_on = timestep_end
    flood_counter = grain_tracker
    field_display = UNIQUE_REGION
  []
  [var_indices]
    type = FeatureFloodCountAux
    variable = var_indices
    execute_on = timestep_end
    flood_counter = grain_tracker
    field_display = VARIABLE_COLORING
  []
  [C1111]
    type = RankFourAux
    variable = C1111
    rank_four_tensor = elasticity_tensor
    index_l = 0
    index_j = 0
    index_k = 0
    index_i = 0
    execute_on = timestep_end
  []
  # [./vonmises_stress]
  #   type = RankTwoScalarAux
  #   variable = vonmises_stress
  #   rank_two_tensor = stress
  #   scalar_type = VonMisesStress
  #   execute_on = timestep_end
  # [../]
  [euler_angle0]
    type = OutputEulerAngles
    variable = euler_angle0
    euler_angle_provider = euler_angle_file
    grain_tracker = grain_tracker
    output_euler_angle = 'phi1'
    execute_on = 'initial timestep_end'
  []
  [euler_angle1]
    type = OutputEulerAngles
    variable = euler_angle1
    euler_angle_provider = euler_angle_file
    grain_tracker = grain_tracker
    output_euler_angle = 'Phi'
    execute_on = 'initial timestep_end'
  []
  [euler_angle2]
    type = OutputEulerAngles
    variable = euler_angle2
    euler_angle_provider = euler_angle_file
    grain_tracker = grain_tracker
    output_euler_angle = 'phi2'
    execute_on = 'initial timestep_end'
  []
  # [./f]
  #   type = MaterialStdVectorAux
  #   index = 0
  #   property = plastic_yield_function
  #   variable = f
  # [../]
  # [./iter]
  #   type = MaterialRealAux
  #   property = plastic_NR_iterations
  #   variable = iter
  # [../]
[]

[BCs]
  [Periodic]
    [All]
      auto_direction = 'x'
      variable = 'gr0 gr1 gr2 gr3 gr4 gr5 gr6 gr7'
    []
  []
  # [./top_displacement]
  #   type = DirichletBC
  #   variable = disp_y
  #   boundary = top
  #   value = 0.1
  # [../]
  # [./y]
  #   type = FunctionDirichletBC
  #   variable = disp_y
  #   boundary = 'front back'
  #   function = '0E-6*y'
  # [../]
  # [./x_anchor]
  #   type = DirichletBC
  #   variable = disp_x
  #   boundary = 'left right'
  #   value = 0.0
  # [../]
  # [./y_anchor]
  #   type = DirichletBC
  #   variable = disp_y
  #   boundary = bottom
  #   value = 0.0
  # [../]
[]

[Materials]
  [Copper]
    type = GBEvolution
    block = 0
    T = 500 # K
    wGB = 15 # nm
    GBmob0 = 2.5e-6 # m^4/(Js) from Schoenfelder 1997
    Q = 0.23 # Migration energy in eV
    GBenergy = 0.708 # GB energy in J/m^2
  []
  [ElasticityTensor]
    type = ComputePolycrystalElasticityTensor
    grain_tracker = grain_tracker
    euler_angle_provider = euler_angle_file
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
  [strain]
    type = ComputeFiniteStrain
    block = 0
    displacements = 'disp_x disp_y'
  []
  # [./stress]
  #   type = ComputeLinearElasticStress
  #   block = 0
  # [../]
  # [./mc]
  #   type = ComputeMultiPlasticityStress
  #   block = 0
  #   ep_plastic_tolerance = 1E-9
  #   plastic_models = j2
  #   debug_fspb = crash
  # [../]
[]

[Postprocessors]
  [ngrains]
    type = FeatureFloodCount
    variable = bnds
    threshold = 0.7
  []
  [dofs]
    type = NumDOFs
  []
  [dt]
    type = TimestepSize
  []
  [run_time]
    type = PerfGraphData
    section_name = "Root"
    data_type = total
  []
  # [./f]
  #   type = PointValue
  #   point = '0 0 0'
  #   variable = f
  # [../]
  # [./iter]
  #   type = PointValue
  #   point = '0 0 0'
  #   variable = iter
  # [../]
[]

[Reporters]
  [updated_ea]
    type = ConstantReporter
    real_vector_names = 'ea0 ea1 ea2'
    real_vector_values = '0; 0; 0' # Dummy value

    dof_id_type_vector_names = 'grain_id'
    dof_id_type_vector_values = '0'
  []
[]

[Preconditioning]
  [SMP]
    type = SMP
    # coupled_groups = 'disp_x,disp_y'
  []
[]

[Executioner]
  type = Transient
  scheme = bdf2
  solve_type = PJFNK
  petsc_options_iname = '-pc_type -pc_hypre_type -ksp_gmres_restart -pc_hypre_boomeramg_strong_threshold'
  petsc_options_value = 'hypre boomeramg 31 0.7'
  l_tol = 1.0e-4
  l_max_its = 30
  nl_max_its = 25
  nl_rel_tol = 1.0e-7
  # start_time = 0.0
  # end_time = 100
  num_steps = 5
  [TimeStepper]
    type = IterationAdaptiveDT
    dt = 1.5
    growth_factor = 1.2
    cutback_factor = 0.8
    optimal_iterations = 8
  []
  # [./Adaptivity]
  #   initial_adaptivity = 2
  #   refine_fraction = 0.8
  #   coarsen_fraction = 0.05
  #   max_h_level = 3
  # [../]
[]

[Outputs]
  # file_base = poly36_grtracker
  exodus = true
[]
