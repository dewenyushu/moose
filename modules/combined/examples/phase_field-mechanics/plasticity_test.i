[Mesh]
  type = GeneratedMesh
  dim = 2
  nx = 40
  ny = 40
  nz = 0
  xmax = 1000
  ymax = 1000
  zmax = 0
  elem_type = QUAD4
  # uniform_refine = 1
[]

# [GlobalParams]
#   op_num = 8
#   var_name_base = gr
#   grain_num = 36
# []

[Variables]
  # [./PolycrystalVariables]
  # [../]
  [./disp_x]
    order = FIRST
    family = LAGRANGE
  [../]
  [./disp_y]
    order = FIRST
    family = LAGRANGE
  [../]
[]

[UserObjects]
  # [./euler_angle_file]
  #   type = EulerAngleFileReader
  #   file_name = grn_36_rand_2D.tex
  # [../]
  # [./voronoi]
  #   type = PolycrystalVoronoi
  #   # coloring_algorithm = bt
  # [../]
  # [./grain_tracker]
  #   type = GrainTrackerElasticity
  #   threshold = 0.2
  #   compute_var_to_feature_map = true
  #   execute_on = 'initial timestep_begin'
  #   flood_entity_type = ELEMENTAL

  #   C_ijkl = '1.27e5 0.708e5 0.708e5 1.27e5 0.708e5 1.27e5 0.7355e5 0.7355e5 0.7355e5'
  #   fill_method = symmetric9
  #   euler_angle_provider = euler_angle_file
  # [../]
  [./str]
    type = TensorMechanicsHardeningConstant
    value = 0.0001
  [../]
  [./j2]
    type = TensorMechanicsPlasticJ2
    yield_strength = str
    yield_function_tolerance = 1E-3
    internal_constraint_tolerance = 1E-9
  [../]
[]

[ICs]
  # [./PolycrystalICs]
  #   [./PolycrystalColoringIC]
  #     polycrystal_ic_uo = voronoi
  #   [../]
  # [../]
[]

[AuxVariables]
  # [./bnds]
  #   order = FIRST
  #   family = LAGRANGE
  # [../]
  [./elastic_strain11]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./elastic_strain22]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./elastic_strain12]
    order = CONSTANT
    family = MONOMIAL
  [../]
  # [./unique_grains]
  #   order = CONSTANT
  #   family = MONOMIAL
  # [../]
  # [./var_indices]
  #   order = CONSTANT
  #   family = MONOMIAL
  # [../]
  # [./vonmises_stress]
  #   order = CONSTANT
  #   family = MONOMIAL
  # [../]
  [./force_x]
    order = FIRST
    family = LAGRANGE
  [../]
  [./force_y]
    order = FIRST
    family = LAGRANGE
  [../]
  [./C1111]
    order = CONSTANT
    family = MONOMIAL
  [../]
  # [./euler_angle]
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

[Modules]
  [./TensorMechanics]
    [./Master]
      displacements = 'disp_x disp_y'
      [./mech]
        add_variables = true
        strain = FINITE
        additional_generate_output = 'stress_yy stress_xy stress_xx strain_xx strain_xy strain_yy strain_zz hydrostatic_stress mid_principal_stress min_principal_stress max_principal_stress vonmises_stress'
        decomposition_method = EigenSolution
        save_in = 'force_x force_y'
        displacements = 'disp_x disp_y'
      [../]
    [../]
  [../]
[]

[AuxKernels]
  # [./BndsCalc]
  #   type = BndsCalcAux
  #   variable = bnds
  #   execute_on = timestep_end
  # [../]
  [./elastic_strain11]
    type = RankTwoAux
    variable = elastic_strain11
    rank_two_tensor = elastic_strain
    index_i = 0
    index_j = 0
    execute_on = timestep_end
  [../]
  [./elastic_strain22]
    type = RankTwoAux
    variable = elastic_strain22
    rank_two_tensor = elastic_strain
    index_i = 1
    index_j = 1
    execute_on = timestep_end
  [../]
  [./elastic_strain12]
    type = RankTwoAux
    variable = elastic_strain12
    rank_two_tensor = elastic_strain
    index_i = 0
    index_j = 1
    execute_on = timestep_end
  [../]
  # [./unique_grains]
  #   type = FeatureFloodCountAux
  #   variable = unique_grains
  #   execute_on = timestep_end
  #   flood_counter = grain_tracker
  #   field_display = UNIQUE_REGION
  # [../]
  # [./var_indices]
  #   type = FeatureFloodCountAux
  #   variable = var_indices
  #   execute_on = timestep_end
  #   flood_counter = grain_tracker
  #   field_display = VARIABLE_COLORING
  # [../]
  [./C1111]
    type = RankFourAux
    variable = C1111
    rank_four_tensor = elasticity_tensor
    index_l = 0
    index_j = 0
    index_k = 0
    index_i = 0
    execute_on = timestep_end
  [../]
  # [./vonmises_stress]
  #   type = RankTwoScalarAux
  #   variable = vonmises_stress
  #   rank_two_tensor = stress
  #   scalar_type = VonMisesStress
  #   execute_on = timestep_end
  # [../]
  # [./euler_angle]
  #   type = OutputEulerAngles
  #   variable = euler_angle
  #   euler_angle_provider = euler_angle_file
  #   grain_tracker = grain_tracker
  #   output_euler_angle = 'phi1'
  #   execute_on = 'initial timestep_end'
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

[BCs]
  # [./Periodic]
  #   [./All]
  #     auto_direction = 'x'
  #     variable = 'gr0 gr1 gr2 gr3 gr4 gr5 gr6 gr7'
  #   [../]
  # [../]
  # [./top_displacement]
  #   type = DirichletBC
  #   variable = disp_y
  #   boundary = top
  #   value = 0.1
  # [../]
    [./y]
      type = FunctionDirichletBC
      variable = disp_y
      boundary = 'top'
      function = '10.0*t'
    [../]
  [./x_anchor]
    type = DirichletBC
    variable = disp_x
    boundary = 'left'
    value = 0.0
  [../]
  [./y_anchor]
    type = DirichletBC
    variable = disp_y
    boundary = bottom
    value = 0.0
  [../]
[]

[Materials]
  # [./Copper]
  #   type = GBEvolution
  #   block = 0
  #   T = 500 # K
  #   wGB = 15 # nm
  #   GBmob0 = 2.5e-6 # m^4/(Js) from Schoenfelder 1997
  #   Q = 0.23 # Migration energy in eV
  #   GBenergy = 0.708 # GB energy in J/m^2
  # [../]
  # [./ElasticityTensor]
  #   type = ComputePolycrystalElasticityTensor
  #   grain_tracker = grain_tracker
  # [../]
  # [./elasticity_tensor]
  #   type = ComputeElasticityTensor
  #   block = 0
  #   fill_method = symmetric9
  #   C_ijkl = '1.27e5 0.708e5 0.708e5 1.27e5 0.708e5 1.27e5 0.7355e5 0.7355e5 0.7355e5'
  #   # C_ijkl = '1.27e6 0.708e6 0.708e6 1.27e6 0.708e6 1.27e6 0.7355e6 0.7355e6 0.7355e6'
  # [../]
  [./elasticity_tensor]
    type = ComputeElasticityTensor
    block = 0
    fill_method = symmetric_isotropic_E_nu
    C_ijkl = '385000 0.23'
  [../]
  # [./strain]
  #   type = ComputeFiniteStrain
  #   block = 0
  #   displacements = 'disp_x disp_y'
  # [../]
  # [./stress]
  #   type = ComputeLinearElasticStress
  #   block = 0
  # [../]
  [./mc]
    type = ComputeMultiPlasticityStress
    block = 0
    ep_plastic_tolerance = 1E-8
    plastic_models = j2
    debug_fspb = crash
  [../]
[]

[Postprocessors]
  # [./ngrains]
  #   type = FeatureFloodCount
  #   variable = bnds
  #   threshold = 0.7
  # [../]
  [./dofs]
    type = NumDOFs
  [../]
  [./dt]
    type = TimestepSize
  [../]
  [./run_time]
    type = PerfGraphData
    section_name = "Root"
    data_type = total
  [../]
  [./f]
    type = PointValue
    point = '100 100 0'
    variable = f
  [../]
  [./iter]
    type = PointValue
    point = '100 100 0'
    variable = iter
  [../]
[]

[Preconditioning]
  [./SMP]
    type = SMP
    # coupled_groups = 'disp_x,disp_y'
    full = true
  [../]
[]

[Executioner]
  type = Transient
  scheme = bdf2
  solve_type = NEWTON
  # line_search = none
  petsc_options_iname = '-pc_type -pc_asm_overlap -sub_pc_type -ksp_type -ksp_gmres_restart'
  petsc_options_value = ' asm      2              lu            gmres     200'
  # petsc_options_iname = '-pc_type -ksp_type -snes_type -pc_factor_shift_type -pc_factor_shift_amount '
  # petsc_options_value = 'lu preonly  vinewtonrsls NONZERO 1e-10'
  # petsc_options_iname = '-pc_type -pc_hypre_type -ksp_gmres_restart -pc_hypre_boomeramg_strong_threshold'
  # petsc_options_value = 'hypre boomeramg 31 0.7'
  # l_tol = 1.0e-3
  # l_max_its = 30
  # nl_max_its = 25
  # nl_rel_tol = 1.0e-6
  # nl_abs_tol = 1e-8
  # start_time = 0.0
  # dt = 1
  # num_steps = 150
  # automatic_scaling = true
  nl_rel_tol = 1e-6  ##nonlinear relative tolerance
  nl_abs_tol = 1e-6
  l_max_its = 10   ##max linear iterations Previous:200
  nl_max_its = 20  ##max nonlinear iterations Previous:50
  start_time=0
  line_search = 'none'
  end_time = 500
  num_steps = 1500
  dt = 1
  dtmin = 1e-15
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
  # file_base = poly36_grtracker
  exodus = true
[]
