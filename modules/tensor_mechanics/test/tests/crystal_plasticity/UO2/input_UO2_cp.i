[GlobalParams]
  displacements = 'disp_x disp_y disp_z'
[]

[Mesh]
  [cube]
    type = GeneratedMeshGenerator
    dim = 3
    elem_type = HEX8
  []
[]

[AuxVariables]
  [fp_zz]
    order = CONSTANT
    family = MONOMIAL
  []
  [slip_increment_1]
    order = CONSTANT
    family = MONOMIAL
  []
  [slip_increment_10]
    order = CONSTANT
    family = MONOMIAL
  []
  [dislocation_density_1]
    order = CONSTANT
    family = MONOMIAL
  []
  [dislocation_density_2]
    order = CONSTANT
    family = MONOMIAL
  []
[]

[Modules/TensorMechanics/Master/all]
  strain = FINITE
  add_variables = true
[]

[AuxKernels]
  [fp_zz]
    type = RankTwoAux
    variable = fp_zz
    rank_two_tensor = plastic_deformation_gradient
    index_j = 2
    index_i = 2
    execute_on = timestep_end
  []
  [slip_increment_1]
   type = MaterialStdVectorAux
   variable = slip_increment_1
   property = slip_increment
   index = 1
   execute_on = timestep_end
  []
  [slip_increment_10]
   type = MaterialStdVectorAux
   variable = slip_increment_10
   property = slip_increment
   index = 10
   execute_on = timestep_end
  []
  [dislocation_density_1]
   type = MaterialStdVectorAux
   variable = dislocation_density_1
   property = forest_dislocation_density
   index = 1
   execute_on = timestep_end
  []
  [dislocation_density_2]
   type = MaterialStdVectorAux
   variable = dislocation_density_2
   property = forest_dislocation_density
   index = 2
   execute_on = timestep_end
  []
[]

[BCs]
  [fix_y]
    type = DirichletBC
    variable = disp_y
    preset = true
    boundary = 'bottom'
    value = 0
  []
  [fix_x]
    type = DirichletBC
    variable = disp_x
    boundary = 'left'
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
    boundary = front
    function = '5.0e-4*t'
  []
[]

[Materials]
  [elasticity_tensor]
    type = ComputeElasticityTensorCP
    C_ijkl = '1.08e5 6.034e4 6.034e4 1.08e5 6.03e4 1.08e5 2.86e4 2.86e4 2.86e4' #Tallon and Wolfenden. J. Phys. Chem. Solids (1979)
    fill_method = symmetric9
    euler_angle_1 = 54.74
    euler_angle_2 = 45.0
    euler_angle_3 = 270.0
  []
  [stress]
    type = ComputeMultipleCrystalPlasticityStress
    crystal_plasticity_models = 'uo2_cp'
    tan_mod_type = exact
  []
  [uo2_cp]
    type = UraniumDioxideCrystalPlasticityUpdate
    number_slip_systems = 24
    slip_sys_file_name = 'uo2_slip_sys.txt'
    parameter_filename = 'uo2_slip_resistance_params.txt'
    stol = 0.01
    print_state_variable_convergence_error_messages = true
  []
[]

[Postprocessors]
  [fp_zz]
    type = ElementAverageValue
    variable = fp_zz
  []
  [slip_increment_1]
    type = ElementAverageValue
    variable = slip_increment_1
  []
  [slip_increment_10]
    type = ElementAverageValue
    variable = slip_increment_10
  []
  [dislocation_density_1]
    type = ElementAverageValue
    variable = dislocation_density_1
  []
  [dislocation_density_2]
    type = ElementAverageValue
    variable = dislocation_density_2
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

  petsc_options_iname = '-pc_type -pc_asm_overlap -sub_pc_type -ksp_type -ksp_gmres_restart'
  petsc_options_value = ' asm      2              lu            gmres     200'
  nl_abs_tol = 1e-10
  nl_rel_tol = 1e-10
  nl_abs_step_tol = 1e-10

  dt = 0.05
  dtmin = 1e-5
  end_time = 1.0
[]

[Outputs]
  csv = true
  perf_graph = true
[]
