neml2_input = 'viscoplasticity_perfect'

[GlobalParams]
  displacements = 'disp_x disp_y'
[]

[Problem]
  type = ReferenceResidualProblem
  extra_tag_vectors = 'ref'
  reference_vector = 'ref'
[]

[Mesh]
  [gmg]
    type = GeneratedMeshGenerator
    dim = 2
    nx = 5
    ny = 60
    xmax = 0.005
    ymax = 0.06
  []
[]

[NEML2]
  input = '${neml2_input}.i'
  model = 'model'
  verbose = true
  mode = PARSE_ONLY
  device = 'cpu'
[]

[UserObjects]
  active = 'model input_strain input_old_strain input_old_stress'
  [input_temperature]
    type = MOOSEVariableToNEML2
    moose_variable = T
    neml2_variable = forces/T
  []
  [input_strain]
    type = MOOSERankTwoTensorMaterialPropertyToNEML2
    moose_material_property = mechanical_strain
    neml2_variable = forces/E
  []
  [input_old_strain]
    type = MOOSEOldRankTwoTensorMaterialPropertyToNEML2
    moose_material_property = mechanical_strain
    neml2_variable = old_forces/E
  []
  [input_old_stress]
    type = MOOSEOldRankTwoTensorMaterialPropertyToNEML2
    moose_material_property = small_stress
    neml2_variable = old_state/S
  []
  [input_old_ep]
    type = MOOSEOldRealMaterialPropertyToNEML2
    moose_material_property = equivalent_plastic_strain
    neml2_variable = old_state/internal/ep
  []
  [input_old_Kp]
    type = MOOSEOldSymmetricRankTwoTensorMaterialPropertyToNEML2
    moose_material_property = kinematic_plastic_strain
    neml2_variable = old_state/internal/Kp
  []
  [input_old_X1]
    type = MOOSEOldSymmetricRankTwoTensorMaterialPropertyToNEML2
    moose_material_property = backstress_1
    neml2_variable = old_state/internal/X1
  []
  [input_old_X2]
    type = MOOSEOldSymmetricRankTwoTensorMaterialPropertyToNEML2
    moose_material_property = backstress_2
    neml2_variable = old_state/internal/X2
  []
  [input_old_Ep]
    type = MOOSEOldSymmetricRankTwoTensorMaterialPropertyToNEML2
    moose_material_property = plastic_strain
    neml2_variable = old_state/internal/Ep
  []
  [input_old_gamma]
    type = MOOSEOldRealMaterialPropertyToNEML2
    moose_material_property = consistency_parameter
    neml2_variable = old_state/internal/gamma
  []

  [model]
    type = ExecuteNEML2Model
    model = model
    # add other gatherers here if needed
    gather_uos = 'input_strain input_old_strain input_old_stress'
  []
[]

[Materials]
  # add other outputs here if needed
  active = 'output_stress_jacobian'
  [output_stress_jacobian]
    type = NEML2StressToMOOSE
    execute_neml2_model_uo = model
    neml2_stress_output = state/S
    neml2_strain_input = forces/E
  []
  [output_ep]
    type = NEML2ToRealMOOSEMaterialProperty
    execute_neml2_model_uo = model
    neml2_variable = state/internal/ep
    moose_material_property = equivalent_plastic_strain
  []
  [output_Kp]
    type = NEML2ToSymmetricRankTwoTensorMOOSEMaterialProperty
    execute_neml2_model_uo = model
    neml2_variable = state/internal/Kp
    moose_material_property = kinematic_plastic_strain
  []
  [output_X1]
    type = NEML2ToSymmetricRankTwoTensorMOOSEMaterialProperty
    execute_neml2_model_uo = model
    neml2_variable = state/internal/X1
    moose_material_property = backstress_1
  []
  [output_X2]
    type = NEML2ToSymmetricRankTwoTensorMOOSEMaterialProperty
    execute_neml2_model_uo = model
    neml2_variable = state/internal/X2
    moose_material_property = backstress_2
  []
  [output_Ep]
    type = NEML2ToSymmetricRankTwoTensorMOOSEMaterialProperty
    execute_neml2_model_uo = model
    neml2_variable = state/internal/Ep
    moose_material_property = plastic_strain
  []
  [output_gamma]
    type = NEML2ToRealMOOSEMaterialProperty
    execute_neml2_model_uo = model
    neml2_variable = state/internal/gamma
    moose_material_property = consistency_parameter
  []
[]

[AuxVariables]
  [T]
  []
  [strain_yy]
    order = CONSTANT
    family = MONOMIAL
  []
[]

[AuxKernels]
  [strain_yy_aux]
    type = MaterialRankTwoTensorAux
    i = 1
    j = 1
    property = mechanical_strain
    variable = strain_yy
  []
[]

[Physics]
  [SolidMechanics]
    [QuasiStatic]
      [all]
        strain = SMALL
        new_system = true
        add_variables = true
        formulation = TOTAL
        volumetric_locking_correction = true
        extra_vector_tags = 'ref'
        generate_output = 'cauchy_stress_xx cauchy_stress_yy cauchy_stress_xy mechanical_strain_xx mechanical_strain_yy mechanical_strain_xy'
      []
    []
  []
[]

[Functions]
  [top_load]
    type = ParsedFunction
    expression = 'if(t<100.0, 0.2539e6*t , 25.39e6)'
  []
[]

[BCs]
  [fix_bottom]
    type = DirichletBC
    variable = disp_y
    boundary = bottom
    value = 0
    extra_vector_tags = 'ref'
  []
  [fix_left]
    type = DirichletBC
    variable = disp_x
    boundary = left
    value = 0
    extra_vector_tags = 'ref'
  []
  [load_top]
    # type = FunctionDirichletBC
    # variable = disp_y
    # boundary = top
    # function = t
    # preset = false
    type = ADFunctionNeumannBC
    use_displaced_mesh = true
    variable = disp_y
    boundary = 'top'
    function = top_load
    extra_vector_tags = 'ref'
  []
[]

[Postprocessors]
  [strain]
    type = ElementAverageValue
    variable = strain_yy
    execute_on = 'initial timestep_end'
  []
  [load]
    type = FunctionValuePostprocessor
    function = top_load
  []
[]

[Executioner]
  type = Transient
  solve_type = NEWTON
  petsc_options_iname = '-pc_type'
  petsc_options_value = 'lu'
  automatic_scaling = true
  line_search = 'none'
  l_max_its = 50
  nl_max_its = 25
  nl_rel_tol = 1e-6
  nl_abs_tol = 5e-8
  start_time = 0
  end_time = 432000.0 ## 120 hours
  # end_time = 2e5
  dtmin = 1e-2
  [TimeStepper]
    type = IterationAdaptiveDT
    dt = 1.0e3
    optimal_iterations = 3
    iteration_window = 1
  []
  residual_and_jacobian_together = true
[]

[Outputs]
  file_base = '${neml2_input}'
  exodus = true
  csv = true
[]
