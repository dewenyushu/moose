[Mesh]
  displacements = 'ux uy'
  [gmg]
    type = GeneratedMeshGenerator
    dim = 2
    nx = 11
    ny = 11
    xmin = -4
    xmax = 4
    ymin = -4
    ymax = 4
  []
[]

[Variables]
  # adjoint
  [ux]
  []
  [uy]
  []
[]

[AuxVariables]
  [dummy]
  []
  [T]
  []
  # displacement variables to be transferred from the forward app
  # we use them to compute stress and stress derivative wrt E
  # let them be 0 for test purpose
  [state_x]
  []
  [state_y]
  []
[]

[Modules]
  [TensorMechanics]
    [Master]
      [all]
        strain = SMALL
        new_system = true
        add_variables = true
        formulation = TOTAL
        incremental = true
        volumetric_locking_correction = false
        generate_output = 'cauchy_stress_xx'
        displacements = 'ux uy'
        # add base name to distinguish between forward and adjoint
        base_name = 'adjoint'
      []
      displacements = 'ux uy'
    []
  []
[]

[BCs]
  [bottom_x]
    type = DirichletBC
    variable = ux
    boundary = bottom
    value = 0.0
  []
  [bottom_y]
    type = DirichletBC
    variable = uy
    boundary = bottom
    value = 0.0
  []
  [top_x]
    type = NeumannBC
    variable = ux
    boundary = top
    value = 1.0
  []
  [top_y]
    type = NeumannBC
    variable = uy
    boundary = top
    value = 1.0
  []
[]

[NEML2]
  # two elasticity models are listed inside "elasticity.i" for forward and adjoint, respectively
  input = 'elasticity.i'
  model = 'adjoint_elasticity_model'
  verbose = false
  temperature = 'T'
  mode = PARSE_ONLY
  device = 'cpu'
[]

[Materials]
  [adjoint_stress]
    type = CauchyStressFromNEML2Receiver
    neml2_uo = adjoint_neml2_stress_UO
    base_name = 'adjoint'
  []
  # forward stress is not used
  # [forward_stress]
  #   type = CauchyStressFromNEML2Receiver
  #   neml2_uo = forward_neml2_stress_UO
  #   base_name = 'forward'
  # []
  [forward_strain]
    type = ComputeSmallStrain
    displacements = 'state_x state_y'
    base_name = 'forward'
  []
  # adjoint and forward use the same young's modulus value
  # need to be separated between adjoint and forward?
  [E_material]
    type = GenericFunctionMaterial
    prop_names = 'E_material'
    prop_values = E
  []
[]

[Functions]
  [E]
    type = NearestReporterCoordinatesFunction
    x_coord_name = parametrization/coordx
    y_coord_name = parametrization/coordy
    value_name = parametrization/youngs_modulus
  []
[]

[Reporters]
  [measure_data]
    type = OptimizationData
    variable = ux
  []
  [parametrization]
    type = ConstantReporter
    real_vector_names = 'coordx coordy youngs_modulus'
    real_vector_values = '0 1 2; 0 1 2; 7.5 7.5 7.5'
  []
[]

[UserObjects]
  # forward stress derivative,to be used in gradient calculation
  [forward_E_batch_material]
    type = BatchPropertyDerivativeRankTwoTensorReal
    material_property = 'E_material'
  []
  [forward_neml2_stress_UO]
    type = CauchyStressFromNEML2UO
    temperature = 'T'
    model = 'forward_elasticity_model'
    scalar_material_property_names = 'E'
    scalar_material_property_values = 'forward_E_batch_material'
    # use forward strain calculated from state_x and state_y
    mechanical_strain = 'forward_mechanical_strain'
  []
  # adjoint stress derivative, not used
  [adjoint_E_batch_material]
    type = BatchPropertyDerivativeRankTwoTensorReal
    material_property = 'E_material'
  []
  [adjoint_neml2_stress_UO]
    type = CauchyStressFromNEML2UO
    temperature = 'T'
    model = 'adjoint_elasticity_model'
    scalar_material_property_names = 'E'
    scalar_material_property_values = 'adjoint_E_batch_material'
    # use adjoint strain calculated tensor mechanics module
    mechanical_strain = 'adjoint_mechanical_strain'
  []
[]

[VectorPostprocessors]
  [grad_youngs_modulus]
    type = AdjointStrainStressGradNEML2InnerProduct
    neml2_uo = 'forward_neml2_stress_UO'
    stress_derivative = 'forward_E_batch_material' # this should come from the forward
    adjoint_strain_name = 'adjoint_mechanical_strain'
    variable = dummy
    function = E
  []
[]

[Executioner]
  type = Steady
  solve_type = NEWTON
  petsc_options_iname = '-pc_type'
  petsc_options_value = 'lu'
  # required for NEML2 material models
  residual_and_jacobian_together = true
[]

[Postprocessors]
  [point1]
    type = PointValue
    point = '-1.0 -1.0 0.0'
    variable = ux
    execute_on = TIMESTEP_END
  []
  [point2]
    type = PointValue
    point = '-1.0 0.0 0.0'
    variable = ux
    execute_on = TIMESTEP_END
  []
  [point3]
    type = PointValue
    point = '-1.0 1.0 0.0'
    variable = ux
    execute_on = TIMESTEP_END
  []
[]

[Outputs]
  console = true
[]
