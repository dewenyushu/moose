[Mesh]
  [gmg]
    type = GeneratedMeshGenerator
    dim = 3
    xmin = -2
    xmax = 2
    ymin = -2
    ymax = 2
    zmin = 0
    zmax = 2
    nx = 16
    ny = 16
    nz = 8
  []
[]

[Variables]
  [adjoint_T]
  []
[]

[Kernels]
  [conduction]
    type = MatDiffusion
    diffusivity = diffusivity
    variable = adjoint_T
  []
[]

[Reporters]
  [misfit]
    type = OptimizationData
  []
  [data]
    type = ConstantReporter
    real_vector_names = 'diffusivity'
    real_vector_values = '10 10'
  []
[]

[DiracKernels]
  [pt]
    type = ReporterPointSource
    variable = adjoint_T
    x_coord_name = misfit/measurement_xcoord
    y_coord_name = misfit/measurement_ycoord
    z_coord_name = misfit/measurement_zcoord
    value_name = misfit/misfit_values
  []
[]

[BCs]
  [bottom]
    type = DirichletBC
    variable = adjoint_T
    boundary = back
    value = 0
  []
[]

[AuxVariables]
  [temperature_forward]
  []
[]

[Functions]
  [diffusivity_function]
    type = ParameterMeshFunction
    exodus_mesh = parameter_mesh_in.e
    family = MONOMIAL
    order = CONSTANT
    parameter_name = data/diffusivity
  []
[]

[Materials] #same material as what was used in the forward model
  [mat]
    type = GenericFunctionMaterial
    prop_names = diffusivity
    prop_values = diffusivity_function
  []
[]

[VectorPostprocessors]
  [gradvec]
    type = ElementOptimizationDiffusionCoefFunctionInnerProduct
    variable = adjoint_T
    forward_variable = temperature_forward
    function = diffusivity_function
  []
[]

[Executioner]
  type = Steady
  solve_type = NEWTON
  petsc_options_iname = '-pc_type'
  petsc_options_value = 'lu'
  nl_forced_its = 1
  line_search = none
  nl_abs_tol = 1e-8
[]

[Outputs]
  console = false
  file_base = 'adjoint'
[]
