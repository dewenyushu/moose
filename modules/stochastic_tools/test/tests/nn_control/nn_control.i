[Mesh]
  [square]
    type = GeneratedMeshGenerator
    nx = 2
    ny = 2
    dim = 2
  []
[]

[Variables]
  [u]
  []
[]

[Kernels]
  [diff]
    type = Diffusion
    variable = u
  []
[]

[BCs]
  [left]
    type = PostprocessorDirichletBC
    variable = u
    boundary = 3
    postprocessor = received_bc
  []

  [right]
    type = DirichletBC
    variable = u
    boundary = 1
    value = 1
  []
[]

[Executioner]
  type = Transient

  solve_type = 'PJFNK'
  start_time = 0.0
  num_steps = 20
  dt = 1

  nl_abs_tol = 1e-10
  l_tol = 1e-05
  line_search = 'none'

  # For picard tests
  # picard_abs_tol = 1e-3
[]

[Postprocessors]
  [integral]
    type = ElementIntegralVariablePostprocessor
    variable = u
    execute_on = 'initial timestep_end'
  []
  [received_bc]
    type = Receiver
    default = 0
  []
[]

[Controls]
  [integral_value]
    type = BasicNNControl
    postprocessor = integral
    target = 1.5
    parameter_pp = 'received_bc'
    K_integral = -1
    K_proportional = -1
    K_derivative = -0.1
    execute_on = 'initial timestep_begin'
  []
[]

[Outputs]
  file_base = out
  exodus = false
  csv = true
[]
