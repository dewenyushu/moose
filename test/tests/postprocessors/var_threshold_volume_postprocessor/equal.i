[Mesh]
  [gmg]
    type = GeneratedMeshGenerator
    dim = 2
    nx = 3
    ny = 3
  []
[]

[Variables]
  [dummy]
    order = FIRST
    family = LAGRANGE
  []
[]

[AuxVariables]
  [u]
    order = FIRST
    family = LAGRANGE
  []
[]

[Kernels]
  [null]
    type = NullKernel
    variable = dummy
  []
[]

[AuxKernels]
  [func_aux]
    type = FunctionAux
    function = func_u
    variable = u
  []
[]

[Functions]
  [func_u]
    type = ParsedFunction
    expression = 'x'
  []
[]

[Postprocessors]
  [var_volume]
    type = VarThresholdVolumePostprocessor
    coupled_var = u
    criterion_type = EQUAL
    threshold = 0.5
    execute_on = TIMESTEP_END
  []
[]

[Executioner]
  type = Transient
  num_steps = 1
  dt = 0.5
  solve_type = PJFNK
  petsc_options_iname = '-pc_type -pc_hypre_type'
  petsc_options_value = 'hypre boomeramg'
[]

[Outputs]
  csv = true
[]
