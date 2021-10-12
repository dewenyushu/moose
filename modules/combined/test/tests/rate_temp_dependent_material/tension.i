[GlobalParams]
  displacements = 'disp_x disp_y disp_z'
[]

[Problem]
  kernel_coverage_check = false
[]

[Mesh]
  [mesh]
    type = GeneratedMeshGenerator
    dim = 3
  []
[]

[AuxVariables]
  [./temp]
    initial_condition = 600
  [../]
[]

[Modules/TensorMechanics/Master/All]
  strain = FINITE
  incremental = true
  add_variables = true
  generate_output = 'stress_xx stress_yy stress_zz stress_xy stress_yz stress_xz strain_yy strain_xx '
                    'strain_zz strain_xy strain_xz strain_yz'
  use_automatic_differentiation = true
[]

[BCs]
  [uz_bottom_fix]
    type = ADDirichletBC
    variable = disp_z
    boundary = back
    value = 0.0
  []

  [uz_top_disp]
    type = FunctionDirichletBC
    variable = disp_z
    boundary = front
    function = '0.05*t'
  []
[]

[Materials]
  [E]
    type = ADPiecewiseLinearInterpolationMaterial
    x = '294.994  1671.48  1721.77'
    y = '210e9 80e9 6.16016e9' # Pa
    property = youngs_modulus
    variable = temp
    extrapolation = false
  []
  [nu]
    type = ADPiecewiseLinearInterpolationMaterial
    x = '294.994 1669.62 1721.77'
    y = '0.246407   0.36961  0.513347'
    property = poissons_ratio
    variable = temp
    extrapolation = false
  []
  [elasticity_tensor]
    type = ADComputeVariableIsotropicElasticityTensor
    youngs_modulus = youngs_modulus
    poissons_ratio = poissons_ratio
  []

  [radial_return_stress]
    type = ADComputeMultipleInelasticStress
    inelastic_models = 'rate_temp_plas'
  []
  [rate_temp_plas]
    type = ADRateTempDependentStressUpdate
    temperature = temp
    Ex = '294.994  1671.48  1721.77'
    Ey = '210e9 80e9 6.16016e9'
    nux = '294.994 1669.62 1721.77'
    nuy = '0.246407   0.36961  0.513347'
    # n2 = ${T_melt}
    absolute_tolerance = 1e-8
    use_substep = false
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

  automatic_scaling = true

  #Preconditioned JFNK (default)
  solve_type = 'NEWTON'

  # # petsc_options = '-snes_ksp'
  # petsc_options_iname = '-pc_type -ksp_type'
  # petsc_options_value = 'lu  preonly'

  petsc_options_iname = '-ksp_type -pc_type -pc_factor_mat_solver_package -pc_factor_shift_type -pc_factor_shift_amount'
  petsc_options_value = 'preonly lu       superlu_dist nonzero 1e-10'

  line_search = 'none'

  l_max_its = 100
  nl_max_its = 20
  nl_rel_tol = 1e-6
  nl_abs_tol = 1e-6
  # l_tol = 1e-5

  start_time = 0.0
  end_time = 20
  dt = 0.5
  dtmin = 1e-6
[]

[Postprocessors]
  [stress_zz]
    type = ElementAverageValue
    variable = stress_zz
  []
  [strain_zz]
    type = ElementAverageValue
    variable = strain_zz
  []
[]

[Outputs]
  [exodus]
    type = Exodus
    interval = 5
  []
  [csv]
    type = CSV
  []
  [console]
    type = Console
    max_rows = 5
  []
[]
