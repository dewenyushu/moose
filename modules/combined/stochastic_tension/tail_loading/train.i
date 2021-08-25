[StochasticTools]
[]

[Distributions]
  [right_circle_y]
    type = Uniform
    lower_bound = -0.05
    upper_bound = 0.05
  []
  [left_circle_y]
    type = Uniform
    lower_bound = -0.05
    upper_bound = 0.05
  []
  [right_circle_r]
    type = Uniform
    lower_bound = 0.725
    upper_bound = 0.775
  []
  [left_circle_r]
    type = Uniform
    lower_bound = 0.725
    upper_bound = 0.775
  []
  [fillet_radius]
    type = Uniform
    lower_bound = 0.25
    upper_bound = 0.35
  []
[]

[GlobalParams]
  sampler = quad
  distributions = 'right_circle_y left_circle_y right_circle_r left_circle_r fillet_radius'
[]

[Samplers]
  [quad]
    type = Quadrature
    order = 4
    sparse_grid = smolyak
    execute_on = PRE_MULTIAPP_SETUP
  []
[]

[MultiApps]
  [sub]
    type = SamplerFullSolveMultiApp
    input_files = sub_tail_loading.i
    mode = batch-reset
  []
[]

[Controls]
  [param]
    type = MultiAppCommandLineControl
    multi_app = sub
    param_names = 'Mesh/gmsh/right_circle_y Mesh/gmsh/left_circle_y Mesh/gmsh/right_circle_r '
                  'Mesh/gmsh/left_circle_r Mesh/gmsh/fillet_radius'
  []
[]

[Transfers]
  [data]
    type = SamplerReporterTransfer
    multi_app = sub
    stochastic_reporter = storage
    from_reporter = 'stress_xx_center/value strain_xx_center/value stress_xx_top/value '
                    'strain_xx_top/value stress_xx_bot/value strain_xx_bot/value'
  []
[]

[Reporters]
  [storage]
    type = StochasticReporter
  []
[]

[Trainers]
  [poly_chaos_stress_xx_center]
    type = PolynomialChaosTrainer
    order = 4
    response = storage/data:stress_xx_center:value
    # response_type = vector_real
    converged_reporter = storage/data:converged
    skip_unconverged_samples = true
  []
  [poly_chaos_strain_xx_center]
    type = PolynomialChaosTrainer
    order = 4
    response = storage/data:strain_xx_center:value
    # response_type = vector_real
    converged_reporter = storage/data:converged
    skip_unconverged_samples = true
  []
  [poly_chaos_stress_xx_top]
    type = PolynomialChaosTrainer
    order = 4
    response = storage/data:stress_xx_top:value
    # response_type = vector_real
    converged_reporter = storage/data:converged
    skip_unconverged_samples = true
  []
  [poly_chaos_strain_xx_top]
    type = PolynomialChaosTrainer
    order = 4
    response = storage/data:strain_xx_top:value
    # response_type = vector_real
    converged_reporter = storage/data:converged
    skip_unconverged_samples = true
  []
  [poly_chaos_stress_xx_bot]
    type = PolynomialChaosTrainer
    order = 4
    response = storage/data:stress_xx_bot:value
    # response_type = vector_real
    converged_reporter = storage/data:converged
    skip_unconverged_samples = true
  []
  [poly_chaos_strain_xx_bot]
    type = PolynomialChaosTrainer
    order = 4
    response = storage/data:strain_xx_bot:value
    # response_type = vector_real
    converged_reporter = storage/data:converged
    skip_unconverged_samples = true
  []
[]

[Outputs]
  execute_on = timestep_end
  [trainer]
    type = SurrogateTrainerOutput
    trainers = 'poly_chaos_stress_xx_center  poly_chaos_strain_xx_center poly_chaos_stress_xx_top '
               'poly_chaos_strain_xx_top poly_chaos_stress_xx_bot poly_chaos_strain_xx_bot'
  []
[]
