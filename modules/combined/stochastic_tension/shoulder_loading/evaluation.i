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
  [right_circle_x]
    type = Uniform
    lower_bound = 5.7
    upper_bound = 5.8
  []
  [left_circle_x]
    type = Uniform
    lower_bound = -5.8
    upper_bound = -5.7
  []
  # [fillet_radius]
  #   type = Uniform
  #   lower_bound = 0.2
  #   upper_bound = 0.4
  # []
[]

[GlobalParams]
  sampler = quad
  distributions = 'right_circle_y left_circle_y right_circle_x left_circle_x'
[]

[Samplers]
  [sample]
    type = MonteCarlo
    num_rows = 200
    execute_on = PRE_MULTIAPP_SETUP
  []
[]

[Surrogates]
  [pc_strain_xx_bot]
    type = PolynomialChaos
    filename = ./5param_surrogates/train_1param_trainer_poly_chaos_strain_xx_bot.rd
  []
  [pc_strain_xx_center]
    type = PolynomialChaos
    filename = ./5param_surrogates/train_1param_trainer_poly_chaos_strain_xx_center.rd
  []
  [pc_strain_xx_top]
    type = PolynomialChaos
    filename = ./5param_surrogates/train_1param_trainer_poly_chaos_strain_xx_top.rd
  []
  [pc_stress_xx_bot]
    type = PolynomialChaos
    filename = ./5param_surrogates/train_1param_trainer_poly_chaos_stress_xx_bot.rd
  []
  [pc_stress_xx_center]
    type = PolynomialChaos
    filename = ./5param_surrogates/train_1param_trainer_poly_chaos_stress_xx_center.rd
  []
  [pc_stress_xx_top]
    type = PolynomialChaos
    filename = ./5param_surrogates/train_1param_trainer_poly_chaos_stress_xx_top.rd
  []
[]

[Reporters]
  [eval_strain_xx_bot]
    type = EvaluateSurrogate
    model = pc_strain_xx_bot
  []
  [eval_strain_xx_center]
    type = EvaluateSurrogate
    model = pc_strain_xx_center
  []
  [eval_strain_xx_top]
    type = EvaluateSurrogate
    model = pc_strain_xx_top
  []
  [eval_strain_xx_left]
    type = EvaluateSurrogate
    model = pc_strain_xx_left
  []
  [eval_strain_xx_right]
    type = EvaluateSurrogate
    model = pc_strain_xx_right
  []
  [eval_stress_xx_bot]
    type = EvaluateSurrogate
    model = pc_stress_xx_bot
  []
  [eval_stress_xx_center]
    type = EvaluateSurrogate
    model = pc_stress_xx_center
  []
  [eval_stress_xx_top]
    type = EvaluateSurrogate
    model = pc_stress_xx_top
  []
  [eval_stress_xx_left]
    type = EvaluateSurrogate
    model = pc_stress_xx_left
  []
  [eval_stress_xx_right]
    type = EvaluateSurrogate
    model = pc_stress_xx_right
  []
[]

[VectorPostprocessors]
  [stats_strain_xx_bot]
    type = PolynomialChaosStatistics
    pc_name = pc_strain_xx_bot
    compute = 'mean stddev'
  []
  [sobol_strain_xx_bot]
    type = PolynomialChaosSobolStatistics
    pc_name = pc_strain_xx_bot
    sensitivity_order = 'first second total'
  []
  [stats_strain_xx_center]
    type = PolynomialChaosStatistics
    pc_name = pc_strain_xx_center
    compute = 'mean stddev'
  []
  [sobol_strain_xx_center]
    type = PolynomialChaosSobolStatistics
    pc_name = pc_strain_xx_center
    sensitivity_order = 'first second total'
  []
  [stats_strain_xx_top]
    type = PolynomialChaosStatistics
    pc_name = pc_strain_xx_top
    compute = 'mean stddev'
  []
  [sobol_strain_xx_top]
    type = PolynomialChaosSobolStatistics
    pc_name = pc_strain_xx_top
    sensitivity_order = 'first second total'
  []
  [stats_strain_xx_left]
    type = PolynomialChaosStatistics
    pc_name = pc_strain_xx_left
    compute = 'mean stddev'
  []
  [sobol_strain_xx_left]
    type = PolynomialChaosSobolStatistics
    pc_name = pc_strain_xx_left
    sensitivity_order = 'first second total'
  []
  [stats_strain_xx_right]
    type = PolynomialChaosStatistics
    pc_name = pc_strain_xx_right
    compute = 'mean stddev'
  []
  [sobol_strain_xx_right]
    type = PolynomialChaosSobolStatistics
    pc_name = pc_strain_xx_right
    sensitivity_order = 'first second total'
  []
  [stats_stress_xx_bot]
    type = PolynomialChaosStatistics
    pc_name = pc_stress_xx_bot
    compute = 'mean stddev'
  []
  [sobol_stress_xx_bot]
    type = PolynomialChaosSobolStatistics
    pc_name = pc_stress_xx_bot
    sensitivity_order = 'first second total'
  []
  [stats_stress_xx_center]
    type = PolynomialChaosStatistics
    pc_name = pc_stress_xx_center
    compute = 'mean stddev'
  []
  [sobol_stress_xx_center]
    type = PolynomialChaosSobolStatistics
    pc_name = pc_stress_xx_center
    sensitivity_order = 'first second total'
  []
  [stats_stress_xx_top]
    type = PolynomialChaosStatistics
    pc_name = pc_stress_xx_top
    compute = 'mean stddev'
  []
  [sobol_stress_xx_top]
    type = PolynomialChaosSobolStatistics
    pc_name = pc_stress_xx_top
    sensitivity_order = 'first second total'
  []
  [stats_stress_xx_left]
    type = PolynomialChaosStatistics
    pc_name = pc_stress_xx_left
    compute = 'mean stddev'
  []
  [sobol_stress_xx_left]
    type = PolynomialChaosSobolStatistics
    pc_name = pc_stress_xx_left
    sensitivity_order = 'first second total'
  []
  [stats_stress_xx_right]
    type = PolynomialChaosStatistics
    pc_name = pc_stress_xx_right
    compute = 'mean stddev'
  []
  [sobol_stress_xx_right]
    type = PolynomialChaosSobolStatistics
    pc_name = pc_stress_xx_right
    sensitivity_order = 'first second total'
  []
[]
[Outputs]
  csv = true
  execute_on = timestep_end
[]
