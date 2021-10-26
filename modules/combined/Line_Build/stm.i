[StochasticTools]
[]

[Distributions]
  [speed_dist]
    type = Uniform
    lower_bound = 1
    upper_bound = 4
  []
  [power_dist]
    type = Uniform
    lower_bound = 300
    upper_bound = 500
  []
[]

[Samplers]
  [pc_sampler]
    type = MonteCarlo
    num_rows = 1000
    distributions = 'speed_dist power_dist'
    execute_on = 'PRE_MULTIAPP_SETUP'
  []
[]

[MultiApps]
  [pc_master]
    type = SamplerFullSolveMultiApp
    input_files = master_app.i
    sampler = pc_sampler
  []
[]

[Controls]
  [pc_cmdline]
    type = MultiAppCommandLineControl
    multi_app = pc_master
    sampler = pc_sampler
    param_names = 'speed power'
  []
[]

[VectorPostprocessors]
  [data]
    type = SamplerData
    sampler = pc_sampler
    execute_on = 'INITIAL TIMESTEP_END'
    parallel_type = DISTRIBUTED
  []
[]

[Outputs]
  file_base = 'output/sampler_data'
  execute_on = 'INITIAL'
  csv = true
[]
