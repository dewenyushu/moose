[StochasticTools]
[]

[Samplers]
  [sample]
    type = CSVSampler
    samples_file = './inputs/single_track_samples.csv'
    column_indices = '0 1 2 3'
    execute_on = 'PRE_MULTIAPP_SETUP'
  []
[]

[MultiApps]
  [sub]
    type = SamplerFullSolveMultiApp
    sampler = sample
    input_files = 'thermal.i'
    mode = batch-reset
  []
[]

[Controls]
  [cmdline]
    type = MultiAppSamplerControl
    multi_app = sub
    sampler = sample
    param_names = 'MW end_time factor r_factor'
  []
[]

[VectorPostprocessors]
  [data]
    type = SamplerData
    sampler = sample
    execute_on = 'initial timestep_end'
  []
[]

[Outputs]
  execute_on = timestep_end
  csv = true
[]

