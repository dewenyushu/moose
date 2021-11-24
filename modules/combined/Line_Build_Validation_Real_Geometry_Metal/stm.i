[StochasticTools]
[]

# [Distributions]
#   [speed_dist]
#     type = Uniform
#     lower_bound = 0.5
#     upper_bound = 4
#   []
#   [power_dist]
#     type = Uniform
#     lower_bound = 200
#     upper_bound = 400
#   []
# []

[Samplers]
  [pc_sampler]
    type = CSVSampler
    samples_file = 'Input_Params_Part.csv'
    column_names = 'power speed r dt'
    execute_on = 'PRE_MULTIAPP_SETUP'
  []
[]

[MultiApps]
  [pc_master]
    type = SamplerFullSolveMultiApp
    input_files = thermal_app_melt_diff_laser_factor.i
    sampler = pc_sampler
    mode = batch-reset
    ignore_solve_not_converge = true
  []
[]

[Controls]
  [pc_cmdline]
    type = MultiAppCommandLineControl
    multi_app = pc_master
    sampler = pc_sampler
    param_names = 'power speed r dt'
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

# [Outputs]
#   file_base = 'output/sampler_data'
#   execute_on = 'INITIAL'
#   csv = true
# []
