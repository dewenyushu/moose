[StochasticTools]
[]

[Postprocessors]
  # [integral]
  #   type = ElementIntegralVariablePostprocessor
  #   variable = u
  #   execute_on = 'initial timestep_end'
  # []
  # [received_bc]
  #   type = Receiver
  #   default = 0
  # []
  [power]
    type = Receiver
    default = 0
  []
  [max_temperature]
    type = Receiver
    default = 300
  []
[]

[Controls]
  [integral_value]
    type = BasicNNControl
    input_pp = 'max_temperature'
    output_pp = 'power'
    target = 1.5
    K_integral = -1
    K_proportional = -1
    K_derivative = -0.1
    execute_on = 'initial timestep_begin'
  []
[]

[MultiApps]
  [thermal_am_app]
    type = TransientMultiApp
    input_files = 'thermal_app.i'
    cli_args = 'Executioner/num_steps=20'
  []
[]

[Transfers]
  [to_sub_power]
    type = MultiAppPostprocessorTransfer
    direction = to_multiapp
    multi_app = thermal_am_app
    to_postprocessor = 'received_power'
    from_postprocessor = 'power'
  []
  [from_sub_max_temperature]
    type = MultiAppPostprocessorTransfer
    direction = from_multiapp
    multi_app = thermal_am_app
    to_postprocessor = 'max_temperature'
    from_postprocessor = 'bead_max_temperature'
    reduction_type = AVERAGE
  []
[]

[Outputs]
  file_base = out
  exodus = false
  csv = true
[]
