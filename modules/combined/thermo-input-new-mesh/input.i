[Mesh]
  [file_mesh]
    type = FileMeshGenerator
    file = 'mesh.e'
  []
[]

[Variables]
  # Temperature variable
  [temp]
    initial_condition = 273
  []
[]

# Adding the thermal problem
[Kernels]
  [heat_time]
    type = ADHeatConductionTimeDerivative
    specific_heat = specific_heat
    density_name = density
    variable = temp
  []
  [heat_conduction]
    type = ADHeatConduction
    variable = temp
    thermal_conductivity = thermal_conductivity
  []
[]

[Materials]
  # Definition of the thermal conductivity
  [thermal_conductivity_w]
    type = ADConstantMaterial
    property_name = thermal_conductivity
    value = 175 # from google
    block = 'W'
  []
  [thermal_conductivity_steel] # Eurofer97 - temperature dependent
    type = ADCoupledValueFunctionMaterial
    prop_name = thermal_conductivity
    function = '-1.158e-10*pow(x, 4) + 2.550e-7*pow(x,3) - 1.621e-4*x*x + 3.575e-2*x + 27.699'
    v = temp
    block = 'Steel'
  []
  [thermal_conductivity_sic]
    type = ADCoupledValueFunctionMaterial
    prop_name = thermal_conductivity
    function = 'if(x>300, 1.0/(-0.0003+1.05e-5*x), 351)'
    v = temp
    block = 'SIC'
  []
  # Definition of the specific heat
  [specific_heat_w]
    type = ADConstantMaterial
    property_name = specific_heat
    value = 130 # J/kg/K, from google
    block = 'W'
  []
  [specific_heat_steel] # Eurofer97
    type = ADConstantMaterial
    property_name = specific_heat
    value = 420 # from google
    block = 'Steel'
  []
  [specific_heat_sic]
    type = ADConstantMaterial
    property_name = specific_heat
    value = 670 # from google
    block = 'SIC'
  []
  # Definition of the density
  [density_w]
    type = ADConstantMaterial
    property_name = density
    value = 19350 # from google
    block = 'W'
  []
  [density_steel] # Eurofer97
    type = ADConstantMaterial
    property_name = density
    value = 7850 # from google
    block = 'Steel'
  []
  [density_sic]
    type = ADConstantMaterial
    property_name = density
    value = 3210 # from google
    block = 'SIC'
  []
[]

[BCs]
  [bottom]  # Adiabatic Wall (or temp sink?)
    type = NeumannBC
    variable = temp
    boundary = bottom
    value = 0
  []
  [top]
    type = NeumannBC
    variable = temp
    boundary = top
    value = 0
  []
  [left] # Fixed heat flux
    type = NeumannBC
    variable = temp
    boundary = left
    value = 0.25e6 #W/m^2
  []
  [right] # Fixed temperature (higher than room temp)
    type = DirichletBC
    variable = temp
    boundary = right
    value = 320
  []
  [he]
    type = ConvectiveHeatFluxBC
    variable = temp
    boundary = he
    T_infinity = 620 #K
    heat_transfer_coefficient = 0 # TBD
  []
  [flow_channel]
    type = ConvectiveHeatFluxBC
    variable = temp
    boundary = flow_channel
    T_infinity = 800 #K TBD
    heat_transfer_coefficient = 0 # TBD
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

  solve_type = NEWTON
  petsc_options_iname = '-pc_type'
  petsc_options_value = 'lu'

  nl_rel_tol = 1e-6
  nl_max_its = 20

  l_tol = 1e-8

  dt = 1
  dtmin = 1e-4

  num_steps = 20
[]

[Outputs]
  file_base = 'thermo_out'
  exodus = true
[]
