T_melt = 1700
T_room = 300

[Problem]
  kernel_coverage_check = false
  # material_coverage_check = false
[]

[Mesh]
  [./mesh]
    type = FileMeshGenerator
    # file = cylinder_product.e
    # file = ./input_mesh/solid_cylinder_product.e
    file = ./input_mesh/solid_cylinder_17_substrate_100x100.e
    # file = ./input_mesh/solid_cylinder_17_substrate_100x100_fine_v0.e
  [../]
  [./add_set1]
    type = SubdomainBoundingBoxGenerator
    input = mesh
    block_id = 2
    bottom_left = '-50 -50 0'
    top_right = '50 50 5'
  [../]
  [./add_set2]
    type = SubdomainBoundingBoxGenerator
    input = add_set1
    block_id = 1
    bottom_left = '50 50 5'
    top_right = '50 50 25'
  [../]
  [./middle]
    type = SideSetsAroundSubdomainGenerator
    input = add_set2
    block = 2
    new_boundary = 'moving_boundary'
    normal = '0 0 1'
  []
[]

[Variables]
  [./temp]
    initial_condition = ${T_room}
    block = '2'
  [../]
[]

[Kernels]
  [./time]
    type = ADHeatConductionTimeDerivative
    variable = temp
    block = '2'
  [../]
  [./heat_conduct]
    type = ADHeatConduction
    variable = temp
    use_displaced_mesh = true
    thermal_conductivity = thermal_conductivity
    block = '2'
  [../]
  [./heatsource]
    type = ADMatHeatSource
    material_property = volumetric_heat
    variable = temp
    scalar = 1
    use_displaced_mesh = true
    block = '2'
  [../]
[]

[Functions]
  [./heat_source_x]
    type = PiecewiseLinear
    data_file = ./input_mat_params/path_eq_t_x.csv
    format = columns
    scale_factor = 1.0
  [../]
  [./heat_source_y]
    type = PiecewiseLinear
    data_file = ./input_mat_params/path_eq_t_y.csv
    format = columns
    scale_factor = 1.0
  [../]
  [./heat_source_z]
    type = PiecewiseLinear
    data_file = ./input_mat_params/path_eq_t_z.csv
    format = columns
    scale_factor = 1.0
  [../]
  [./specific_heat]
    type = PiecewiseLinear
    data_file = ./input_mat_params/Specific_Heat.csv
    format = columns
    scale_factor = 1.0
  [../]
  [./thermal_conductivity]
    type = PiecewiseLinear
    data_file = ./input_mat_params/Thermal_Conductivity.csv
    format = columns
    scale_factor = 0.05e-3
  [../]
[]

[BCs]
  [./temp_bottom_fix]
    type = ADDirichletBC
    variable = temp
    boundary = 1
    value = ${T_room}
  [../]


  # [./dirichlet]
  #   type = ADDirichletBC
  #   variable = temp
  #   boundary = 'moving_boundary'
  #   value = 500
  # [../]
  [./convective]
    type = ConvectiveFluxFunction # Convective flux, e.g. q'' = h*(Tw - Tf)
    variable = temp
    boundary = 'moving_boundary'
    coefficient = 2e-5
    T_infinity = ${T_room}
  [../]
  [./convective_substrate]
    type = ConvectiveFluxFunction # Convective flux, e.g. q'' = h*(Tw - Tf)
    variable = temp
    boundary = 2
    coefficient = 2e-5
    T_infinity = ${T_room}
  [../]
[]

[Materials]
  [./thermal_expansion_strain]
    type = ADComputeThermalExpansionEigenstrain
    stress_free_temperature = ${T_melt}
    thermal_expansion_coeff = 0.5e-5
    temperature = temp
    eigenstrain_name = thermal_eigenstrain
    # activated_elem_aux = activated_elem
    block = '2'
  [../]
  [./volumetric_heat]
    type = FunctionPathEllipsoidHeatSource
    a = 4.0
    b = 4.0
    c = 4.0
    power = 450
    efficienty = 0.36
    factor = 1.0
    function_x= heat_source_x
    function_y= heat_source_y
    function_z= heat_source_z
  [../]
  [./density]
    type = ADDensity
    density = 7.609e-6
    block = '2'
  [../]
  [./heat]
    type = ADHeatConductionMaterial
    specific_heat_temperature_function = specific_heat
    thermal_conductivity_temperature_function = thermal_conductivity
    temp = temp
    block = '2'
  [../]
[]

[Preconditioning]
  [./smp]
    type = SMP
    full = true
  [../]
[]


[Executioner]
  type = Transient

  automatic_scaling = true

  #Preconditioned JFNK (default)
  solve_type = 'NEWTON'

  # petsc_options = '-snes_ksp'
  petsc_options_iname = '-pc_type -ksp_type'
  petsc_options_value = 'lu  preonly'

  line_search = 'none'

  l_max_its = 10
  nl_max_its = 20
  nl_rel_tol = 1e-4
  nl_abs_tol = 1e-12
  # l_tol = 1e-5

  start_time = 0.0
  end_time = 2.0
  dt = 0.1
  dtmin = 1e-4
[]

[AuxVariables]
  [./temp_aux]
    order = FIRST
    family = LAGRANGE
  [../]
[]

[ICs]
  [./temp_aux_ic]
    variable = temp_aux
    value = 290
    type = ConstantIC
  [../]
[]

[UserObjects]
  [./activated_elem_uo]
    type = ActivateElementTemp
    execute_on = timestep_begin
    function_x= heat_source_x
    function_y= heat_source_y
    function_z= heat_source_z
    # activate_tol=1e-2
    active_subdomain_id = 2
    variable_activation = false
    expand_boundary_name = 'moving_boundary'
  [../]
[]

# [Postprocessors]
#   [./temperature]
#     type = ElementAverageValue
#     variable = temp
#     block = '2'
#   [../]
# []

[Outputs]
  file_base = './output/EA_thermal_out'
  [./exodus]
    type = Exodus
  [../]
  [./csv]
    type = CSV
    # file_base = 'tension_ramp_T${T_init}_out'
    execute_on = 'INITIAL TIMESTEP_END'
  [../]
[]
