T_melt = 1000
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
    # file = ./input_mesh/solid_cylinder_17_substrate_100x100.e
    file = ./input_mesh/solid_cylinder_17_substrate_100x100_fine_v0.e
  [../]
  [./add_set1]
    type = SubdomainBoundingBoxGenerator
    input = mesh
    block_id = 3
    bottom_left = '-50 -50 0'
    top_right = '50 50 5'
  [../]
  [./add_set2]
    type = SubdomainBoundingBoxGenerator
    input = add_set1
    block_id = 1
    bottom_left = '-50 -50 5'
    top_right = '50 50 25'
  [../]
  [./add_set3]
    type = SubdomainBoundingBoxGenerator
    input = add_set2
    block_id = 2
    bottom_left = '-0.5 -0.5 5'
    top_right = '0.5 0.5 6'
  [../]
  [./moving_boundary]
    type = SideSetsAroundSubdomainGenerator
    input = add_set3
    block = 2
    new_boundary = 'moving_boundary'
  []
  [./middle]
    type = SideSetsAroundSubdomainGenerator
    input = moving_boundary
    block = 3
    new_boundary = 'middle'
    normal = '0 0 1'
  []
[]

[Variables]
  [./temp]
    block = '2 3'
  [../]
[]

[ICs]
  [./temp_substrate]
    type = ConstantIC
    variable = temp
    value = ${T_room}
    block = '3'
  [../]
  [./temp_product]
    type = ConstantIC
    variable = temp
    value = ${T_melt}
    block = '2'
  [../]
[]

[AuxVariables]
  [./density]
    order = CONSTANT
    family = MONOMIAL
    block = '2 3'
  [../]
  [./ad_prop1]
    order = CONSTANT
    family = MONOMIAL
    block = '2 3'
  [../]
[]

[Kernels]
  # [./null_kernel]
  #   type = NullKernel
  #   variable = temp
  #   block = '2 3'
  #   jacobian_fill = 1e-5
  # [../]
  [./time]
    type = ADHeatConductionTimeDerivative
    variable = temp
    block = '2 3'
  [../]
  [./heat_conduct]
    type = ADHeatConduction
    variable = temp
    use_displaced_mesh = true
    thermal_conductivity = 5e-10
    block = '2 3'
  [../]
  [./heatsource]
    type = ADMatHeatSource
    material_property = volumetric_heat
    variable = temp
    scalar = 1
    use_displaced_mesh = true
    block = '2 3'
  [../]
[]

[AuxKernels]
  [./density_output]
    type = ADMaterialRealAux
    variable = density
    property = density
    block = '2 3'
  [../]
  [./ad_prop1_output]
    type = ADMaterialRealAux
    variable = ad_prop1
    property = ad_test_mat_property
    block = '2 3'
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
  #
  #
  # [./convective]
  #   # type = ConvectiveFluxFunction # Convective flux, e.g. q'' = h*(Tw - Tf)
  #   type = ADConvectiveHeatFluxBC
  #   variable = temp
  #   boundary = 'moving_boundary'
  #   # coefficient = 2e-5
  #   heat_transfer_coefficient = 2e-4
  #   T_infinity = ${T_room}
  # [../]
  # [./convective_substrate]
  #   # type = ConvectiveFluxFunction # Convective flux, e.g. q'' = h*(Tw - Tf)
  #   type = ADConvectiveHeatFluxBC
  #   variable = temp
  #   boundary = 2
  #   # coefficient = 2e-5
  #   heat_transfer_coefficient = 2e-4
  #   T_infinity = ${T_room}
  # [../]
  # [./convective_middle]
  #   # type = ConvectiveFluxFunction # Convective flux, e.g. q'' = h*(Tw - Tf)
  #   type = ADConvectiveHeatFluxBC
  #   variable = temp
  #   boundary = 'middle'
  #   # coefficient = 2e-5
  #   heat_transfer_coefficient = 2e-4
  #   T_infinity = ${T_room}
  # [../]
[]

[Materials]
  [./thermal_expansion_strain_product]
    type = ADComputeThermalExpansionEigenstrain
    stress_free_temperature = ${T_melt}
    thermal_expansion_coeff = 0.5e-5
    temperature = temp
    eigenstrain_name = thermal_eigenstrain_product
    block = '2'
  [../]
  [./thermal_expansion_strain_substrate]
    type = ADComputeThermalExpansionEigenstrain
    stress_free_temperature = ${T_room}
    thermal_expansion_coeff = 0.5e-5
    temperature = temp
    eigenstrain_name = thermal_eigenstrain_substrate
    block = '3'
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
    block = '2 3'
  [../]
  [./heat]
    type = ADHeatConductionMaterial
    specific_heat_temperature_function = specific_heat
    thermal_conductivity_temperature_function = thermal_conductivity
    temp = temp
    block = '2 3'
  [../]
  [./ad_stateful]
    type = ADStatefulTest
    prop_names = ad_test_mat_property
    prop_values = 1.0
    block = '2 3'
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
  end_time = 0.5
  dt = 0.1
  dtmin = 1e-4
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
    inactive_subdomain_id = 3
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
