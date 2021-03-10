[GlobalParams]
  displacements = 'disp_x disp_y disp_z'
[]

[Problem]
  kernel_coverage_check = false
[]

[Mesh]
  [./mesh]
    type = FileMeshGenerator
    file = hollow.e
  [../]
  # Substrate subdomain
  [./add_set1]
    type = SubdomainBoundingBoxGenerator
    input = mesh
    block_id = 3
    bottom_left = '-50 -50 0'
    top_right = '50 50 5'
  [../]
  # Product subdomain (entire product)
  # Parameters:
  #       bottom_left, top_right bounding box coordinates based on the geometry from your mesh
  [./add_set2]
    type = SubdomainBoundingBoxGenerator
    input = add_set1
    block_id = 1
    bottom_left = '50 50 5'
    top_right = '50 50 20'
  [../]
  # Active product subdomain
  # This is the subdomain where:
  #      more alements are added to during simulation
  #      thermal-mechanical properties are defined and solved
  # Note: Initially, this subdomain contains several elements from the product subdomain (bllock 1)
  #       This is not ideal, but is necessary because MOOSE does not allow an empty subdomain
  # Parameters:
  #       bottom_left, top_right bounding box coordinates based on the geometry from your mesh
  [./add_set3]
    type = SubdomainBoundingBoxGenerator
    input = add_set2
    block_id = 2
    bottom_left = '10 -1 5'
    top_right = '8 1 6'
  [../]
  # The interface that gets updated as more material is desposited
  [./moving_boundary]
    type = SideSetsAroundSubdomainGenerator
    input = add_set3
    block = 2
    new_boundary = 'moving_boundary'
  [../]
  displacements='disp_x disp_y disp_z'
[]

# Only define variables on the substrate + active product subdomains
[Variables]
  [./disp_x]
    block = '2 3'
  [../]
  [./disp_y]
    block = '2 3'
  [../]
  [./disp_z]
    block = '2 3'
  [../]
  [./temp]
    block = '2 3'
  [../]
[]

[ICs]
  # The substrate has environment temperature as initial temperature
  # Parameters:
  #       value = environment/ambient temperature [K]
  [./temp_substrate]
    type = ConstantIC
    variable = temp
    value = 300
    block = '3'
  [../]
  # The product has melt temperature as initial temperature
  # Because we activate the deposited material when it begins to solidify
  # Parameters:
  #       value = melt temperature [K]
  [./temp_product]
    type = ConstantIC
    variable = temp
    value = 800
    block = '2'
  [../]
[]

# Additional variables that will be output in the exodus file
[AuxVariables]
  [./von_mises]
    order = CONSTANT
    family = MONOMIAL
    block = '2 3'
  [../]
[]

# Mechanical kernel
[Modules/TensorMechanics/Master]
    strain = FINITE
    incremental = true
    add_variables = true
    generate_output = 'stress_xx stress_yy stress_zz stress_xy stress_yz stress_xz strain_yy strain_xx strain_zz strain_xy strain_xz strain_yz'
    use_automatic_differentiation = true
  [product]
    block = '2'
    eigenstrain_names = 'thermal_eigenstrain_product'
    use_automatic_differentiation = true
  []
  [substrate]
    block = '3'
    eigenstrain_names = 'thermal_eigenstrain_substrate'
    use_automatic_differentiation = true
  []
[]

# Thermal kernel
[Kernels]
  [./time]
    type = ADHeatConductionTimeDerivative
    variable = temp
    block = '2 3'
  [../]
  [./heat_conduct]
    type = ADHeatConduction
    variable = temp
    use_displaced_mesh = true
    thermal_conductivity = thermal_conductivity
    block = '2 3'
  [../]
  # The heat source mimics the moving laser during additive manufacturing
  # It is implemented as a path dependent heat source, 'volumetric_heat', defined in the Materials block
  [./heat_source]
    type = ADMatHeatSource
    material_property = volumetric_heat
    variable = temp
    scalar = 1
    use_displaced_mesh = true
    block = '2 3'
  [../]
[]

# AuxKernels corresponding to the aux variables
[AuxKernels]
  [./von_mises_kernel]
    type = ADRankTwoScalarAux
    variable = von_mises
    rank_two_tensor = stress
    execute_on = timestep_end
    scalar_type = VonMisesStress
    block = '2 3'
  [../]
[]

[BCs]
  # Fix the bottom displacements
  [./bottom_fix_x]
    type = ADDirichletBC
    variable = disp_x
    boundary = 1
    value = 0.0
  [../]
  [./bottom_fix_y]
    type = ADDirichletBC
    variable = disp_y
    boundary = 1
    value = 0.0
  [../]
  [./bottom_fix_z]
    type = ADDirichletBC
    variable = disp_z
    boundary = 1
    value = 0.0
  [../]
  # Fix the bottom temperature
  [./temp_bottom_fix]
    type = ADDirichletBC
    variable = temp
    boundary = 1
    value = 300
  [../]
  # Convective BC for the moving interface
  # Parameters:
  #       Add convective BCs for the substrate if necessary
  [./convective]
    type = ADConvectiveHeatFluxBC # Convective flux, e.g. q'' = h*(Tw - Tf)
    variable = temp
    boundary = 'moving_boundary'
    heat_transfer_coefficient = 0.1
    T_infinity = 300
  [../]
[]

[Materials]
  # Temperature dependent Young's modulus
  # The values are linearly interpolated from y-data based on the temperature
  # Parameters:
  #       xy_data, first column - temperature, second column - Young's modulus
  [./youngs_modulus]
    type = ADPiecewiseLinearInterpolationMaterial
    xy_data = '295        10e+8
               1000       10e+6
               99900      10e+5'
    property = youngs_modulus
    variable = temp
    block = '2 3'
  [../]
  # Parameters:
  #       poissons_ratio - can be defined with temperature dependency too, just like youngs_modulus
  [./elasticity_tensor]
    type = ADComputeVariableIsotropicElasticityTensor
    youngs_modulus = youngs_modulus
    poissons_ratio = 0.3
    block = '2 3'
  [../]

  # The stress calculation here is for demonstration
  # Parameters:
  #         Based on your problem
  [./stress]
    type = ADComputeFiniteStrainElasticStress
    block = '2 3'
  [../]

  # Eigenstrain material properties consider the deformation caused by thermal expansion
  # Parameters:
  #         stress_free_temperature and thermal_expansion_coeff values
  [./thermal_expansion_strain_product]
    type = ADComputeThermalExpansionEigenstrain
    stress_free_temperature = 800
    thermal_expansion_coeff = 2e-8
    temperature = temp
    eigenstrain_name = thermal_eigenstrain_product
    block = '2'
  [../]
  [./thermal_expansion_strain_substrate]
    type = ADComputeThermalExpansionEigenstrain
    stress_free_temperature = 300
    thermal_expansion_coeff = 2e-8
    temperature = temp
    eigenstrain_name = thermal_eigenstrain_substrate
    block = '3'
  [../]
  # Heat conduction material properties
  # Parameters:
  #       specific_heat and thermal_conductivity values
  [./heat]
    type = ADHeatConductionMaterial
    specific_heat = 603
    thermal_conductivity = 10e-2
    block = '2 3'
  [../]
  # Density
  # Parameters:
  #       density value
  [./density]
    type = ADDensity
    density = 4.43e-6
    block = '2 3'
  [../]
  # Heat source that mimics the moving laser during additive manufacturing
  # Parameters:
  #       a - transverse ellipsoid axe
  #       b - longitudinal ellipsoid axe
  #       c - depth ellipsoid axe
  #       power - laser power, P
  #       efficiency - process efficiency, eta
  #       factor - scaling factor, f
  #       function_x - The x component heating spot travel path, function_x(t) = xt
  #       function_y - The y component heating spot travel path, function_y(t) = yt
  #       function_z - The z component heating spot travel path, function_z(t) = zt
  #       path_x, y, z are defined in the `Functions` block
  # The heat source point (x, y, z) at time t:
  #     Q(x, y, z, t) = 6.0 * sqrt(3.0) * P * eta * f/(a * b * c * pi^1.5) *
  #             exp(-(3.0 * (x - xt)^2.0 / a^2.0 +
  #                   3.0 * (y - yt)^2.0 / b^2.0 +
  #                   3.0 * (z - zt)^2.0 / c^2.0));
  [./volumetric_heat]
    type = FunctionPathEllipsoidHeatSource
    a = 5
    b = 5
    c = 3
    power = 1000
    efficiency = 0.5
    factor = 2
    function_x= path_x
    function_y= path_y
    function_z= path_z
  [../]
[]

# Functions that mimics the scanning path as a function of time
[Functions]
  # x and y follows a circular path
  # Parameters:
  #         c - determines how long it takes to finish scanning a circular shape (4s in this case)
  #         For complex geometry, you may consider using tabulated data from experiments, using type = PiecewiseLinear
  [./path_x]
    type = ParsedFunction
    value = 9.0*cos(2.0*pi/c*t)
    vars = 'c'
    vals = '4'
  [../]
  [./path_y]
    type = ParsedFunction
    value = 9.0*sin(2.0*pi/c*t)
    vars = 'c'
    vals = '4'
  [../]
  # z increases for each layer
  # Parameters:
  #      first column: time
  #      second column: z_coordinate (height)
  [./path_z]
    type = PiecewiseConstant
    xy_data = '0 5
               4 6
               8 7
               12 8'
  [../]
[]

# This is where element activation is included
# Parameters:
#      function_x, y, z - same as that defined above
#      active_subdomain_id - the subdomain that we add elements to
#      inactive_subdomain_id - the substrate subdomain id (note, the 'inactive' here means we do not add elements to the substrate)
#      expand_boundary_name - the boundary on which we update boundary info
#      activate_distance - we activate the the element only if the distance between an element and the laser center is below the activate_distance
[UserObjects]
  [./activated_elem_uo]
    type = ActivateElementsByPath
    execute_on = timestep_begin
    function_x= path_x
    function_y= path_y
    function_z= path_z
    active_subdomain_id = 2
    inactive_subdomain_id = 3
    expand_boundary_name = 'moving_boundary'
    activate_distance = 0.7
  [../]
[]

[Preconditioning]
  [./full]
    type = SMP
    full = true
  [../]
[]

[Executioner]
  type = Transient
  solve_type = PJFNK

  nl_rel_tol = 1e-6
  nl_abs_tol = 1e-6

  petsc_options_iname = '-ksp_type -pc_type -pc_factor_mat_solver_package'
  petsc_options_value = 'preonly lu       superlu_dist'

  l_max_its = 100

  # How long to simulate, time step size, and minimum timestep size
  end_time = 8
  dt = 0.05
  dtmin = 1e-4
[]

[Outputs]
  exodus = true
[]
