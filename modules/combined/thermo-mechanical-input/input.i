[Mesh]
  [file_mesh]
    type = FileMeshGenerator
    file = 'mesh.e'
  []
[]

[GlobalParams]
  displacements = 'disp_x disp_y'
[]

[Variables]
  # Displacement variables
  [disp_x]
  []
  [disp_y]
  []
  # Temperature variable
  [temp]
    initial_condition = 300
  []
[]

# Adding the mechanical problem
[Physics]
  [SolidMechanics]
    [QuasiStatic]
      [all]
        strain = FINITE
        add_variables = true
        use_finite_deform_jacobian = true
        volumetric_locking_correction = true
        use_automatic_differentiation = true
        generate_output = 'stress_xx stress_yy stress_xy strain_xx strain_yy strain_xy'
      []
    []
  []
[]

# Adding the thermal problem
# Properties are set as constant for simplicity, except for the thermal_conductivity, which varies based on the specific material subdomains in the `Materials` block.
[Kernels]
  [heat_time]
    type = ADHeatConductionTimeDerivative
    specific_heat = 1e-3
    density_name = 1e1
    variable = temp
  []
  [heat_conduction]
    type = ADHeatConduction
    variable = temp
    thermal_conductivity = thermal_conductivity
  []
  [heat_source]
    type = HeatSource
    value = 10
    variable = temp
  []
[]

[Materials]
  # Definition of the mechanical model for each component
  # Subdomain 1, W
  [elasticity_tensor_w]
    type = ADComputeIsotropicElasticityTensor
    youngs_modulus = 1e1
    poissons_ratio = 0.3
    block = 'W'
  []
  [stress_w]
    type = ADComputeFiniteStrainElasticStress
    block = 'W'
  []
  # Subdomain 2, Steel
  [elasticity_tensor_steel]
    type = ADComputeIsotropicElasticityTensor
    youngs_modulus = 1e2
    poissons_ratio = 0.3
    block = 'Steel'
  []
  [stress_steel]
    type = ADComputeFiniteStrainElasticStress
    block = 'Steel'
  []
  # Subdomain 3, SIC
  [elasticity_tensor_sic]
    type = ADComputeIsotropicElasticityTensor
    youngs_modulus = 1e3
    poissons_ratio = 0.3
    block = 'SIC'
  []
  [stress_sic]
    type = ADComputeFiniteStrainElasticStress
    block = 'SIC'
  []
  # Definition of the thermal conductivity
  [thermal_cond_w]
    type = ADConstantMaterial
    property_name = thermal_conductivity
    value = 1e-1
    block = 'W'
  []
  [thermal_cond_steel]
    type = ADConstantMaterial
    property_name = thermal_conductivity
    value = 1e-2
    block = 'Steel'
  []
  [thermal_cond_sic]
    type = ADConstantMaterial
    property_name = thermal_conductivity
    value = 1e-3
    block = 'SIC'
  []
[]

[BCs]
  # mechanical boundary conditions
  [move_top]
    type = FunctionDirichletBC
    variable = disp_y
    boundary = top
    function = '1e-4*t'
  []
  [move_bottom]
    type = DirichletBC
    variable = disp_y
    boundary = bottom
    value = 0.0
  []
  # thermal boundary conditions
  [flow_channel]
    type = NeumannBC
    variable = temp
    boundary = flow_channel
    value = 0
  []
  [he]
    type = NeumannBC
    variable = temp
    boundary = he
    value = 0
  []
  [left]
    type = DirichletBC
    variable = temp
    boundary = left
    value = 300
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
  file_base = 'thermomechanical_out'
  exodus = true
[]
