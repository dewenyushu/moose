T_room = 300
T_ambient = 300
T_melt = 1000

speed = 10e-3
power = 300e-3
r = 200e-3
dt = 1

[GlobalParams]
  displacements = 'disp_x disp_y disp_z'
[]

[Mesh]
  [mesh]
    type = GeneratedMeshGenerator
    dim = 3
    xmax = 1
    ymax = 5
    zmax = 1
    nx = 10
    ny = 50
    nz = 10
  []
  [add_set1]
    type = SubdomainBoundingBoxGenerator
    input = mesh
    block_id = 3
    bottom_left = '0 0 0'
    top_right = '1 5 0.5'
  []
  [add_set2]
    type = SubdomainBoundingBoxGenerator
    input = add_set1
    block_id = 1
    bottom_left = '0 0 0.5'
    top_right = '1 5 1'
  []
  [add_set3]
    type = SubdomainBoundingBoxGenerator
    input = add_set2
    block_id = 2
    bottom_left = '0.5 0.5 0.5'
    top_right = '0.6 0.6 0.6'
  []
  [moving_boundary]
    type = SideSetsAroundSubdomainGenerator
    input = add_set3
    block = 2
    new_boundary = 'moving_boundary'
  []
  [middle]
    type = SideSetsAroundSubdomainGenerator
    input = moving_boundary
    block = 3
    new_boundary = 'middle'
    normal = '0 0 1'
  []
  # displacements = 'disp_x disp_y disp_z'
[]

[Problem]
  kernel_coverage_check = false
  material_coverage_check = false
[]

[Variables]
  [disp_x]
    block = '2 3'
  []
  [disp_y]
    block = '2 3'
  []
  [disp_z]
    block = '2 3'
  []
[]

[AuxVariables]
  [temp_aux]
    order = FIRST
    family = LAGRANGE
  []
[]

[Modules/TensorMechanics/Master]
  # [./all]
  strain = FINITE
  incremental = true
  add_variables = true
  generate_output = 'stress_xx stress_yy stress_zz stress_xy stress_yz stress_xz strain_yy strain_xx '
                    'strain_zz strain_xy strain_xz strain_yz'
  use_automatic_differentiation = true
  # eigenstrain_names = 'thermal_eigenstrain'
  # block = '2 3'
  # [../]
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

[BCs]
  [ux_bottom_fix]
    type = ADDirichletBC
    variable = disp_x
    boundary = 1
    value = 0.0
  []
  [uy_bottom_fix]
    type = ADDirichletBC
    variable = disp_y
    boundary = 1
    value = 0.0
  []
  [uz_bottom_fix]
    type = ADDirichletBC
    variable = disp_z
    boundary = 1
    value = 0.0
  []
[]

[Materials]
  [E]
    type = ADPiecewiseLinearInterpolationMaterial
    x = '294.994  1671.48  1721.77'
    y = '201.232 80.0821 6.16016' # 10^9 Pa = 10^9 kg/m/s^2 = kg/mm/ms^2
    property = youngs_modulus
    variable = temp_aux
    extrapolation = false
    block = '2 3'
  []
  [nu]
    type = ADPiecewiseLinearInterpolationMaterial
    x = '294.994 1669.62 1721.77'
    y = '0.246407   0.36961  0.513347'
    property = poissons_ratio
    variable = temp_aux
    extrapolation = false
    block = '2 3'
  []
  [elasticity_tensor]
    type = ADComputeVariableIsotropicElasticityTensor
    youngs_modulus = youngs_modulus
    poissons_ratio = poissons_ratio
    block = '2 3'
  []

  [stress]
    type = ADComputeFiniteStrainElasticStress
    block = '2 3'
  []

  [thermal_expansion_strain_product]
    type = ADComputeThermalExpansionEigenstrain
    stress_free_temperature = ${T_room}
    # thermal_expansion_coeff = 1.72e-5
    thermal_expansion_coeff = 6.72e-6 #1.72e-5 /K
    temperature = temp_aux
    eigenstrain_name = thermal_eigenstrain_product
    block = '2'
  []
  [thermal_expansion_strain_substrate]
    type = ADComputeThermalExpansionEigenstrain
    stress_free_temperature = ${T_room}
    # thermal_expansion_coeff = 1.72e-5
    thermal_expansion_coeff = 6.72e-6 #1.72e-5 /K
    temperature = temp_aux
    eigenstrain_name = thermal_eigenstrain_substrate
    block = '3'
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
  solve_type = 'PJFNK'

  # automatic_scaling = true

  petsc_options_iname = '-ksp_type -pc_type -pc_factor_mat_solver_package'
  petsc_options_value = 'preonly lu       superlu_dist'

  line_search = 'none'

  l_max_its = 100
  nl_max_its = 20
  nl_rel_tol = 1e-5
  nl_abs_tol = 1e-6

  start_time = 0.0
  end_time = '${fparse 5/speed}'
  dt = ${dt}
  num_steps = 40
  dtmin = 1e-6
[]

[UserObjects]
  [activated_elem_uo]
    type = CoupledVarThresholdElementSubdomainModifier
    execute_on = 'TIMESTEP_BEGIN'
    coupled_var = temp_aux
    block = 1
    subdomain_id = 2
    criterion_type = ABOVE
    threshold = ${T_melt}
    moving_boundary_name = 'moving_boundary'
    apply_initial_conditions = false
  []
[]

# [Adaptivity]
#   steps = 1
#   marker = marker
#   initial_marker = marker
#   max_h_level = 2
#   [Indicators/indicator]
#     type = GradientJumpIndicator
#     variable = temp_aux
#     scale_by_flux_faces = true
#   []
#   [Markers/marker]
#     type = ErrorFractionMarker
#     indicator = indicator
#     coarsen = 0
#     refine = 0.7
#   []
# []

[Outputs]
  file_base = 'output/thermal_activation_amr/sub_out'
  [exodus]
    type = Exodus
    # execute_on = 'INITIAL TIMESTEP_END'
    interval = 1
  []
[]
