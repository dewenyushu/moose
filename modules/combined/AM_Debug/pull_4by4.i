T_melt = 1000
# T_room = 300

[GlobalParams]
  displacements = 'disp_x disp_y disp_z'
[]

[Problem]
  kernel_coverage_check = false
[]

[Mesh]
  [./gen]
    type = GeneratedMeshGenerator
    dim = 3
    xmin =0
    xmax =1
    ymin =0
    ymax =1
    zmin =0
    zmax =1
    nx=4
    ny=4
    nz=4
  [../]
  [./rename_bnd]
    type = RenameBoundaryGenerator
    input = gen
    # bottom = 0, front = 1, right = 2, back = 3, left = 4, top = 5
    old_boundary_id='0 1 2 3 4 5'
    new_boundary_name='bottom front right back left top'
  [../]
  displacements='disp_x disp_y disp_z'
[]

[Variables]
  [./disp_x]
  [../]
  [./disp_y]
  [../]
  [./disp_z]
  [../]
[]

[AuxVariables]
  [./temp]
    initial_condition = ${T_melt}
  [../]
[]

[Modules/TensorMechanics/Master]
  [./all]
    strain = FINITE
    incremental = true
    add_variables = true
    generate_output = 'stress_xx stress_yy stress_zz stress_xy stress_yz stress_xz strain_yy strain_xx strain_zz strain_xy strain_xz strain_yz'
    use_automatic_differentiation = true
  [../]
[]

[Functions]
  [./pull_top]
    # type = ParsedFunction
    # # value = 'if(t>0.5, -0.1, -0.2*t)'
    # value = 0.2*t

    type = PiecewiseLinear
    x = '0 1'
    y = '0 -0.4'
  [../]
[]

[BCs]
  [./ux_bottom_fix]
    type = ADDirichletBC
    variable = disp_x
    boundary = 'bottom'
    value = 0.0
  [../]
  [./uy_bottom_fix]
    type = ADDirichletBC
    variable = disp_y
    boundary = 'bottom'
    value = 0.0
  [../]
  [./uz_bottom_fix]
    type = ADDirichletBC
    variable = disp_z
    boundary = 'bottom'
    value = 0.0
  [../]

  # [./u_left_fix]
  #   type = ADDirichletBC
  #   variable = disp_x
  #   boundary = 'left'
  #   value = 0.0
  # [../]
  # [./u_back_fix]
  #   type = ADDirichletBC
  #   variable = disp_y
  #   boundary = 'back'
  #   value = 0.0
  # [../]

  [./uz_top_pull]
    type = ADFunctionDirichletBC
    variable = disp_z
    boundary = 'top'
    function = pull_top
  [../]
[]

[Materials]
  [./E]
    type = ADPiecewiseLinearInterpolationMaterial
    x = '294.994  1671.48  1721.77'
    y = '201.232e3 80.0821e3 6.16016e3' # [N/mm^3], i.e., MPa
    property = youngs_modulus
    variable = temp
    extrapolation = false
  [../]
  [./nu]
    type = ADPiecewiseLinearInterpolationMaterial
    x = '294.994 1669.62 1721.77'
    y = '0.246407   0.36961  0.513347'
    property = poissons_ratio
    variable = temp
    extrapolation = false
  [../]
  [./elasticity_tensor]
    type = ADComputeVariableIsotropicElasticityTensor
    youngs_modulus = youngs_modulus
    poissons_ratio = poissons_ratio
  [../]

  [./radial_return_stress]
    type = ADComputeMultipleInelasticStress
    inelastic_models = 'rate_temp_plas'
    # block = '2 3'
  [../]
  [./rate_temp_plas]
    type = ADRateTempDependentStressUpdate
    temperature = temp
    Y0 = 5.264e03 # [MPa]
    Rd1 = 8.565e-4 # [MPa]
    hxi = 1.670e-06 # [mm/(s*MPa)]
    Ex = '294.994  1671.48  1721.77'
    Ey = '201.232e3 80.0821e3 6.16016e3' # [N/mm^3], i.e., MPa
    nux = '294.994 1669.62 1721.77'
    nuy = '0.246407   0.36961  0.513347'
    absolute_tolerance = 1e-8

    internal_solve_full_iteration_history = true
    internal_solve_output_on = "on_error"
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

  petsc_options_iname = '-ksp_type -pc_type -pc_factor_mat_solver_package'
  petsc_options_value = 'preonly lu       superlu_dist'

  line_search = 'none'

  l_max_its = 10
  nl_max_its = 20
  nl_rel_tol = 1e-4
  nl_abs_tol = 1e-12
  # l_tol = 1e-5

  start_time = 0.0
  end_time = 1.0
  dt = 0.1
  dtmin = 1e-4
[]



[Outputs]
  file_base = 'pull_4by4/out'
  [./exodus]
    type = Exodus
  [../]
  # [checkpoint]
  #   type = Checkpoint
  #   num_files = 2
  # []
[]
