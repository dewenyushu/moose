[GlobalParams]
  displacements = 'disp_x disp_y disp_z'
[]

# specimen length = 10mm
# maximum displacement 0.5mm (5% deformation)
# total time 100s, timestep size 1s
# displacement per step, 5e-3 mm
[Mesh]
  [fmg]
    type = FileMeshGenerator
    file = 16mm_gauge_gyroid_3mm_cell_600um_wall_1mm_shell_trunc.e
  []
[]

[Modules/TensorMechanics/Master]
  [all]
    strain = FINITE
    incremental = true
    add_variables = true
    generate_output = 'stress_xx stress_yy stress_xy stress_yz stress_zz strain_xx strain_yy '
                      'strain_xy strain_xz strain_yz strain_zz elastic_strain_yy vonmises_stress '
                      'plastic_strain_yy'
    use_automatic_differentiation = false
  []
[]

[BCs]
  # fix the bottom
  [bottom_x]
    type = DirichletBC
    variable = disp_x
    boundary = traction_bottom
    value = 0
  []
  [bottom_y]
    type = DirichletBC
    variable = disp_y
    boundary = traction_bottom
    value = 0
  []
  [bottom_z]
    type = DirichletBC
    variable = disp_z
    boundary = traction_bottom
    value = 0
  []
  # [top_x] # free in the x and z directions
  #   type = DirichletBC
  #   variable = disp_x
  #   boundary = traction_top
  #   value = 0
  # []
  [top_y]
    type = FunctionDirichletBC
    variable = disp_y
    boundary = traction_top
    function = 5e-3*t # tension & displacement control
  []
  # [top_z] # free in the x and z directions
  #   type = DirichletBC
  #   variable = disp_z
  #   boundary = traction_top
  #   value = 0
  # []
[]

[Materials]
  [elasticity]
    type = ComputeIsotropicElasticityTensor
    youngs_modulus = 99e3 # MPa
    poissons_ratio = 0.35
  []
  # power law hardening
  [radial_return_stress]
    type = ComputeMultipleInelasticStress
    inelastic_models = 'power_law_hardening'
  []
  [power_law_hardening]
    # material parameters for SS316 or SS304, need to double check
    type = IsotropicPowerLawHardeningStressUpdate
    strength_coefficient = 847 #K
    strain_hardening_exponent = 0.06 #n
    relative_tolerance = 1e-6
    absolute_tolerance = 1e-6
  []
[]

# consider all off-diagonal Jacobians for preconditioning
[Preconditioning]
  [SMP]
    type = SMP
    full = true
  []
[]

[Executioner]
  type = Transient

  solve_type = 'NEWTON'

  line_search = none
  petsc_options_iname = '-pc_type -pc_factor_mat_solver_package -snes_atol -snes_rtol '
  petsc_options_value = 'lu mumps 0 1E-3'

  [TimeStepper]
    type = IterationAdaptiveDT
    dt = 1
    optimal_iterations = 8
    iteration_window = 3
  []
  dtmax = 5
  dtmin = 0.001
  end_time = 100

  nl_max_its = 30
  nl_rel_tol = 1e-6
  nl_abs_tol = 1e-6

  [Predictor]
    type = SimplePredictor
    scale = 1
  []
[]

[Outputs]
  file_base = './output/input_disp_control_plastic_out'
  exodus = true
  [check_point]
    type = Checkpoint
    interval = 5
  []
[]
