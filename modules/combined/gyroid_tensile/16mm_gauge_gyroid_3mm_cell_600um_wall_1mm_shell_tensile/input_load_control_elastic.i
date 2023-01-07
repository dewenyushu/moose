[GlobalParams]
  displacements = 'disp_x disp_y disp_z'
[]

[Mesh]
  [./fmg]
    type = FileMeshGenerator
    file = 16mm_gauge_gyroid_3mm_cell_600um_wall_1mm_shell_tensile.e
  []
[]

[Modules/TensorMechanics/Master]
  [all]
    add_variables = true
    generate_output = 'strain_xx strain_yy strain_zz stress_xx stress_yy stress_zz'
  []
[]

#
# Added boundary/loading conditions
# https://mooseframework.inl.gov/modules/tensor_mechanics/tutorials/introduction/step02.html
#
[BCs]
  [bottom_x]
    type = DirichletBC
    variable = disp_x
    boundary = traction_bottom
    # boundary = bottom
    value = 0
  []
  [bottom_y]
    type = DirichletBC
    variable = disp_y
    boundary = traction_bottom
    # boundary = bottom
    value = 0
  []
  [bottom_z]
    type = DirichletBC
    variable = disp_z
    boundary = traction_bottom
    # boundary = bottom
    value = 0
  []
  [Pressure]
    [top]
      boundary = traction_top
      # boundary =top
      function = -1e7*t # -> should we do displacement control instead?
    []
  []
[]

[Materials]
  [elasticity]
    type = ComputeIsotropicElasticityTensor
    youngs_modulus = 1e9
    poissons_ratio = 0.3
  []
  [stress]
    type = ComputeLinearElasticStress
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
  # scheme = crank-nicolson
  line_search = none
  petsc_options_iname = '-pc_type -pc_factor_mat_solver_package -snes_atol -snes_rtol '
  petsc_options_value = 'lu mumps 0 1E-3'

  nl_abs_tol = 1e-8
  nl_rel_tol = 1e-5

  end_time = 5
  dt = 0.5

  [Predictor]
    type = SimplePredictor
    scale = 1
  []
[]

[Outputs]
  file_base = './output/input_load_control_elastic_out'
  exodus = true
  [check_point]
    type = Checkpoint
    interval = 20
  []
[]
