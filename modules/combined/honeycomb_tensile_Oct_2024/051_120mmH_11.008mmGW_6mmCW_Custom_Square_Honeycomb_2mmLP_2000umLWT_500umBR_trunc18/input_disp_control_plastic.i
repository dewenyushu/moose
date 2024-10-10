[GlobalParams]
  displacements = 'disp_x disp_y disp_z'
  volumetric_locking_correction = true
[]

# specimen length = 12mm
# maximum displacement 0.6mm (5% deformation)
# total time 100s, timestep size 1s
# displacement per step, 6e-3 mm
[Mesh]
  [fmg]
    type = FileMeshGenerator
    file = convert_in.e
  []
  [cnode]
    type = ExtraNodesetGenerator
    nodes = '256501 257274'
    new_boundary = 100
    input = fmg
  []
  patch_update_strategy = auto
[]

[Modules]

  [TensorMechanics]

    [Master]
      [all]
        strain = FINITE
        incremental = true
        add_variables = true
        generate_output = 'stress_zz strain_zz vonmises_stress'
        use_automatic_differentiation = false
        # saves residuals
        save_in = 'saved_x saved_y saved_z'
      []
    []
  []
[]

[AuxVariables]
  [saved_x]
  []
  [saved_y]
  []
  [saved_z]
  []
[]

[BCs]
  [node_x]
    type = DirichletBC
    variable = disp_x
    boundary = 100
    value = 0
  []
  [node_y]
    type = DirichletBC
    variable = disp_y
    boundary = 100
    value = 0
  []
  [bottom_z]
    type = FunctionDirichletBC
    variable = disp_z
    boundary = traction_bottom
    function = -6.0e-3*t # tension & displacement control
  []
  [top_z]
    type = FunctionDirichletBC
    variable = disp_z
    boundary = traction_top
    function = 6.0e-3*t # tension & displacement control
  []
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
    # relative_tolerance = 1e-6
    # absolute_tolerance = 1e-6
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

  petsc_options_iname = '-pc_type -pc_hypre_type'
  petsc_options_value = 'hypre boomeramg'

  # [TimeStepper]
  #   type = IterationAdaptiveDT
  #   dt = 0.1
  #   optimal_iterations = 8
  #   iteration_window = 3
  # []
  # dtmax = 0.2
  # dtmin = 0.001
  # end_time = 100
  # # end_time = 7

  dt = 0.1
  dtmin = 0.001
  end_time = 20

  nl_max_its = 50
  nl_rel_tol = 1e-5
  nl_abs_tol = 1e-5
  l_abs_tol = 1e-8

  [Predictor]
    type = SimplePredictor
    scale = 1.0
  []
[]

[Outputs]
  exodus = true
  csv = true

  [checkpoint]
    type = Checkpoint
    interval = 5
    num_files = 2
  []
[]

[Postprocessors]
  # sum of nodal reaction force
  [force_z]
    type = NodalSum
    variable = saved_z
    boundary = traction_top
  []
  # use the sideset reaction kernel
  [reaction_z]
    type = SidesetReaction
    direction = '0 0 1'
    stress_tensor = stress
    boundary = traction_top
  []
  # use integration along the side
  [integral_stress_zz]
    type = SideIntegralVariablePostprocessor
    variable = stress_zz
    boundary = traction_top
  []
  [max_von_mises]
    type = ElementExtremeValue
    variable = vonmises_stress
    value_type = max
  []
  [max_stress_zz]
    type = ElementExtremeValue
    variable = stress_zz
    value_type = max
  []
  [integral_sigma_zz]
    type = ElementIntegralVariablePostprocessor
    variable = stress_zz
  []
  [integral_epsilon_zz]
    type = ElementIntegralVariablePostprocessor
    variable = strain_zz
  []
  [volume]
    type = VolumePostprocessor
    execute_on = 'INITIAL TIMESTEP_END'
  []
[]

[VectorPostprocessors]
  [hist_vm_stress]
    type = VariableValueVolumeHistogram
    variable = 'vonmises_stress'
    min_value = 0
    max_value = 1e3
    bin_number = 100
  []
[]
