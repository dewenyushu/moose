E = 1.7e11
D0 = 100.0

[GlobalParams]
  order = FIRST
  family = LAGRANGE
  displacements = 'disp_x disp_y disp_z'
  volumetric_locking_correction = true
[]

[Mesh]
  [block]
    type = GeneratedMeshGenerator
	dim = 3
	nx = 5
	ny = 5
	nz = 50
	xmax = 1
	ymax = 1
	zmax = 10
  []
[]

[Problem]
  type = ReferenceResidualProblem
  extra_tag_vectors = 'ref'
  reference_vector = 'ref'
  group_variables = 'disp_x disp_y disp_z'
  acceptable_multiplier = 10
  acceptable_iterations = 10
[]

[Variables]
  [disp_x]
  []
  [disp_y]
  []
  [disp_z]
  []
  [dose]
  []
[]

[Physics/SolidMechanics/QuasiStatic]
  [block]
    strain = FINITE
    volumetric_locking_correction = true
    #eigenstrain_names = 'swelling_strain'
    decomposition_method = EigenSolution
    generate_output = 'vonmises_stress stress_zz creep_strain_xx creep_strain_yy creep_strain_zz elastic_strain_zz'
    use_finite_deform_jacobian = true
    extra_vector_tags = 'ref'
	incremental = true
  []
[]

[Kernels]
  [value]
    type = MaterialPropertyValue
	prop_name = dpa
	variable = dose
  []
[]

#[AuxVariables]
#  [dose]
#  []
#[]

#[AuxKernels]
#  [damage_dose]
#    type = FunctionAux
#    function = damage_dose
#    variable = dose
#  []
#[]

[Functions]
  [damage_dose]
    type = ParsedFunction
	expression = 't*D0'		## damage does as a function of t, from 0, 0.1, ... 1.0
	symbol_names = 'D0'
	symbol_values = '${D0}'
  []
[]

[BCs]
  [outerPressure]
    type = Pressure
    boundary = front
    variable = disp_z
    factor = 1e6			## 1MPa pressure for creep strain
  []
  [no_x_block]
    type = DirichletBC
    variable = disp_x
    boundary = 'left'
    value = 0.0
  []
  [no_y_block]
    type = DirichletBC
    variable = disp_y
    boundary = 'bottom'
    value = 0.0
  []
  [no_z_block]
    type = DirichletBC
    variable = disp_z
    boundary = 'back'
    value = 0.0
  []
[]

[Materials]
  [elastic_tensor]
    type = ComputeIsotropicElasticityTensor
    youngs_modulus = ${E}
    poissons_ratio = 0.3
  []
#  [stress]
#    type = ComputeFiniteStrainElasticStress
#  []
#  [volumetric_swelling]
#    type = ParsedMaterial
#    property_name = volumetric_swelling
#    coupled_variables = 'dose'
#    expression = '0.03*dose/100'
#  []
#  [irradiation_swelling_strain]
#    type = ComputeVolumetricEigenstrain
#    volumetric_materials = volumetric_swelling
#    eigenstrain_name = swelling_strain
#    args = ''
#  []
#  [damage_dose]
#    type = ParsedMaterial
#    property_name = dpa
#    coupled_variables = 'dose'
#    expression = 't*dose'
#  []
  [damage_dose]
    type = GenericFunctionMaterial
    prop_names = dpa
    prop_values = damage_dose
  []
  [irradiation_creep_strain]
    type = VP5iaeaCreepMat
    dose = dose
    youngs_modulus = ${E}
  []
  [creep_plas]
    type = ComputeMultipleInelasticStress
    inelastic_models = 'irradiation_creep_strain'
  []
[]

[Preconditioning]
  [SMP]
    type = SMP
    full = true
  []
[]

[Executioner]
  automatic_scaling = true
  off_diagonals_in_auto_scaling = true
  solve_type = 'PJFNK'

  type = Transient
  petsc_options = '-ksp_snes_ew'
  petsc_options_iname = '-pc_type -pc_factor_mat_solver_package'
  petsc_options_value = 'lu       superlu_dist'
  line_search = basic

  l_max_its = 4
  nl_max_its = 20
  nl_rel_tol = 1e-6
  nl_abs_tol = 5e-10

  start_time = 0.0
  dt = 0.1
  end_time = 1.0

  dtmax = 1
  dtmin = 0.005

  [Predictor]
    type = SimplePredictor
    scale = 1.0
  []
[]

[Postprocessors]
  [pdata]
    type = PerfGraphData
    data_type = total
    section_name = "Root"
    execute_on = timestep_end
  []
  [block_disp_x]
    type = SideAverageValue
    variable = disp_x
    boundary = front
  []
  [block_disp_y]
    type = SideAverageValue
    variable = disp_y
    boundary = front
  []
  [block_disp_z]
    type = SideAverageValue
    variable = disp_z
    boundary = front
  []
[]

[Outputs]
  exodus = true
  perf_graph = true
  [block_displace]
    type = CSV
    file_base = block_displace
    execute_on = timestep_end
    sort_columns = True
    show = 'block_disp_x block_disp_y block_disp_z'
  []
  [runtime]
    type = CSV
    file_base = runtime
    show = pdata
    execute_on = timestep_end
  []
[]
