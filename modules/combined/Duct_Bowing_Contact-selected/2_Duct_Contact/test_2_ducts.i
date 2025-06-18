fbase = './outputs/'

[GlobalParams]
  displacements = 'disp_x disp_y disp_z'
[]

# units are cm - do not forget to convert to meter
duct_length = 400.00 		# [cm] including nozzle, minus handling socket at top not expanded
duct_outer_ftf = 13.38129 	# [cm]
duct_inner_ftf = 12.77718 	# [cm]
duct_gap = 0.60410			# [cm]


aclp_location = 300.00 	# [cm] not expanded
tlp_location = 400.00 	# [cm] not expanded
lp_outer_ftf = 13.93505 # [cm]
lp_length = 5.00 		# [cm]

# discretization
duct_n_ax = 320
lp_n_ax = 4

ns = 4
duct_intervals_perishperic =  '1 1' # '2 2' #

[Mesh]
  [dummy]
    type = PolygonConcentricCircleMeshGenerator
    num_sides = 6
    num_sectors_per_side = '${ns} ${ns} ${ns} ${ns} ${ns} ${ns}'
    background_intervals = 1
    background_block_ids = '100'
	polygon_size_style = 'apothem'
    polygon_size = ${fparse duct_outer_ftf /2/100 + duct_gap /2 /100}
	duct_sizes_style = apothem
    duct_sizes = '${fparse duct_inner_ftf /2 /100} ${fparse duct_outer_ftf /2 /100}'
    duct_intervals = ${duct_intervals_perishperic}
    duct_block_ids = '1500 100'
	preserve_volumes = on
	quad_center_elements = true
    outward_interface_boundary_names = 'dummy duct'
    interface_boundary_id_shift = 10
  []
  [duct1]
    type = PolygonConcentricCircleMeshGenerator
    num_sides = 6
    num_sectors_per_side = '${ns} ${ns} ${ns} ${ns} ${ns} ${ns}'
    background_intervals = 1
    background_block_ids = '100'
	polygon_size_style = 'apothem'
    polygon_size = ${fparse duct_outer_ftf /2/100 + duct_gap /2 /100}
	duct_sizes_style = apothem
    duct_sizes = '${fparse duct_inner_ftf /2 /100} ${fparse duct_outer_ftf /2 /100}'
    duct_intervals = ${duct_intervals_perishperic}
    duct_block_ids = '1 100'
	preserve_volumes = on
	quad_center_elements = true
    outward_interface_boundary_names = 'dummy duct'
    interface_boundary_id_shift = 10
  []
  [duct2]
    type = PolygonConcentricCircleMeshGenerator
    num_sides = 6
    num_sectors_per_side = '${ns} ${ns} ${ns} ${ns} ${ns} ${ns}'
    background_intervals = 1
    background_block_ids = '100'
	polygon_size_style = 'apothem'
    polygon_size = ${fparse duct_outer_ftf /2/100 + duct_gap /2 /100}
	duct_sizes_style = apothem
    duct_sizes = '${fparse duct_inner_ftf /2 /100} ${fparse duct_outer_ftf /2 /100}'
    duct_intervals = ${duct_intervals_perishperic}
    duct_block_ids = '2 100'
	preserve_volumes = on
	quad_center_elements = true
    outward_interface_boundary_names = 'dummy duct'
    interface_boundary_id_shift = 10
  []
  [duct_pattern]
    type = PatternedHexMeshGenerator
    inputs = 'dummy duct1 duct2'
    pattern =
			  '0 0;
              0 1 2;
               0 0'
    pattern_boundary = none
  []
  [duct_center_removal]
    type = BlockDeletionGenerator
	input = duct_pattern
    block = '100 1500'
  []
  [duct_extrude]
    type = AdvancedExtruderGenerator
	input = duct_center_removal
    direction = '0 0 1'
    heights = '${fparse duct_length/100}'
    num_layers = '${duct_n_ax}'
  []
  [duct_boundary]
    type = RenameBoundaryGenerator
	input = duct_extrude
    old_boundary = '14  15'
    new_boundary = 'duct_bottom duct_top'
  []

  [aclp1]
    type = PolygonConcentricCircleMeshGenerator
    num_sides = 6
    num_sectors_per_side = '${ns} ${ns} ${ns} ${ns} ${ns} ${ns}'
    background_intervals = 1
    background_block_ids = '100'
	polygon_size_style = 'apothem'
    polygon_size = ${fparse duct_outer_ftf /2/100 + duct_gap /2 /100}
	duct_sizes_style = apothem
    duct_sizes = '${fparse duct_outer_ftf /2 /100} ${fparse lp_outer_ftf /2 /100}'
    duct_intervals = ${duct_intervals_perishperic}
    duct_block_ids = '1 100'
	preserve_volumes = on
	quad_center_elements = true
    create_inward_interface_boundaries = True
    interface_boundary_id_shift = 20
  []
  [aclp2]
    type = PolygonConcentricCircleMeshGenerator
    num_sides = 6
    num_sectors_per_side = '${ns} ${ns} ${ns} ${ns} ${ns} ${ns}'
    background_intervals = 1
    background_block_ids = '100'
	polygon_size_style = 'apothem'
    polygon_size = ${fparse duct_outer_ftf /2/100 + duct_gap /2 /100}
	duct_sizes_style = apothem
    duct_sizes = '${fparse duct_outer_ftf /2 /100} ${fparse lp_outer_ftf /2 /100}'
    duct_intervals = ${duct_intervals_perishperic}
    duct_block_ids = '2 100'
	preserve_volumes = on
	quad_center_elements = true
    create_inward_interface_boundaries = True
    interface_boundary_id_shift = 20
  []
  [aclp_pattern]
    type = PatternedHexMeshGenerator
    inputs = 'dummy aclp1 aclp2'
    pattern =
			  '0 0;
              0 1 2;
               0 0'
	interface_boundary_id_shift_pattern =
			     '300 300;
                300 100 200;
                  300 300'
    pattern_boundary = none
  []
  [aclp_center_removal]
    type = BlockDeletionGenerator
	input = aclp_pattern
    block = '100 1500'
  []
  [aclp_translate]
    type = TransformGenerator
	input = aclp_center_removal
	transform = translate
    vector_value = '0 0 ${fparse aclp_location/100 - lp_length/2/100}'
  []
  [aclp_extrude]
    type = AdvancedExtruderGenerator
	input = aclp_translate
    direction = '0 0 1'
    heights = '${fparse lp_length/100}'
    num_layers = '${lp_n_ax}'
  []
  [aclp_boundary]
    type = RenameBoundaryGenerator
	input = aclp_extrude
    old_boundary = '122 222 '
				   '123 223 '
				   '224 225'
    new_boundary = 'aclp_inside aclp_inside '
				   'aclp_1 aclp_2 '
				   'aclp_bottom aclp_top'
  []

  [tlp1]
    type = PolygonConcentricCircleMeshGenerator
    num_sides = 6
    num_sectors_per_side = '${ns} ${ns} ${ns} ${ns} ${ns} ${ns}'
    background_intervals = 1
    background_block_ids = '100'
	polygon_size_style = 'apothem'
    polygon_size = ${fparse duct_outer_ftf /2/100 + duct_gap /2 /100}
	duct_sizes_style = apothem
    duct_sizes = '${fparse duct_outer_ftf /2 /100} ${fparse lp_outer_ftf /2 /100}'
    duct_intervals = ${duct_intervals_perishperic}
    duct_block_ids = '1 100'
	preserve_volumes = on
	quad_center_elements = true
    create_inward_interface_boundaries = True
    interface_boundary_id_shift = 30
  []
  [tlp2]
    type = PolygonConcentricCircleMeshGenerator
    num_sides = 6
    num_sectors_per_side = '${ns} ${ns} ${ns} ${ns} ${ns} ${ns}'
    background_intervals = 1
    background_block_ids = '100'
	polygon_size_style = 'apothem'
    polygon_size = ${fparse duct_outer_ftf /2/100 + duct_gap /2 /100}
	duct_sizes_style = apothem
    duct_sizes = '${fparse duct_outer_ftf /2 /100} ${fparse lp_outer_ftf /2 /100}'
    duct_intervals = ${duct_intervals_perishperic}
    duct_block_ids = '2 100'
	preserve_volumes = on
	quad_center_elements = true
    create_inward_interface_boundaries = True
    interface_boundary_id_shift = 30
  []
  [tlp_pattern]
    type = PatternedHexMeshGenerator
    inputs = 'dummy tlp1 tlp2'
    pattern =
			  '0 0;
              0 1 2;
               0 0'
	interface_boundary_id_shift_pattern =
			     '300 300;
                300 100 200;
                  300 300'
    pattern_boundary = none
  []
  [tlp_center_removal]
    type = BlockDeletionGenerator
	input = tlp_pattern
    block = '100 1500'
  []
  [tlp_translate]
    type = TransformGenerator
	input = tlp_center_removal
	transform = translate
    vector_value = '0 0 ${fparse tlp_location/100 - lp_length/2/100}'
  []
  [tlp_extrude]
    type = AdvancedExtruderGenerator
	input = tlp_translate
    direction = '0 0 1'
    heights = '${fparse lp_length/100/2}'
    num_layers = '${fparse lp_n_ax/2}'
  []
  [tlp_boundary]
    type = RenameBoundaryGenerator
	input = tlp_extrude
    old_boundary = '132 232 '
	               '133 233 '
				   '234 235'
    new_boundary = 'tlp_inside tlp_inside '
				   'tlp_1 tlp_2 '
				   'tlp_bottom tlp_top'
  []

  ## Stitching with Sidesets
  [tlp_stitching]
    type = StitchedMeshGenerator
	inputs = 'duct_boundary tlp_boundary'
	stitch_boundaries_pairs = 'duct tlp_inside'
	clear_stitched_boundary_ids = False
	prevent_boundary_ids_overlap = False
  []
  [aclp_stitching]
    type = StitchedMeshGenerator
	inputs = 'tlp_stitching aclp_boundary'
	stitch_boundaries_pairs = 'duct aclp_inside'
	clear_stitched_boundary_ids = False
	prevent_boundary_ids_overlap = False
  []
  patch_update_strategy = auto
  patch_size = 20
  partitioner = centroid
  centroid_partitioner_direction = z
[]

[Problem]
   type = ReferenceResidualProblem
   extra_tag_vectors = 'ref'
   reference_vector = 'ref'
   group_variables = 'disp_x disp_y disp_z'
   acceptable_multiplier = 10
[]

[Variables]
  [disp_x]
  []
  [disp_y]
  []
  [disp_z]
  []
[]

[AuxVariables]
  [temp]
    initial_condition = 400.0
  []
  [aclp1_2_contact]
  []
  [tlp1_2_contact]
  []
[]

[Physics/SolidMechanics/QuasiStatic]
  [duct_1]
    strain = FINITE
	volumetric_locking_correction = true
    add_variables = true
    eigenstrain_names = 'thermal_expansion'
	decomposition_method = EigenSolution
	generate_output = 'vonmises_stress'
	temperature = temp
	use_finite_deform_jacobian = true
	extra_vector_tags = 'ref'
	block = '1'
  []
  [duct_2]
    strain = FINITE
	volumetric_locking_correction = true
    add_variables = true
	decomposition_method = EigenSolution
	generate_output = 'vonmises_stress'
	temperature = temp
	use_finite_deform_jacobian = true
	extra_vector_tags = 'ref'
	block = '2'
  []
[]

[AuxKernels]
  [tfunc]
    type = FunctionAux
	variable = temp
	function = temp_func
	block = '1'
  []
  [aclp1_2]
    type = PenetrationAux
    variable = aclp1_2_contact
    boundary = 'aclp_1'
    paired_boundary = 'aclp_2'
    quantity = normal_force_magnitude
  []
  [tlp1_2]
    type = PenetrationAux
    variable = tlp1_2_contact
    boundary = 'tlp_1'
    paired_boundary = 'tlp_2'
    quantity = normal_force_magnitude
  []
[]

[Functions]
# The duct temperatures are defined at the corners and linearly vary in the axial direction
# and along the face of the duct.

  [temp_func]
    type = ParsedFunction
	#At center of wall, y=+-0.075m
	#T varies across the cross-section from 500C to 550C, ramps up to that from 400C at z=1.5m to 2.5m
	expression = 'if(t>1.0,400+if(z>2.5,1*(125-25/.075513*(y)*1),if(z>1.5,1*(z-1.5)/1.0*(125-25/.075513*(y)*1),0)), 400+if(z>2.5,t*(125-25/.075513*(y)*t),if(z>1.5,t*(z-1.5)/1.0*(125-25/.075513*(y)*t),0)))'
  []
[]

[BCs]
  [no_x]
    type = DirichletBC
    variable = disp_x
    boundary = duct_bottom
    value = 0.0
  []
  [no_y]
    type = DirichletBC
    variable = disp_y
    boundary = duct_bottom
    value = 0.0
  []
  [no_z]
    type = DirichletBC
    variable = disp_z
    boundary = duct_bottom
    value = 0.0
  []
[]

[Materials]
  [elasticity_tensor]
    type = ComputeIsotropicElasticityTensor
    youngs_modulus = 1.7e11
    poissons_ratio = 0.3
	block = '1 2'
  []
  [stress]
    type = ComputeFiniteStrainElasticStress
	block = '1 2'
  []
  [thermal_expansion_strain]
    type = ComputeThermalExpansionEigenstrain
    stress_free_temperature = 400.0
    thermal_expansion_coeff = 18.0e-6
    temperature = temp
    eigenstrain_name = thermal_expansion
	block = '1'
  []
[]

[Contact]
  [aclp]
    primary =   'aclp_2'
    secondary = 'aclp_1'
    model = frictionless
    formulation = penalty
    normalize_penalty = true
    penalty = 1e10
    tangential_tolerance = 0.0
    normal_smoothing_distance = 0.0
  []
  [tlp]
    primary =   'tlp_2'
    secondary = 'tlp_1'
    model = frictionless
    formulation = penalty
    normalize_penalty = true
    penalty = 5e10
    tangential_tolerance = 0.0
    normal_smoothing_distance = 0.0
  []
[]

[Preconditioning]
  active = 'smp1'
  [smp1]
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
  end_time = 1.0
  dt = 0.1

  dtmax = 0.5
  dtmin = 0.01

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
  [aclp_force]
    type = NodalSum
    variable = aclp1_2_contact
    boundary = 'aclp_1'
  []
  [tlp_force]
    type = NodalSum
    variable = tlp1_2_contact
    boundary = 'tlp_1'
  []
[]

[VectorPostprocessors]
  [duct_1_average]
    type = AverageSectionValueSampler
    axis_direction = '0 0 1'
    block = '1'
    variables = 'disp_x disp_y disp_z'
    reference_point = '0.0 0.0 0.0'
	require_equal_node_counts = false
  []
  [duct_2_average]
    type = AverageSectionValueSampler
    axis_direction = '0 0 1'
    block = '2'
    variables = 'disp_x disp_y disp_z'
    reference_point = '0.0 0.136833 0.0'
	require_equal_node_counts = false
  []
[]

[Outputs]
  exodus = true
  perf_graph = true
  [duct_displace]
    type = CSV
    file_base = ${fbase}average_section_disp
    execute_on = timestep_end
    show = 'duct_1_average duct_2_average'
  []
  [force_plots]
    type = CSV
    file_base = ${fbase}contact_force
    execute_on = timestep_end
    show = 'aclp_force tlp_force'
  []
  [runtime]
    type = CSV
    file_base = ${fbase}runtime
    show = pdata
    execute_on = timestep_end
  []
  [cpt]
    type = Checkpoint
    num_files = 2
  []
[]
