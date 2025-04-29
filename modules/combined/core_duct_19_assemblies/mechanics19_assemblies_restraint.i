# Structural Mechanics mesh for 19 ABR assemblies
# Dimensions are pre-expanded for initial thermal expansion heat up from room temp (20 degC) to 
# hot nominal condition (380 degC)
# 
# Units are [cm] - do not forget to convert to meter when using fparse below

hot_nominal_temp = 380      # [deg C]

duct_length = 480.6143 		# [cm] including nozzle, minus handling socket at top, assumes fixed support at the bottom
duct_outer_ftf = 15.8118 	# [cm]
duct_inner_ftf = 15.0187 	# [cm]
duct_gap = 0.4348 			# [cm]

aclp_location = 253.0895 	# [cm]
tlp_location = 442.2674 	# [cm]
lp_outer_ftf = 16.2043 		# [cm]
lp_length = 10.2258			# [cm]

active_fuel_bot = 161.0569  # [cm]
active_fuel_top = 242.8636  # [cm]

# discretization
duct_n_ax = 188
lp_n_ax = 4

ns = 4
duct_intervals_perishperic =  '1 1' # '2 2' #

#theta2 = 1.570796326794897
theta3 = 2.617993877991494
theta4 = 3.665191429188092
#theta5 = 4.712388980384690
theta6 = 5.759586531581287
theta7 = 0.523598775598299

[GlobalParams]
  displacements = 'disp_x disp_y disp_z'
[]

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
  [duct3]
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
    duct_block_ids = '3 100'
	preserve_volumes = on
	quad_center_elements = true
    outward_interface_boundary_names = 'dummy duct'
    interface_boundary_id_shift = 10
  []
  [duct4]
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
    duct_block_ids = '4 100'
	preserve_volumes = on
	quad_center_elements = true
    outward_interface_boundary_names = 'dummy duct'
    interface_boundary_id_shift = 10
  []
  [duct5]
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
    duct_block_ids = '5 100'
	preserve_volumes = on
	quad_center_elements = true
    outward_interface_boundary_names = 'dummy duct'
    interface_boundary_id_shift = 10
  []
  [duct6]
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
    duct_block_ids = '6 100'
	preserve_volumes = on
	quad_center_elements = true
    outward_interface_boundary_names = 'dummy duct'
    interface_boundary_id_shift = 10
  []
  [duct7]
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
    duct_block_ids = '7 100'
	preserve_volumes = on
	quad_center_elements = true
    outward_interface_boundary_names = 'dummy duct'
    interface_boundary_id_shift = 10
  []
  [duct8]
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
    duct_block_ids = '8 100'
	preserve_volumes = on
	quad_center_elements = true
    outward_interface_boundary_names = 'dummy duct'
    interface_boundary_id_shift = 10
  []
  [duct9]
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
    duct_block_ids = '9 100'
	preserve_volumes = on
	quad_center_elements = true
    outward_interface_boundary_names = 'dummy duct'
    interface_boundary_id_shift = 10
  []
  [duct10]
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
    duct_block_ids = '10 100'
	preserve_volumes = on
	quad_center_elements = true
    outward_interface_boundary_names = 'dummy duct'
    interface_boundary_id_shift = 10
  []
  [duct11]
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
    duct_block_ids = '11 100'
	preserve_volumes = on
	quad_center_elements = true
    outward_interface_boundary_names = 'dummy duct'
    interface_boundary_id_shift = 10
  []
  [duct12]
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
    duct_block_ids = '12 100'
	preserve_volumes = on
	quad_center_elements = true
    outward_interface_boundary_names = 'dummy duct'
    interface_boundary_id_shift = 10
  []
  [duct13]
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
    duct_block_ids = '13 100'
	preserve_volumes = on
	quad_center_elements = true
    outward_interface_boundary_names = 'dummy duct'
    interface_boundary_id_shift = 10
  []
  [duct14]
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
    duct_block_ids = '14 100'
	preserve_volumes = on
	quad_center_elements = true
    outward_interface_boundary_names = 'dummy duct'
    interface_boundary_id_shift = 10
  []
  [duct15]
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
    duct_block_ids = '15 100'
	preserve_volumes = on
	quad_center_elements = true
    outward_interface_boundary_names = 'dummy duct'
    interface_boundary_id_shift = 10
  []
  [duct16]
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
    duct_block_ids = '16 100'
	preserve_volumes = on
	quad_center_elements = true
    outward_interface_boundary_names = 'dummy duct'
    interface_boundary_id_shift = 10
  []
  [duct17]
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
    duct_block_ids = '17 100'
	preserve_volumes = on
	quad_center_elements = true
    outward_interface_boundary_names = 'dummy duct'
    interface_boundary_id_shift = 10
  []
  [duct18]
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
    duct_block_ids = '18 100'
	preserve_volumes = on
	quad_center_elements = true
    outward_interface_boundary_names = 'dummy duct'
    interface_boundary_id_shift = 10
  []
  [duct19]
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
    duct_block_ids = '19 100'
	preserve_volumes = on
	quad_center_elements = true
    outward_interface_boundary_names = 'dummy duct'
    interface_boundary_id_shift = 10
  []
  [duct_pattern]
    type = PatternedHexMeshGenerator
    inputs = 'dummy duct1 duct2 duct3 duct4 duct5 duct6 duct7 duct8 duct9 duct10 duct11 duct12 duct13 duct14 duct15 duct16 duct17 duct18 duct19'
    pattern =
             '12 11 10;
			  13 4 3 9;
             14 5 1 2 8;
              15 6 7 19;
		       16 17 18'
    pattern_boundary = none
  []
  [duct_center_removal]
    type = BlockDeletionGenerator
	input = duct_pattern
    block = '100'
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
  [aclp3]
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
    duct_block_ids = '3 100'
	preserve_volumes = on
	quad_center_elements = true
    create_inward_interface_boundaries = True
    interface_boundary_id_shift = 20
  []
  [aclp4]
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
    duct_block_ids = '4 100'
	preserve_volumes = on
	quad_center_elements = true
    create_inward_interface_boundaries = True
    interface_boundary_id_shift = 20
  []
  [aclp5]
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
    duct_block_ids = '5 100'
	preserve_volumes = on
	quad_center_elements = true
    create_inward_interface_boundaries = True
    interface_boundary_id_shift = 20
  []
  [aclp6]
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
    duct_block_ids = '6 100'
	preserve_volumes = on
	quad_center_elements = true
    create_inward_interface_boundaries = True
    interface_boundary_id_shift = 20
  []
  [aclp7]
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
    duct_block_ids = '7 100'
	preserve_volumes = on
	quad_center_elements = true
    create_inward_interface_boundaries = True
    interface_boundary_id_shift = 20
  []
  [aclp8]
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
    duct_block_ids = '8 100'
	preserve_volumes = on
	quad_center_elements = true
    create_inward_interface_boundaries = True
    interface_boundary_id_shift = 20
  []
  [aclp9]
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
    duct_block_ids = '9 100'
	preserve_volumes = on
	quad_center_elements = true
    create_inward_interface_boundaries = True
    interface_boundary_id_shift = 20
  []
  [aclp10]
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
    duct_block_ids = '10 100'
	preserve_volumes = on
	quad_center_elements = true
    create_inward_interface_boundaries = True
    interface_boundary_id_shift = 20
  []
  [aclp11]
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
    duct_block_ids = '11 100'
	preserve_volumes = on
	quad_center_elements = true
    create_inward_interface_boundaries = True
    interface_boundary_id_shift = 20
  []
  [aclp12]
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
    duct_block_ids = '12 100'
	preserve_volumes = on
	quad_center_elements = true
    create_inward_interface_boundaries = True
    interface_boundary_id_shift = 20
  []
  [aclp13]
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
    duct_block_ids = '13 100'
	preserve_volumes = on
	quad_center_elements = true
    create_inward_interface_boundaries = True
    interface_boundary_id_shift = 20
  []
  [aclp14]
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
    duct_block_ids = '14 100'
	preserve_volumes = on
	quad_center_elements = true
    create_inward_interface_boundaries = True
    interface_boundary_id_shift = 20
  []
  [aclp15]
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
    duct_block_ids = '15 100'
	preserve_volumes = on
	quad_center_elements = true
    create_inward_interface_boundaries = True
    interface_boundary_id_shift = 20
  []
  [aclp16]
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
    duct_block_ids = '16 100'
	preserve_volumes = on
	quad_center_elements = true
    create_inward_interface_boundaries = True
    interface_boundary_id_shift = 20
  []
  [aclp17]
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
    duct_block_ids = '17 100'
	preserve_volumes = on
	quad_center_elements = true
    create_inward_interface_boundaries = True
    interface_boundary_id_shift = 20
  []
  [aclp18]
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
    duct_block_ids = '18 100'
	preserve_volumes = on
	quad_center_elements = true
    create_inward_interface_boundaries = True
    interface_boundary_id_shift = 20
  []
  [aclp19]
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
    duct_block_ids = '19 100'
	preserve_volumes = on
	quad_center_elements = true
    create_inward_interface_boundaries = True
    interface_boundary_id_shift = 20
  []
  [aclp_pattern]
    type = PatternedHexMeshGenerator
    inputs = 'dummy aclp1 aclp2 aclp3 aclp4 aclp5 aclp6 aclp7 aclp8 aclp9 aclp10 aclp11 aclp12 aclp13 aclp14 aclp15 aclp16 aclp17 aclp18 aclp19'
    pattern =
             '12 11 10;
			  13 4 3 9;
             14 5 1 2 8;
              15 6 7 19;
		       16 17 18'
	interface_boundary_id_shift_pattern =
               '1200 1100 1000; 
			   1300 400 300 900;
              1400 500 100 200 800;
               1500 600 700 1900;
			    1600 1700 1800'
    pattern_boundary = none
  []
  [aclp_center_removal]
    type = BlockDeletionGenerator
	input = aclp_pattern
    block = '100'
  []
  [aclp_translate]
    type = TransformGenerator
    #input = aclp_rename_boundary
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
    old_boundary = '122 222 322 422 522 622 722 822 922 1022 1122 1222 1322 1422 1522 1622 1722 1822 1922 '
				   '123 223 323 423 523 623 723 823 923 1023 1123 1223 1323 1423 1523 1623 1723 1823 1923 '
				   '1924 1925'
    new_boundary = 'aclp_inside aclp_inside aclp_inside aclp_inside aclp_inside aclp_inside aclp_inside '
	               'aclp_inside aclp_inside aclp_inside aclp_inside aclp_inside aclp_inside aclp_inside '
				   'aclp_inside aclp_inside aclp_inside aclp_inside aclp_inside '
				   'aclp_1 aclp_2 aclp_3 aclp_4 aclp_5 aclp_6 aclp_7 aclp_8 aclp_9 aclp_10 aclp_11 aclp_12 aclp_13 aclp_14 '
				   'aclp_15 aclp_16 aclp_17 aclp_18 aclp_19 '
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
  [tlp3]
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
    duct_block_ids = '3 100'
	preserve_volumes = on
	quad_center_elements = true
    create_inward_interface_boundaries = True
    interface_boundary_id_shift = 30
  []
  [tlp4]
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
    duct_block_ids = '4 100'
	preserve_volumes = on
	quad_center_elements = true
    create_inward_interface_boundaries = True
    interface_boundary_id_shift = 30
  []
  [tlp5]
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
    duct_block_ids = '5 100'
	preserve_volumes = on
	quad_center_elements = true
    create_inward_interface_boundaries = True
    interface_boundary_id_shift = 30
  []
  [tlp6]
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
    duct_block_ids = '6 100'
	preserve_volumes = on
	quad_center_elements = true
    create_inward_interface_boundaries = True
    interface_boundary_id_shift = 30
  []
  [tlp7]
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
    duct_block_ids = '7 100'
	preserve_volumes = on
	quad_center_elements = true
    create_inward_interface_boundaries = True
    interface_boundary_id_shift = 30
  []
  [tlp8]
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
    duct_block_ids = '8 100'
	preserve_volumes = on
	quad_center_elements = true
    create_inward_interface_boundaries = True
    interface_boundary_id_shift = 30
  []
  [tlp9]
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
    duct_block_ids = '9 100'
	preserve_volumes = on
	quad_center_elements = true
    create_inward_interface_boundaries = True
    interface_boundary_id_shift = 30
  []
  [tlp10]
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
    duct_block_ids = '10 100'
	preserve_volumes = on
	quad_center_elements = true
    create_inward_interface_boundaries = True
    interface_boundary_id_shift = 30
  []
  [tlp11]
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
    duct_block_ids = '11 100'
	preserve_volumes = on
	quad_center_elements = true
    create_inward_interface_boundaries = True
    interface_boundary_id_shift = 30
  []
  [tlp12]
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
    duct_block_ids = '12 100'
	preserve_volumes = on
	quad_center_elements = true
    create_inward_interface_boundaries = True
    interface_boundary_id_shift = 30
  []
  [tlp13]
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
    duct_block_ids = '13 100'
	preserve_volumes = on
	quad_center_elements = true
    create_inward_interface_boundaries = True
    interface_boundary_id_shift = 30
  []
  [tlp14]
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
    duct_block_ids = '14 100'
	preserve_volumes = on
	quad_center_elements = true
    create_inward_interface_boundaries = True
    interface_boundary_id_shift = 30
  []
  [tlp15]
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
    duct_block_ids = '15 100'
	preserve_volumes = on
	quad_center_elements = true
    create_inward_interface_boundaries = True
    interface_boundary_id_shift = 30
  []
  [tlp16]
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
    duct_block_ids = '16 100'
	preserve_volumes = on
	quad_center_elements = true
    create_inward_interface_boundaries = True
    interface_boundary_id_shift = 30
  []
  [tlp17]
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
    duct_block_ids = '17 100'
	preserve_volumes = on
	quad_center_elements = true
    create_inward_interface_boundaries = True
    interface_boundary_id_shift = 30
  []
  [tlp18]
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
    duct_block_ids = '18 100'
	preserve_volumes = on
	quad_center_elements = true
    create_inward_interface_boundaries = True
    interface_boundary_id_shift = 30
  []
  [tlp19]
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
    duct_block_ids = '19 100'
	preserve_volumes = on
	quad_center_elements = true
    create_inward_interface_boundaries = True
    interface_boundary_id_shift = 30
  []
  [tlp_pattern]
    type = PatternedHexMeshGenerator
    inputs = 'dummy tlp1 tlp2 tlp3 tlp4 tlp5 tlp6 tlp7 tlp8 tlp9 tlp10 tlp11 tlp12 tlp13 tlp14 tlp15 tlp16 tlp17 tlp18 tlp19'
    pattern =
             '12 11 10;
			  13 4 3 9;
             14 5 1 2 8;
              15 6 7 19;
		       16 17 18'
	interface_boundary_id_shift_pattern =
               '1200 1100 1000; 
			   1300 400 300 900;
              1400 500 100 200 800;
               1500 600 700 1900;
			    1600 1700 1800'
    pattern_boundary = none
  []
  [tlp_center_removal]
    type = BlockDeletionGenerator
	input = tlp_pattern
    block = '100'
  []
  [tlp_translate]
    type = TransformGenerator
    #input = tlp_rename_boundary
	input = tlp_center_removal
	transform = translate
    vector_value = '0 0 ${fparse tlp_location/100 - lp_length/2/100}'
  [] 
  [tlp_extrude]
    type = AdvancedExtruderGenerator
	input = tlp_translate
    direction = '0 0 1'
    heights = '${fparse lp_length/100}'
    num_layers = '${lp_n_ax}'
  []
  [tlp_boundary]
    type = RenameBoundaryGenerator
	input = tlp_extrude
    old_boundary = '132 232 332 432 532 632 732 832 932 1032 1132 1232 1332 1432 1532 1632 1732 1832 1932 '
	               '133 233 333 433 533 633 733 833 933 1033 1133 1233 1333 1433 1533 1633 1733 1833 1933 '
				   '1934 1935'
    new_boundary = 'tlp_inside tlp_inside tlp_inside tlp_inside tlp_inside tlp_inside tlp_inside tlp_inside tlp_inside tlp_inside tlp_inside tlp_inside tlp_inside tlp_inside '
	               'tlp_inside tlp_inside tlp_inside tlp_inside tlp_inside '
				   'tlp_1 tlp_2 tlp_3 tlp_4 tlp_5 tlp_6 tlp_7 tlp_8 tlp_9 tlp_10 tlp_11 tlp_12 tlp_13 tlp_14 '
				   'tlp_15 tlp_16 tlp_17 tlp_18 tlp_19 '
				   'tlp_bottom tlp_top'
  []

  [restraint_ring]
    type = PolygonConcentricCircleMeshGenerator
    num_sides = 6
    num_sectors_per_side = '${ns} ${ns} ${ns} ${ns} ${ns} ${ns}'
    background_intervals = 1
    background_block_ids = '101'
	polygon_size_style = 'apothem'
    polygon_size = ${fparse duct_outer_ftf /2/100 + duct_gap /2 /100}
	duct_sizes_style = apothem
    duct_sizes = '${fparse lp_outer_ftf /2 /100}'
    duct_intervals = 1
    duct_block_ids = '21'
	preserve_volumes = on
	quad_center_elements = true
    create_inward_interface_boundaries = True
    interface_boundary_id_shift = 30
  []
  [restraint_pattern]
    type = PatternedHexMeshGenerator
    inputs = 'dummy restraint_ring'
    pattern =
             '1 1 1 1;
			 1 0 0 0 1;
			1 0 0 0 0 1;
           1 0 0 0 0 0 1;
            1 0 0 0 0 1;
		     1 0 0 0 1;
			  1 1 1 1'
    pattern_boundary = none
  []
  [remove_dummy]
    type = BlockDeletionGenerator
	input = restraint_pattern
    block = '21 1500 100'
  []
  [restraint_translate]
    type = TransformGenerator
	input = remove_dummy
	transform = translate
    vector_value = '0 0 ${fparse tlp_location/100 - lp_length/2/100}'
  [] 
  [restraint_extrude]
    type = AdvancedExtruderGenerator
	input = restraint_translate
    direction = '0 0 1'
    heights = '${fparse lp_length/100}'
    num_layers = '${lp_n_ax}'
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
  [cmbn]
    type = CombinerGenerator
	inputs = 'aclp_stitching restraint_extrude'
  []
  #final_generator = restraint_extrude
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

[Physics/SolidMechanics/QuasiStatic]
  [ducts]
    strain = FINITE
	add_variables = true
    volumetric_locking_correction = true
    eigenstrain_names = thermal_expansion
    decomposition_method = EigenSolution
    generate_output = 'vonmises_stress'
    temperature = temp
    use_finite_deform_jacobian = true
    extra_vector_tags = 'ref'
	block = '1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 101'
	displacements = 'disp_x disp_y disp_z'
  []
[]

[AuxVariables]
  [temp]
	order = FIRST
	family = LAGRANGE
	initial_condition = ${hot_nominal_temp}
  []
  [aclp1_2_force]
  []
  [aclp1_3_force]
  []
  [aclp1_4_force]
  []
  [aclp1_5_force]
  []
  [aclp1_6_force]
  []
  [aclp1_7_force]
  []
  [aclp2_3_force]
  []
  [aclp2_7_force]
  []
  [aclp2_8_force]
  []
  [aclp2_9_force]
  []
  [aclp2_19_force]
  []
  [aclp3_4_force]
  []
  [aclp3_9_force]
  []
  [aclp3_10_force]
  []
  [aclp3_11_force]
  []
  [aclp4_5_force]
  []
  [aclp4_11_force]
  []
  [aclp4_12_force]
  []
  [aclp4_13_force]
  []
  [aclp5_6_force]
  []
  [aclp5_13_force]
  []
  [aclp5_14_force]
  []
  [aclp5_15_force]
  []
  [aclp6_7_force]
  []
  [aclp6_15_force]
  []
  [aclp6_16_force]
  []
  [aclp6_17_force]
  []
  [aclp7_17_force]
  []
  [aclp7_18_force]
  []
  [aclp7_19_force]
  []
  [aclp8_9_force]
  []
  [aclp8_19_force]
  []
  [aclp9_10_force]
  []
  [aclp10_11_force]
  []
  [aclp11_12_force]
  []
  [aclp12_13_force]
  []
  [aclp13_14_force]
  []
  [aclp14_15_force]
  []
  [aclp15_16_force]
  []
  [aclp16_17_force]
  []
  [aclp17_18_force]
  []
  [aclp18_19_force]
  []
  
  [tlp1_2_force]
  []
  [tlp1_3_force]
  []
  [tlp1_4_force]
  []
  [tlp1_5_force]
  []
  [tlp1_6_force]
  []
  [tlp1_7_force]
  []
  [tlp2_3_force]
  []
  [tlp2_7_force]
  []
  [tlp2_8_force]
  []
  [tlp2_9_force]
  []
  [tlp2_19_force]
  []
  [tlp3_4_force]
  []
  [tlp3_9_force]
  []
  [tlp3_10_force]
  []
  [tlp3_11_force]
  []
  [tlp4_5_force]
  []
  [tlp4_11_force]
  []
  [tlp4_12_force]
  []
  [tlp4_13_force]
  []
  [tlp5_6_force]
  []
  [tlp5_13_force]
  []
  [tlp5_14_force]
  []
  [tlp5_15_force]
  []
  [tlp6_7_force]
  []
  [tlp6_15_force]
  []
  [tlp6_16_force]
  []
  [tlp6_17_force]
  []
  [tlp7_17_force]
  []
  [tlp7_18_force]
  []
  [tlp7_19_force]
  []
  [tlp8_9_force]
  []
  [tlp8_19_force]
  []
  [tlp9_10_force]
  []
  [tlp10_11_force]
  []
  [tlp11_12_force]
  []
  [tlp12_13_force]
  []
  [tlp13_14_force]
  []
  [tlp14_15_force]
  []
  [tlp15_16_force]
  []
  [tlp16_17_force]
  []
  [tlp17_18_force]
  []
  [tlp18_19_force]
  []
  
  [tlp8_restraint_force]
  []
  [tlp9_restraint_force]
  []
  [tlp10_restraint_force]
  []
  [tlp11_restraint_force]
  []
  [tlp12_restraint_force]
  []
  [tlp13_restraint_force]
  []
  [tlp14_restraint_force]
  []
  [tlp15_restraint_force]
  []
  [tlp16_restraint_force]
  []
  [tlp17_restraint_force]
  []
  [tlp18_restraint_force]
  []
  [tlp19_restraint_force]
  []
[]

[AuxKernels]
  [tfunc_duct2]
    type = FunctionAux
    function = temp_duct2
    block = 2
    variable = temp
  []
  [tfunc_duct3]
    type = FunctionAux
    function = temp_duct3
    block = 3
    variable = temp
  []
  [tfunc_duct4]
    type = FunctionAux
    function = temp_duct4
    block = 4
    variable = temp
  []
  [tfunc_duct5]
    type = FunctionAux
    function = temp_duct5
    block = 5
    variable = temp
  []
  [tfunc_duct6]
    type = FunctionAux
    function = temp_duct6
    block = 6
    variable = temp
  []
  [tfunc_duct7]
    type = FunctionAux
    function = temp_duct7
    block = 7
    variable = temp
  []
  [aclp1_2]
    type = PenetrationAux
    variable = aclp1_2_force
    boundary = 'aclp_2'
    paired_boundary = 'aclp_1'
    quantity = normal_force_magnitude
  []
  [aclp1_3]
    type = PenetrationAux
    variable = aclp1_3_force
    boundary = 'aclp_3'
    paired_boundary = 'aclp_1'
    quantity = normal_force_magnitude
  []
  [aclp1_4]
    type = PenetrationAux
    variable = aclp1_4_force
    boundary = 'aclp_4'
    paired_boundary = 'aclp_1'
    quantity = normal_force_magnitude
  []
  [aclp1_5]
    type = PenetrationAux
    variable = aclp1_5_force
    boundary = 'aclp_5'
    paired_boundary = 'aclp_1'
    quantity = normal_force_magnitude
  []
  [aclp1_6]
    type = PenetrationAux
    variable = aclp1_6_force
    boundary = 'aclp_6'
    paired_boundary = 'aclp_1'
    quantity = normal_force_magnitude
  []
  [aclp1_7]
    type = PenetrationAux
    variable = aclp1_7_force
    boundary = 'aclp_7'
    paired_boundary = 'aclp_1'
    quantity = normal_force_magnitude
  []
  [aclp2_3]
    type = PenetrationAux
    variable = aclp2_3_force
    boundary = 'aclp_3'
    paired_boundary = 'aclp_2'
    quantity = normal_force_magnitude
  []
  [aclp2_7]
    type = PenetrationAux
    variable = aclp2_7_force
    boundary = 'aclp_7'
    paired_boundary = 'aclp_2'
    quantity = normal_force_magnitude
  []
  [aclp2_8]
    type = PenetrationAux
    variable = aclp2_8_force
    boundary = 'aclp_8'
    paired_boundary = 'aclp_2'
    quantity = normal_force_magnitude
  []
  [aclp2_9]
    type = PenetrationAux
    variable = aclp2_9_force
    boundary = 'aclp_9'
    paired_boundary = 'aclp_2'
    quantity = normal_force_magnitude
  []
  [aclp2_19]
    type = PenetrationAux
    variable = aclp2_19_force
    boundary = 'aclp_19'
    paired_boundary = 'aclp_2'
    quantity = normal_force_magnitude
  []
  [aclp3_4]
    type = PenetrationAux
    variable = aclp3_4_force
    boundary = 'aclp_4'
    paired_boundary = 'aclp_3'
    quantity = normal_force_magnitude
  []
  [aclp3_9]
    type = PenetrationAux
    variable = aclp3_9_force
    boundary = 'aclp_9'
    paired_boundary = 'aclp_3'
    quantity = normal_force_magnitude
  []
  [aclp3_10]
    type = PenetrationAux
    variable = aclp3_10_force
    boundary = 'aclp_10'
    paired_boundary = 'aclp_3'
    quantity = normal_force_magnitude
  []
  [aclp3_11]
    type = PenetrationAux
    variable = aclp3_11_force
    boundary = 'aclp_11'
    paired_boundary = 'aclp_3'
    quantity = normal_force_magnitude
  []
  [aclp4_5]
    type = PenetrationAux
    variable = aclp4_5_force
    boundary = 'aclp_5'
    paired_boundary = 'aclp_4'
    quantity = normal_force_magnitude
  []
  [aclp4_11]
    type = PenetrationAux
    variable = aclp4_11_force
    boundary = 'aclp_11'
    paired_boundary = 'aclp_4'
    quantity = normal_force_magnitude
  []
  [aclp4_12]
    type = PenetrationAux
    variable = aclp4_12_force
    boundary = 'aclp_12'
    paired_boundary = 'aclp_4'
    quantity = normal_force_magnitude
  []
  [aclp4_13]
    type = PenetrationAux
    variable = aclp4_13_force
    boundary = 'aclp_13'
    paired_boundary = 'aclp_4'
    quantity = normal_force_magnitude
  []
  [aclp5_6]
    type = PenetrationAux
    variable = aclp5_6_force
    boundary = 'aclp_6'
    paired_boundary = 'aclp_5'
    quantity = normal_force_magnitude
  []
  [aclp5_13]
    type = PenetrationAux
    variable = aclp5_13_force
    boundary = 'aclp_13'
    paired_boundary = 'aclp_5'
    quantity = normal_force_magnitude
  []
  [aclp5_14]
    type = PenetrationAux
    variable = aclp5_14_force
    boundary = 'aclp_14'
    paired_boundary = 'aclp_5'
    quantity = normal_force_magnitude
  []
  [aclp5_15]
    type = PenetrationAux
    variable = aclp5_15_force
    boundary = 'aclp_15'
    paired_boundary = 'aclp_5'
    quantity = normal_force_magnitude
  []
  [aclp6_7]
    type = PenetrationAux
    variable = aclp6_7_force
    boundary = 'aclp_7'
    paired_boundary = 'aclp_6'
    quantity = normal_force_magnitude
  []
  [aclp6_15]
    type = PenetrationAux
    variable = aclp6_15_force
    boundary = 'aclp_15'
    paired_boundary = 'aclp_6'
    quantity = normal_force_magnitude
  []
  [aclp6_16]
    type = PenetrationAux
    variable = aclp6_16_force
    boundary = 'aclp_16'
    paired_boundary = 'aclp_6'
    quantity = normal_force_magnitude
  []
  [aclp6_17]
    type = PenetrationAux
    variable = aclp6_17_force
    boundary = 'aclp_17'
    paired_boundary = 'aclp_6'
    quantity = normal_force_magnitude
  []
  [aclp7_17]
    type = PenetrationAux
    variable = aclp7_17_force
    boundary = 'aclp_17'
    paired_boundary = 'aclp_7'
    quantity = normal_force_magnitude
  []
  [aclp7_18]
    type = PenetrationAux
    variable = aclp7_18_force
    boundary = 'aclp_18'
    paired_boundary = 'aclp_7'
    quantity = normal_force_magnitude
  []
  [aclp7_19]
    type = PenetrationAux
    variable = aclp7_19_force
    boundary = 'aclp_19'
    paired_boundary = 'aclp_7'
    quantity = normal_force_magnitude
  []
  [aclp8_9]
    type = PenetrationAux
    variable = aclp8_9_force
    boundary = 'aclp_9'
    paired_boundary = 'aclp_8'
    quantity = normal_force_magnitude
  []
  [aclp8_19]
    type = PenetrationAux
    variable = aclp8_19_force
    boundary = 'aclp_19'
    paired_boundary = 'aclp_8'
    quantity = normal_force_magnitude
  []
  [aclp9_10]
    type = PenetrationAux
    variable = aclp9_10_force
    boundary = 'aclp_10'
    paired_boundary = 'aclp_9'
    quantity = normal_force_magnitude
  []
  [aclp10_11]
    type = PenetrationAux
    variable = aclp10_11_force
    boundary = 'aclp_11'
    paired_boundary = 'aclp_10'
    quantity = normal_force_magnitude
  []
  [aclp11_12]
    type = PenetrationAux
    variable = aclp11_12_force
    boundary = 'aclp_12'
    paired_boundary = 'aclp_11'
    quantity = normal_force_magnitude
  []
  [aclp12_13]
    type = PenetrationAux
    variable = aclp12_13_force
    boundary = 'aclp_13'
    paired_boundary = 'aclp_12'
    quantity = normal_force_magnitude
  []
  [aclp13_14]
    type = PenetrationAux
    variable = aclp13_14_force
    boundary = 'aclp_14'
    paired_boundary = 'aclp_13'
    quantity = normal_force_magnitude
  []
  [aclp14_15]
    type = PenetrationAux
    variable = aclp14_15_force
    boundary = 'aclp_15'
    paired_boundary = 'aclp_14'
    quantity = normal_force_magnitude
  []
  [aclp15_16]
    type = PenetrationAux
    variable = aclp15_16_force
    boundary = 'aclp_16'
    paired_boundary = 'aclp_15'
    quantity = normal_force_magnitude
  []
  [aclp16_17]
    type = PenetrationAux
    variable = aclp16_17_force
    boundary = 'aclp_17'
    paired_boundary = 'aclp_16'
    quantity = normal_force_magnitude
  []
  [aclp17_18]
    type = PenetrationAux
    variable = aclp17_18_force
    boundary = 'aclp_18'
    paired_boundary = 'aclp_17'
    quantity = normal_force_magnitude
  []
  [aclp18_19]
    type = PenetrationAux
    variable = aclp18_19_force
    boundary = 'aclp_19'
    paired_boundary = 'aclp_18'
    quantity = normal_force_magnitude
  []

  [tlp1_2]
    type = PenetrationAux
    variable = tlp1_2_force
    boundary = 'tlp_2'
    paired_boundary = 'tlp_1'
    quantity = normal_force_magnitude
  []
  [tlp1_3]
    type = PenetrationAux
    variable = tlp1_3_force
    boundary = 'tlp_3'
    paired_boundary = 'tlp_1'
    quantity = normal_force_magnitude
  []
  [tlp1_4]
    type = PenetrationAux
    variable = tlp1_4_force
    boundary = 'tlp_4'
    paired_boundary = 'tlp_1'
    quantity = normal_force_magnitude
  []
  [tlp1_5]
    type = PenetrationAux
    variable = tlp1_5_force
    boundary = 'tlp_5'
    paired_boundary = 'tlp_1'
    quantity = normal_force_magnitude
  []
  [tlp1_6]
    type = PenetrationAux
    variable = tlp1_6_force
    boundary = 'tlp_6'
    paired_boundary = 'tlp_1'
    quantity = normal_force_magnitude
  []
  [tlp1_7]
    type = PenetrationAux
    variable = tlp1_7_force
    boundary = 'tlp_7'
    paired_boundary = 'tlp_1'
    quantity = normal_force_magnitude
  []
  [tlp2_3]
    type = PenetrationAux
    variable = tlp2_3_force
    boundary = 'tlp_3'
    paired_boundary = 'tlp_2'
    quantity = normal_force_magnitude
  []
  [tlp2_7]
    type = PenetrationAux
    variable = tlp2_7_force
    boundary = 'tlp_7'
    paired_boundary = 'tlp_2'
    quantity = normal_force_magnitude
  []
  [tlp2_8]
    type = PenetrationAux
    variable = tlp2_8_force
    boundary = 'tlp_8'
    paired_boundary = 'tlp_2'
    quantity = normal_force_magnitude
  []
  [tlp2_9]
    type = PenetrationAux
    variable = tlp2_9_force
    boundary = 'tlp_9'
    paired_boundary = 'tlp_2'
    quantity = normal_force_magnitude
  []
  [tlp2_19]
    type = PenetrationAux
    variable = tlp2_19_force
    boundary = 'tlp_19'
    paired_boundary = 'tlp_2'
    quantity = normal_force_magnitude
  []
  [tlp3_4]
    type = PenetrationAux
    variable = tlp3_4_force
    boundary = 'tlp_4'
    paired_boundary = 'tlp_3'
    quantity = normal_force_magnitude
  []
  [tlp3_9]
    type = PenetrationAux
    variable = tlp3_9_force
    boundary = 'tlp_9'
    paired_boundary = 'tlp_3'
    quantity = normal_force_magnitude
  []
  [tlp3_10]
    type = PenetrationAux
    variable = tlp3_10_force
    boundary = 'tlp_10'
    paired_boundary = 'tlp_3'
    quantity = normal_force_magnitude
  []
  [tlp3_11]
    type = PenetrationAux
    variable = tlp3_11_force
    boundary = 'tlp_11'
    paired_boundary = 'tlp_3'
    quantity = normal_force_magnitude
  []
  [tlp4_5]
    type = PenetrationAux
    variable = tlp4_5_force
    boundary = 'tlp_5'
    paired_boundary = 'tlp_4'
    quantity = normal_force_magnitude
  []
  [tlp4_11]
    type = PenetrationAux
    variable = tlp4_11_force
    boundary = 'tlp_11'
    paired_boundary = 'tlp_4'
    quantity = normal_force_magnitude
  []
  [tlp4_12]
    type = PenetrationAux
    variable = tlp4_12_force
    boundary = 'tlp_12'
    paired_boundary = 'tlp_4'
    quantity = normal_force_magnitude
  []
  [tlp4_13]
    type = PenetrationAux
    variable = tlp4_13_force
    boundary = 'tlp_13'
    paired_boundary = 'tlp_4'
    quantity = normal_force_magnitude
  []
  [tlp5_6]
    type = PenetrationAux
    variable = tlp5_6_force
    boundary = 'tlp_6'
    paired_boundary = 'tlp_5'
    quantity = normal_force_magnitude
  []
  [tlp5_13]
    type = PenetrationAux
    variable = tlp5_13_force
    boundary = 'tlp_13'
    paired_boundary = 'tlp_5'
    quantity = normal_force_magnitude
  []
  [tlp5_14]
    type = PenetrationAux
    variable = tlp5_14_force
    boundary = 'tlp_14'
    paired_boundary = 'tlp_5'
    quantity = normal_force_magnitude
  []
  [tlp5_15]
    type = PenetrationAux
    variable = tlp5_15_force
    boundary = 'tlp_15'
    paired_boundary = 'tlp_5'
    quantity = normal_force_magnitude
  []
  [tlp6_7]
    type = PenetrationAux
    variable = tlp6_7_force
    boundary = 'tlp_7'
    paired_boundary = 'tlp_6'
    quantity = normal_force_magnitude
  []
  [tlp6_15]
    type = PenetrationAux
    variable = tlp6_15_force
    boundary = 'tlp_15'
    paired_boundary = 'tlp_6'
    quantity = normal_force_magnitude
  []
  [tlp6_16]
    type = PenetrationAux
    variable = tlp6_16_force
    boundary = 'tlp_16'
    paired_boundary = 'tlp_6'
    quantity = normal_force_magnitude
  []
  [tlp6_17]
    type = PenetrationAux
    variable = tlp6_17_force
    boundary = 'tlp_17'
    paired_boundary = 'tlp_6'
    quantity = normal_force_magnitude
  []
  [tlp7_17]
    type = PenetrationAux
    variable = tlp7_17_force
    boundary = 'tlp_17'
    paired_boundary = 'tlp_7'
    quantity = normal_force_magnitude
  []
  [tlp7_18]
    type = PenetrationAux
    variable = tlp7_18_force
    boundary = 'tlp_18'
    paired_boundary = 'tlp_7'
    quantity = normal_force_magnitude
  []
  [tlp7_19]
    type = PenetrationAux
    variable = tlp7_19_force
    boundary = 'tlp_19'
    paired_boundary = 'tlp_7'
    quantity = normal_force_magnitude
  []
  [tlp8_9]
    type = PenetrationAux
    variable = tlp8_9_force
    boundary = 'tlp_9'
    paired_boundary = 'tlp_8'
    quantity = normal_force_magnitude
  []
  [tlp8_19]
    type = PenetrationAux
    variable = tlp8_19_force
    boundary = 'tlp_19'
    paired_boundary = 'tlp_8'
    quantity = normal_force_magnitude
  []
  [tlp9_10]
    type = PenetrationAux
    variable = tlp9_10_force
    boundary = 'tlp_10'
    paired_boundary = 'tlp_9'
    quantity = normal_force_magnitude
  []
  [tlp10_11]
    type = PenetrationAux
    variable = tlp10_11_force
    boundary = 'tlp_11'
    paired_boundary = 'tlp_10'
    quantity = normal_force_magnitude
  []
  [tlp11_12]
    type = PenetrationAux
    variable = tlp11_12_force
    boundary = 'tlp_12'
    paired_boundary = 'tlp_11'
    quantity = normal_force_magnitude
  []
  [tlp12_13]
    type = PenetrationAux
    variable = tlp12_13_force
    boundary = 'tlp_13'
    paired_boundary = 'tlp_12'
    quantity = normal_force_magnitude
  []
  [tlp13_14]
    type = PenetrationAux
    variable = tlp13_14_force
    boundary = 'tlp_14'
    paired_boundary = 'tlp_13'
    quantity = normal_force_magnitude
  []
  [tlp14_15]
    type = PenetrationAux
    variable = tlp14_15_force
    boundary = 'tlp_15'
    paired_boundary = 'tlp_14'
    quantity = normal_force_magnitude
  []
  [tlp15_16]
    type = PenetrationAux
    variable = tlp15_16_force
    boundary = 'tlp_16'
    paired_boundary = 'tlp_15'
    quantity = normal_force_magnitude
  []
  [tlp16_17]
    type = PenetrationAux
    variable = tlp16_17_force
    boundary = 'tlp_17'
    paired_boundary = 'tlp_16'
    quantity = normal_force_magnitude
  []
  [tlp17_18]
    type = PenetrationAux
    variable = tlp17_18_force
    boundary = 'tlp_18'
    paired_boundary = 'tlp_17'
    quantity = normal_force_magnitude
  []
  [tlp18_19]
    type = PenetrationAux
    variable = tlp18_19_force
    boundary = 'tlp_19'
    paired_boundary = 'tlp_18'
    quantity = normal_force_magnitude
  []
  
  [tlp8_restraint]
    type = PenetrationAux
    variable = tlp8_restraint_force
    boundary = 31
    paired_boundary = 'tlp_8'
    quantity = normal_force_magnitude
  []
  
  [tlp9_restraint]
    type = PenetrationAux
    variable = tlp9_restraint_force
    boundary = 31
    paired_boundary = 'tlp_9'
    quantity = normal_force_magnitude
  []
  
  [tlp10_restraint]
    type = PenetrationAux
    variable = tlp10_restraint_force
    boundary = 31
    paired_boundary = 'tlp_10'
    quantity = normal_force_magnitude
  []
  
  [tlp11_restraint]
    type = PenetrationAux
    variable = tlp11_restraint_force
    boundary = 31
    paired_boundary = 'tlp_11'
    quantity = normal_force_magnitude
  []
  
  [tlp12_restraint]
    type = PenetrationAux
    variable = tlp12_restraint_force
    boundary = 31
    paired_boundary = 'tlp_12'
    quantity = normal_force_magnitude
  []
  
  [tlp13_restraint]
    type = PenetrationAux
    variable = tlp13_restraint_force
    boundary = 31
    paired_boundary = 'tlp_13'
    quantity = normal_force_magnitude
  []
  
  [tlp14_restraint]
    type = PenetrationAux
    variable = tlp14_restraint_force
    boundary = 31
    paired_boundary = 'tlp_14'
    quantity = normal_force_magnitude
  []
  
  [tlp15_restraint]
    type = PenetrationAux
    variable = tlp15_restraint_force
    boundary = 31
    paired_boundary = 'tlp_15'
    quantity = normal_force_magnitude
  []
  
  [tlp16_restraint]
    type = PenetrationAux
    variable = tlp16_restraint_force
    boundary = 31
    paired_boundary = 'tlp_16'
    quantity = normal_force_magnitude
  []
  
  [tlp17_restraint]
    type = PenetrationAux
    variable = tlp17_restraint_force
    boundary = 31
    paired_boundary = 'tlp_17'
    quantity = normal_force_magnitude
  []
  
  [tlp18_restraint]
    type = PenetrationAux
    variable = tlp18_restraint_force
    boundary = 31
    paired_boundary = 'tlp_18'
    quantity = normal_force_magnitude
  []
  
  [tlp19_restraint]
    type = PenetrationAux
    variable = tlp19_restraint_force
    boundary = 31
    paired_boundary = 'tlp_19'
    quantity = normal_force_magnitude
  []
[]

[Functions]
  [temp_duct2]
    type = ParsedFunction
	expression = '${hot_nominal_temp} + if(z > ${fparse active_fuel_top/100}, t*(125 - 25/${fparse duct_outer_ftf/2/100} * ((y - ${fparse (duct_outer_ftf + duct_gap)/100}))*t), '
									   'if(z > ${fparse active_fuel_bot/100}, t*(z - ${fparse active_fuel_bot/100})/(${fparse active_fuel_top/100 - active_fuel_bot/100})*(125 - 25/${fparse duct_outer_ftf/2/100} * ((y - ${fparse (duct_outer_ftf + duct_gap)/100}))*t), 0))'
  []
  [temp_duct3]
    type = ParsedFunction
	expression = '${hot_nominal_temp} + if(z > ${fparse active_fuel_top/100}, t*(125 - 25/${fparse duct_outer_ftf/2/100} * ( (x - ${fparse (duct_outer_ftf + duct_gap)/100} * cos(${theta3}))*sin(-pi/3) + (y - ${fparse (duct_outer_ftf + duct_gap)/100} * sin(${theta3}))*cos(-pi/3))*t), '
									   'if(z > ${fparse active_fuel_bot/100}, t*(z - ${fparse active_fuel_bot/100})/(${fparse active_fuel_top/100 - active_fuel_bot/100})*(125 - 25/${fparse duct_outer_ftf/2/100} * ( (x - ${fparse (duct_outer_ftf + duct_gap)/100} * cos(${theta3}))*sin(-pi/3) + (y - ${fparse (duct_outer_ftf + duct_gap)/100} * sin(${theta3}))*cos(-pi/3))*t), 0))'
  []
  [temp_duct4]
    type = ParsedFunction
	expression = '${hot_nominal_temp} + if(z > ${fparse active_fuel_top/100}, t*(125 + 25/${fparse duct_outer_ftf/2/100} * ( (x - ${fparse (duct_outer_ftf + duct_gap)/100} * cos(${theta4}))*sin(pi/3) + (y - ${fparse (duct_outer_ftf + duct_gap)/100} * sin(${theta4}))*cos(pi/3))*t), '
									   'if(z > ${fparse active_fuel_bot/100}, t*(z - ${fparse active_fuel_bot/100})/(${fparse active_fuel_top/100 - active_fuel_bot/100})*(125 + 25/${fparse duct_outer_ftf/2/100} * ( (x - ${fparse (duct_outer_ftf + duct_gap)/100} * cos(${theta4}))*sin(pi/3) + (y - ${fparse (duct_outer_ftf + duct_gap)/100} * sin(${theta4}))*cos(pi/3))*t), 0))'
  []
  [temp_duct5]
    type = ParsedFunction
	expression = '${hot_nominal_temp} + if(z > ${fparse active_fuel_top/100}, t*(125 + 25/${fparse duct_outer_ftf/2/100} * ((y + ${fparse (duct_outer_ftf + duct_gap)/100}))*t), '
									   'if(z > ${fparse active_fuel_bot/100}, t*(z - ${fparse active_fuel_bot/100})/(${fparse active_fuel_top/100 - active_fuel_bot/100})*(125 + 25/${fparse duct_outer_ftf/2/100} * ((y + ${fparse (duct_outer_ftf + duct_gap)/100}))*t), 0))'
  []
  [temp_duct6]
    type = ParsedFunction
	expression = '${hot_nominal_temp} + if(z > ${fparse active_fuel_top/100}, t*(125 + 25/${fparse duct_outer_ftf/2/100} * ( (x - ${fparse (duct_outer_ftf + duct_gap)/100} * cos(${theta6}))*sin(-pi/3) + (y - ${fparse (duct_outer_ftf + duct_gap)/100} * sin(${theta6}))*cos(-pi/3))*t), '
									   'if(z > ${fparse active_fuel_bot/100}, t*(z - ${fparse active_fuel_bot/100})/(${fparse active_fuel_top/100 - active_fuel_bot/100})*(125 + 25/${fparse duct_outer_ftf/2/100} * ( (x - ${fparse (duct_outer_ftf + duct_gap)/100} * cos(${theta6}))*sin(-pi/3) + (y - ${fparse (duct_outer_ftf + duct_gap)/100} * sin(${theta6}))*cos(-pi/3))*t), 0))'
  []
  [temp_duct7]
    type = ParsedFunction
	expression = '${hot_nominal_temp} + if(z > ${fparse active_fuel_top/100}, t*(125 - 25/${fparse duct_outer_ftf/2/100} * ( (x - ${fparse (duct_outer_ftf + duct_gap)/100} * cos(${theta7}))*sin(pi/3) + (y - ${fparse (duct_outer_ftf + duct_gap)/100} * sin(${theta7}))*cos(pi/3))*t), '
									   'if(z > ${fparse active_fuel_bot/100}, t*(z - ${fparse active_fuel_bot/100})/(${fparse active_fuel_top/100 - active_fuel_bot/100})*(125 - 25/${fparse duct_outer_ftf/2/100} * ( (x - ${fparse (duct_outer_ftf + duct_gap)/100} * cos(${theta7}))*sin(pi/3) + (y - ${fparse (duct_outer_ftf + duct_gap)/100} * sin(${theta7}))*cos(pi/3))*t), 0))'
  []
[]

[BCs]
  [fixed_x]
    type = DirichletBC
    variable = disp_x
    boundary = 'duct_bottom'
    value = 0.0
  []
  [fixed_y]
    type = DirichletBC
    variable = disp_y
    boundary = 'duct_bottom'
    value = 0.0
  []
  [fixed_z]
    type = DirichletBC
    variable = disp_z
    boundary = 'duct_bottom'
    value = 0.0
  []
  [fix_restraint_x]
    type = DirichletBC
	variable = disp_x
	boundary = 31
	value = 0.0
  []
  [fix_restraint_y]
    type = DirichletBC
	variable = disp_y
	boundary = 31
	value = 0.0
  []
  [fix_restraint_z]
    type = DirichletBC
	variable = disp_z
	boundary = 31
	value = 0.0
  []
[]

[Materials]
  [elasticity_tensor]
    type = ComputeIsotropicElasticityTensor
    youngs_modulus = 1.7e11
    poissons_ratio = 0.3
	block = '1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 101'
  []
  [small_stress]
    type = ComputeFiniteStrainElasticStress
	block = '1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 101'
  []
  [thermal_expansion_strain]
    type = ComputeThermalExpansionEigenstrain
    stress_free_temperature = ${fparse hot_nominal_temp}
    thermal_expansion_coeff = 18.0e-6
    temperature = temp
    eigenstrain_name = thermal_expansion
	block = '1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 101'
  []
[]

[Contact]
  [aclp]
    primary =   'aclp_1  aclp_1  aclp_1  aclp_1  aclp_1  aclp_1 '
	            'aclp_2  aclp_2  aclp_2  aclp_2  aclp_2  aclp_3 '
				'aclp_3  aclp_3  aclp_3  aclp_4  aclp_4  aclp_4 '
				'aclp_4  aclp_5  aclp_5  aclp_5  aclp_5  aclp_6 '
				'aclp_6  aclp_6  aclp_6  aclp_7  aclp_7  aclp_7 '
				'aclp_8  aclp_8  aclp_9  aclp_10 aclp_11 aclp_12 '
				'aclp_13 aclp_14 aclp_15 aclp_16 aclp_17 aclp_18'
	secondary = 'aclp_2  aclp_3  aclp_4  aclp_5  aclp_6  aclp_7 '
	            'aclp_3  aclp_7  aclp_8  aclp_9  aclp_19 aclp_4 '
				'aclp_9  aclp_10 aclp_11 aclp_5  aclp_11 aclp_12 '
				'aclp_13 aclp_6  aclp_13 aclp_14 aclp_15 aclp_7 '
				'aclp_15 aclp_16 aclp_17 aclp_17 aclp_18 aclp_19 '
				'aclp_9  aclp_19 aclp_10 aclp_11 aclp_12 aclp_13 '
				'aclp_14 aclp_15 aclp_16 aclp_17 aclp_18 aclp_19'
    model = frictionless
    formulation = penalty
    normalize_penalty = true
    penalty = 1e10
    tangential_tolerance = 0.0
    normal_smoothing_distance = 0.0
  []
  [tlp]
    primary =   'tlp_1  tlp_1  tlp_1  tlp_1  tlp_1  tlp_1 '
	            'tlp_2  tlp_2  tlp_2  tlp_2  tlp_2  tlp_3 '
				'tlp_3  tlp_3  tlp_3  tlp_4  tlp_4  tlp_4 '
				'tlp_4  tlp_5  tlp_5  tlp_5  tlp_5  tlp_6 '
				'tlp_6  tlp_6  tlp_6  tlp_7  tlp_7  tlp_7 '
				'tlp_8  tlp_8  tlp_9  tlp_10 tlp_11 tlp_12 '
				'tlp_13 tlp_14 tlp_15 tlp_16 tlp_17 tlp_18'
	secondary = 'tlp_2  tlp_3  tlp_4  tlp_5  tlp_6  tlp_7 '
	            'tlp_3  tlp_7  tlp_8  tlp_9  tlp_19 tlp_4 '
				'tlp_9  tlp_10 tlp_11 tlp_5  tlp_11 tlp_12 '
				'tlp_13 tlp_6  tlp_13 tlp_14 tlp_15 tlp_7 '
				'tlp_15 tlp_16 tlp_17 tlp_17 tlp_18 tlp_19 '
				'tlp_9  tlp_19 tlp_10 tlp_11 tlp_12 tlp_13 '
				'tlp_14 tlp_15 tlp_16 tlp_17 tlp_18 tlp_19'
    model = frictionless
    formulation = penalty
    normalize_penalty = true
    penalty = 5e10
    tangential_tolerance = 0.0
    normal_smoothing_distance = 0.0
  []
  [restraint]
    primary =   'tlp_8  tlp_9  tlp_10 tlp_11 tlp_12 tlp_13 '
				'tlp_14 tlp_15 tlp_16 tlp_17 tlp_18 tlp_19'
	secondary = '   31     31      31     31     31     31 '
				'   31     31      31     31     31     31 '
    model = frictionless
    formulation = penalty
    normalize_penalty = true
    penalty = 5e10
    tangential_tolerance = 0.0
    normal_smoothing_distance = 0.0
  []
[]

[Postprocessors]
  [aclp_force_1_2]
    type = NodalSum
    variable = aclp1_2_force
    boundary = 'aclp_2'
  []
  [aclp_force_1_3]
    type = NodalSum
    variable = aclp1_3_force
    boundary = 'aclp_3'
  []
  [aclp_force_1_4]
    type = NodalSum
    variable = aclp1_4_force
    boundary = 'aclp_4'
  []
  [aclp_force_1_5]
    type = NodalSum
    variable = aclp1_5_force
    boundary = 'aclp_5'
  []
  [aclp_force_1_6]
    type = NodalSum
    variable = aclp1_6_force
    boundary = 'aclp_6'
  []
  [aclp_force_1_7]
    type = NodalSum
    variable = aclp1_7_force
    boundary = 'aclp_7'
  []
  [aclp_force_2_3]
    type = NodalSum
    variable = aclp2_3_force
    boundary = 'aclp_3'
  []
  [aclp_force_2_7]
    type = NodalSum
    variable = aclp2_7_force
    boundary = 'aclp_7'
  []
  [aclp_force_2_8]
    type = NodalSum
    variable = aclp2_8_force
    boundary = 'aclp_8'
  []
  [aclp_force_2_9]
    type = NodalSum
    variable = aclp2_9_force
    boundary = 'aclp_9'
  []
  [aclp_force_2_19]
    type = NodalSum
    variable = aclp2_19_force
    boundary = 'aclp_19'
  []
  [aclp_force_3_4]
    type = NodalSum
    variable = aclp3_4_force
    boundary = 'aclp_4'
  []
  [aclp_force_3_9]
    type = NodalSum
    variable = aclp3_9_force
    boundary = 'aclp_9'
  []
  [aclp_force_3_10]
    type = NodalSum
    variable = aclp3_10_force
    boundary = 'aclp_10'
  []
  [aclp_force_3_11]
    type = NodalSum
    variable = aclp3_11_force
    boundary = 'aclp_11'
  []
  [aclp_force_4_5]
    type = NodalSum
    variable = aclp4_5_force
    boundary = 'aclp_5'
  []
  [aclp_force_4_11]
    type = NodalSum
    variable = aclp4_11_force
    boundary = 'aclp_11'
  []
  [aclp_force_4_12]
    type = NodalSum
    variable = aclp4_12_force
    boundary = 'aclp_12'
  []
  [aclp_force_4_13]
    type = NodalSum
    variable = aclp4_13_force
    boundary = 'aclp_13'
  []
  [aclp_force_5_6]
    type = NodalSum
    variable = aclp5_6_force
    boundary = 'aclp_6'
  []
  [aclp_force_5_13]
    type = NodalSum
    variable = aclp5_13_force
    boundary = 'aclp_13'
  []
  [aclp_force_5_14]
    type = NodalSum
    variable = aclp5_14_force
    boundary = 'aclp_14'
  []
  [aclp_force_5_15]
    type = NodalSum
    variable = aclp5_15_force
    boundary = 'aclp_15'
  []
  [aclp_force_6_7]
    type = NodalSum
    variable = aclp6_7_force
    boundary = 'aclp_7'
  []
  [aclp_force_6_15]
    type = NodalSum
    variable = aclp6_15_force
    boundary = 'aclp_15'
  []
  [aclp_force_6_16]
    type = NodalSum
    variable = aclp6_16_force
    boundary = 'aclp_16'
  []
  [aclp_force_6_17]
    type = NodalSum
    variable = aclp6_17_force
    boundary = 'aclp_17'
  []
  [aclp_force_7_17]
    type = NodalSum
    variable = aclp7_17_force
    boundary = 'aclp_17'
  []
  [aclp_force_7_18]
    type = NodalSum
    variable = aclp7_18_force
    boundary = 'aclp_18'
  []
  [aclp_force_7_19]
    type = NodalSum
    variable = aclp7_19_force
    boundary = 'aclp_19'
  []
  [aclp_force_8_9]
    type = NodalSum
    variable = aclp8_9_force
    boundary = 'aclp_9'
  []
  [aclp_force_8_19]
    type = NodalSum
    variable = aclp8_19_force
    boundary = 'aclp_19'
  []
  [aclp_force_9_10]
    type = NodalSum
    variable = aclp9_10_force
    boundary = 'aclp_10'
  []
  [aclp_force_10_11]
    type = NodalSum
    variable = aclp10_11_force
    boundary = 'aclp_11'
  []
  [aclp_force_11_12]
    type = NodalSum
    variable = aclp11_12_force
    boundary = 'aclp_12'
  []
  [aclp_force_12_13]
    type = NodalSum
    variable = aclp12_13_force
    boundary = 'aclp_13'
  []
  [aclp_force_13_14]
    type = NodalSum
    variable = aclp13_14_force
    boundary = 'aclp_14'
  []
  [aclp_force_14_15]
    type = NodalSum
    variable = aclp14_15_force
    boundary = 'aclp_15'
  []
  [aclp_force_15_16]
    type = NodalSum
    variable = aclp15_16_force
    boundary = 'aclp_16'
  []
  [aclp_force_16_17]
    type = NodalSum
    variable = aclp16_17_force
    boundary = 'aclp_17'
  []
  [aclp_force_17_18]
    type = NodalSum
    variable = aclp17_18_force
    boundary = 'aclp_18'
  []
  [aclp_force_18_19]
    type = NodalSum
    variable = aclp18_19_force
    boundary = 'aclp_19'
  []


  [tlp_force_1_2]
    type = NodalSum
    variable = tlp1_2_force
    boundary = 'tlp_2'
  []
  [tlp_force_1_3]
    type = NodalSum
    variable = tlp1_3_force
    boundary = 'tlp_3'
  []
  [tlp_force_1_4]
    type = NodalSum
    variable = tlp1_4_force
    boundary = 'tlp_4'
  []
  [tlp_force_1_5]
    type = NodalSum
    variable = tlp1_5_force
    boundary = 'tlp_5'
  []
  [tlp_force_1_6]
    type = NodalSum
    variable = tlp1_6_force
    boundary = 'tlp_6'
  []
  [tlp_force_1_7]
    type = NodalSum
    variable = tlp1_7_force
    boundary = 'tlp_7'
  []
  [tlp_force_2_3]
    type = NodalSum
    variable = tlp2_3_force
    boundary = 'tlp_3'
  []
  [tlp_force_2_7]
    type = NodalSum
    variable = tlp2_7_force
    boundary = 'tlp_7'
  []
  [tlp_force_2_8]
    type = NodalSum
    variable = tlp2_8_force
    boundary = 'tlp_8'
  []
  [tlp_force_2_9]
    type = NodalSum
    variable = tlp2_9_force
    boundary = 'tlp_9'
  []
  [tlp_force_2_19]
    type = NodalSum
    variable = tlp2_19_force
    boundary = 'tlp_19'
  []
  [tlp_force_3_4]
    type = NodalSum
    variable = tlp3_4_force
    boundary = 'tlp_4'
  []
  [tlp_force_3_9]
    type = NodalSum
    variable = tlp3_9_force
    boundary = 'tlp_9'
  []
  [tlp_force_3_10]
    type = NodalSum
    variable = tlp3_10_force
    boundary = 'tlp_10'
  []
  [tlp_force_3_11]
    type = NodalSum
    variable = tlp3_11_force
    boundary = 'tlp_11'
  []
  [tlp_force_4_5]
    type = NodalSum
    variable = tlp4_5_force
    boundary = 'tlp_5'
  []
  [tlp_force_4_11]
    type = NodalSum
    variable = tlp4_11_force
    boundary = 'tlp_11'
  []
  [tlp_force_4_12]
    type = NodalSum
    variable = tlp4_12_force
    boundary = 'tlp_12'
  []
  [tlp_force_4_13]
    type = NodalSum
    variable = tlp4_13_force
    boundary = 'tlp_13'
  []
  [tlp_force_5_6]
    type = NodalSum
    variable = tlp5_6_force
    boundary = 'tlp_6'
  []
  [tlp_force_5_13]
    type = NodalSum
    variable = tlp5_13_force
    boundary = 'tlp_13'
  []
  [tlp_force_5_14]
    type = NodalSum
    variable = tlp5_14_force
    boundary = 'tlp_14'
  []
  [tlp_force_5_15]
    type = NodalSum
    variable = tlp5_15_force
    boundary = 'tlp_15'
  []
  [tlp_force_6_7]
    type = NodalSum
    variable = tlp6_7_force
    boundary = 'tlp_7'
  []
  [tlp_force_6_15]
    type = NodalSum
    variable = tlp6_15_force
    boundary = 'tlp_15'
  []
  [tlp_force_6_16]
    type = NodalSum
    variable = tlp6_16_force
    boundary = 'tlp_16'
  []
  [tlp_force_6_17]
    type = NodalSum
    variable = tlp6_17_force
    boundary = 'tlp_17'
  []
  [tlp_force_7_17]
    type = NodalSum
    variable = tlp7_17_force
    boundary = 'tlp_17'
  []
  [tlp_force_7_18]
    type = NodalSum
    variable = tlp7_18_force
    boundary = 'tlp_18'
  []
  [tlp_force_7_19]
    type = NodalSum
    variable = tlp7_19_force
    boundary = 'tlp_19'
  []
  [tlp_force_8_9]
    type = NodalSum
    variable = tlp8_9_force
    boundary = 'tlp_9'
  []
  [tlp_force_8_19]
    type = NodalSum
    variable = tlp8_19_force
    boundary = 'tlp_19'
  []
  [tlp_force_9_10]
    type = NodalSum
    variable = tlp9_10_force
    boundary = 'tlp_10'
  []
  [tlp_force_10_11]
    type = NodalSum
    variable = tlp10_11_force
    boundary = 'tlp_11'
  []
  [tlp_force_11_12]
    type = NodalSum
    variable = tlp11_12_force
    boundary = 'tlp_12'
  []
  [tlp_force_12_13]
    type = NodalSum
    variable = tlp12_13_force
    boundary = 'tlp_13'
  []
  [tlp_force_13_14]
    type = NodalSum
    variable = tlp13_14_force
    boundary = 'tlp_14'
  []
  [tlp_force_14_15]
    type = NodalSum
    variable = tlp14_15_force
    boundary = 'tlp_15'
  []
  [tlp_force_15_16]
    type = NodalSum
    variable = tlp15_16_force
    boundary = 'tlp_16'
  []
  [tlp_force_16_17]
    type = NodalSum
    variable = tlp16_17_force
    boundary = 'tlp_17'
  []
  [tlp_force_17_18]
    type = NodalSum
    variable = tlp17_18_force
    boundary = 'tlp_18'
  []
  [tlp_force_18_19]
    type = NodalSum
    variable = tlp18_19_force
    boundary = 'tlp_19'
  []
  
  [tlp_force_8_restraint]
    type = NodalSum
    variable = tlp8_restraint_force
    boundary = 31
  []
  [tlp_force_9_restraint]
    type = NodalSum
    variable = tlp9_restraint_force
    boundary = 31
  []
  [tlp_force_10_restraint]
    type = NodalSum
    variable = tlp10_restraint_force
    boundary = 31
  []
  [tlp_force_11_restraint]
    type = NodalSum
    variable = tlp11_restraint_force
    boundary = 31
  []
  [tlp_force_12_restraint]
    type = NodalSum
    variable = tlp12_restraint_force
    boundary = 31
  []
  [tlp_force_13_restraint]
    type = NodalSum
    variable = tlp13_restraint_force
    boundary = 31
  []
  [tlp_force_14_restraint]
    type = NodalSum
    variable = tlp14_restraint_force
    boundary = 31
  []
  [tlp_force_15_restraint]
    type = NodalSum
    variable = tlp15_restraint_force
    boundary = 31
  []
  [tlp_force_16_restraint]
    type = NodalSum
    variable = tlp16_restraint_force
    boundary = 31
  []
  [tlp_force_17_restraint]
    type = NodalSum
    variable = tlp17_restraint_force
    boundary = 31
  []
  [tlp_force_18_restraint]
    type = NodalSum
    variable = tlp18_restraint_force
    boundary = 31
  []
  [tlp_force_19_restraint]
    type = NodalSum
    variable = tlp19_restraint_force
    boundary = 31
  []
[]

[VectorPostprocessors]
  [duct1_disp_average]
    type = AverageSectionValueSampler
    axis_direction = '0 0 1'
    block = '1'
    variables = 'disp_x disp_y disp_z'
    reference_point = '0.0 0.0 0.0'
	require_equal_node_counts = false
  []
  [duct2_disp_average]
    type = AverageSectionValueSampler
    axis_direction = '0 0 1'
    block = '2'
    variables = 'disp_x disp_y disp_z'
    reference_point = '0.0 0.0 0.0'
	require_equal_node_counts = false
  []
  [duct3_disp_average]
    type = AverageSectionValueSampler
    axis_direction = '0 0 1'
    block = '3'
    variables = 'disp_x disp_y disp_z'
    reference_point = '0.0 0.0 0.0'
	require_equal_node_counts = false
  []
  [duct4_disp_average]
    type = AverageSectionValueSampler
    axis_direction = '0 0 1'
    block = '4'
    variables = 'disp_x disp_y disp_z'
    reference_point = '0.0 0.0 0.0'
	require_equal_node_counts = false
  []
  [duct5_disp_average]
    type = AverageSectionValueSampler
    axis_direction = '0 0 1'
    block = '5'
    variables = 'disp_x disp_y disp_z'
    reference_point = '0.0 0.0 0.0'
	require_equal_node_counts = false
  []
  [duct6_disp_average]
    type = AverageSectionValueSampler
    axis_direction = '0 0 1'
    block = '6'
    variables = 'disp_x disp_y disp_z'
    reference_point = '0.0 0.0 0.0'
	require_equal_node_counts = false
  []
  [duct7_disp_average]
    type = AverageSectionValueSampler
    axis_direction = '0 0 1'
    block = '7'
    variables = 'disp_x disp_y disp_z'
    reference_point = '0.0 0.0 0.0'
	require_equal_node_counts = false
  []
  [duct8_disp_average]
    type = AverageSectionValueSampler
    axis_direction = '0 0 1'
    block = '8'
    variables = 'disp_x disp_y disp_z'
    reference_point = '0.0 0.0 0.0'
	require_equal_node_counts = false
  []
  [duct9_disp_average]
    type = AverageSectionValueSampler
    axis_direction = '0 0 1'
    block = '9'
    variables = 'disp_x disp_y disp_z'
    reference_point = '0.0 0.0 0.0'
	require_equal_node_counts = false
  []
  [duct10_disp_average]
    type = AverageSectionValueSampler
    axis_direction = '0 0 1'
    block = '10'
    variables = 'disp_x disp_y disp_z'
    reference_point = '0.0 0.0 0.0'
	require_equal_node_counts = false
  []
  [duct11_disp_average]
    type = AverageSectionValueSampler
    axis_direction = '0 0 1'
    block = '11'
    variables = 'disp_x disp_y disp_z'
    reference_point = '0.0 0.0 0.0'
	require_equal_node_counts = false
  []
  [duct12_disp_average]
    type = AverageSectionValueSampler
    axis_direction = '0 0 1'
    block = '12'
    variables = 'disp_x disp_y disp_z'
    reference_point = '0.0 0.0 0.0'
	require_equal_node_counts = false
  []
  [duct13_disp_average]
    type = AverageSectionValueSampler
    axis_direction = '0 0 1'
    block = '13'
    variables = 'disp_x disp_y disp_z'
    reference_point = '0.0 0.0 0.0'
	require_equal_node_counts = false
  []
  [duct14_disp_average]
    type = AverageSectionValueSampler
    axis_direction = '0 0 1'
    block = '14'
    variables = 'disp_x disp_y disp_z'
    reference_point = '0.0 0.0 0.0'
	require_equal_node_counts = false
  []
  [duct15_disp_average]
    type = AverageSectionValueSampler
    axis_direction = '0 0 1'
    block = '15'
    variables = 'disp_x disp_y disp_z'
    reference_point = '0.0 0.0 0.0'
	require_equal_node_counts = false
  []
  [duct16_disp_average]
    type = AverageSectionValueSampler
    axis_direction = '0 0 1'
    block = '16'
    variables = 'disp_x disp_y disp_z'
    reference_point = '0.0 0.0 0.0'
	require_equal_node_counts = false
  []
  [duct17_disp_average]
    type = AverageSectionValueSampler
    axis_direction = '0 0 1'
    block = '17'
    variables = 'disp_x disp_y disp_z'
    reference_point = '0.0 0.0 0.0'
	require_equal_node_counts = false
  []
  [duct18_disp_average]
    type = AverageSectionValueSampler
    axis_direction = '0 0 1'
    block = '18'
    variables = 'disp_x disp_y disp_z'
    reference_point = '0.0 0.0 0.0'
	require_equal_node_counts = false
  []
  [duct19_disp_average]
    type = AverageSectionValueSampler
    axis_direction = '0 0 1'
    block = '19'
    variables = 'disp_x disp_y disp_z'
    reference_point = '0.0 0.0 0.0'
	require_equal_node_counts = false
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
  dt = 0.1
  end_time = 1.0

  dtmax = 1
  dtmin = 0.01

  [Predictor]
    type = SimplePredictor
    scale = 1.0
  []
[]

[Outputs]
  exodus = true
  perf_graph = true
  [test_average_section]
    type = CSV
    file_base = average_section_disp
    execute_on = final
    show = 'duct1_disp_average  duct2_disp_average  duct3_disp_average  duct4_disp_average  duct5_disp_average ' 
	       'duct6_disp_average  duct7_disp_average  duct8_disp_average  duct9_disp_average  duct10_disp_average '
		   'duct11_disp_average duct12_disp_average duct13_disp_average duct14_disp_average duct15_disp_average '
		   'duct16_disp_average duct17_disp_average duct18_disp_average duct19_disp_average'
  []
  [aclp_force_plots]
    type = CSV
    file_base = aclp_force_plots
    execute_on = timestep_end
    show = 'aclp_force_1_2   aclp_force_1_3   aclp_force_1_4   aclp_force_1_5   aclp_force_1_6   aclp_force_1_7 '
	       'aclp_force_2_3   aclp_force_2_7   aclp_force_2_8   aclp_force_2_9   aclp_force_2_19  aclp_force_3_4 '
		   'aclp_force_3_9   aclp_force_3_10  aclp_force_3_11  aclp_force_4_5   aclp_force_4_11  aclp_force_4_12 '
		   'aclp_force_4_13  aclp_force_5_6   aclp_force_5_13  aclp_force_5_14  aclp_force_5_15  aclp_force_6_7 '
		   'aclp_force_6_15  aclp_force_6_16  aclp_force_6_17  aclp_force_7_17  aclp_force_7_18  aclp_force_7_19 '
		   'aclp_force_8_9   aclp_force_8_19  aclp_force_9_10  aclp_force_10_11 aclp_force_11_12 aclp_force_12_13 '
		   'aclp_force_13_14 aclp_force_14_15 aclp_force_15_16 aclp_force_16_17 aclp_force_17_18 aclp_force_18_19'
  []
  [tlp_force_plots]
    type = CSV
    file_base = tlp_force_plots
    execute_on = timestep_end
    show = 'tlp_force_1_2   tlp_force_1_3   tlp_force_1_4   tlp_force_1_5   tlp_force_1_6   tlp_force_1_7 '
	       'tlp_force_2_3   tlp_force_2_7   tlp_force_2_8   tlp_force_2_9   tlp_force_2_19  tlp_force_3_4 '
		   'tlp_force_3_9   tlp_force_3_10  tlp_force_3_11  tlp_force_4_5   tlp_force_4_11  tlp_force_4_12 '
		   'tlp_force_4_13  tlp_force_5_6   tlp_force_5_13  tlp_force_5_14  tlp_force_5_15  tlp_force_6_7 '
		   'tlp_force_6_15  tlp_force_6_16  tlp_force_6_17  tlp_force_7_17  tlp_force_7_18  tlp_force_7_19 '
		   'tlp_force_8_9   tlp_force_8_19  tlp_force_9_10  tlp_force_10_11 tlp_force_11_12 tlp_force_12_13 '
		   'tlp_force_13_14 tlp_force_14_15 tlp_force_15_16 tlp_force_16_17 tlp_force_17_18 tlp_force_18_19 '
		   'tlp_force_8_restraint  tlp_force_9_restraint  tlp_force_10_restraint tlp_force_11_restraint '
		   'tlp_force_12_restraint tlp_force_13_restraint tlp_force_14_restraint tlp_force_15_restraint '
		   'tlp_force_16_restraint tlp_force_17_restraint tlp_force_18_restraint tlp_force_19_restraint'
  []
[]



