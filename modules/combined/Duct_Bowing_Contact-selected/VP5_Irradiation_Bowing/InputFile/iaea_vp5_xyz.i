E = 1.7e11
D0 = 100.0
a = 0.697
L = 1.0

[GlobalParams]
  order = FIRST
  family = LAGRANGE
  displacements = 'disp_x disp_y disp_z'
  volumetric_locking_correction = true
[]

[Mesh]
## File mesh generator
  [fmg]
    type = FileMeshGenerator
    file = './../Mesh/vp5_mesh.e'
  []
  patch_update_strategy = iteration
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
  [rad_ducts]
    strain = FINITE
    volumetric_locking_correction = true
    eigenstrain_names = 'swelling_strain'
    decomposition_method = EigenSolution
    generate_output = 'vonmises_stress'
    temperature = temp
    use_finite_deform_jacobian = true
    extra_vector_tags = 'ref'
	block = '42 43 44'
  []
  [non_rad]
    strain = FINITE
    volumetric_locking_correction = true
    decomposition_method = EigenSolution
    generate_output = 'vonmises_stress'
    temperature = temp
    use_finite_deform_jacobian = true
    extra_vector_tags = 'ref'
	block = '1 3 10 11 23 24 67 68 69 98 99 100 101 200 201 300 1000'
  []
[]

[Kernels]
  [value]
    type = MaterialPropertyValue
	prop_name = dpa
	variable = dose
  []
[]

[Functions]
  [damage_dose]
    type = ParsedFunction
	expression = 'if(z<1.5,0,if(z<2.5,if(r<a,t*D0*(1.0-r*r/(a*a))*sin(pi*(z-1.5)/L),0),0))'
	symbol_names = 'r D0 a L'
	symbol_values = 'radius ${D0} ${a} ${L}'
  []
  [radius]
    type = ParsedFunction
	expression = 'sqrt(x*x + y*y)'
  []
[]

	[AuxVariables]
		[temp]
			initial_condition = 400
		[]
		[r]
		[]
		[aclp1_3_contact_x]
		[]
		[aclp1_3_contact_y]
		[]
		[aclp1_3_contact_z]
		[]
		[aclp3_10_contact_x]
		[]
		[aclp3_10_contact_y]
		[]
		[aclp3_10_contact_z]
		[]
		[aclp3_11_contact_x]
		[]
		[aclp3_11_contact_y]
		[]
		[aclp3_11_contact_z]
		[]
		[aclp10_11_contact_x]
		[]
		[aclp10_11_contact_y]
		[]
		[aclp10_11_contact_z]
		[]
		[aclp10_23_contact_x]
		[]
		[aclp10_23_contact_y]
		[]
		[aclp10_23_contact_z]
		[]
		[aclp10_24_contact_x]
		[]
		[aclp10_24_contact_y]
		[]
		[aclp10_24_contact_z]
		[]
		[aclp11_24_contact_x]
		[]
		[aclp11_24_contact_y]
		[]
		[aclp11_24_contact_z]
		[]
		[aclp23_24_contact_x]
		[]
		[aclp23_24_contact_y]
		[]
		[aclp23_24_contact_z]
		[]
		[aclp23_42_contact_x]
		[]
		[aclp23_42_contact_y]
		[]
		[aclp23_42_contact_z]
		[]
		[aclp23_43_contact_x]
		[]
		[aclp23_43_contact_y]
		[]
		[aclp23_43_contact_z]
		[]
		[aclp24_43_contact_x]
		[]
		[aclp24_43_contact_y]
		[]
		[aclp24_43_contact_z]
		[]
		[aclp24_44_contact_x]
		[]
		[aclp24_44_contact_y]
		[]
		[aclp24_44_contact_z]
		[]
		[aclp42_43_contact_x]
		[]
		[aclp42_43_contact_y]
		[]
		[aclp42_43_contact_z]
		[]
		[aclp42_67_contact_x]
		[]
		[aclp42_67_contact_y]
		[]
		[aclp42_67_contact_z]
		[]
		[aclp42_68_contact_x]
		[]
		[aclp42_68_contact_y]
		[]
		[aclp42_68_contact_z]
		[]
		[aclp43_44_contact_x]
		[]
		[aclp43_44_contact_y]
		[]
		[aclp43_44_contact_z]
		[]
		[aclp43_68_contact_x]
		[]
		[aclp43_68_contact_y]
		[]
		[aclp43_68_contact_z]
		[]
		[aclp43_69_contact_x]
		[]
		[aclp43_69_contact_y]
		[]
		[aclp43_69_contact_z]
		[]
		[aclp44_69_contact_x]
		[]
		[aclp44_69_contact_y]
		[]
		[aclp44_69_contact_z]
		[]
		[aclp67_68_contact_x]
		[]
		[aclp67_68_contact_y]
		[]
		[aclp67_68_contact_z]
		[]
		[aclp67_98_contact_x]
		[]
		[aclp67_98_contact_y]
		[]
		[aclp67_98_contact_z]
		[]
		[aclp67_99_contact_x]
		[]
		[aclp67_99_contact_y]
		[]
		[aclp67_99_contact_z]
		[]
		[aclp68_69_contact_x]
		[]
		[aclp68_69_contact_y]
		[]
		[aclp68_69_contact_z]
		[]
		[aclp68_99_contact_x]
		[]
		[aclp68_99_contact_y]
		[]
		[aclp68_99_contact_z]
		[]
		[aclp68_100_contact_x]
		[]
		[aclp68_100_contact_y]
		[]
		[aclp68_100_contact_z]
		[]
		[aclp69_100_contact_x]
		[]
		[aclp69_100_contact_y]
		[]
		[aclp69_100_contact_z]
		[]
		[aclp69_101_contact_x]
		[]
		[aclp69_101_contact_y]
		[]
		[aclp69_101_contact_z]
		[]
		[aclp98_99_contact_x]
		[]
		[aclp98_99_contact_y]
		[]
		[aclp98_99_contact_z]
		[]
		[aclp99_100_contact_x]
		[]
		[aclp99_100_contact_y]
		[]
		[aclp99_100_contact_z]
		[]
		[aclp100_101_contact_x]
		[]
		[aclp100_101_contact_y]
		[]
		[aclp100_101_contact_z]
		[]
		[tlp1_3_contact_x]
		[]
		[tlp1_3_contact_y]
		[]
		[tlp1_3_contact_z]
		[]
		[tlp3_10_contact_x]
		[]
		[tlp3_10_contact_y]
		[]
		[tlp3_10_contact_z]
		[]
		[tlp3_11_contact_x]
		[]
		[tlp3_11_contact_y]
		[]
		[tlp3_11_contact_z]
		[]
		[tlp10_11_contact_x]
		[]
		[tlp10_11_contact_y]
		[]
		[tlp10_11_contact_z]
		[]
		[tlp10_23_contact_x]
		[]
		[tlp10_23_contact_y]
		[]
		[tlp10_23_contact_z]
		[]
		[tlp10_24_contact_x]
		[]
		[tlp10_24_contact_y]
		[]
		[tlp10_24_contact_z]
		[]
		[tlp11_24_contact_x]
		[]
		[tlp11_24_contact_y]
		[]
		[tlp11_24_contact_z]
		[]
		[tlp23_24_contact_x]
		[]
		[tlp23_24_contact_y]
		[]
		[tlp23_24_contact_z]
		[]
		[tlp23_42_contact_x]
		[]
		[tlp23_42_contact_y]
		[]
		[tlp23_42_contact_z]
		[]
		[tlp23_43_contact_x]
		[]
		[tlp23_43_contact_y]
		[]
		[tlp23_43_contact_z]
		[]
		[tlp24_43_contact_x]
		[]
		[tlp24_43_contact_y]
		[]
		[tlp24_43_contact_z]
		[]
		[tlp24_44_contact_x]
		[]
		[tlp24_44_contact_y]
		[]
		[tlp24_44_contact_z]
		[]
		[tlp42_43_contact_x]
		[]
		[tlp42_43_contact_y]
		[]
		[tlp42_43_contact_z]
		[]
		[tlp42_67_contact_x]
		[]
		[tlp42_67_contact_y]
		[]
		[tlp42_67_contact_z]
		[]
		[tlp42_68_contact_x]
		[]
		[tlp42_68_contact_y]
		[]
		[tlp42_68_contact_z]
		[]
		[tlp43_44_contact_x]
		[]
		[tlp43_44_contact_y]
		[]
		[tlp43_44_contact_z]
		[]
		[tlp43_68_contact_x]
		[]
		[tlp43_68_contact_y]
		[]
		[tlp43_68_contact_z]
		[]
		[tlp43_69_contact_x]
		[]
		[tlp43_69_contact_y]
		[]
		[tlp43_69_contact_z]
		[]
		[tlp44_69_contact_x]
		[]
		[tlp44_69_contact_y]
		[]
		[tlp44_69_contact_z]
		[]
		[tlp67_68_contact_x]
		[]
		[tlp67_68_contact_y]
		[]
		[tlp67_68_contact_z]
		[]
		[tlp67_98_contact_x]
		[]
		[tlp67_98_contact_y]
		[]
		[tlp67_98_contact_z]
		[]
		[tlp67_99_contact_x]
		[]
		[tlp67_99_contact_y]
		[]
		[tlp67_99_contact_z]
		[]
		[tlp68_69_contact_x]
		[]
		[tlp68_69_contact_y]
		[]
		[tlp68_69_contact_z]
		[]
		[tlp68_99_contact_x]
		[]
		[tlp68_99_contact_y]
		[]
		[tlp68_99_contact_z]
		[]
		[tlp68_100_contact_x]
		[]
		[tlp68_100_contact_y]
		[]
		[tlp68_100_contact_z]
		[]
		[tlp69_100_contact_x]
		[]
		[tlp69_100_contact_y]
		[]
		[tlp69_100_contact_z]
		[]
		[tlp69_101_contact_x]
		[]
		[tlp69_101_contact_y]
		[]
		[tlp69_101_contact_z]
		[]
		[tlp98_99_contact_x]
		[]
		[tlp98_99_contact_y]
		[]
		[tlp98_99_contact_z]
		[]
		[tlp99_100_contact_x]
		[]
		[tlp99_100_contact_y]
		[]
		[tlp99_100_contact_z]
		[]
		[tlp100_101_contact_x]
		[]
		[tlp100_101_contact_y]
		[]
		[tlp100_101_contact_z]
		[]
		[block1_tlp_contact_x]
		[]
		[block1_tlp_contact_y]
		[]
		[block1_tlp_contact_z]
		[]
		[block2_tlp_contact_x]
		[]
		[block2_tlp_contact_y]
		[]
		[block2_tlp_contact_z]
		[]
		[block3_tlp_contact_x]
		[]
		[block3_tlp_contact_y]
		[]
		[block3_tlp_contact_z]
		[]
		[block1_aclp_contact_x]
		[]
		[block1_aclp_contact_y]
		[]
		[block1_aclp_contact_z]
		[]
		[block2_aclp_contact_x]
		[]
		[block2_aclp_contact_y]
		[]
		[block2_aclp_contact_z]
		[]
		[block3_aclp_contact_x]
		[]
		[block3_aclp_contact_y]
		[]
		[block3_aclp_contact_z]
		[]
		[tlp101_restraint_contact_x]
		[]
		[tlp101_restraint_contact_y]
		[]
		[tlp101_restraint_contact_z]
		[]
		[aclp101_restraint_contact_x]
		[]
		[aclp101_restraint_contact_y]
		[]
		[aclp101_restraint_contact_z]
		[]
	[]


	[AuxKernels]
		[tfunc]
			type = ConstantAux
			value = 400
			variable = temp
		[]
		[radius]
			type = FunctionAux
			function = radius
			variable = r
		[]
		[aclp1_3_x]
			type = PenetrationAux
			variable = aclp1_3_contact_x
			boundary = 'ACLP1_3'
			paired_boundary = 'ACLP3_1'
			quantity = normal_force_x
		[]
		[aclp1_3_y]
			type = PenetrationAux
			variable = aclp1_3_contact_y
			boundary = 'ACLP1_3'
			paired_boundary = 'ACLP3_1'
			quantity = normal_force_y
		[]
		[aclp1_3_z]
			type = PenetrationAux
			variable = aclp1_3_contact_z
			boundary = 'ACLP1_3'
			paired_boundary = 'ACLP3_1'
			quantity = normal_force_z
		[]
		[aclp3_10_x]
			type = PenetrationAux
			variable = aclp3_10_contact_x
			boundary = 'ACLP3_10'
			paired_boundary = 'ACLP10_3'
			quantity = normal_force_x
		[]
		[aclp3_10_y]
			type = PenetrationAux
			variable = aclp3_10_contact_y
			boundary = 'ACLP3_10'
			paired_boundary = 'ACLP10_3'
			quantity = normal_force_y
		[]
		[aclp3_10_z]
			type = PenetrationAux
			variable = aclp3_10_contact_z
			boundary = 'ACLP3_10'
			paired_boundary = 'ACLP10_3'
			quantity = normal_force_z
		[]
		[aclp3_11_x]
			type = PenetrationAux
			variable = aclp3_11_contact_x
			boundary = 'ACLP3_11'
			paired_boundary = 'ACLP11_3'
			quantity = normal_force_x
		[]
		[aclp3_11_y]
			type = PenetrationAux
			variable = aclp3_11_contact_y
			boundary = 'ACLP3_11'
			paired_boundary = 'ACLP11_3'
			quantity = normal_force_y
		[]
		[aclp3_11_z]
			type = PenetrationAux
			variable = aclp3_11_contact_z
			boundary = 'ACLP3_11'
			paired_boundary = 'ACLP11_3'
			quantity = normal_force_z
		[]
		[aclp10_11_x]
			type = PenetrationAux
			variable = aclp10_11_contact_x
			boundary = 'ACLP10_11'
			paired_boundary = 'ACLP11_10'
			quantity = normal_force_x
		[]
		[aclp10_11_y]
			type = PenetrationAux
			variable = aclp10_11_contact_y
			boundary = 'ACLP10_11'
			paired_boundary = 'ACLP11_10'
			quantity = normal_force_y
		[]
		[aclp10_11_z]
			type = PenetrationAux
			variable = aclp10_11_contact_z
			boundary = 'ACLP10_11'
			paired_boundary = 'ACLP11_10'
			quantity = normal_force_z
		[]
		[aclp10_23_x]
			type = PenetrationAux
			variable = aclp10_23_contact_x
			boundary = 'ACLP10_23'
			paired_boundary = 'ACLP23_10'
			quantity = normal_force_x
		[]
		[aclp10_23_y]
			type = PenetrationAux
			variable = aclp10_23_contact_y
			boundary = 'ACLP10_23'
			paired_boundary = 'ACLP23_10'
			quantity = normal_force_y
		[]
		[aclp10_23_z]
			type = PenetrationAux
			variable = aclp10_23_contact_z
			boundary = 'ACLP10_23'
			paired_boundary = 'ACLP23_10'
			quantity = normal_force_z
		[]
		[aclp10_24_x]
			type = PenetrationAux
			variable = aclp10_24_contact_x
			boundary = 'ACLP10_24'
			paired_boundary = 'ACLP24_10'
			quantity = normal_force_x
		[]
		[aclp10_24_y]
			type = PenetrationAux
			variable = aclp10_24_contact_y
			boundary = 'ACLP10_24'
			paired_boundary = 'ACLP24_10'
			quantity = normal_force_y
		[]
		[aclp10_24_z]
			type = PenetrationAux
			variable = aclp10_24_contact_z
			boundary = 'ACLP10_24'
			paired_boundary = 'ACLP24_10'
			quantity = normal_force_z
		[]
		[aclp11_24_x]
			type = PenetrationAux
			variable = aclp11_24_contact_x
			boundary = 'ACLP11_24'
			paired_boundary = 'ACLP24_11'
			quantity = normal_force_x
		[]
		[aclp11_24_y]
			type = PenetrationAux
			variable = aclp11_24_contact_y
			boundary = 'ACLP11_24'
			paired_boundary = 'ACLP24_11'
			quantity = normal_force_y
		[]
		[aclp11_24_z]
			type = PenetrationAux
			variable = aclp11_24_contact_z
			boundary = 'ACLP11_24'
			paired_boundary = 'ACLP24_11'
			quantity = normal_force_z
		[]
		[aclp23_24_x]
			type = PenetrationAux
			variable = aclp23_24_contact_x
			boundary = 'ACLP23_24'
			paired_boundary = 'ACLP24_23'
			quantity = normal_force_x
		[]
		[aclp23_24_y]
			type = PenetrationAux
			variable = aclp23_24_contact_y
			boundary = 'ACLP23_24'
			paired_boundary = 'ACLP24_23'
			quantity = normal_force_y
		[]
		[aclp23_24_z]
			type = PenetrationAux
			variable = aclp23_24_contact_z
			boundary = 'ACLP23_24'
			paired_boundary = 'ACLP24_23'
			quantity = normal_force_z
		[]
		[aclp23_42_x]
			type = PenetrationAux
			variable = aclp23_42_contact_x
			boundary = 'ACLP23_42'
			paired_boundary = 'ACLP42_23'
			quantity = normal_force_x
		[]
		[aclp23_42_y]
			type = PenetrationAux
			variable = aclp23_42_contact_y
			boundary = 'ACLP23_42'
			paired_boundary = 'ACLP42_23'
			quantity = normal_force_y
		[]
		[aclp23_42_z]
			type = PenetrationAux
			variable = aclp23_42_contact_z
			boundary = 'ACLP23_42'
			paired_boundary = 'ACLP42_23'
			quantity = normal_force_z
		[]
		[aclp23_43_x]
			type = PenetrationAux
			variable = aclp23_43_contact_x
			boundary = 'ACLP23_43'
			paired_boundary = 'ACLP43_23'
			quantity = normal_force_x
		[]
		[aclp23_43_y]
			type = PenetrationAux
			variable = aclp23_43_contact_y
			boundary = 'ACLP23_43'
			paired_boundary = 'ACLP43_23'
			quantity = normal_force_y
		[]
		[aclp23_43_z]
			type = PenetrationAux
			variable = aclp23_43_contact_z
			boundary = 'ACLP23_43'
			paired_boundary = 'ACLP43_23'
			quantity = normal_force_z
		[]
		[aclp24_43_x]
			type = PenetrationAux
			variable = aclp24_43_contact_x
			boundary = 'ACLP24_43'
			paired_boundary = 'ACLP43_24'
			quantity = normal_force_x
		[]
		[aclp24_43_y]
			type = PenetrationAux
			variable = aclp24_43_contact_y
			boundary = 'ACLP24_43'
			paired_boundary = 'ACLP43_24'
			quantity = normal_force_y
		[]
		[aclp24_43_z]
			type = PenetrationAux
			variable = aclp24_43_contact_z
			boundary = 'ACLP24_43'
			paired_boundary = 'ACLP43_24'
			quantity = normal_force_z
		[]
		[aclp24_44_x]
			type = PenetrationAux
			variable = aclp24_44_contact_x
			boundary = 'ACLP24_44'
			paired_boundary = 'ACLP44_24'
			quantity = normal_force_x
		[]
		[aclp24_44_y]
			type = PenetrationAux
			variable = aclp24_44_contact_y
			boundary = 'ACLP24_44'
			paired_boundary = 'ACLP44_24'
			quantity = normal_force_y
		[]
		[aclp24_44_z]
			type = PenetrationAux
			variable = aclp24_44_contact_z
			boundary = 'ACLP24_44'
			paired_boundary = 'ACLP44_24'
			quantity = normal_force_z
		[]
		[aclp42_43_x]
			type = PenetrationAux
			variable = aclp42_43_contact_x
			boundary = 'ACLP42_43'
			paired_boundary = 'ACLP43_42'
			quantity = normal_force_x
		[]
		[aclp42_43_y]
			type = PenetrationAux
			variable = aclp42_43_contact_y
			boundary = 'ACLP42_43'
			paired_boundary = 'ACLP43_42'
			quantity = normal_force_y
		[]
		[aclp42_43_z]
			type = PenetrationAux
			variable = aclp42_43_contact_z
			boundary = 'ACLP42_43'
			paired_boundary = 'ACLP43_42'
			quantity = normal_force_z
		[]
		[aclp42_67_x]
			type = PenetrationAux
			variable = aclp42_67_contact_x
			boundary = 'ACLP42_67'
			paired_boundary = 'ACLP67_42'
			quantity = normal_force_x
		[]
		[aclp42_67_y]
			type = PenetrationAux
			variable = aclp42_67_contact_y
			boundary = 'ACLP42_67'
			paired_boundary = 'ACLP67_42'
			quantity = normal_force_y
		[]
		[aclp42_67_z]
			type = PenetrationAux
			variable = aclp42_67_contact_z
			boundary = 'ACLP42_67'
			paired_boundary = 'ACLP67_42'
			quantity = normal_force_z
		[]
		[aclp42_68_x]
			type = PenetrationAux
			variable = aclp42_68_contact_x
			boundary = 'ACLP42_68'
			paired_boundary = 'ACLP68_42'
			quantity = normal_force_x
		[]
		[aclp42_68_y]
			type = PenetrationAux
			variable = aclp42_68_contact_y
			boundary = 'ACLP42_68'
			paired_boundary = 'ACLP68_42'
			quantity = normal_force_y
		[]
		[aclp42_68_z]
			type = PenetrationAux
			variable = aclp42_68_contact_z
			boundary = 'ACLP42_68'
			paired_boundary = 'ACLP68_42'
			quantity = normal_force_z
		[]
		[aclp43_44_x]
			type = PenetrationAux
			variable = aclp43_44_contact_x
			boundary = 'ACLP43_44'
			paired_boundary = 'ACLP44_43'
			quantity = normal_force_x
		[]
		[aclp43_44_y]
			type = PenetrationAux
			variable = aclp43_44_contact_y
			boundary = 'ACLP43_44'
			paired_boundary = 'ACLP44_43'
			quantity = normal_force_y
		[]
		[aclp43_44_z]
			type = PenetrationAux
			variable = aclp43_44_contact_z
			boundary = 'ACLP43_44'
			paired_boundary = 'ACLP44_43'
			quantity = normal_force_z
		[]
		[aclp43_68_x]
			type = PenetrationAux
			variable = aclp43_68_contact_x
			boundary = 'ACLP43_68'
			paired_boundary = 'ACLP68_43'
			quantity = normal_force_x
		[]
		[aclp43_68_y]
			type = PenetrationAux
			variable = aclp43_68_contact_y
			boundary = 'ACLP43_68'
			paired_boundary = 'ACLP68_43'
			quantity = normal_force_y
		[]
		[aclp43_68_z]
			type = PenetrationAux
			variable = aclp43_68_contact_z
			boundary = 'ACLP43_68'
			paired_boundary = 'ACLP68_43'
			quantity = normal_force_z
		[]
		[aclp43_69_x]
			type = PenetrationAux
			variable = aclp43_69_contact_x
			boundary = 'ACLP43_69'
			paired_boundary = 'ACLP69_43'
			quantity = normal_force_x
		[]
		[aclp43_69_y]
			type = PenetrationAux
			variable = aclp43_69_contact_y
			boundary = 'ACLP43_69'
			paired_boundary = 'ACLP69_43'
			quantity = normal_force_y
		[]
		[aclp43_69_z]
			type = PenetrationAux
			variable = aclp43_69_contact_z
			boundary = 'ACLP43_69'
			paired_boundary = 'ACLP69_43'
			quantity = normal_force_z
		[]
		[aclp44_69_x]
			type = PenetrationAux
			variable = aclp44_69_contact_x
			boundary = 'ACLP44_69'
			paired_boundary = 'ACLP69_44'
			quantity = normal_force_x
		[]
		[aclp44_69_y]
			type = PenetrationAux
			variable = aclp44_69_contact_y
			boundary = 'ACLP44_69'
			paired_boundary = 'ACLP69_44'
			quantity = normal_force_y
		[]
		[aclp44_69_z]
			type = PenetrationAux
			variable = aclp44_69_contact_z
			boundary = 'ACLP44_69'
			paired_boundary = 'ACLP69_44'
			quantity = normal_force_z
		[]
		[aclp67_68_x]
			type = PenetrationAux
			variable = aclp67_68_contact_x
			boundary = 'ACLP67_68'
			paired_boundary = 'ACLP68_67'
			quantity = normal_force_x
		[]
		[aclp67_68_y]
			type = PenetrationAux
			variable = aclp67_68_contact_y
			boundary = 'ACLP67_68'
			paired_boundary = 'ACLP68_67'
			quantity = normal_force_y
		[]
		[aclp67_68_z]
			type = PenetrationAux
			variable = aclp67_68_contact_z
			boundary = 'ACLP67_68'
			paired_boundary = 'ACLP68_67'
			quantity = normal_force_z
		[]
		[aclp67_98_x]
			type = PenetrationAux
			variable = aclp67_98_contact_x
			boundary = 'ACLP67_98'
			paired_boundary = 'ACLP98_67'
			quantity = normal_force_x
		[]
		[aclp67_98_y]
			type = PenetrationAux
			variable = aclp67_98_contact_y
			boundary = 'ACLP67_98'
			paired_boundary = 'ACLP98_67'
			quantity = normal_force_y
		[]
		[aclp67_98_z]
			type = PenetrationAux
			variable = aclp67_98_contact_z
			boundary = 'ACLP67_98'
			paired_boundary = 'ACLP98_67'
			quantity = normal_force_z
		[]
		[aclp67_99_x]
			type = PenetrationAux
			variable = aclp67_99_contact_x
			boundary = 'ACLP67_99'
			paired_boundary = 'ACLP99_67'
			quantity = normal_force_x
		[]
		[aclp67_99_y]
			type = PenetrationAux
			variable = aclp67_99_contact_y
			boundary = 'ACLP67_99'
			paired_boundary = 'ACLP99_67'
			quantity = normal_force_y
		[]
		[aclp67_99_z]
			type = PenetrationAux
			variable = aclp67_99_contact_z
			boundary = 'ACLP67_99'
			paired_boundary = 'ACLP99_67'
			quantity = normal_force_z
		[]
		[aclp68_69_x]
			type = PenetrationAux
			variable = aclp68_69_contact_x
			boundary = 'ACLP68_69'
			paired_boundary = 'ACLP69_68'
			quantity = normal_force_x
		[]
		[aclp68_69_y]
			type = PenetrationAux
			variable = aclp68_69_contact_y
			boundary = 'ACLP68_69'
			paired_boundary = 'ACLP69_68'
			quantity = normal_force_y
		[]
		[aclp68_69_z]
			type = PenetrationAux
			variable = aclp68_69_contact_z
			boundary = 'ACLP68_69'
			paired_boundary = 'ACLP69_68'
			quantity = normal_force_z
		[]
		[aclp68_99_x]
			type = PenetrationAux
			variable = aclp68_99_contact_x
			boundary = 'ACLP68_99'
			paired_boundary = 'ACLP99_68'
			quantity = normal_force_x
		[]
		[aclp68_99_y]
			type = PenetrationAux
			variable = aclp68_99_contact_y
			boundary = 'ACLP68_99'
			paired_boundary = 'ACLP99_68'
			quantity = normal_force_y
		[]
		[aclp68_99_z]
			type = PenetrationAux
			variable = aclp68_99_contact_z
			boundary = 'ACLP68_99'
			paired_boundary = 'ACLP99_68'
			quantity = normal_force_z
		[]
		[aclp68_100_x]
			type = PenetrationAux
			variable = aclp68_100_contact_x
			boundary = 'ACLP68_100'
			paired_boundary = 'ACLP100_68'
			quantity = normal_force_x
		[]
		[aclp68_100_y]
			type = PenetrationAux
			variable = aclp68_100_contact_y
			boundary = 'ACLP68_100'
			paired_boundary = 'ACLP100_68'
			quantity = normal_force_y
		[]
		[aclp68_100_z]
			type = PenetrationAux
			variable = aclp68_100_contact_z
			boundary = 'ACLP68_100'
			paired_boundary = 'ACLP100_68'
			quantity = normal_force_z
		[]
		[aclp69_100_x]
			type = PenetrationAux
			variable = aclp69_100_contact_x
			boundary = 'ACLP69_100'
			paired_boundary = 'ACLP100_69'
			quantity = normal_force_x
		[]
		[aclp69_100_y]
			type = PenetrationAux
			variable = aclp69_100_contact_y
			boundary = 'ACLP69_100'
			paired_boundary = 'ACLP100_69'
			quantity = normal_force_y
		[]
		[aclp69_100_z]
			type = PenetrationAux
			variable = aclp69_100_contact_z
			boundary = 'ACLP69_100'
			paired_boundary = 'ACLP100_69'
			quantity = normal_force_z
		[]
		[aclp69_101_x]
			type = PenetrationAux
			variable = aclp69_101_contact_x
			boundary = 'ACLP69_101'
			paired_boundary = 'ACLP101_69'
			quantity = normal_force_x
		[]
		[aclp69_101_y]
			type = PenetrationAux
			variable = aclp69_101_contact_y
			boundary = 'ACLP69_101'
			paired_boundary = 'ACLP101_69'
			quantity = normal_force_y
		[]
		[aclp69_101_z]
			type = PenetrationAux
			variable = aclp69_101_contact_z
			boundary = 'ACLP69_101'
			paired_boundary = 'ACLP101_69'
			quantity = normal_force_z
		[]
		[aclp98_99_x]
			type = PenetrationAux
			variable = aclp98_99_contact_x
			boundary = 'ACLP98_99'
			paired_boundary = 'ACLP99_98'
			quantity = normal_force_x
		[]
		[aclp98_99_y]
			type = PenetrationAux
			variable = aclp98_99_contact_y
			boundary = 'ACLP98_99'
			paired_boundary = 'ACLP99_98'
			quantity = normal_force_y
		[]
		[aclp98_99_z]
			type = PenetrationAux
			variable = aclp98_99_contact_z
			boundary = 'ACLP98_99'
			paired_boundary = 'ACLP99_98'
			quantity = normal_force_z
		[]
		[aclp99_100_x]
			type = PenetrationAux
			variable = aclp99_100_contact_x
			boundary = 'ACLP99_100'
			paired_boundary = 'ACLP100_99'
			quantity = normal_force_x
		[]
		[aclp99_100_y]
			type = PenetrationAux
			variable = aclp99_100_contact_y
			boundary = 'ACLP99_100'
			paired_boundary = 'ACLP100_99'
			quantity = normal_force_y
		[]
		[aclp99_100_z]
			type = PenetrationAux
			variable = aclp99_100_contact_z
			boundary = 'ACLP99_100'
			paired_boundary = 'ACLP100_99'
			quantity = normal_force_z
		[]
		[aclp100_101_x]
			type = PenetrationAux
			variable = aclp100_101_contact_x
			boundary = 'ACLP100_101'
			paired_boundary = 'ACLP101_100'
			quantity = normal_force_x
		[]
		[aclp100_101_y]
			type = PenetrationAux
			variable = aclp100_101_contact_y
			boundary = 'ACLP100_101'
			paired_boundary = 'ACLP101_100'
			quantity = normal_force_y
		[]
		[aclp100_101_z]
			type = PenetrationAux
			variable = aclp100_101_contact_z
			boundary = 'ACLP100_101'
			paired_boundary = 'ACLP101_100'
			quantity = normal_force_z
		[]
		[tlp1_3_x]
			type = PenetrationAux
			variable = tlp1_3_contact_x
			boundary = 'TLP1_3'
			paired_boundary = 'TLP3_1'
			quantity = normal_force_x
		[]
		[tlp1_3_y]
			type = PenetrationAux
			variable = tlp1_3_contact_y
			boundary = 'TLP1_3'
			paired_boundary = 'TLP3_1'
			quantity = normal_force_y
		[]
		[tlp1_3_z]
			type = PenetrationAux
			variable = tlp1_3_contact_z
			boundary = 'TLP1_3'
			paired_boundary = 'TLP3_1'
			quantity = normal_force_z
		[]
		[tlp3_10_x]
			type = PenetrationAux
			variable = tlp3_10_contact_x
			boundary = 'TLP3_10'
			paired_boundary = 'TLP10_3'
			quantity = normal_force_x
		[]
		[tlp3_10_y]
			type = PenetrationAux
			variable = tlp3_10_contact_y
			boundary = 'TLP3_10'
			paired_boundary = 'TLP10_3'
			quantity = normal_force_y
		[]
		[tlp3_10_z]
			type = PenetrationAux
			variable = tlp3_10_contact_z
			boundary = 'TLP3_10'
			paired_boundary = 'TLP10_3'
			quantity = normal_force_z
		[]
		[tlp3_11_x]
			type = PenetrationAux
			variable = tlp3_11_contact_x
			boundary = 'TLP3_11'
			paired_boundary = 'TLP11_3'
			quantity = normal_force_x
		[]
		[tlp3_11_y]
			type = PenetrationAux
			variable = tlp3_11_contact_y
			boundary = 'TLP3_11'
			paired_boundary = 'TLP11_3'
			quantity = normal_force_y
		[]
		[tlp3_11_z]
			type = PenetrationAux
			variable = tlp3_11_contact_z
			boundary = 'TLP3_11'
			paired_boundary = 'TLP11_3'
			quantity = normal_force_z
		[]
		[tlp10_11_x]
			type = PenetrationAux
			variable = tlp10_11_contact_x
			boundary = 'TLP10_11'
			paired_boundary = 'TLP11_10'
			quantity = normal_force_x
		[]
		[tlp10_11_y]
			type = PenetrationAux
			variable = tlp10_11_contact_y
			boundary = 'TLP10_11'
			paired_boundary = 'TLP11_10'
			quantity = normal_force_y
		[]
		[tlp10_11_z]
			type = PenetrationAux
			variable = tlp10_11_contact_z
			boundary = 'TLP10_11'
			paired_boundary = 'TLP11_10'
			quantity = normal_force_z
		[]
		[tlp10_23_x]
			type = PenetrationAux
			variable = tlp10_23_contact_x
			boundary = 'TLP10_23'
			paired_boundary = 'TLP23_10'
			quantity = normal_force_x
		[]
		[tlp10_23_y]
			type = PenetrationAux
			variable = tlp10_23_contact_y
			boundary = 'TLP10_23'
			paired_boundary = 'TLP23_10'
			quantity = normal_force_y
		[]
		[tlp10_23_z]
			type = PenetrationAux
			variable = tlp10_23_contact_z
			boundary = 'TLP10_23'
			paired_boundary = 'TLP23_10'
			quantity = normal_force_z
		[]
		[tlp10_24_x]
			type = PenetrationAux
			variable = tlp10_24_contact_x
			boundary = 'TLP10_24'
			paired_boundary = 'TLP24_10'
			quantity = normal_force_x
		[]
		[tlp10_24_y]
			type = PenetrationAux
			variable = tlp10_24_contact_y
			boundary = 'TLP10_24'
			paired_boundary = 'TLP24_10'
			quantity = normal_force_y
		[]
		[tlp10_24_z]
			type = PenetrationAux
			variable = tlp10_24_contact_z
			boundary = 'TLP10_24'
			paired_boundary = 'TLP24_10'
			quantity = normal_force_z
		[]
		[tlp11_24_x]
			type = PenetrationAux
			variable = tlp11_24_contact_x
			boundary = 'TLP11_24'
			paired_boundary = 'TLP24_11'
			quantity = normal_force_x
		[]
		[tlp11_24_y]
			type = PenetrationAux
			variable = tlp11_24_contact_y
			boundary = 'TLP11_24'
			paired_boundary = 'TLP24_11'
			quantity = normal_force_y
		[]
		[tlp11_24_z]
			type = PenetrationAux
			variable = tlp11_24_contact_z
			boundary = 'TLP11_24'
			paired_boundary = 'TLP24_11'
			quantity = normal_force_z
		[]
		[tlp23_24_x]
			type = PenetrationAux
			variable = tlp23_24_contact_x
			boundary = 'TLP23_24'
			paired_boundary = 'TLP24_23'
			quantity = normal_force_x
		[]
		[tlp23_24_y]
			type = PenetrationAux
			variable = tlp23_24_contact_y
			boundary = 'TLP23_24'
			paired_boundary = 'TLP24_23'
			quantity = normal_force_y
		[]
		[tlp23_24_z]
			type = PenetrationAux
			variable = tlp23_24_contact_z
			boundary = 'TLP23_24'
			paired_boundary = 'TLP24_23'
			quantity = normal_force_z
		[]
		[tlp23_42_x]
			type = PenetrationAux
			variable = tlp23_42_contact_x
			boundary = 'TLP23_42'
			paired_boundary = 'TLP42_23'
			quantity = normal_force_x
		[]
		[tlp23_42_y]
			type = PenetrationAux
			variable = tlp23_42_contact_y
			boundary = 'TLP23_42'
			paired_boundary = 'TLP42_23'
			quantity = normal_force_y
		[]
		[tlp23_42_z]
			type = PenetrationAux
			variable = tlp23_42_contact_z
			boundary = 'TLP23_42'
			paired_boundary = 'TLP42_23'
			quantity = normal_force_z
		[]
		[tlp23_43_x]
			type = PenetrationAux
			variable = tlp23_43_contact_x
			boundary = 'TLP23_43'
			paired_boundary = 'TLP43_23'
			quantity = normal_force_x
		[]
		[tlp23_43_y]
			type = PenetrationAux
			variable = tlp23_43_contact_y
			boundary = 'TLP23_43'
			paired_boundary = 'TLP43_23'
			quantity = normal_force_y
		[]
		[tlp23_43_z]
			type = PenetrationAux
			variable = tlp23_43_contact_z
			boundary = 'TLP23_43'
			paired_boundary = 'TLP43_23'
			quantity = normal_force_z
		[]
		[tlp24_43_x]
			type = PenetrationAux
			variable = tlp24_43_contact_x
			boundary = 'TLP24_43'
			paired_boundary = 'TLP43_24'
			quantity = normal_force_x
		[]
		[tlp24_43_y]
			type = PenetrationAux
			variable = tlp24_43_contact_y
			boundary = 'TLP24_43'
			paired_boundary = 'TLP43_24'
			quantity = normal_force_y
		[]
		[tlp24_43_z]
			type = PenetrationAux
			variable = tlp24_43_contact_z
			boundary = 'TLP24_43'
			paired_boundary = 'TLP43_24'
			quantity = normal_force_z
		[]
		[tlp24_44_x]
			type = PenetrationAux
			variable = tlp24_44_contact_x
			boundary = 'TLP24_44'
			paired_boundary = 'TLP44_24'
			quantity = normal_force_x
		[]
		[tlp24_44_y]
			type = PenetrationAux
			variable = tlp24_44_contact_y
			boundary = 'TLP24_44'
			paired_boundary = 'TLP44_24'
			quantity = normal_force_y
		[]
		[tlp24_44_z]
			type = PenetrationAux
			variable = tlp24_44_contact_z
			boundary = 'TLP24_44'
			paired_boundary = 'TLP44_24'
			quantity = normal_force_z
		[]
		[tlp42_43_x]
			type = PenetrationAux
			variable = tlp42_43_contact_x
			boundary = 'TLP42_43'
			paired_boundary = 'TLP43_42'
			quantity = normal_force_x
		[]
		[tlp42_43_y]
			type = PenetrationAux
			variable = tlp42_43_contact_y
			boundary = 'TLP42_43'
			paired_boundary = 'TLP43_42'
			quantity = normal_force_y
		[]
		[tlp42_43_z]
			type = PenetrationAux
			variable = tlp42_43_contact_z
			boundary = 'TLP42_43'
			paired_boundary = 'TLP43_42'
			quantity = normal_force_z
		[]
		[tlp42_67_x]
			type = PenetrationAux
			variable = tlp42_67_contact_x
			boundary = 'TLP42_67'
			paired_boundary = 'TLP67_42'
			quantity = normal_force_x
		[]
		[tlp42_67_y]
			type = PenetrationAux
			variable = tlp42_67_contact_y
			boundary = 'TLP42_67'
			paired_boundary = 'TLP67_42'
			quantity = normal_force_y
		[]
		[tlp42_67_z]
			type = PenetrationAux
			variable = tlp42_67_contact_z
			boundary = 'TLP42_67'
			paired_boundary = 'TLP67_42'
			quantity = normal_force_z
		[]
		[tlp42_68_x]
			type = PenetrationAux
			variable = tlp42_68_contact_x
			boundary = 'TLP42_68'
			paired_boundary = 'TLP68_42'
			quantity = normal_force_x
		[]
		[tlp42_68_y]
			type = PenetrationAux
			variable = tlp42_68_contact_y
			boundary = 'TLP42_68'
			paired_boundary = 'TLP68_42'
			quantity = normal_force_y
		[]
		[tlp42_68_z]
			type = PenetrationAux
			variable = tlp42_68_contact_z
			boundary = 'TLP42_68'
			paired_boundary = 'TLP68_42'
			quantity = normal_force_z
		[]
		[tlp43_44_x]
			type = PenetrationAux
			variable = tlp43_44_contact_x
			boundary = 'TLP43_44'
			paired_boundary = 'TLP44_43'
			quantity = normal_force_x
		[]
		[tlp43_44_y]
			type = PenetrationAux
			variable = tlp43_44_contact_y
			boundary = 'TLP43_44'
			paired_boundary = 'TLP44_43'
			quantity = normal_force_y
		[]
		[tlp43_44_z]
			type = PenetrationAux
			variable = tlp43_44_contact_z
			boundary = 'TLP43_44'
			paired_boundary = 'TLP44_43'
			quantity = normal_force_z
		[]
		[tlp43_68_x]
			type = PenetrationAux
			variable = tlp43_68_contact_x
			boundary = 'TLP43_68'
			paired_boundary = 'TLP68_43'
			quantity = normal_force_x
		[]
		[tlp43_68_y]
			type = PenetrationAux
			variable = tlp43_68_contact_y
			boundary = 'TLP43_68'
			paired_boundary = 'TLP68_43'
			quantity = normal_force_y
		[]
		[tlp43_68_z]
			type = PenetrationAux
			variable = tlp43_68_contact_z
			boundary = 'TLP43_68'
			paired_boundary = 'TLP68_43'
			quantity = normal_force_z
		[]
		[tlp43_69_x]
			type = PenetrationAux
			variable = tlp43_69_contact_x
			boundary = 'TLP43_69'
			paired_boundary = 'TLP69_43'
			quantity = normal_force_x
		[]
		[tlp43_69_y]
			type = PenetrationAux
			variable = tlp43_69_contact_y
			boundary = 'TLP43_69'
			paired_boundary = 'TLP69_43'
			quantity = normal_force_y
		[]
		[tlp43_69_z]
			type = PenetrationAux
			variable = tlp43_69_contact_z
			boundary = 'TLP43_69'
			paired_boundary = 'TLP69_43'
			quantity = normal_force_z
		[]
		[tlp44_69_x]
			type = PenetrationAux
			variable = tlp44_69_contact_x
			boundary = 'TLP44_69'
			paired_boundary = 'TLP69_44'
			quantity = normal_force_x
		[]
		[tlp44_69_y]
			type = PenetrationAux
			variable = tlp44_69_contact_y
			boundary = 'TLP44_69'
			paired_boundary = 'TLP69_44'
			quantity = normal_force_y
		[]
		[tlp44_69_z]
			type = PenetrationAux
			variable = tlp44_69_contact_z
			boundary = 'TLP44_69'
			paired_boundary = 'TLP69_44'
			quantity = normal_force_z
		[]
		[tlp67_68_x]
			type = PenetrationAux
			variable = tlp67_68_contact_x
			boundary = 'TLP67_68'
			paired_boundary = 'TLP68_67'
			quantity = normal_force_x
		[]
		[tlp67_68_y]
			type = PenetrationAux
			variable = tlp67_68_contact_y
			boundary = 'TLP67_68'
			paired_boundary = 'TLP68_67'
			quantity = normal_force_y
		[]
		[tlp67_68_z]
			type = PenetrationAux
			variable = tlp67_68_contact_z
			boundary = 'TLP67_68'
			paired_boundary = 'TLP68_67'
			quantity = normal_force_z
		[]
		[tlp67_98_x]
			type = PenetrationAux
			variable = tlp67_98_contact_x
			boundary = 'TLP67_98'
			paired_boundary = 'TLP98_67'
			quantity = normal_force_x
		[]
		[tlp67_98_y]
			type = PenetrationAux
			variable = tlp67_98_contact_y
			boundary = 'TLP67_98'
			paired_boundary = 'TLP98_67'
			quantity = normal_force_y
		[]
		[tlp67_98_z]
			type = PenetrationAux
			variable = tlp67_98_contact_z
			boundary = 'TLP67_98'
			paired_boundary = 'TLP98_67'
			quantity = normal_force_z
		[]
		[tlp67_99_x]
			type = PenetrationAux
			variable = tlp67_99_contact_x
			boundary = 'TLP67_99'
			paired_boundary = 'TLP99_67'
			quantity = normal_force_x
		[]
		[tlp67_99_y]
			type = PenetrationAux
			variable = tlp67_99_contact_y
			boundary = 'TLP67_99'
			paired_boundary = 'TLP99_67'
			quantity = normal_force_y
		[]
		[tlp67_99_z]
			type = PenetrationAux
			variable = tlp67_99_contact_z
			boundary = 'TLP67_99'
			paired_boundary = 'TLP99_67'
			quantity = normal_force_z
		[]
		[tlp68_69_x]
			type = PenetrationAux
			variable = tlp68_69_contact_x
			boundary = 'TLP68_69'
			paired_boundary = 'TLP69_68'
			quantity = normal_force_x
		[]
		[tlp68_69_y]
			type = PenetrationAux
			variable = tlp68_69_contact_y
			boundary = 'TLP68_69'
			paired_boundary = 'TLP69_68'
			quantity = normal_force_y
		[]
		[tlp68_69_z]
			type = PenetrationAux
			variable = tlp68_69_contact_z
			boundary = 'TLP68_69'
			paired_boundary = 'TLP69_68'
			quantity = normal_force_z
		[]
		[tlp68_99_x]
			type = PenetrationAux
			variable = tlp68_99_contact_x
			boundary = 'TLP68_99'
			paired_boundary = 'TLP99_68'
			quantity = normal_force_x
		[]
		[tlp68_99_y]
			type = PenetrationAux
			variable = tlp68_99_contact_y
			boundary = 'TLP68_99'
			paired_boundary = 'TLP99_68'
			quantity = normal_force_y
		[]
		[tlp68_99_z]
			type = PenetrationAux
			variable = tlp68_99_contact_z
			boundary = 'TLP68_99'
			paired_boundary = 'TLP99_68'
			quantity = normal_force_z
		[]
		[tlp68_100_x]
			type = PenetrationAux
			variable = tlp68_100_contact_x
			boundary = 'TLP68_100'
			paired_boundary = 'TLP100_68'
			quantity = normal_force_x
		[]
		[tlp68_100_y]
			type = PenetrationAux
			variable = tlp68_100_contact_y
			boundary = 'TLP68_100'
			paired_boundary = 'TLP100_68'
			quantity = normal_force_y
		[]
		[tlp68_100_z]
			type = PenetrationAux
			variable = tlp68_100_contact_z
			boundary = 'TLP68_100'
			paired_boundary = 'TLP100_68'
			quantity = normal_force_z
		[]
		[tlp69_100_x]
			type = PenetrationAux
			variable = tlp69_100_contact_x
			boundary = 'TLP69_100'
			paired_boundary = 'TLP100_69'
			quantity = normal_force_x
		[]
		[tlp69_100_y]
			type = PenetrationAux
			variable = tlp69_100_contact_y
			boundary = 'TLP69_100'
			paired_boundary = 'TLP100_69'
			quantity = normal_force_y
		[]
		[tlp69_100_z]
			type = PenetrationAux
			variable = tlp69_100_contact_z
			boundary = 'TLP69_100'
			paired_boundary = 'TLP100_69'
			quantity = normal_force_z
		[]
		[tlp69_101_x]
			type = PenetrationAux
			variable = tlp69_101_contact_x
			boundary = 'TLP69_101'
			paired_boundary = 'TLP101_69'
			quantity = normal_force_x
		[]
		[tlp69_101_y]
			type = PenetrationAux
			variable = tlp69_101_contact_y
			boundary = 'TLP69_101'
			paired_boundary = 'TLP101_69'
			quantity = normal_force_y
		[]
		[tlp69_101_z]
			type = PenetrationAux
			variable = tlp69_101_contact_z
			boundary = 'TLP69_101'
			paired_boundary = 'TLP101_69'
			quantity = normal_force_z
		[]
		[tlp98_99_x]
			type = PenetrationAux
			variable = tlp98_99_contact_x
			boundary = 'TLP98_99'
			paired_boundary = 'TLP99_98'
			quantity = normal_force_x
		[]
		[tlp98_99_y]
			type = PenetrationAux
			variable = tlp98_99_contact_y
			boundary = 'TLP98_99'
			paired_boundary = 'TLP99_98'
			quantity = normal_force_y
		[]
		[tlp98_99_z]
			type = PenetrationAux
			variable = tlp98_99_contact_z
			boundary = 'TLP98_99'
			paired_boundary = 'TLP99_98'
			quantity = normal_force_z
		[]
		[tlp99_100_x]
			type = PenetrationAux
			variable = tlp99_100_contact_x
			boundary = 'TLP99_100'
			paired_boundary = 'TLP100_99'
			quantity = normal_force_x
		[]
		[tlp99_100_y]
			type = PenetrationAux
			variable = tlp99_100_contact_y
			boundary = 'TLP99_100'
			paired_boundary = 'TLP100_99'
			quantity = normal_force_y
		[]
		[tlp99_100_z]
			type = PenetrationAux
			variable = tlp99_100_contact_z
			boundary = 'TLP99_100'
			paired_boundary = 'TLP100_99'
			quantity = normal_force_z
		[]
		[tlp100_101_x]
			type = PenetrationAux
			variable = tlp100_101_contact_x
			boundary = 'TLP100_101'
			paired_boundary = 'TLP101_100'
			quantity = normal_force_x
		[]
		[tlp100_101_y]
			type = PenetrationAux
			variable = tlp100_101_contact_y
			boundary = 'TLP100_101'
			paired_boundary = 'TLP101_100'
			quantity = normal_force_y
		[]
		[tlp100_101_z]
			type = PenetrationAux
			variable = tlp100_101_contact_z
			boundary = 'TLP100_101'
			paired_boundary = 'TLP101_100'
			quantity = normal_force_z
		[]
		[tlp_block1_3_x]
			type = PenetrationAux
			variable = block1_tlp_contact_x
			paired_boundary = 'block1_tlp'
			boundary = 'TLP3_4'
			quantity = normal_force_x
		[]
		[tlp_block1_3_y]
			type = PenetrationAux
			variable = block1_tlp_contact_y
			paired_boundary = 'block1_tlp'
			boundary = 'TLP3_4'
			quantity = normal_force_y
		[]
		[tlp_block1_3_z]
			type = PenetrationAux
			variable = block1_tlp_contact_z
			paired_boundary = 'block1_tlp'
			boundary = 'TLP3_4'
			quantity = normal_force_z
		[]
		[tlp_block2_24_x]
			type = PenetrationAux
			variable = block2_tlp_contact_x
			paired_boundary = 'block2_tlp'
			boundary = 'TLP24_25'
			quantity = normal_force_x
		[]
		[tlp_block2_24_y]
			type = PenetrationAux
			variable = block2_tlp_contact_y
			paired_boundary = 'block2_tlp'
			boundary = 'TLP24_25'
			quantity = normal_force_y
		[]
		[tlp_block2_24_z]
			type = PenetrationAux
			variable = block2_tlp_contact_z
			paired_boundary = 'block2_tlp'
			boundary = 'TLP24_25'
			quantity = normal_force_z
		[]
		[tlp_block3_69_x]
			type = PenetrationAux
			variable = block3_tlp_contact_x
			paired_boundary = 'block3_tlp'
			boundary = 'TLP69_70'
			quantity = normal_force_x
		[]
		[tlp_block3_69_y]
			type = PenetrationAux
			variable = block3_tlp_contact_y
			paired_boundary = 'block3_tlp'
			boundary = 'TLP69_70'
			quantity = normal_force_y
		[]
		[tlp_block3_69_z]
			type = PenetrationAux
			variable = block3_tlp_contact_z
			paired_boundary = 'block3_tlp'
			boundary = 'TLP69_70'
			quantity = normal_force_z
		[]
		[aclp_block1_3_x]
			type = PenetrationAux
			variable = block1_aclp_contact_x
			paired_boundary = 'block1_aclp'
			boundary = 'ACLP3_4'
			quantity = normal_force_x
		[]
		[aclp_block1_3_y]
			type = PenetrationAux
			variable = block1_aclp_contact_y
			paired_boundary = 'block1_aclp'
			boundary = 'ACLP3_4'
			quantity = normal_force_y
		[]
		[aclp_block1_3_z]
			type = PenetrationAux
			variable = block1_aclp_contact_z
			paired_boundary = 'block1_aclp'
			boundary = 'ACLP3_4'
			quantity = normal_force_z
		[]
		[aclp_block2_24_x]
			type = PenetrationAux
			variable = block2_aclp_contact_x
			paired_boundary = 'block2_aclp'
			boundary = 'ACLP24_25'
			quantity = normal_force_x
		[]
		[aclp_block2_24_y]
			type = PenetrationAux
			variable = block2_aclp_contact_y
			paired_boundary = 'block2_aclp'
			boundary = 'ACLP24_25'
			quantity = normal_force_y
		[]
		[aclp_block2_24_z]
			type = PenetrationAux
			variable = block2_aclp_contact_z
			paired_boundary = 'block2_aclp'
			boundary = 'ACLP24_25'
			quantity = normal_force_z
		[]
		[aclp_block3_69_x]
			type = PenetrationAux
			variable = block3_aclp_contact_x
			paired_boundary = 'block3_aclp'
			boundary = 'ACLP69_70'
			quantity = normal_force_x
		[]
		[aclp_block3_69_y]
			type = PenetrationAux
			variable = block3_aclp_contact_y
			paired_boundary = 'block3_aclp'
			boundary = 'ACLP69_70'
			quantity = normal_force_y
		[]
		[aclp_block3_69_z]
			type = PenetrationAux
			variable = block3_aclp_contact_z
			paired_boundary = 'block3_aclp'
			boundary = 'ACLP69_70'
			quantity = normal_force_z
		[]
		[tlp101_restraint_x]
			type = PenetrationAux
			variable = tlp101_restraint_contact_x
			paired_boundary = 'TLP101_TLPRR4_2'
			boundary = 'TLPRR4_2'
			quantity = normal_force_x
		[]
		[tlp101_restraint_y]
			type = PenetrationAux
			variable = tlp101_restraint_contact_y
			paired_boundary = 'TLP101_TLPRR4_2'
			boundary = 'TLPRR4_2'
			quantity = normal_force_y
		[]
		[tlp101_restraint_z]
			type = PenetrationAux
			variable = tlp101_restraint_contact_z
			paired_boundary = 'TLP101_TLPRR4_2'
			boundary = 'TLPRR4_2'
			quantity = normal_force_z
		[]
		[aclp101_restraint_x]
			type = PenetrationAux
			variable = aclp101_restraint_contact_x
			paired_boundary = 'ACLP101_ACLPRR4_2'
			boundary = 'ACLPRR4_2'
			quantity = normal_force_x
		[]
		[aclp101_restraint_y]
			type = PenetrationAux
			variable = aclp101_restraint_contact_y
			paired_boundary = 'ACLP101_ACLPRR4_2'
			boundary = 'ACLPRR4_2'
			quantity = normal_force_y
		[]
		[aclp101_restraint_z]
			type = PenetrationAux
			variable = aclp101_restraint_contact_z
			paired_boundary = 'ACLP101_ACLPRR4_2'
			boundary = 'ACLPRR4_2'
			quantity = normal_force_z
		[]
	[]


[BCs]
  [restraint_x]
    type = DirichletBC
    variable = disp_x
    boundary = 'fix_restraint'
    value = 0.0
  []
  [restraint_y]
    type = DirichletBC
    variable = disp_y
    boundary = 'fix_restraint'
    value = 0.0
  []
  [restraint_z]
    type = DirichletBC
    variable = disp_z
    boundary = 'fix_restraint'
    value = 0.0
  []
  [InclinedNoDisplacementBC]
    [penalty_symmetry_duct]
      displacements = 'disp_x disp_y disp_z'
      boundary = 'right_symm'
      penalty = 1.0e17
    []
    [penalty_symmetry_nozzle]
      displacements = 'disp_x disp_y disp_z'
      boundary = 'right_symm_nozzle'
      penalty = 1.0e17
    []
  []
  [symm_left_x_duct]
    type = DirichletBC
    variable = disp_x
    boundary = 'left_symm'
    value = 0.0
  []
  [symm_left_x_nozzle]
    type = DirichletBC
    variable = disp_x
    boundary = 'left_symm_nozzle'
    value = 0.0
  []
  [no_x_block]
    type = DirichletBC
    variable = disp_x
    boundary = 'block_left block_front block_back block_top block_bottom'
    value = 0.0
  []
  [no_y_block]
    type = DirichletBC
    variable = disp_y
    boundary = 'block_left block_front block_back block_top block_bottom'
    value = 0.0
  []
  [no_z_block]
    type = DirichletBC
    variable = disp_z
    boundary = 'block_left block_front block_back block_top block_bottom'
    value = 0.0
  []
  [nozzle_bot_x]
    type = DirichletBC
    variable = disp_x
    boundary = 'nozzle_bottom'
    value = 0.0
  []
  [nozzle_bot_y]
    type = DirichletBC
    variable = disp_y
    boundary = 'nozzle_bottom'
    value = 0.0
  []
  [nozzle_bot_z]
    type = DirichletBC
    variable = disp_z
    boundary = 'nozzle_bottom'
    value = 0.0
  []
  [nozzle_top_x]
    type = DirichletBC
    variable = disp_x
    boundary = 'nozzle_top'
    value = 0.0
  []
  [nozzle_top_y]
    type = DirichletBC
    variable = disp_y
    boundary = 'nozzle_top'
    value = 0.0
  []
[]

[Materials]
  [elastic_tensor]
    type = ComputeIsotropicElasticityTensor
    youngs_modulus = ${E}
    poissons_ratio = 0.3
    block = '1 3 10 11 23 24 42 43 44 67 68 69 98 99 100 101 200 201'
  []
  # Fully constrained block. No need to increase stiffness
  [elastic_tensor2]
    type = ComputeIsotropicElasticityTensor
    youngs_modulus = ${E}
    poissons_ratio = 0.0
    block = '300'
  []
  [elastic_tensor_nozzle]
    type = ComputeIsotropicElasticityTensor
    youngs_modulus = 1.7e10
    poissons_ratio = 0.3
    block = '1000'
  []
  [stress]
    type = ComputeFiniteStrainElasticStress
    block = '1 3 10 11 23 24 67 68 69 98 99 100 101 200 201 300 1000'
  []
  [volumetric_swelling]
    type = ParsedMaterial
	property_name = volumetric_swelling
	coupled_variables = 'dose'
	expression = '0.03*dose/100'
	block = '42 43 44'
  []
  [irradiation_swelling_strain]
    type = ComputeVolumetricEigenstrain
	volumetric_materials = volumetric_swelling
	eigenstrain_name = swelling_strain
	args = ''
	block = '42 43 44'
  []
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
	block = '42 43 44'
  []
[]

[Contact]
  [aclp]
    primary =   'ACLP3_1 ACLP10_3 ACLP11_3 ACLP11_10 ACLP23_10 ACLP24_10 ACLP24_11 ACLP24_23 ACLP42_23 '
	            'ACLP43_23 ACLP43_24 ACLP44_24 ACLP43_42 ACLP67_42 ACLP68_42 ACLP44_43 ACLP68_43 '
				'ACLP69_43 ACLP69_44 ACLP68_67 ACLP98_67 ACLP99_67 ACLP69_68 ACLP99_68 ACLP100_68 '
				'ACLP100_69 ACLP101_69 ACLP99_98 ACLP100_99 ACLP101_100'
    secondary = 'ACLP1_3 ACLP3_10 ACLP3_11 ACLP10_11 ACLP10_23 ACLP10_24 ACLP11_24 ACLP23_24 '
	            'ACLP23_42 ACLP23_43 ACLP24_43 ACLP24_44 ACLP42_43 ACLP42_67 ACLP42_68 ACLP43_44 '
				'ACLP43_68 ACLP43_69 ACLP44_69 ACLP67_68 ACLP67_98 ACLP67_99 ACLP68_69 ACLP68_99 '
				'ACLP68_100 ACLP69_100 ACLP69_101 ACLP98_99 ACLP99_100 ACLP100_101'
    model = frictionless
    formulation = penalty
    normalize_penalty = true
    penalty = 5e10
    tangential_tolerance = 0.0
    normal_smoothing_distance = 0.0
  []
  [tlp]
    primary =   'TLP3_1 TLP10_3 TLP11_3 TLP11_10 TLP23_10 TLP24_10 TLP24_11 TLP24_23 TLP42_23 TLP43_23 '
	            'TLP43_24 TLP44_24 TLP43_42 TLP67_42 TLP68_42 TLP44_43 TLP68_43 TLP69_43 TLP69_44 '
				'TLP68_67 TLP98_67 TLP99_67 TLP69_68 TLP99_68 TLP100_68 TLP100_69 TLP101_69 TLP99_98 '
				'TLP100_99 TLP101_100'
    secondary = 'TLP1_3 TLP3_10 TLP3_11 TLP10_11 TLP10_23 TLP10_24 TLP11_24 TLP23_24 TLP23_42 '
	            'TLP23_43 TLP24_43 TLP24_44 TLP42_43 TLP42_67 TLP42_68 TLP43_44 TLP43_68 TLP43_69 '
				'TLP44_69 TLP67_68 TLP67_98 TLP67_99 TLP68_69 TLP68_99 TLP68_100 TLP69_100 TLP69_101 '
				'TLP98_99 TLP99_100 TLP100_101'
    model = frictionless
    formulation = penalty
    normalize_penalty = true
    penalty = 5e10
    tangential_tolerance = 0.0
    normal_smoothing_distance = 0.0
  []
  [aclp_rr]
    primary = 'ACLP98_ACLPRR1 ACLP98_ACLPRR2_1 ACLP99_ACLPRR2_2 ACLP99_ACLPRR3_1 ACLP100_ACLPRR3_2 '
              'ACLP100_ACLPRR4_1 ACLP101_ACLPRR4_2'
    secondary = 'ACLPRR1 ACLPRR2_1 ACLPRR2_2 ACLPRR3_1 ACLPRR3_2 ACLPRR4_1 ACLPRR4_2'
    model = frictionless
    formulation = penalty
    normalize_penalty = true
    penalty = 5e10
    tangential_tolerance = 0.0
    normal_smoothing_distance = 0.0
  []
  [tlp_rr]
    primary = 'TLP98_TLPRR1 TLP98_TLPRR2_1 TLP99_TLPRR2_2 TLP99_TLPRR3_1 TLP100_TLPRR3_2 '
              'TLP100_TLPRR4_1 TLP101_TLPRR4_2'
    secondary = 'TLPRR1 TLPRR2_1 TLPRR2_2 TLPRR3_1 TLPRR3_2 TLPRR4_1 TLPRR4_2'
    model = frictionless
    formulation = penalty
    normalize_penalty = true
    penalty = 5e10
    tangential_tolerance = 0.0
    normal_smoothing_distance = 0.0
  []
  [aclp_block]
    primary = 'ACLP3_4 ACLP24_25 ACLP69_70'
    secondary = 'block1_aclp block2_aclp block3_aclp'
    model = frictionless
    formulation = penalty
    normalize_penalty = true
    penalty = 5e10
    tangential_tolerance = 0.0
    normal_smoothing_distance = 0.0
  []
  [tlp_block]
    primary = 'TLP3_4 TLP24_25 TLP69_70'
    secondary = 'block1_tlp block2_tlp block3_tlp'
    model = frictionless
    formulation = penalty
    normalize_penalty = true
    penalty = 5e10
    tangential_tolerance = 0.0
    normal_smoothing_distance = 0.0
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
  #dt = 0.1
  end_time = 1.0
   [TimeStepper]
     type = TimeSequenceStepper
	 time_sequence = '0.05 0.10 0.12 0.14 0.16 0.19 0.23 0.27 0.32 0.37 0.44 0.50 0.61 0.72 0.85 1.0'
   []

  dtmax = 1
  dtmin = 0.005

  [Predictor]
    type = SimplePredictor
    scale = 0.8
  []
[]

	[Postprocessors]
		[pdata]
			type = PerfGraphData
			data_type = total
			section_name = 'Root'
			execute_on = timestep_end
		[]
		[aclp_force_1_2_half_force_x]
			type = NodalSum
			variable = aclp1_3_contact_x
			boundary = 'ACLP1_3'
		[]
		[aclp_force_1_2_half_force_y]
			type = NodalSum
			variable = aclp1_3_contact_y
			boundary = 'ACLP1_3'
		[]
		[aclp_force_1_2_half_force_z]
			type = NodalSum
			variable = aclp1_3_contact_z
			boundary = 'ACLP1_3'
		[]
		[aclp_force_1_2]
			type = ParsedPostprocessor
			expression = '2* sqrt(aclp_force_1_2_half_force_x*aclp_force_1_2_half_force_x+aclp_force_1_2_half_force_y*aclp_force_1_2_half_force_y+aclp_force_1_2_half_force_z*aclp_force_1_2_half_force_z)'
			pp_names = 'aclp_force_1_2_half_force_x aclp_force_1_2_half_force_y aclp_force_1_2_half_force_z'
		[]

		[aclp_force_11_1_x]
			type = NodalSum
			variable = aclp10_11_contact_x
			boundary = 'ACLP10_11'
		[]
		[aclp_force_11_1_y]
			type = NodalSum
			variable = aclp10_11_contact_y
			boundary = 'ACLP10_11'
		[]
		[aclp_force_11_1_z]
			type = NodalSum
			variable = aclp10_11_contact_z
			boundary = 'ACLP10_11'
		[]
		[aclp_force_11_1]
			type = ParsedPostprocessor
			expression = 'sqrt(aclp_force_11_1_x*aclp_force_11_1_x+aclp_force_11_1_y*aclp_force_11_1_y+aclp_force_11_1_z*aclp_force_11_1_z)'
			pp_names = 'aclp_force_11_1_x aclp_force_11_1_y aclp_force_11_1_z'
		[]

		[aclp_force_11_2_x]
			type = NodalSum
			variable = aclp11_24_contact_x
			boundary = 'ACLP11_24'
		[]
		[aclp_force_11_2_y]
			type = NodalSum
			variable = aclp11_24_contact_y
			boundary = 'ACLP11_24'
		[]
		[aclp_force_11_2_z]
			type = NodalSum
			variable = aclp11_24_contact_z
			boundary = 'ACLP11_24'
		[]
		[aclp_force_11_2]
			type = ParsedPostprocessor
			expression = 'sqrt(aclp_force_11_2_x*aclp_force_11_2_x+aclp_force_11_2_y*aclp_force_11_2_y+aclp_force_11_2_z*aclp_force_11_2_z)'
			pp_names = 'aclp_force_11_2_x aclp_force_11_2_y aclp_force_11_2_z'
		[]

		[aclp_force_11_6_x]
			type = NodalSum
			variable = aclp3_11_contact_x
			boundary = 'ACLP3_11'
		[]
		[aclp_force_11_6_y]
			type = NodalSum
			variable = aclp3_11_contact_y
			boundary = 'ACLP3_11'
		[]
		[aclp_force_11_6_z]
			type = NodalSum
			variable = aclp3_11_contact_z
			boundary = 'ACLP3_11'
		[]
		[aclp_force_11_6]
			type = ParsedPostprocessor
			expression = 'sqrt(aclp_force_11_6_x*aclp_force_11_6_x+aclp_force_11_6_y*aclp_force_11_6_y+aclp_force_11_6_z*aclp_force_11_6_z)'
			pp_names = 'aclp_force_11_6_x aclp_force_11_6_y aclp_force_11_6_z'
		[]

		[aclp_force_42_2_half_force_x]
			type = NodalSum
			variable = aclp42_67_contact_x
			boundary = 'ACLP42_67'
		[]
		[aclp_force_42_2_half_force_y]
			type = NodalSum
			variable = aclp42_67_contact_y
			boundary = 'ACLP42_67'
		[]
		[aclp_force_42_2_half_force_z]
			type = NodalSum
			variable = aclp42_67_contact_z
			boundary = 'ACLP42_67'
		[]
		[aclp_force_42_2]
			type = ParsedPostprocessor
			expression = '2* sqrt(aclp_force_42_2_half_force_x*aclp_force_42_2_half_force_x+aclp_force_42_2_half_force_y*aclp_force_42_2_half_force_y+aclp_force_42_2_half_force_z*aclp_force_42_2_half_force_z)'
			pp_names = 'aclp_force_42_2_half_force_x aclp_force_42_2_half_force_y aclp_force_42_2_half_force_z'
		[]

		[aclp_force_42_3_x]
			type = NodalSum
			variable = aclp42_68_contact_x
			boundary = 'ACLP42_68'
		[]
		[aclp_force_42_3_y]
			type = NodalSum
			variable = aclp42_68_contact_y
			boundary = 'ACLP42_68'
		[]
		[aclp_force_42_3_z]
			type = NodalSum
			variable = aclp42_68_contact_z
			boundary = 'ACLP42_68'
		[]
		[aclp_force_42_3]
			type = ParsedPostprocessor
			expression = 'sqrt(aclp_force_42_3_x*aclp_force_42_3_x+aclp_force_42_3_y*aclp_force_42_3_y+aclp_force_42_3_z*aclp_force_42_3_z)'
			pp_names = 'aclp_force_42_3_x aclp_force_42_3_y aclp_force_42_3_z'
		[]

		[aclp_force_42_5_x]
			type = NodalSum
			variable = aclp23_42_contact_x
			boundary = 'ACLP23_42'
		[]
		[aclp_force_42_5_y]
			type = NodalSum
			variable = aclp23_42_contact_y
			boundary = 'ACLP23_42'
		[]
		[aclp_force_42_5_z]
			type = NodalSum
			variable = aclp23_42_contact_z
			boundary = 'ACLP23_42'
		[]
		[aclp_force_42_5]
			type = ParsedPostprocessor
			expression = 'sqrt(aclp_force_42_5_x*aclp_force_42_5_x+aclp_force_42_5_y*aclp_force_42_5_y+aclp_force_42_5_z*aclp_force_42_5_z)'
			pp_names = 'aclp_force_42_5_x aclp_force_42_5_y aclp_force_42_5_z'
		[]

		[aclp_force_43_1_x]
			type = NodalSum
			variable = aclp42_43_contact_x
			boundary = 'ACLP42_43'
		[]
		[aclp_force_43_1_y]
			type = NodalSum
			variable = aclp42_43_contact_y
			boundary = 'ACLP42_43'
		[]
		[aclp_force_43_1_z]
			type = NodalSum
			variable = aclp42_43_contact_z
			boundary = 'ACLP42_43'
		[]
		[aclp_force_43_1]
			type = ParsedPostprocessor
			expression = 'sqrt(aclp_force_43_1_x*aclp_force_43_1_x+aclp_force_43_1_y*aclp_force_43_1_y+aclp_force_43_1_z*aclp_force_43_1_z)'
			pp_names = 'aclp_force_43_1_x aclp_force_43_1_y aclp_force_43_1_z'
		[]

		[aclp_force_43_2_x]
			type = NodalSum
			variable = aclp43_68_contact_x
			boundary = 'ACLP43_68'
		[]
		[aclp_force_43_2_y]
			type = NodalSum
			variable = aclp43_68_contact_y
			boundary = 'ACLP43_68'
		[]
		[aclp_force_43_2_z]
			type = NodalSum
			variable = aclp43_68_contact_z
			boundary = 'ACLP43_68'
		[]
		[aclp_force_43_2]
			type = ParsedPostprocessor
			expression = 'sqrt(aclp_force_43_2_x*aclp_force_43_2_x+aclp_force_43_2_y*aclp_force_43_2_y+aclp_force_43_2_z*aclp_force_43_2_z)'
			pp_names = 'aclp_force_43_2_x aclp_force_43_2_y aclp_force_43_2_z'
		[]

		[aclp_force_43_3_x]
			type = NodalSum
			variable = aclp43_69_contact_x
			boundary = 'ACLP43_69'
		[]
		[aclp_force_43_3_y]
			type = NodalSum
			variable = aclp43_69_contact_y
			boundary = 'ACLP43_69'
		[]
		[aclp_force_43_3_z]
			type = NodalSum
			variable = aclp43_69_contact_z
			boundary = 'ACLP43_69'
		[]
		[aclp_force_43_3]
			type = ParsedPostprocessor
			expression = 'sqrt(aclp_force_43_3_x*aclp_force_43_3_x+aclp_force_43_3_y*aclp_force_43_3_y+aclp_force_43_3_z*aclp_force_43_3_z)'
			pp_names = 'aclp_force_43_3_x aclp_force_43_3_y aclp_force_43_3_z'
		[]

		[aclp_force_43_5_x]
			type = NodalSum
			variable = aclp24_43_contact_x
			boundary = 'ACLP24_43'
		[]
		[aclp_force_43_5_y]
			type = NodalSum
			variable = aclp24_43_contact_y
			boundary = 'ACLP24_43'
		[]
		[aclp_force_43_5_z]
			type = NodalSum
			variable = aclp24_43_contact_z
			boundary = 'ACLP24_43'
		[]
		[aclp_force_43_5]
			type = ParsedPostprocessor
			expression = 'sqrt(aclp_force_43_5_x*aclp_force_43_5_x+aclp_force_43_5_y*aclp_force_43_5_y+aclp_force_43_5_z*aclp_force_43_5_z)'
			pp_names = 'aclp_force_43_5_x aclp_force_43_5_y aclp_force_43_5_z'
		[]

		[aclp_force_43_6_x]
			type = NodalSum
			variable = aclp23_43_contact_x
			boundary = 'ACLP23_43'
		[]
		[aclp_force_43_6_y]
			type = NodalSum
			variable = aclp23_43_contact_y
			boundary = 'ACLP23_43'
		[]
		[aclp_force_43_6_z]
			type = NodalSum
			variable = aclp23_43_contact_z
			boundary = 'ACLP23_43'
		[]
		[aclp_force_43_6]
			type = ParsedPostprocessor
			expression = 'sqrt(aclp_force_43_6_x*aclp_force_43_6_x+aclp_force_43_6_y*aclp_force_43_6_y+aclp_force_43_6_z*aclp_force_43_6_z)'
			pp_names = 'aclp_force_43_6_x aclp_force_43_6_y aclp_force_43_6_z'
		[]

		[aclp_force_44_1_x]
			type = NodalSum
			variable = aclp43_44_contact_x
			boundary = 'ACLP43_44'
		[]
		[aclp_force_44_1_y]
			type = NodalSum
			variable = aclp43_44_contact_y
			boundary = 'ACLP43_44'
		[]
		[aclp_force_44_1_z]
			type = NodalSum
			variable = aclp43_44_contact_z
			boundary = 'ACLP43_44'
		[]
		[aclp_force_44_1]
			type = ParsedPostprocessor
			expression = 'sqrt(aclp_force_44_1_x*aclp_force_44_1_x+aclp_force_44_1_y*aclp_force_44_1_y+aclp_force_44_1_z*aclp_force_44_1_z)'
			pp_names = 'aclp_force_44_1_x aclp_force_44_1_y aclp_force_44_1_z'
		[]

		[aclp_force_44_2_x]
			type = NodalSum
			variable = aclp44_69_contact_x
			boundary = 'ACLP44_69'
		[]
		[aclp_force_44_2_y]
			type = NodalSum
			variable = aclp44_69_contact_y
			boundary = 'ACLP44_69'
		[]
		[aclp_force_44_2_z]
			type = NodalSum
			variable = aclp44_69_contact_z
			boundary = 'ACLP44_69'
		[]
		[aclp_force_44_2]
			type = ParsedPostprocessor
			expression = 'sqrt(aclp_force_44_2_x*aclp_force_44_2_x+aclp_force_44_2_y*aclp_force_44_2_y+aclp_force_44_2_z*aclp_force_44_2_z)'
			pp_names = 'aclp_force_44_2_x aclp_force_44_2_y aclp_force_44_2_z'
		[]

		[aclp_force_44_6_x]
			type = NodalSum
			variable = aclp24_44_contact_x
			boundary = 'ACLP24_44'
		[]
		[aclp_force_44_6_y]
			type = NodalSum
			variable = aclp24_44_contact_y
			boundary = 'ACLP24_44'
		[]
		[aclp_force_44_6_z]
			type = NodalSum
			variable = aclp24_44_contact_z
			boundary = 'ACLP24_44'
		[]
		[aclp_force_44_6]
			type = ParsedPostprocessor
			expression = 'sqrt(aclp_force_44_6_x*aclp_force_44_6_x+aclp_force_44_6_y*aclp_force_44_6_y+aclp_force_44_6_z*aclp_force_44_6_z)'
			pp_names = 'aclp_force_44_6_x aclp_force_44_6_y aclp_force_44_6_z'
		[]

		[tlp_force_1_2_half_force_x]
			type = NodalSum
			variable = tlp1_3_contact_x
			boundary = 'TLP1_3'
		[]
		[tlp_force_1_2_half_force_y]
			type = NodalSum
			variable = tlp1_3_contact_y
			boundary = 'TLP1_3'
		[]
		[tlp_force_1_2_half_force_z]
			type = NodalSum
			variable = tlp1_3_contact_z
			boundary = 'TLP1_3'
		[]
		[tlp_force_1_2]
			type = ParsedPostprocessor
			expression = '2* sqrt(tlp_force_1_2_half_force_x*tlp_force_1_2_half_force_x+tlp_force_1_2_half_force_y*tlp_force_1_2_half_force_y+tlp_force_1_2_half_force_z*tlp_force_1_2_half_force_z)'
			pp_names = 'tlp_force_1_2_half_force_x tlp_force_1_2_half_force_y tlp_force_1_2_half_force_z'
		[]

		[tlp_force_11_1_x]
			type = NodalSum
			variable = tlp10_11_contact_x
			boundary = 'TLP10_11'
		[]
		[tlp_force_11_1_y]
			type = NodalSum
			variable = tlp10_11_contact_y
			boundary = 'TLP10_11'
		[]
		[tlp_force_11_1_z]
			type = NodalSum
			variable = tlp10_11_contact_z
			boundary = 'TLP10_11'
		[]
		[tlp_force_11_1]
			type = ParsedPostprocessor
			expression = 'sqrt(tlp_force_11_1_x*tlp_force_11_1_x+tlp_force_11_1_y*tlp_force_11_1_y+tlp_force_11_1_z*tlp_force_11_1_z)'
			pp_names = 'tlp_force_11_1_x tlp_force_11_1_y tlp_force_11_1_z'
		[]

		[tlp_force_11_2_x]
			type = NodalSum
			variable = tlp11_24_contact_x
			boundary = 'TLP11_24'
		[]
		[tlp_force_11_2_y]
			type = NodalSum
			variable = tlp11_24_contact_y
			boundary = 'TLP11_24'
		[]
		[tlp_force_11_2_z]
			type = NodalSum
			variable = tlp11_24_contact_z
			boundary = 'TLP11_24'
		[]
		[tlp_force_11_2]
			type = ParsedPostprocessor
			expression = 'sqrt(tlp_force_11_2_x*tlp_force_11_2_x+tlp_force_11_2_y*tlp_force_11_2_y+tlp_force_11_2_z*tlp_force_11_2_z)'
			pp_names = 'tlp_force_11_2_x tlp_force_11_2_y tlp_force_11_2_z'
		[]

		[tlp_force_11_6_x]
			type = NodalSum
			variable = tlp3_11_contact_x
			boundary = 'TLP3_11'
		[]
		[tlp_force_11_6_y]
			type = NodalSum
			variable = tlp3_11_contact_y
			boundary = 'TLP3_11'
		[]
		[tlp_force_11_6_z]
			type = NodalSum
			variable = tlp3_11_contact_z
			boundary = 'TLP3_11'
		[]
		[tlp_force_11_6]
			type = ParsedPostprocessor
			expression = 'sqrt(tlp_force_11_6_x*tlp_force_11_6_x+tlp_force_11_6_y*tlp_force_11_6_y+tlp_force_11_6_z*tlp_force_11_6_z)'
			pp_names = 'tlp_force_11_6_x tlp_force_11_6_y tlp_force_11_6_z'
		[]

		[tlp_force_42_2_half_force_x]
			type = NodalSum
			variable = tlp42_67_contact_x
			boundary = 'TLP42_67'
		[]
		[tlp_force_42_2_half_force_y]
			type = NodalSum
			variable = tlp42_67_contact_y
			boundary = 'TLP42_67'
		[]
		[tlp_force_42_2_half_force_z]
			type = NodalSum
			variable = tlp42_67_contact_z
			boundary = 'TLP42_67'
		[]
		[tlp_force_42_2]
			type = ParsedPostprocessor
			expression = '2* sqrt(tlp_force_42_2_half_force_x*tlp_force_42_2_half_force_x+tlp_force_42_2_half_force_y*tlp_force_42_2_half_force_y+tlp_force_42_2_half_force_z*tlp_force_42_2_half_force_z)'
			pp_names = 'tlp_force_42_2_half_force_x tlp_force_42_2_half_force_y tlp_force_42_2_half_force_z'
		[]

		[tlp_force_42_3_x]
			type = NodalSum
			variable = tlp42_68_contact_x
			boundary = 'TLP42_68'
		[]
		[tlp_force_42_3_y]
			type = NodalSum
			variable = tlp42_68_contact_y
			boundary = 'TLP42_68'
		[]
		[tlp_force_42_3_z]
			type = NodalSum
			variable = tlp42_68_contact_z
			boundary = 'TLP42_68'
		[]
		[tlp_force_42_3]
			type = ParsedPostprocessor
			expression = 'sqrt(tlp_force_42_3_x*tlp_force_42_3_x+tlp_force_42_3_y*tlp_force_42_3_y+tlp_force_42_3_z*tlp_force_42_3_z)'
			pp_names = 'tlp_force_42_3_x tlp_force_42_3_y tlp_force_42_3_z'
		[]

		[tlp_force_42_5_x]
			type = NodalSum
			variable = tlp23_42_contact_x
			boundary = 'TLP23_42'
		[]
		[tlp_force_42_5_y]
			type = NodalSum
			variable = tlp23_42_contact_y
			boundary = 'TLP23_42'
		[]
		[tlp_force_42_5_z]
			type = NodalSum
			variable = tlp23_42_contact_z
			boundary = 'TLP23_42'
		[]
		[tlp_force_42_5]
			type = ParsedPostprocessor
			expression = 'sqrt(tlp_force_42_5_x*tlp_force_42_5_x+tlp_force_42_5_y*tlp_force_42_5_y+tlp_force_42_5_z*tlp_force_42_5_z)'
			pp_names = 'tlp_force_42_5_x tlp_force_42_5_y tlp_force_42_5_z'
		[]

		[tlp_force_43_1_x]
			type = NodalSum
			variable = tlp42_43_contact_x
			boundary = 'TLP42_43'
		[]
		[tlp_force_43_1_y]
			type = NodalSum
			variable = tlp42_43_contact_y
			boundary = 'TLP42_43'
		[]
		[tlp_force_43_1_z]
			type = NodalSum
			variable = tlp42_43_contact_z
			boundary = 'TLP42_43'
		[]
		[tlp_force_43_1]
			type = ParsedPostprocessor
			expression = 'sqrt(tlp_force_43_1_x*tlp_force_43_1_x+tlp_force_43_1_y*tlp_force_43_1_y+tlp_force_43_1_z*tlp_force_43_1_z)'
			pp_names = 'tlp_force_43_1_x tlp_force_43_1_y tlp_force_43_1_z'
		[]

		[tlp_force_43_2_x]
			type = NodalSum
			variable = tlp43_68_contact_x
			boundary = 'TLP43_68'
		[]
		[tlp_force_43_2_y]
			type = NodalSum
			variable = tlp43_68_contact_y
			boundary = 'TLP43_68'
		[]
		[tlp_force_43_2_z]
			type = NodalSum
			variable = tlp43_68_contact_z
			boundary = 'TLP43_68'
		[]
		[tlp_force_43_2]
			type = ParsedPostprocessor
			expression = 'sqrt(tlp_force_43_2_x*tlp_force_43_2_x+tlp_force_43_2_y*tlp_force_43_2_y+tlp_force_43_2_z*tlp_force_43_2_z)'
			pp_names = 'tlp_force_43_2_x tlp_force_43_2_y tlp_force_43_2_z'
		[]

		[tlp_force_43_3_x]
			type = NodalSum
			variable = tlp43_69_contact_x
			boundary = 'TLP43_69'
		[]
		[tlp_force_43_3_y]
			type = NodalSum
			variable = tlp43_69_contact_y
			boundary = 'TLP43_69'
		[]
		[tlp_force_43_3_z]
			type = NodalSum
			variable = tlp43_69_contact_z
			boundary = 'TLP43_69'
		[]
		[tlp_force_43_3]
			type = ParsedPostprocessor
			expression = 'sqrt(tlp_force_43_3_x*tlp_force_43_3_x+tlp_force_43_3_y*tlp_force_43_3_y+tlp_force_43_3_z*tlp_force_43_3_z)'
			pp_names = 'tlp_force_43_3_x tlp_force_43_3_y tlp_force_43_3_z'
		[]

		[tlp_force_43_5_x]
			type = NodalSum
			variable = tlp24_43_contact_x
			boundary = 'TLP24_43'
		[]
		[tlp_force_43_5_y]
			type = NodalSum
			variable = tlp24_43_contact_y
			boundary = 'TLP24_43'
		[]
		[tlp_force_43_5_z]
			type = NodalSum
			variable = tlp24_43_contact_z
			boundary = 'TLP24_43'
		[]
		[tlp_force_43_5]
			type = ParsedPostprocessor
			expression = 'sqrt(tlp_force_43_5_x*tlp_force_43_5_x+tlp_force_43_5_y*tlp_force_43_5_y+tlp_force_43_5_z*tlp_force_43_5_z)'
			pp_names = 'tlp_force_43_5_x tlp_force_43_5_y tlp_force_43_5_z'
		[]

		[tlp_force_43_6_x]
			type = NodalSum
			variable = tlp23_43_contact_x
			boundary = 'TLP23_43'
		[]
		[tlp_force_43_6_y]
			type = NodalSum
			variable = tlp23_43_contact_y
			boundary = 'TLP23_43'
		[]
		[tlp_force_43_6_z]
			type = NodalSum
			variable = tlp23_43_contact_z
			boundary = 'TLP23_43'
		[]
		[tlp_force_43_6]
			type = ParsedPostprocessor
			expression = 'sqrt(tlp_force_43_6_x*tlp_force_43_6_x+tlp_force_43_6_y*tlp_force_43_6_y+tlp_force_43_6_z*tlp_force_43_6_z)'
			pp_names = 'tlp_force_43_6_x tlp_force_43_6_y tlp_force_43_6_z'
		[]

		[tlp_force_44_1_x]
			type = NodalSum
			variable = tlp43_44_contact_x
			boundary = 'TLP43_44'
		[]
		[tlp_force_44_1_y]
			type = NodalSum
			variable = tlp43_44_contact_y
			boundary = 'TLP43_44'
		[]
		[tlp_force_44_1_z]
			type = NodalSum
			variable = tlp43_44_contact_z
			boundary = 'TLP43_44'
		[]
		[tlp_force_44_1]
			type = ParsedPostprocessor
			expression = 'sqrt(tlp_force_44_1_x*tlp_force_44_1_x+tlp_force_44_1_y*tlp_force_44_1_y+tlp_force_44_1_z*tlp_force_44_1_z)'
			pp_names = 'tlp_force_44_1_x tlp_force_44_1_y tlp_force_44_1_z'
		[]

		[tlp_force_44_2_x]
			type = NodalSum
			variable = tlp44_69_contact_x
			boundary = 'TLP44_69'
		[]
		[tlp_force_44_2_y]
			type = NodalSum
			variable = tlp44_69_contact_y
			boundary = 'TLP44_69'
		[]
		[tlp_force_44_2_z]
			type = NodalSum
			variable = tlp44_69_contact_z
			boundary = 'TLP44_69'
		[]
		[tlp_force_44_2]
			type = ParsedPostprocessor
			expression = 'sqrt(tlp_force_44_2_x*tlp_force_44_2_x+tlp_force_44_2_y*tlp_force_44_2_y+tlp_force_44_2_z*tlp_force_44_2_z)'
			pp_names = 'tlp_force_44_2_x tlp_force_44_2_y tlp_force_44_2_z'
		[]

		[tlp_force_44_6_x]
			type = NodalSum
			variable = tlp24_44_contact_x
			boundary = 'TLP24_44'
		[]
		[tlp_force_44_6_y]
			type = NodalSum
			variable = tlp24_44_contact_y
			boundary = 'TLP24_44'
		[]
		[tlp_force_44_6_z]
			type = NodalSum
			variable = tlp24_44_contact_z
			boundary = 'TLP24_44'
		[]
		[tlp_force_44_6]
			type = ParsedPostprocessor
			expression = 'sqrt(tlp_force_44_6_x*tlp_force_44_6_x+tlp_force_44_6_y*tlp_force_44_6_y+tlp_force_44_6_z*tlp_force_44_6_z)'
			pp_names = 'tlp_force_44_6_x tlp_force_44_6_y tlp_force_44_6_z'
		[]

		[tlp_force_101_1_x]
			type = NodalSum
			variable = tlp100_101_contact_x
			boundary = 'TLP100_101'
		[]
		[tlp_force_101_1_y]
			type = NodalSum
			variable = tlp100_101_contact_y
			boundary = 'TLP100_101'
		[]
		[tlp_force_101_1_z]
			type = NodalSum
			variable = tlp100_101_contact_z
			boundary = 'TLP100_101'
		[]
		[tlp_force_101_1]
			type = ParsedPostprocessor
			expression = 'sqrt(tlp_force_101_1_x*tlp_force_101_1_x+tlp_force_101_1_y*tlp_force_101_1_y+tlp_force_101_1_z*tlp_force_101_1_z)'
			pp_names = 'tlp_force_101_1_x tlp_force_101_1_y tlp_force_101_1_z'
		[]

		[tlp_force_101_2_x]
			type = NodalSum
			variable = tlp101_restraint_contact_x
			boundary = 'TLPRR4_2'
		[]
		[tlp_force_101_2_y]
			type = NodalSum
			variable = tlp101_restraint_contact_y
			boundary = 'TLPRR4_2'
		[]
		[tlp_force_101_2_z]
			type = NodalSum
			variable = tlp101_restraint_contact_z
			boundary = 'TLPRR4_2'
		[]
		[tlp_force_101_2]
			type = ParsedPostprocessor
			expression = 'sqrt(tlp_force_101_2_x*tlp_force_101_2_x+tlp_force_101_2_y*tlp_force_101_2_y+tlp_force_101_2_z*tlp_force_101_2_z)'
			pp_names = 'tlp_force_101_2_x tlp_force_101_2_y tlp_force_101_2_z'
		[]

		[tlp_force_101_6_x]
			type = NodalSum
			variable = tlp69_101_contact_x
			boundary = 'TLP69_101'
		[]
		[tlp_force_101_6_y]
			type = NodalSum
			variable = tlp69_101_contact_y
			boundary = 'TLP69_101'
		[]
		[tlp_force_101_6_z]
			type = NodalSum
			variable = tlp69_101_contact_z
			boundary = 'TLP69_101'
		[]
		[tlp_force_101_6]
			type = ParsedPostprocessor
			expression = 'sqrt(tlp_force_101_6_x*tlp_force_101_6_x+tlp_force_101_6_y*tlp_force_101_6_y+tlp_force_101_6_z*tlp_force_101_6_z)'
			pp_names = 'tlp_force_101_6_x tlp_force_101_6_y tlp_force_101_6_z'
		[]

	[]



[VectorPostprocessors]
  [duct_11_average]
    type = AverageSectionValueSampler
    axis_direction = '0 0 1'
    block = '11'
    variables = 'disp_x disp_y disp_z'
    reference_point = '0.0 0.242234 0.0'
	require_equal_node_counts = false
	tolerance = 1e-5
	symmetry_plane = '-1 0 0'
  []
  [duct_42_average]
    type = AverageSectionValueSampler
    axis_direction = '0 0 1'
    block = '42'
    variables = 'disp_x disp_y disp_z'
    reference_point = '0.279708 0.484468 0.0'
	require_equal_node_counts = false
	tolerance = 1e-5
	symmetry_plane = '0.866 -0.5 0'
  []
  [duct_43_average]
    type = AverageSectionValueSampler
    axis_direction = '0 0 1'
    block = '43'
    variables = 'disp_x disp_y disp_z'
    reference_point = '0.139854 0.484468 0.0'
	require_equal_node_counts = false
	tolerance = 1e-5
  []
  [duct_44_average]
    type = AverageSectionValueSampler
    axis_direction = '0 0 1'
    block = '44'
    variables = 'disp_x disp_y disp_z'
    reference_point = '0.0 0.484468 0.0'
	require_equal_node_counts = false
	tolerance = 1e-5
	symmetry_plane = '-1 0 0'
  []
  [duct_101_average]
    type = AverageSectionValueSampler
    axis_direction = '0 0 1'
    block = '101'
    variables = 'disp_x disp_y disp_z'
    reference_point = '0.0 0.726702 0.0'
	require_equal_node_counts = false
	tolerance = 1e-5
	symmetry_plane = '-1 0 0'
  []
[]

[Outputs]
  exodus = true
  perf_graph = true
  [test_average_section]
    type = CSV
    file_base = average_section_disp
    execute_on = final
    show = 'duct_11_average duct_42_average duct_43_average duct_44_average duct_101_average'
  []
  [tlp_force_plots]
    type = CSV
    file_base = tlp_force_plots
    execute_on = timestep_end
    show = 'tlp_force_1_2  tlp_force_11_1 tlp_force_11_2 tlp_force_11_6 tlp_force_42_2 tlp_force_42_3 '
	       'tlp_force_42_5 tlp_force_43_1 tlp_force_43_2 tlp_force_43_3 tlp_force_43_5 tlp_force_43_6 '
		   'tlp_force_44_1 tlp_force_44_2 tlp_force_44_6 tlp_force_101_1 tlp_force_101_2 tlp_force_101_6'
  []
  [aclp_force_plots]
    type = CSV
    file_base = aclp_force_plots
    execute_on = timestep_end
    show = 'aclp_force_1_2  aclp_force_11_1 aclp_force_11_2 aclp_force_11_6 aclp_force_42_2 aclp_force_42_3 '
	       'aclp_force_42_5 aclp_force_43_1 aclp_force_43_2 aclp_force_43_3 aclp_force_43_5 aclp_force_43_6 '
		   'aclp_force_44_1 aclp_force_44_2 aclp_force_44_6'
  []
  [runtime]
    type = CSV
    file_base = runtime
    show = pdata
    execute_on = timestep_end
  []
[]
