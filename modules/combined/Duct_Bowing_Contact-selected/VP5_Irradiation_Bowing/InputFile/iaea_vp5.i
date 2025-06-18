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
    file = vp5_mesh.e
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
  [aclp1_3_contact]
  []
  [aclp3_10_contact]
  []
  [aclp3_11_contact]
  []
  [aclp10_11_contact]
  []
  [aclp10_23_contact]
  []
  [aclp10_24_contact]
  []
  [aclp11_24_contact]
  []
  [aclp23_24_contact]
  []
  [aclp23_42_contact]
  []
  [aclp23_43_contact]
  []
  [aclp24_43_contact]
  []
  [aclp24_44_contact]
  []
  [aclp42_43_contact]
  []
  [aclp42_67_contact]
  []
  [aclp42_68_contact]
  []
  [aclp43_44_contact]
  []
  [aclp43_68_contact]
  []
  [aclp43_69_contact]
  []
  [aclp44_69_contact]
  []
  [aclp67_68_contact]
  []
  [aclp67_98_contact]
  []
  [aclp67_99_contact]
  []
  [aclp68_69_contact]
  []
  [aclp68_99_contact]
  []
  [aclp68_100_contact]
  []
  [aclp69_100_contact]
  []
  [aclp69_101_contact]
  []
  [aclp98_99_contact]
  []
  [aclp99_100_contact]
  []
  [aclp100_101_contact]
  []
  [tlp1_3_contact]
  []
  [tlp3_10_contact]
  []
  [tlp3_11_contact]
  []
  [tlp10_11_contact]
  []
  [tlp10_23_contact]
  []
  [tlp10_24_contact]
  []
  [tlp11_24_contact]
  []
  [tlp23_24_contact]
  []
  [tlp23_42_contact]
  []
  [tlp23_43_contact]
  []
  [tlp24_43_contact]
  []
  [tlp24_44_contact]
  []
  [tlp42_43_contact]
  []
  [tlp42_67_contact]
  []
  [tlp42_68_contact]
  []
  [tlp43_44_contact]
  []
  [tlp43_68_contact]
  []
  [tlp43_69_contact]
  []
  [tlp44_69_contact]
  []
  [tlp67_68_contact]
  []
  [tlp67_98_contact]
  []
  [tlp67_99_contact]
  []
  [tlp68_69_contact]
  []
  [tlp68_99_contact]
  []
  [tlp68_100_contact]
  []
  [tlp69_100_contact]
  []
  [tlp69_101_contact]
  []
  [tlp98_99_contact]
  []
  [tlp99_100_contact]
  []
  [tlp100_101_contact]
  []
  [block1_tlp_contact]
  []
  [block2_tlp_contact]
  []
  [block3_tlp_contact]
  []
  [block1_aclp_contact]
  []
  [block2_aclp_contact]
  []
  [block3_aclp_contact]
  []
  [aclp101_restraint_contact]
  []
  [tlp101_restraint_contact]
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
  [aclp1_3]
    type = PenetrationAux
    variable = aclp1_3_contact
    boundary = 'ACLP1_3'
    paired_boundary = 'ACLP3_1'
    quantity = normal_force_magnitude
  []
  [aclp3_10]
    type = PenetrationAux
    variable = aclp3_10_contact
    boundary = 'ACLP3_10'
    paired_boundary = 'ACLP10_3'
    quantity = normal_force_magnitude
  []
  [aclp3_11]
    type = PenetrationAux
    variable = aclp3_11_contact
    boundary = 'ACLP3_11'
    paired_boundary = 'ACLP11_3'
    quantity = normal_force_magnitude
  []
  [aclp10_11]
    type = PenetrationAux
    variable = aclp10_11_contact
    boundary = 'ACLP10_11'
    paired_boundary = 'ACLP11_10'
    quantity = normal_force_magnitude
  []
  [aclp10_23]
    type = PenetrationAux
    variable = aclp10_23_contact
    boundary = 'ACLP10_23'
    paired_boundary = 'ACLP23_10'
    quantity = normal_force_magnitude
  []
  [aclp10_24]
    type = PenetrationAux
    variable = aclp10_24_contact
    boundary = 'ACLP10_24'
    paired_boundary = 'ACLP24_10'
    quantity = normal_force_magnitude
  []
  [aclp11_24]
    type = PenetrationAux
    variable = aclp11_24_contact
    boundary = 'ACLP11_24'
    paired_boundary = 'ACLP24_11'
    quantity = normal_force_magnitude
  []
  [aclp23_24]
    type = PenetrationAux
    variable = aclp23_24_contact
    boundary = 'ACLP23_24'
    paired_boundary = 'ACLP24_23'
    quantity = normal_force_magnitude
  []
  [aclp23_42]
    type = PenetrationAux
    variable = aclp23_42_contact
    boundary = 'ACLP23_42'
    paired_boundary = 'ACLP42_23'
    quantity = normal_force_magnitude
  []
  [aclp23_43]
    type = PenetrationAux
    variable = aclp23_43_contact
    boundary = 'ACLP23_43'
    paired_boundary = 'ACLP43_23'
    quantity = normal_force_magnitude
  []
  [aclp24_43]
    type = PenetrationAux
    variable = aclp24_43_contact
    boundary = 'ACLP24_43'
    paired_boundary = 'ACLP43_24'
    quantity = normal_force_magnitude
  []
  [aclp24_44]
    type = PenetrationAux
    variable = aclp24_44_contact
    boundary = 'ACLP24_44'
    paired_boundary = 'ACLP44_24'
    quantity = normal_force_magnitude
  []
  [aclp42_43]
    type = PenetrationAux
    variable = aclp42_43_contact
    boundary = 'ACLP42_43'
    paired_boundary = 'ACLP43_42'
    quantity = normal_force_magnitude
  []
  [aclp42_67]
    type = PenetrationAux
    variable = aclp42_67_contact
    boundary = 'ACLP42_67'
    paired_boundary = 'ACLP67_42'
    quantity = normal_force_magnitude
  []
  [aclp42_68]
    type = PenetrationAux
    variable = aclp42_68_contact
    boundary = 'ACLP42_68'
    paired_boundary = 'ACLP68_42'
    quantity = normal_force_magnitude
  []
  [aclp43_44]
    type = PenetrationAux
    variable = aclp43_44_contact
    boundary = 'ACLP43_44'
    paired_boundary = 'ACLP44_43'
    quantity = normal_force_magnitude
  []
  [aclp43_68]
    type = PenetrationAux
    variable = aclp43_68_contact
    boundary = 'ACLP43_68'
    paired_boundary = 'ACLP68_43'
    quantity = normal_force_magnitude
  []
  [aclp43_69]
    type = PenetrationAux
    variable = aclp43_69_contact
    boundary = 'ACLP43_69'
    paired_boundary = 'ACLP69_43'
    quantity = normal_force_magnitude
  []
  [aclp44_69]
    type = PenetrationAux
    variable = aclp44_69_contact
    boundary = 'ACLP44_69'
    paired_boundary = 'ACLP69_44'
    quantity = normal_force_magnitude
  []
  [aclp67_68]
    type = PenetrationAux
    variable = aclp67_68_contact
    boundary = 'ACLP67_68'
    paired_boundary = 'ACLP68_67'
    quantity = normal_force_magnitude
  []
  [aclp67_98]
    type = PenetrationAux
    variable = aclp67_98_contact
    boundary = 'ACLP67_98'
    paired_boundary = 'ACLP98_67'
    quantity = normal_force_magnitude
  []
  [aclp67_99]
    type = PenetrationAux
    variable = aclp67_99_contact
    boundary = 'ACLP67_99'
    paired_boundary = 'ACLP99_67'
    quantity = normal_force_magnitude
  []
  [aclp68_69]
    type = PenetrationAux
    variable = aclp68_69_contact
    boundary = 'ACLP68_69'
    paired_boundary = 'ACLP69_68'
    quantity = normal_force_magnitude
  []
  [aclp68_99]
    type = PenetrationAux
    variable = aclp68_99_contact
    boundary = 'ACLP68_99'
    paired_boundary = 'ACLP99_68'
    quantity = normal_force_magnitude
  []
  [aclp68_100]
    type = PenetrationAux
    variable = aclp68_100_contact
    boundary = 'ACLP68_100'
    paired_boundary = 'ACLP100_68'
    quantity = normal_force_magnitude
  []
  [aclp69_100]
    type = PenetrationAux
    variable = aclp69_100_contact
    boundary = 'ACLP69_100'
    paired_boundary = 'ACLP100_69'
    quantity = normal_force_magnitude
  []
  [aclp69_101]
    type = PenetrationAux
    variable = aclp69_101_contact
    boundary = 'ACLP69_101'
    paired_boundary = 'ACLP101_69'
    quantity = normal_force_magnitude
  []
  [aclp98_99]
    type = PenetrationAux
    variable = aclp98_99_contact
    boundary = 'ACLP98_99'
    paired_boundary = 'ACLP99_98'
    quantity = normal_force_magnitude
  []
  [aclp99_100]
    type = PenetrationAux
    variable = aclp99_100_contact
    boundary = 'ACLP99_100'
    paired_boundary = 'ACLP100_99'
    quantity = normal_force_magnitude
  []
  [aclp100_101]
    type = PenetrationAux
    variable = aclp100_101_contact
    boundary = 'ACLP100_101'
    paired_boundary = 'ACLP101_100'
    quantity = normal_force_magnitude
  []
  [tlp1_3]
    type = PenetrationAux
    variable = tlp1_3_contact
    boundary = 'TLP1_3'
    paired_boundary = 'TLP3_1'
    quantity = normal_force_magnitude
  []
  [tlp3_10]
    type = PenetrationAux
    variable = tlp3_10_contact
    boundary = 'TLP3_10'
    paired_boundary = 'TLP10_3'
    quantity = normal_force_magnitude
  []
  [tlp3_11]
    type = PenetrationAux
    variable = tlp3_11_contact
    boundary = 'TLP3_11'
    paired_boundary = 'TLP11_3'
    quantity = normal_force_magnitude
  []
  [tlp10_11]
    type = PenetrationAux
    variable = tlp10_11_contact
    boundary = 'TLP10_11'
    paired_boundary = 'TLP11_10'
    quantity = normal_force_magnitude
  []
  [tlp10_23]
    type = PenetrationAux
    variable = tlp10_23_contact
    boundary = 'TLP10_23'
    paired_boundary = 'TLP23_10'
    quantity = normal_force_magnitude
  []
  [tlp10_24]
    type = PenetrationAux
    variable = tlp10_24_contact
    boundary = 'TLP10_24'
    paired_boundary = 'TLP24_10'
    quantity = normal_force_magnitude
  []
  [tlp11_24]
    type = PenetrationAux
    variable = tlp11_24_contact
    boundary = 'TLP11_24'
    paired_boundary = 'TLP24_11'
    quantity = normal_force_magnitude
  []
  [tlp23_24]
    type = PenetrationAux
    variable = tlp23_24_contact
    boundary = 'TLP23_24'
    paired_boundary = 'TLP24_23'
    quantity = normal_force_magnitude
  []
  [tlp23_42]
    type = PenetrationAux
    variable = tlp23_42_contact
    boundary = 'TLP23_42'
    paired_boundary = 'TLP42_23'
    quantity = normal_force_magnitude
  []
  [tlp23_43]
    type = PenetrationAux
    variable = tlp23_43_contact
    boundary = 'TLP23_43'
    paired_boundary = 'TLP43_23'
    quantity = normal_force_magnitude
  []
  [tlp24_43]
    type = PenetrationAux
    variable = tlp24_43_contact
    boundary = 'TLP24_43'
    paired_boundary = 'TLP43_24'
    quantity = normal_force_magnitude
  []
  [tlp24_44]
    type = PenetrationAux
    variable = tlp24_44_contact
    boundary = 'TLP24_44'
    paired_boundary = 'TLP44_24'
    quantity = normal_force_magnitude
  []
  [tlp42_43]
    type = PenetrationAux
    variable = tlp42_43_contact
    boundary = 'TLP42_43'
    paired_boundary = 'TLP43_42'
    quantity = normal_force_magnitude
  []
  [tlp42_67]
    type = PenetrationAux
    variable = tlp42_67_contact
    boundary = 'TLP42_67'
    paired_boundary = 'TLP67_42'
    quantity = normal_force_magnitude
  []
  [tlp42_68]
    type = PenetrationAux
    variable = tlp42_68_contact
    boundary = 'TLP42_68'
    paired_boundary = 'TLP68_42'
    quantity = normal_force_magnitude
  []
  [tlp43_44]
    type = PenetrationAux
    variable = tlp43_44_contact
    boundary = 'TLP43_44'
    paired_boundary = 'TLP44_43'
    quantity = normal_force_magnitude
  []
  [tlp43_68]
    type = PenetrationAux
    variable = tlp43_68_contact
    boundary = 'TLP43_68'
    paired_boundary = 'TLP68_43'
    quantity = normal_force_magnitude
  []
  [tlp43_69]
    type = PenetrationAux
    variable = tlp43_69_contact
    boundary = 'TLP43_69'
    paired_boundary = 'TLP69_43'
    quantity = normal_force_magnitude
  []
  [tlp44_69]
    type = PenetrationAux
    variable = tlp44_69_contact
    boundary = 'TLP44_69'
    paired_boundary = 'TLP69_44'
    quantity = normal_force_magnitude
  []
  [tlp67_68]
    type = PenetrationAux
    variable = tlp67_68_contact
    boundary = 'TLP67_68'
    paired_boundary = 'TLP68_67'
    quantity = normal_force_magnitude
  []
  [tlp67_98]
    type = PenetrationAux
    variable = tlp67_98_contact
    boundary = 'TLP67_98'
    paired_boundary = 'TLP98_67'
    quantity = normal_force_magnitude
  []
  [tlp67_99]
    type = PenetrationAux
    variable = tlp67_99_contact
    boundary = 'TLP67_99'
    paired_boundary = 'TLP99_67'
    quantity = normal_force_magnitude
  []
  [tlp68_69]
    type = PenetrationAux
    variable = tlp68_69_contact
    boundary = 'TLP68_69'
    paired_boundary = 'TLP69_68'
    quantity = normal_force_magnitude
  []
  [tlp68_99]
    type = PenetrationAux
    variable = tlp68_99_contact
    boundary = 'TLP68_99'
    paired_boundary = 'TLP99_68'
    quantity = normal_force_magnitude
  []
  [tlp68_100]
    type = PenetrationAux
    variable = tlp68_100_contact
    boundary = 'TLP68_100'
    paired_boundary = 'TLP100_68'
    quantity = normal_force_magnitude
  []
  [tlp69_100]
    type = PenetrationAux
    variable = tlp69_100_contact
    boundary = 'TLP69_100'
    paired_boundary = 'TLP100_69'
    quantity = normal_force_magnitude
  []
  [tlp69_101]
    type = PenetrationAux
    variable = tlp69_101_contact
    boundary = 'TLP69_101'
    paired_boundary = 'TLP101_69'
    quantity = normal_force_magnitude
  []
  [tlp98_99]
    type = PenetrationAux
    variable = tlp98_99_contact
    boundary = 'TLP98_99'
    paired_boundary = 'TLP99_98'
    quantity = normal_force_magnitude
  []
  [tlp99_100]
    type = PenetrationAux
    variable = tlp99_100_contact
    boundary = 'TLP99_100'
    paired_boundary = 'TLP100_99'
    quantity = normal_force_magnitude
  []
  [tlp100_101]
    type = PenetrationAux
    variable = tlp100_101_contact
    boundary = 'TLP100_101'
    paired_boundary = 'TLP101_100'
    quantity = normal_force_magnitude
  []
  [tlp_block1_3]
    type = PenetrationAux
    variable = block1_tlp_contact
    paired_boundary = 'block1_tlp'
    boundary = 'TLP3_4'
    quantity = normal_force_magnitude
  []
  [tlp_block2_24]
    type = PenetrationAux
    variable = block2_tlp_contact
    paired_boundary = 'block2_tlp'
    boundary = 'TLP24_25'
    quantity = normal_force_magnitude
  []
  [tlp_block3_69]
    type = PenetrationAux
    variable = block3_tlp_contact
    paired_boundary = 'block3_tlp'
    boundary = 'TLP69_70'
    quantity = normal_force_magnitude
  []
  [aclp_block1_3]
    type = PenetrationAux
    variable = block1_aclp_contact
    paired_boundary = 'block1_aclp'
    boundary = 'ACLP3_4'
    quantity = normal_force_magnitude
  []
  [aclp_block2_24]
    type = PenetrationAux
    variable = block2_aclp_contact
    paired_boundary = 'block2_aclp'
    boundary = 'ACLP24_25'
    quantity = normal_force_magnitude
  []
  [aclp_block3_69]
    type = PenetrationAux
    variable = block3_aclp_contact
    paired_boundary = 'block3_aclp'
    boundary = 'ACLP69_70'
    quantity = normal_force_magnitude
  []
  [aclp_101_restraint]
    type = PenetrationAux
    variable = aclp101_restraint_contact
    paired_boundary = 'ACLPRR4_2'
    boundary = 'ACLP101_ACLPRR4_2'
    quantity = normal_force_magnitude
  []
  [tlp_101_restraint]
    type = PenetrationAux
    variable = tlp101_restraint_contact
    paired_boundary = 'TLP101_TLPRR4_2'
    boundary = 'TLPRR4_2'
    quantity = normal_force_magnitude
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
    section_name = "Root"
    execute_on = timestep_end
  []
  [aclp_force_1_2_half_force]
    type = NodalSum
	variable = aclp1_3_contact
	boundary = 'ACLP1_3'
  []
  [aclp_force_1_2]
    type = ParsedPostprocessor
	expression = '2*aclp_force_1_2_half_force'
	pp_names = 'aclp_force_1_2_half_force'
  []
  [aclp_force_11_1]
    type = NodalSum
    variable = aclp10_11_contact
    boundary = 'ACLP10_11'
  []
  [aclp_force_11_2]
    type = NodalSum
    variable = aclp11_24_contact
    boundary = 'ACLP11_24'
  []
  [aclp_force_11_6]
    type = NodalSum
    variable = aclp3_11_contact
    boundary = 'ACLP3_11'
  []
  [aclp_force_42_2_half_force]
    type = NodalSum
    variable = aclp42_67_contact
    boundary = 'ACLP42_67'
  []
  [aclp_force_42_2]
    type = ParsedPostprocessor
	expression = '2*aclp_force_42_2_half_force'
	pp_names = 'aclp_force_42_2_half_force'
  []
  [aclp_force_42_3]
    type = NodalSum
    variable = aclp42_68_contact
    boundary = 'ACLP42_68'
  []
  [aclp_force_42_5]
    type = NodalSum
    variable = aclp23_42_contact
    boundary = 'ACLP23_42'
  []
  [aclp_force_43_1]
    type = NodalSum
    variable = aclp42_43_contact
    boundary = 'ACLP42_43'
  []
  [aclp_force_43_2]
    type = NodalSum
    variable = aclp43_68_contact
    boundary = 'ACLP43_68'
  []
  [aclp_force_43_3]
    type = NodalSum
    variable = aclp43_69_contact
    boundary = 'ACLP43_69'
  []
  [aclp_force_43_5]
    type = NodalSum
    variable = aclp24_43_contact
    boundary = 'ACLP24_43'
  []
  [aclp_force_43_6]
    type = NodalSum
    variable = aclp23_43_contact
    boundary = 'ACLP23_43'
  []
  [aclp_force_44_1]
    type = NodalSum
    variable = aclp43_44_contact
    boundary = 'ACLP43_44'
  []
  [aclp_force_44_2]
    type = NodalSum
    variable = aclp44_69_contact
    boundary = 'ACLP44_69'
  []
  [aclp_force_44_6]
    type = NodalSum
    variable = aclp24_44_contact
    boundary = 'ACLP24_44'
  []

  [tlp_force_1_2_half_force]
    type = NodalSum
	variable = tlp1_3_contact
	boundary = 'TLP1_3'
  []
  [tlp_force_1_2]
    type = ParsedPostprocessor
	expression = '2*tlp_force_1_2_half_force'
	pp_names = 'tlp_force_1_2_half_force'
  []
  [tlp_force_11_1]
    type = NodalSum
    variable = tlp10_11_contact
    boundary = 'TLP10_11'
  []
  [tlp_force_11_2]
    type = NodalSum
    variable = tlp11_24_contact
    boundary = 'TLP11_24'
  []
  [tlp_force_11_6]
    type = NodalSum
    variable = tlp3_11_contact
    boundary = 'TLP3_11'
  []
  [tlp_force_42_2_half_force]
    type = NodalSum
    variable = tlp42_67_contact
    boundary = 'TLP42_67'
  []
  [tlp_force_42_2]
    type = ParsedPostprocessor
	expression = '2*tlp_force_42_2_half_force'
	pp_names = 'tlp_force_42_2_half_force'
  []
  [tlp_force_42_3]
    type = NodalSum
    variable = tlp42_68_contact
    boundary = 'TLP42_68'
  []
  [tlp_force_42_5]
    type = NodalSum
    variable = tlp23_42_contact
    boundary = 'TLP23_42'
  []
  [tlp_force_43_1]
    type = NodalSum
    variable = tlp42_43_contact
    boundary = 'TLP42_43'
  []
  [tlp_force_43_2]
    type = NodalSum
    variable = tlp43_68_contact
    boundary = 'TLP43_68'
  []
  [tlp_force_43_3]
    type = NodalSum
    variable = tlp43_69_contact
    boundary = 'TLP43_69'
  []
  [tlp_force_43_5]
    type = NodalSum
    variable = tlp24_43_contact
    boundary = 'TLP24_43'
  []
  [tlp_force_43_6]
    type = NodalSum
    variable = tlp23_43_contact
    boundary = 'TLP23_43'
  []
  [tlp_force_44_1]
    type = NodalSum
    variable = tlp43_44_contact
    boundary = 'TLP43_44'
  []
  [tlp_force_44_2]
    type = NodalSum
    variable = tlp44_69_contact
    boundary = 'TLP44_69'
  []
  [tlp_force_44_6]
    type = NodalSum
    variable = tlp24_44_contact
    boundary = 'TLP24_44'
  []
  [tlp_force_101_1]
    type = NodalSum
    variable = tlp100_101_contact
    boundary = 'TLP100_101'
  []
  [tlp_force_101_2]
    type = NodalSum
    variable = tlp101_restraint_contact
    boundary = 'TLPRR4_2'
  []
  [tlp_force_101_6]
    type = NodalSum
    variable = tlp69_101_contact
    boundary = 'TLP69_101'
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
