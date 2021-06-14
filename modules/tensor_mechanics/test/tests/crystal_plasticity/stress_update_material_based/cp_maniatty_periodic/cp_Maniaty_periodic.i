[Mesh]
   type = GeneratedMesh
   dim = 3
   elem_type = Hex8
   displacements = 'disp_x disp_y disp_z'
   nx = 1
   ny = 1
   nz = 1
[]

[Variables]
  [./disp_x]
  [../]
  [./disp_y]
  [../]
  [./disp_z]
  [../]
[]

[GlobalParams]
  volumetric_locking_correction = true
[]

[AuxVariables]

  [./pk2_zz]
    order = CONSTANT
    family = MONOMIAL
  [../]

  [./vonmises]
  order = CONSTANT
  family = MONOMIAL
  [../]

  [./stress_yy]
    order = CONSTANT
    family = MONOMIAL
  [../]

  [./e_yy]
    order = CONSTANT
    family = MONOMIAL
  [../]

  [./fp_00]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./fp_01]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./fp_02]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./fp_10]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./fp_11]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./fp_12]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./fp_20]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./fp_21]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./fp_22]
    order = CONSTANT
    family = MONOMIAL
  [../]

#  [./fe_00]
#    order = CONSTANT
#    family = MONOMIAL
#  [../]
#  [./fe_01]
#    order = CONSTANT
#    family = MONOMIAL
#  [../]
#  [./fe_02]
#    order = CONSTANT
#    family = MONOMIAL
#  [../]
#  [./fe_10]
#    order = CONSTANT
#    family = MONOMIAL
#  [../]
#  [./fe_11]
#    order = CONSTANT
#    family = MONOMIAL
#  [../]
#  [./fe_12]
#    order = CONSTANT
#    family = MONOMIAL
#  [../]
#  [./fe_20]
#    order = CONSTANT
#    family = MONOMIAL
#  [../]
#  [./fe_21]
#    order = CONSTANT
#    family = MONOMIAL
#  [../]
#  [./fe_22]
#    order = CONSTANT
#    family = MONOMIAL
#  [../]

  [./gss0]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./gss1]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./gss2]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./gss3]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./gss4]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./gss5]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./gss6]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./gss7]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./gss8]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./gss9]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./gss10]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./gss11]
    order = CONSTANT
    family = MONOMIAL
  [../]

  [./euler1]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./euler2]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./euler3]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./stress_zz]
    order = CONSTANT
    family = MONOMIAL 
  [../]
  [./e_zz]
    order = CONSTANT
    family = MONOMIAL    
  [../]

  [./rot00]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./rot01]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./rot02]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./rot10]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./rot11]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./rot12]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./rot20]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./rot21]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./rot22]
    order = CONSTANT
    family = MONOMIAL
  [../]

#  [./eqv_slip_incr00]
#    order = CONSTANT
#    family = MONOMIAL
#  [../]
#  [./eqv_slip_incr01]
#    order = CONSTANT
#    family = MONOMIAL
#  [../]
#  [./eqv_slip_incr02]
#    order = CONSTANT
#    family = MONOMIAL
#  [../]
#  [./eqv_slip_incr10]
#    order = CONSTANT
#    family = MONOMIAL
#  [../]
#  [./eqv_slip_incr11]
#    order = CONSTANT
#    family = MONOMIAL
#  [../]
#  [./eqv_slip_incr12]
#    order = CONSTANT
#    family = MONOMIAL
#  [../]
#  [./eqv_slip_incr20]
#    order = CONSTANT
#    family = MONOMIAL
#  [../]
#  [./eqv_slip_incr21]
#    order = CONSTANT
#    family = MONOMIAL
#  [../]
#  [./eqv_slip_incr22]
#    order = CONSTANT
#    family = MONOMIAL
#  [../]

  [./tau0]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./tau1]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./tau2]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./tau3]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./tau4]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./tau5]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./tau6]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./tau7]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./tau8]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./tau9]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./tau10]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./tau11]
    order = CONSTANT
    family = MONOMIAL
  [../]

#  [./hab00]
#    order = CONSTANT
#    family = MONOMIAL
#  [../]

  [./hb0]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./hb1]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./hb2]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./hb3]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./hb4]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./hb5]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./hb6]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./hb7]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./hb8]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./hb9]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./hb10]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./hb11]
    order = CONSTANT
    family = MONOMIAL
  [../]

  [./val0]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./val1]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./val2]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./val3]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./val4]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./val5]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./val6]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./val7]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./val8]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./val9]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./val10]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./val11]
    order = CONSTANT
    family = MONOMIAL
  [../]

  [./crysrot00]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./crysrot01]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./crysrot02]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./crysrot10]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./crysrot11]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./crysrot12]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./crysrot20]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./crysrot21]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./crysrot22]
    order = CONSTANT
    family = MONOMIAL
  [../]
[]

[Functions]
  [./tdisp]
    type = ParsedFunction
    value = 0.0001*t
  [../]
[]

[UserObjects]
  [./prop_read]
    type = ElementPropertyReadFile
    prop_file_name = 'euler_112_1.txt'
    # Enter file data as prop#1, prop#2, .., prop#nprop
    nprop = 3
    execute_on = timestep_end
    read_type = element
  [../]
[]

[AuxKernels]

  [./vonmises]
  type = RankTwoScalarAux
  rank_two_tensor = stress
  variable = vonmises
  scalar_type = VonMisesStress
  execute_on = timestep_end
  [../]

  [./stress_yy]
    type = RankTwoAux
    variable = stress_yy
    rank_two_tensor = stress
    index_j = 1
    index_i = 1
    execute_on = timestep_end
  [../]
  [./e_yy]
    type = RankTwoAux
    variable = e_yy
    rank_two_tensor = lage
    index_j = 1
    index_i = 1
    execute_on = timestep_end
  [../]
  [./stress_zz]
    type = RankTwoAux
    variable =stress_zz
    rank_two_tensor = stress
    index_j = 2
    index_i = 2
    execute_on = timestep_end    
  [../]
 
  [./e_zz]
    type = RankTwoAux
    variable = e_zz
    rank_two_tensor = lage
    index_j = 2
    index_i = 2
    execute_on = timestep_end
    
  [../]

  [./fp_00]
    type = RankTwoAux
    variable = fp_00
    rank_two_tensor = fp
    index_j = 0
    index_i = 0
    execute_on = timestep_end
  [../]

  [./fp_01]
    type = RankTwoAux
    variable = fp_01
    rank_two_tensor = fp
    index_j = 1
    index_i = 0
    execute_on = timestep_end
  [../]

  [./fp_02]
    type = RankTwoAux
    variable = fp_02
    rank_two_tensor = fp
    index_j = 2
    index_i = 0
    execute_on = timestep_end
  [../]

  [./fp_10]
    type = RankTwoAux
    variable = fp_10
    rank_two_tensor = fp
    index_j = 0
    index_i = 1
    execute_on = timestep_end
  [../]

  [./fp_11]
    type = RankTwoAux
    variable = fp_11
    rank_two_tensor = fp
    index_j = 1
    index_i = 1
    execute_on = timestep_end
  [../]

  [./fp_12]
    type = RankTwoAux
    variable = fp_12
    rank_two_tensor = fp
    index_j = 2
    index_i = 1
    execute_on = timestep_end
  [../]

  [./fp_20]
    type = RankTwoAux
    variable = fp_20
    rank_two_tensor = fp
    index_j = 0
    index_i = 2
    execute_on = timestep_end
  [../]

  [./fp_21]
    type = RankTwoAux
    variable = fp_21
    rank_two_tensor = fp
    index_j = 1
    index_i = 2
    execute_on = timestep_end
  [../]

  [./fp_22]
    type = RankTwoAux
    variable = fp_22
    rank_two_tensor = fp
    index_j = 2
    index_i = 2
    execute_on = timestep_end
  [../]

#  [./eqv_slip_incr00]
#    type = RankTwoAux
#    variable = eqv_slip_incr00
#    rank_two_tensor = eqv_slip_incr
#    index_j = 0
#    index_i = 0
#    execute_on = timestep_end
#  [../]
#  [./eqv_slip_incr01]
#    type = RankTwoAux
#    variable = eqv_slip_incr01
#    rank_two_tensor = eqv_slip_incr
#    index_j = 1
#    index_i = 0
#    execute_on = timestep_end
#  [../]
#  [./eqv_slip_incr02]
#    type = RankTwoAux
#    variable = eqv_slip_incr02
#    rank_two_tensor = eqv_slip_incr
#    index_j = 2
#    index_i = 0
#    execute_on = timestep_end
#  [../]
#  [./eqv_slip_incr10]
#    type = RankTwoAux
#    variable = eqv_slip_incr10
#    rank_two_tensor = eqv_slip_incr
#    index_j = 0
#    index_i = 1
#    execute_on = timestep_end
#  [../]
#  [./eqv_slip_incr11]
#    type = RankTwoAux
#    variable = eqv_slip_incr11
#    rank_two_tensor = eqv_slip_incr
#    index_j = 1
#    index_i = 1
#    execute_on = timestep_end
#  [../]
#  [./eqv_slip_incr12]
#    type = RankTwoAux
#    variable = eqv_slip_incr12
#    rank_two_tensor = eqv_slip_incr
#    index_j = 2
#    index_i = 1
#    execute_on = timestep_end
#  [../]
#  [./eqv_slip_incr20]
#    type = RankTwoAux
#    variable = eqv_slip_incr20
#    rank_two_tensor = eqv_slip_incr
#    index_j = 0
#    index_i = 2
#    execute_on = timestep_end
#  [../]
#  [./eqv_slip_incr21]
#    type = RankTwoAux
#    variable = eqv_slip_incr21
#    rank_two_tensor = eqv_slip_incr
#    index_j = 1
#    index_i = 2
#    execute_on = timestep_end
#  [../]
#  [./eqv_slip_incr22]
#    type = RankTwoAux
#    variable = eqv_slip_incr22
#    rank_two_tensor = eqv_slip_incr
#    index_j = 2
#    index_i = 2
#    execute_on = timestep_end
#  [../]

#  [./fe_00]
#    type = RankTwoAux
#    variable = fe_00
#    rank_two_tensor = fe
#    index_j = 0
#    index_i = 0
#    execute_on = timestep_end
#  [../]

#  [./fe_01]
#    type = RankTwoAux
#    variable = fe_01
#    rank_two_tensor = fe
#    index_j = 1
#    index_i = 0
#    execute_on = timestep_end
#  [../]

#  [./fe_02]
#    type = RankTwoAux
#    variable = fe_02
#    rank_two_tensor = fe
#    index_j = 2
#    index_i = 0
#    execute_on = timestep_end
#  [../]

#  [./fe_10]
#    type = RankTwoAux
#    variable = fe_10
#    rank_two_tensor = fe
#    index_j = 0
#    index_i = 1
#    execute_on = timestep_end
#  [../]

#  [./fe_11]
#    type = RankTwoAux
#    variable = fe_11
#    rank_two_tensor = fe
#    index_j = 1
#    index_i = 1
#    execute_on = timestep_end
#  [../]

#  [./fe_12]
#    type = RankTwoAux
#    variable = fe_12
#    rank_two_tensor = fe
#    index_j = 2
#    index_i = 1
#    execute_on = timestep_end
#  [../]

#  [./fe_20]
#    type = RankTwoAux
#    variable = fe_20
#    rank_two_tensor = fe
#    index_j = 0
#    index_i = 2
#    execute_on = timestep_end
#  [../]

#  [./fe_21]
#    type = RankTwoAux
#    variable = fe_21
#    rank_two_tensor = fe
#    index_j = 1
#    index_i = 2
#    execute_on = timestep_end
#  [../]

#  [./fe_22]
#    type = RankTwoAux
#    variable = fe_22
#    rank_two_tensor = fe
#    index_j = 2
#    index_i = 2
#    execute_on = timestep_end
#  [../]

  [./tau0]
    type = MaterialStdVectorAux
    variable = tau0
    property = slip_rate_gss
    index = 0
    execute_on = timestep_end
  [../]
  [./tau1]
    type = MaterialStdVectorAux
    variable = tau1
    property = slip_rate_gss
    index = 1
    execute_on = timestep_end
  [../]
  [./tau2]
    type = MaterialStdVectorAux
    variable = tau2
    property = slip_rate_gss
    index = 2
    execute_on = timestep_end
  [../]
  [./tau3]
    type = MaterialStdVectorAux
    variable = tau3
    property = slip_rate_gss
    index = 3
    execute_on = timestep_end
  [../]
  [./tau4]
    type = MaterialStdVectorAux
    variable = tau4
    property = slip_rate_gss
    index = 4
    execute_on = timestep_end
  [../]
  [./tau5]
    type = MaterialStdVectorAux
    variable = tau5
    property = slip_rate_gss
    index = 5
    execute_on = timestep_end
  [../]
  [./tau6]
    type = MaterialStdVectorAux
    variable = tau6
    property = slip_rate_gss
    index = 6
    execute_on = timestep_end
  [../]
  [./tau7]
    type = MaterialStdVectorAux
    variable = tau7
    property = slip_rate_gss
    index = 7
    execute_on = timestep_end
  [../]
  [./tau8]
    type = MaterialStdVectorAux
    variable = tau8
    property = slip_rate_gss
    index = 8
    execute_on = timestep_end
  [../]
  [./tau9]
    type = MaterialStdVectorAux
    variable = tau9
    property = slip_rate_gss
    index = 9
    execute_on = timestep_end
  [../]
  [./tau10]
    type = MaterialStdVectorAux
    variable = tau10
    property = slip_rate_gss
    index = 10
    execute_on = timestep_end
  [../]
  [./tau11]
    type = MaterialStdVectorAux
    variable = tau11
    property = slip_rate_gss
    index = 11
    execute_on = timestep_end
  [../]

#  [./hab00]
#    type = RankTwoAux
#    variable = hab00
#    rank_two_tensor = hab
#    index_j = 0
#    index_i = 0
#    execute_on = timestep_end
#  [../]

  [./hb0]
    type = MaterialStdVectorAux
    variable = hb0
    property = state_var_evol_rate_comp_gss
    index = 0
    execute_on = timestep_end
  [../]
  [./hb1]
    type = MaterialStdVectorAux
    variable = hb1
    property = state_var_evol_rate_comp_gss
    index = 1
    execute_on = timestep_end
  [../]
  [./hb2]
    type = MaterialStdVectorAux
    variable = hb2
    property = state_var_evol_rate_comp_gss
    index = 2
    execute_on = timestep_end
  [../]
  [./hb3]
    type = MaterialStdVectorAux
    variable = hb3
    property = state_var_evol_rate_comp_gss
    index = 3
    execute_on = timestep_end
  [../]
  [./hb4]
    type = MaterialStdVectorAux
    variable = hb4
    property = state_var_evol_rate_comp_gss
    index = 4
    execute_on = timestep_end
  [../]
  [./hb5]
    type = MaterialStdVectorAux
    variable = hb5
    property = state_var_evol_rate_comp_gss
    index = 5
    execute_on = timestep_end
  [../]
  [./hb6]
    type = MaterialStdVectorAux
    variable = hb6
    property = state_var_evol_rate_comp_gss
    index = 6
    execute_on = timestep_end
  [../]
  [./hb7]
    type = MaterialStdVectorAux
    variable = hb7
    property = state_var_evol_rate_comp_gss
    index = 7
    execute_on = timestep_end
  [../]
  [./hb8]
    type = MaterialStdVectorAux
    variable = hb8
    property = state_var_evol_rate_comp_gss
    index = 8
    execute_on = timestep_end
  [../]
  [./hb9]
    type = MaterialStdVectorAux
    variable = hb9
    property = state_var_evol_rate_comp_gss
    index = 9
    execute_on = timestep_end
  [../]
  [./hb10]
    type = MaterialStdVectorAux
    variable = hb10
    property = state_var_evol_rate_comp_gss
    index = 10
    execute_on = timestep_end
  [../]
  [./hb11]
    type = MaterialStdVectorAux
    variable = hb11
    property = state_var_evol_rate_comp_gss
    index = 11
    execute_on = timestep_end
  [../]

  [./val0]
    type = MaterialStdVectorAux
    variable = val0
    property = slip_rate_gss
    index = 0
    execute_on = timestep_end
  [../]
  [./val1]
    type = MaterialStdVectorAux
    variable = val1
    property = slip_rate_gss
    index = 1
    execute_on = timestep_end
  [../]
  [./val2]
    type = MaterialStdVectorAux
    variable = val2
    property = slip_rate_gss
    index = 2
    execute_on = timestep_end
  [../]
  [./val3]
    type = MaterialStdVectorAux
    variable = val3
    property = slip_rate_gss
    index = 3
    execute_on = timestep_end
  [../]
  [./val4]
    type = MaterialStdVectorAux
    variable = val4
    property = slip_rate_gss
    index = 4
    execute_on = timestep_end
  [../]
  [./val5]
    type = MaterialStdVectorAux
    variable = val5
    property = slip_rate_gss
    index = 5
    execute_on = timestep_end
  [../]
  [./val6]
    type = MaterialStdVectorAux
    variable = val6
    property = slip_rate_gss
    index = 6
    execute_on = timestep_end
  [../]
  [./val7]
    type = MaterialStdVectorAux
    variable = val7
    property = slip_rate_gss
    index = 7
    execute_on = timestep_end
  [../]
  [./val8]
    type = MaterialStdVectorAux
    variable = val8
    property = slip_rate_gss
    index = 8
    execute_on = timestep_end
  [../]
  [./val9]
    type = MaterialStdVectorAux
    variable = val9
    property = slip_rate_gss
    index = 9
    execute_on = timestep_end
  [../]
  [./val10]
    type = MaterialStdVectorAux
    variable = val10
    property = slip_rate_gss
    index = 10
    execute_on = timestep_end
  [../]
  [./val11]
    type = MaterialStdVectorAux
    variable = val11
    property = slip_rate_gss
    index = 11
    execute_on = timestep_end
  [../]


  [./gss0]
    type = MaterialStdVectorAux
    variable = gss0
    property = state_var_gss
    index = 0
    execute_on = timestep_end
  [../]
  [./gss1]
    type = MaterialStdVectorAux
    variable = gss1
    property = state_var_gss
    index = 1
    execute_on = timestep_end
  [../]
  [./gss2]
    type = MaterialStdVectorAux
    variable = gss2
    property = state_var_gss
    index = 2
    execute_on = timestep_end
  [../]
  [./gss3]
    type = MaterialStdVectorAux
    variable = gss3
    property = state_var_gss
    index = 3
    execute_on = timestep_end
  [../]
  [./gss4]
    type = MaterialStdVectorAux
    variable = gss4
    property = state_var_gss
    index = 4
    execute_on = timestep_end
  [../]
  [./gss5]
    type = MaterialStdVectorAux
    variable = gss5
    property = state_var_gss
    index = 5
    execute_on = timestep_end
  [../]
  [./gss6]
    type = MaterialStdVectorAux
    variable = gss6
    property = state_var_gss
    index = 6
    execute_on = timestep_end
  [../]
  [./gss7]
    type = MaterialStdVectorAux
    variable = gss7
    property = state_var_gss
    index = 7
    execute_on = timestep_end
  [../]
  [./gss8]
    type = MaterialStdVectorAux
    variable = gss8
    property = state_var_gss
    index = 8
    execute_on = timestep_end
  [../]
  [./gss9]
    type = MaterialStdVectorAux
    variable = gss9
    property = state_var_gss
    index = 9
    execute_on = timestep_end
  [../]
  [./gss10]
    type = MaterialStdVectorAux
    variable = gss10
    property = state_var_gss
    index = 10
    execute_on = timestep_end
  [../]
  [./gss11]
    type = MaterialStdVectorAux
    variable = gss11
    property = state_var_gss
    index = 11
    execute_on = timestep_end
  [../]

  [./euler1]
    type = MaterialRealVectorValueAux
    variable = euler1
    property = Euler_angles
    component = 0
    execute_on = timestep_end
  [../]
  [./euler2]
    type = MaterialRealVectorValueAux
    variable = euler2
    property = Euler_angles
    component = 1
    execute_on = timestep_end
  [../]
  [./euler3]
    type = MaterialRealVectorValueAux
    variable = euler3
    property = Euler_angles
    component = 2
    execute_on = timestep_end
  [../]
  [./pk2_zz]
    type = RankTwoAux
    variable = pk2_zz
    rank_two_tensor = pk2
    index_j = 2
    index_i = 2
    execute_on = timestep_end
  [../]

  [./rot00]
    type = RankTwoAux
    variable = rot00
    rank_two_tensor = update_rot
    index_j = 0
    index_i = 0
    execute_on = timestep_end
  [../]
  [./rot01]
    type = RankTwoAux
    variable = rot01
    rank_two_tensor = update_rot
    index_j = 1
    index_i = 0
    execute_on = timestep_end
  [../]
  [./rot02]
    type = RankTwoAux
    variable = rot02
    rank_two_tensor = update_rot
    index_j = 2
    index_i = 0
    execute_on = timestep_end
  [../]
  [./rot10]
    type = RankTwoAux
    variable = rot10
    rank_two_tensor = update_rot
    index_j = 0
    index_i = 1
    execute_on = timestep_end
  [../]
  [./rot11]
    type = RankTwoAux
    variable = rot11
    rank_two_tensor = update_rot
    index_j = 1
    index_i = 1
    execute_on = timestep_end
  [../]
  [./rot12]
    type = RankTwoAux
    variable = rot12
    rank_two_tensor = update_rot
    index_j = 2
    index_i = 1
    execute_on = timestep_end
  [../]
  [./rot20]
    type = RankTwoAux
    variable = rot20
    rank_two_tensor = update_rot
    index_j = 0
    index_i = 2
    execute_on = timestep_end
  [../]
  [./rot21]
    type = RankTwoAux
    variable = rot21
    rank_two_tensor = update_rot
    index_j = 1
    index_i = 2
    execute_on = timestep_end
  [../]
  [./rot22]
    type = RankTwoAux
    variable = rot22
    rank_two_tensor = update_rot
    index_j = 2
    index_i = 2
    execute_on = timestep_end
  [../]

  [./crysrot00]
    type = RankTwoAux
    variable = crysrot00
    rank_two_tensor = crysrot
    index_j = 0
    index_i = 0
    execute_on = timestep_end
  [../]
  [./crysrot01]
    type = RankTwoAux
    variable = crysrot01
    rank_two_tensor = crysrot
    index_j = 1
    index_i = 0
    execute_on = timestep_end
  [../]
  [./crysrot02]
    type = RankTwoAux
    variable = crysrot02
    rank_two_tensor = crysrot
    index_j = 2
    index_i = 0
    execute_on = timestep_end
  [../]
  [./crysrot10]
    type = RankTwoAux
    variable = crysrot10
    rank_two_tensor = crysrot
    index_j = 0
    index_i = 1
    execute_on = timestep_end
  [../]
  [./crysrot11]
    type = RankTwoAux
    variable = crysrot11
    rank_two_tensor = crysrot
    index_j = 1
    index_i = 1
    execute_on = timestep_end
  [../]
  [./crysrot12]
    type = RankTwoAux
    variable = crysrot12
    rank_two_tensor = crysrot
    index_j = 2
    index_i = 1
    execute_on = timestep_end
  [../]
  [./crysrot20]
    type = RankTwoAux
    variable = crysrot20
    rank_two_tensor = crysrot
    index_j = 0
    index_i = 2
    execute_on = timestep_end
  [../]
  [./crysrot21]
    type = RankTwoAux
    variable = crysrot21
    rank_two_tensor = crysrot
    index_j = 1
    index_i = 2
    execute_on = timestep_end
  [../]
  [./crysrot22]
    type = RankTwoAux
    variable = crysrot22
    rank_two_tensor = crysrot
    index_j = 2
    index_i = 2
    execute_on = timestep_end
  [../]
[]

[BCs]

  [./Periodic]
    [./all]
      variable = 'disp_x'
      auto_direction = 'z'
    [../]
  [../]

  [./fix_x]
    type = DirichletBC
    variable = disp_x
    boundary = 'left'
    value = 0
  [../]

  [./fix_y]
    type = DirichletBC
    variable = disp_y
    boundary = 'bottom'
    value = 0
  [../]

  [./fix_z]
    type = DirichletBC
    variable = disp_z
    boundary = 'back'
    value = 0
  [../]

  [./tdisp]
    type = FunctionDirichletBC
    variable = disp_z
    boundary = 'front'
    function = tdisp
  [../]
[]

[UserObjects]
  [./slip_rate_gss]
    type = CrystalPlasticitySlipRateGSS
    variable_size = 12
    slip_sys_file_name = slip_sys_fcc.txt
    num_slip_sys_flowrate_props = 2
    flowprops = '1 12 0.001 0.05'
    uo_state_var_name = state_var_gss
  [../]
  [./slip_resistance_gss]
    type = CrystalPlasticitySlipResistanceGSS
    variable_size = 12
    uo_state_var_name = state_var_gss
  [../]
  [./state_var_gss]
    type = CrystalPlasticityStateVariable
    variable_size = 12
    groups = '0 12'
    group_values = '31'
    uo_state_var_evol_rate_comp_name = state_var_evol_rate_comp_gss
    scale_factor = 1.0
  [../]
  [./state_var_evol_rate_comp_gss]
    type = CrystalPlasticityStateVarRateComponentGSS
    variable_size = 12
    hprops = '1.2 75 63 2.25'
    uo_slip_rate_name = slip_rate_gss
    uo_state_var_name = state_var_gss
  [../]
[]

[Materials]
  [./crysp]
    type = FiniteStrainUObasedCP
    stol = 1e-2
    tan_mod_type = exact
    uo_slip_rates = 'slip_rate_gss'
    uo_slip_resistances = 'slip_resistance_gss'
    uo_state_vars = 'state_var_gss'
    uo_state_var_evol_rate_comps = 'state_var_evol_rate_comp_gss'
#    read_prop_user_object = prop_read   #RC 8th September, 2020
  [../]
  [./strain]
    type = ComputeFiniteStrain
    displacements = 'disp_x disp_y disp_z'
  [../]
  [./elasticity_tensor]
    type = ComputeElasticityTensorCP
    C_ijkl = '1.0675e5 0.6041e5 0.6041e5 1.0675e5 0.6041e5 1.0675e5 0.2834e5 0.2834e5 0.2834e5'
    fill_method = symmetric9
    read_prop_user_object = prop_read
  [../]
[]

[Postprocessors]

  [./stress_zz]
    type = ElementAverageValue
    variable = stress_zz
  [../]
  [./e_zz]
    type = ElementAverageValue
    variable = e_zz
  [../]
  [./vonMises]
    type = ElementAverageValue
    variable = vonmises
  [../]
  [./stress_yy]
    type = ElementAverageValue
    variable = stress_yy
  [../]
  [./e_yy]
    type = ElementAverageValue
    variable = e_yy
  [../]

  [./fp_00]
    type = ElementAverageValue
    variable = fp_00
  [../]

  [./fp_01]
    type = ElementAverageValue
    variable = fp_01
  [../]

  [./fp_02]
    type = ElementAverageValue
    variable = fp_02
  [../]

  [./fp_10]
    type = ElementAverageValue
    variable = fp_10
  [../]

  [./fp_11]
    type = ElementAverageValue
    variable = fp_11
  [../]

  [./fp_12]
    type = ElementAverageValue
    variable = fp_12
  [../]

  [./fp_20]
    type = ElementAverageValue
    variable = fp_20
  [../]

  [./fp_21]
    type = ElementAverageValue
    variable = fp_21
  [../]

  [./fp_22]
    type = ElementAverageValue
    variable = fp_22
  [../]

#  [./eqv_slip_incr00]
#    type = ElementAverageValue
#    variable = eqv_slip_incr00
#  [../]
#  [./eqv_slip_incr01]
#    type = ElementAverageValue
#    variable = eqv_slip_incr01
#  [../]
#  [./eqv_slip_incr02]
#    type = ElementAverageValue
#    variable = eqv_slip_incr02
#  [../]
#  [./eqv_slip_incr10]
#    type = ElementAverageValue
#    variable = eqv_slip_incr10
#  [../]
#  [./eqv_slip_incr11]
#    type = ElementAverageValue
#    variable = eqv_slip_incr11
#  [../]
#  [./eqv_slip_incr12]
#    type = ElementAverageValue
#    variable = eqv_slip_incr12
#  [../]
#  [./eqv_slip_incr20]
#    type = ElementAverageValue
#    variable = eqv_slip_incr20
#  [../]
#  [./eqv_slip_incr21]
#    type = ElementAverageValue
#    variable = eqv_slip_incr21
#  [../]
#  [./eqv_slip_incr22]
#    type = ElementAverageValue
#    variable = eqv_slip_incr22
#  [../]

  [./tau0]
    type = ElementAverageValue
    variable = tau0
  [../]
  [./tau1]
    type = ElementAverageValue
    variable = tau1
  [../]
  [./tau2]
    type = ElementAverageValue
    variable = tau2
  [../]
  [./tau3]
    type = ElementAverageValue
    variable = tau3
  [../]
  [./tau4]
    type = ElementAverageValue
    variable = tau4
  [../]
  [./tau5]
    type = ElementAverageValue
    variable = tau5
  [../]
  [./tau6]
    type = ElementAverageValue
    variable = tau6
  [../]
  [./tau7]
    type = ElementAverageValue
    variable = tau7
  [../]
  [./tau8]
    type = ElementAverageValue
    variable = tau8
  [../]
  [./tau9]
    type = ElementAverageValue
    variable = tau9
  [../]
  [./tau10]
    type = ElementAverageValue
    variable = tau10
  [../]
  [./tau11]
    type = ElementAverageValue
    variable = tau11
  [../]

#  [./hab00]
#    type = ElementAverageValue
#    variable = hab00
#  [../]

  [./hb0]
    type = ElementAverageValue
    variable = hb0
  [../]
  [./hb1]
    type = ElementAverageValue
    variable = hb1
  [../]
  [./hb2]
    type = ElementAverageValue
    variable = hb2
  [../]
  [./hb3]
    type = ElementAverageValue
    variable = hb3
  [../]
  [./hb4]
    type = ElementAverageValue
    variable = hb4
  [../]
  [./hb5]
    type = ElementAverageValue
    variable = hb5
  [../]
  [./hb6]
    type = ElementAverageValue
    variable = hb6
  [../]
  [./hb7]
    type = ElementAverageValue
    variable = hb7
  [../]
  [./hb8]
    type = ElementAverageValue
    variable = hb8
  [../]
  [./hb9]
    type = ElementAverageValue
    variable = hb9
  [../]
  [./hb10]
    type = ElementAverageValue
    variable = hb10
  [../]
  [./hb11]
    type = ElementAverageValue
    variable = hb11
  [../]

  [./val0]
    type = ElementAverageValue
    variable = val0
  [../]
  [./val1]
    type = ElementAverageValue
    variable = val1
  [../]
  [./val2]
    type = ElementAverageValue
    variable = val2
  [../]
  [./val3]
    type = ElementAverageValue
    variable = val3
  [../]
  [./val4]
    type = ElementAverageValue
    variable = val4
  [../]
  [./val5]
    type = ElementAverageValue
    variable = val5
  [../]
  [./val6]
    type = ElementAverageValue
    variable = val6
  [../]
  [./val7]
    type = ElementAverageValue
    variable = val7
  [../]
  [./val8]
    type = ElementAverageValue
    variable = val8
  [../]
  [./val9]
    type = ElementAverageValue
    variable = val9
  [../]
  [./val10]
    type = ElementAverageValue
    variable = val10
  [../]
  [./val11]
    type = ElementAverageValue
    variable = val11
  [../]

#  [./fe_00]
#    type = ElementAverageValue
#    variable = fe_00
#  [../]

#  [./fe_01]
#    type = ElementAverageValue
#    variable = fe_01
#  [../]

#  [./fe_02]
#    type = ElementAverageValue
#    variable = fe_02
#  [../]

#  [./fe_10]
#    type = ElementAverageValue
#    variable = fe_10
#  [../]

#  [./fe_11]
#    type = ElementAverageValue
#    variable = fe_11
#  [../]

#  [./fe_12]
#    type = ElementAverageValue
#    variable = fe_12
#  [../]

#  [./fe_20]
#    type = ElementAverageValue
#    variable = fe_20
#  [../]

#  [./fe_21]
#    type = ElementAverageValue
#    variable = fe_21
#  [../]

#  [./fe_22]
#    type = ElementAverageValue
#    variable = fe_22
#  [../]

  [./gss0]
    type = ElementAverageValue
    variable = gss0
  [../]
  [./gss1]
    type = ElementAverageValue
    variable = gss1
  [../]
  [./gss2]
    type = ElementAverageValue
    variable = gss2
  [../]
  [./gss3]
    type = ElementAverageValue
    variable = gss3
  [../]
  [./gss4]
    type = ElementAverageValue
    variable = gss4
  [../]
  [./gss5]
    type = ElementAverageValue
    variable = gss5
  [../]
  [./gss6]
    type = ElementAverageValue
    variable = gss6
  [../]
  [./gss7]
    type = ElementAverageValue
    variable = gss7
  [../]
  [./gss8]
    type = ElementAverageValue
    variable = gss8
  [../]
  [./gss9]
    type = ElementAverageValue
    variable = gss9
  [../]
  [./gss10]
    type = ElementAverageValue
    variable = gss10
  [../]
  [./gss11]
    type = ElementAverageValue
    variable = gss11
  [../]

  [./pk2_zz]
    type = ElementAverageValue
    variable = pk2_zz
  [../]

  [./rot00]
    type = ElementAverageValue
    variable = rot00
  [../]
  [./rot01]
    type = ElementAverageValue
    variable = rot01
  [../]
  [./rot02]
    type = ElementAverageValue
    variable = rot02
  [../]
  [./rot10]
    type = ElementAverageValue
    variable = rot10
  [../]
  [./rot11]
    type = ElementAverageValue
    variable = rot11
  [../]
  [./rot12]
    type = ElementAverageValue
    variable = rot12
  [../]
  [./rot20]
    type = ElementAverageValue
    variable = rot20
  [../]
  [./rot21]
    type = ElementAverageValue
    variable = rot21
  [../]
  [./rot22]
    type = ElementAverageValue
    variable = rot22
  [../]

  [./crysrot00]
    type = ElementAverageValue
    variable = crysrot00
  [../]
  [./crysrot01]
    type = ElementAverageValue
    variable = crysrot01
  [../]
  [./crysrot02]
    type = ElementAverageValue
    variable = crysrot02
  [../]
  [./crysrot10]
    type = ElementAverageValue
    variable = crysrot10
  [../]
  [./crysrot11]
    type = ElementAverageValue
    variable = crysrot11
  [../]
  [./crysrot12]
    type = ElementAverageValue
    variable = crysrot12
  [../]
  [./crysrot20]
    type = ElementAverageValue
    variable = crysrot20
  [../]
  [./crysrot21]
    type = ElementAverageValue
    variable = crysrot21
  [../]
  [./crysrot22]
    type = ElementAverageValue
    variable = crysrot22
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
  dt = 0.5
  solve_type = 'PJFNK'

  petsc_options_iname = -pc_hypre_type
  petsc_options_value = boomerang
  nl_abs_tol = 1e-10
  nl_rel_step_tol = 1e-10
  dtmax = 10.0
  nl_rel_tol = 1e-10
#  end_time = 1
  dtmin = 0.0000001
#  num_steps = 12260
  num_steps = 1000
  nl_abs_step_tol = 1e-10
[]

[Outputs]
  exodus = true
  csv = true
  interval = 10
[]

[Kernels]
  [./TensorMechanics]
    displacements = 'disp_x disp_y disp_z'
    use_displaced_mesh = true
  [../]
[]
