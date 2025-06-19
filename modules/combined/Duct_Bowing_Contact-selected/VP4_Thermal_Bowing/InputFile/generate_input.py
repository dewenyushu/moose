import numpy as np
import matplotlib.pyplot as plt

# direction substrate
directions = ["_x", "_y", "_z"]

locations = ["tlp", "aclp"]

# block boundaries
block_bnds = ["3_4", "24_25", "69_70"]

# aclp and tlp contact names
contact_names = ["1_3", "3_10", "3_11", "10_11", "10_23", "10_24", "11_24", "23_24", "23_42", "23_43", "24_43", "24_44", "42_43", "42_67", "42_68", "43_44", "43_68", "43_69", "44_69", "67_68", "67_98", "67_99", "68_69", "68_99", "68_100", "69_100", "69_101", "98_99", "99_100", "100_101"]

# for creating postprocessors, contact boundary name : pp name
pp_names = {"1_3":"1_2", "10_11" : "11_1", "11_24":"11_2", "3_11":"11_6", "42_67":"42_2", "42_68":"42_3", "23_42":"42_5", "42_43":"43_1", "43_68":"43_2", "43_69":"43_3", "24_43":"43_5", "23_43":"43_6", "43_44":"44_1", "44_69":"44_2", "24_44":"44_6", "100_101":"101_1", "4_2":"101_2", "69_101":"101_6"}
half_force_set = {}

# create AuxVariables
with open("AuxVariables.txt", "w") as file:
    file.write("\t[AuxVariables]\n")
    file.write("\t\t[temp]\n")
    file.write("\t\t\tinitial_condition = 400\n")
    file.write("\t\t[]\n")
    file.write("\t\t[r]\n")
    file.write("\t\t[]\n")
    # aclp
    for name in contact_names:
        for dir in directions:
          file.write("\t\t[aclp"+name+"_contact"+dir+"]\n")
          file.write("\t\t[]\n")
    # tlp
    for name in contact_names:
        for dir in directions:
          file.write("\t\t[tlp"+name+"_contact"+dir+"]\n")
          file.write("\t\t[]\n")
    # block
    for loc in locations:
       for n in np.arange(1,4):
          for dir in directions:
             file.write("\t\t[block"+str(n)+"_"+loc+"_contact"+dir+"]\n")
             file.write("\t\t[]\n")

    # restraint
    for loc in locations:
       for dir in directions:
        file.write("\t\t["+loc+"101_restraint_contact"+dir+"]\n")
        file.write("\t\t[]\n")

    # closing bracket
    file.write("\t[]\n")

# create AuxKernels
with open("AuxKernels.txt", "w") as file:
  file.write("\t[AuxKernels]\n")
  #### skip 'tfunc_angled_right' and 'tfunc_vertical'
  # file.write("\t\t[tfunc]\n")
  # file.write("\t\t\ttype = ConstantAux\n")
  # file.write("\t\t\tvalue = 400\n")
  # file.write("\t\t\tvariable = temp\n")
  # file.write("\t\t[]\n")
  # file.write("\t\t[radius]\n")
  # file.write("\t\t\ttype = FunctionAux\n")
  # file.write("\t\t\tfunction = radius\n")
  # file.write("\t\t\tvariable = r\n")
  # file.write("\t\t[]\n")
  # aclp
  for name in contact_names:
    for dir in directions:
      file.write("\t\t[aclp"+name+dir+"]\n")
      file.write("\t\t\ttype = PenetrationAux\n")
      file.write("\t\t\tvariable = aclp"+name+"_contact"+dir+"\n")
      file.write("\t\t\tboundary = 'ACLP"+name+"'\n")
      # get paired boundary name
      bnds = name.split("_")
      file.write("\t\t\tpaired_boundary = 'ACLP"+bnds[1]+"_"+bnds[0]+"'\n")
      file.write("\t\t\tquantity = normal_force"+dir+"\n")
      file.write("\t\t[]\n")

  # tlp
  for name in contact_names:
    for dir in directions:
      file.write("\t\t[tlp"+name+dir+"]\n")
      file.write("\t\t\ttype = PenetrationAux\n")
      file.write("\t\t\tvariable = tlp"+name+"_contact"+dir+"\n")
      file.write("\t\t\tboundary = 'TLP"+name+"'\n")
      # get paired boundary name
      bnds = name.split("_")
      file.write("\t\t\tpaired_boundary = 'TLP"+bnds[1]+"_"+bnds[0]+"'\n")
      file.write("\t\t\tquantity = normal_force"+dir+"\n")
      file.write("\t\t[]\n")

  # block
  for loc in locations:
    for n in np.arange(1,4):
      for dir in directions:
        # get paired boundary name
        bnds = block_bnds[n-1].split("_")
        file.write("\t\t["+loc+"_block"+str(n)+"_"+str(bnds[0])+dir+"]\n")
        file.write("\t\t\ttype = PenetrationAux\n")
        file.write("\t\t\tvariable = block"+str(n)+"_"+loc+"_contact"+dir+"\n")
        file.write("\t\t\tpaired_boundary = 'block"+str(n)+"_"+loc+"'\n")
        file.write("\t\t\tboundary = '"+loc.upper()+block_bnds[n-1]+"'\n")
        file.write("\t\t\tquantity = normal_force"+dir+"\n")
        file.write("\t\t[]\n")

  # restraint
  for loc in locations:
    for dir in directions:
      file.write("\t\t["+loc+"101_restraint"+dir+"]\n")
      file.write("\t\t\ttype = PenetrationAux\n")
      file.write("\t\t\tvariable = "+loc+"101_restraint_contact"+dir+"\n")
      file.write("\t\t\tpaired_boundary = '"+loc.upper()+"101_"+loc.upper()+"RR4_2'\n")
      file.write("\t\t\tboundary = '"+loc.upper()+"RR4_2'\n")
      file.write("\t\t\tquantity = normal_force"+dir+"\n")
      file.write("\t\t[]\n")

  # closing bracket
  file.write("\t[]\n")

# create Postprocessors
with open("Postprocessors.txt", "w") as file:
  file.write("\t[Postprocessors]\n")
  file.write("\t\t[pdata]\n")
  file.write("\t\t\ttype = PerfGraphData\n")
  file.write("\t\t\tdata_type = total\n")
  file.write("\t\t\tsection_name = 'Root'\n")
  file.write("\t\t\texecute_on = timestep_end\n")
  file.write("\t\t[]\n")

  # aclp
  for c_name, pc_name in pp_names.items():
    pp_composed_name = "aclp_force_"+pc_name
    pp_name = pp_composed_name
    if ((pc_name in half_force_set)):
      pp_name = pp_name+"_half_force"
    for dir in directions:
      file.write("\t\t["+pp_name+dir+"]\n")
      file.write("\t\t\ttype = NodalSum\n")
      file.write("\t\t\tvariable = aclp"+c_name+"_contact"+dir+"\n")
      file.write("\t\t\tboundary = 'ACLP"+c_name+"'\n")
      file.write("\t\t[]\n")
    # get the magnitude
    file.write("\t\t["+ pp_composed_name +"]\n")
    file.write("\t\t\ttype = ParsedPostprocessor\n")
    xx = pp_name+"_x"
    yy = pp_name+"_y"
    zz = pp_name+"_z"
    if ((pc_name in half_force_set)):
      file.write("\t\t\texpression = '2* sqrt("+xx+"*"+xx+"+"+yy+"*"+yy+"+"+zz+"*"+zz +")'\n")
    else:
      file.write("\t\t\texpression = 'sqrt("+xx+"*"+xx+"+"+yy+"*"+yy+"+"+zz+"*"+zz +")'\n")
    file.write("\t\t\tpp_names = '"+xx + " "+yy + " "+zz+"' \n")
    file.write("\t\t[]\n\n")

  #tlp
  for c_name, pc_name in pp_names.items():
    pp_composed_name = "tlp_force_"+pc_name
    pp_name = pp_composed_name
    if ((pc_name in half_force_set)):
      pp_name = pp_name+"_half_force"
    for dir in directions:
      file.write("\t\t["+pp_name+dir+"]\n")
      file.write("\t\t\ttype = NodalSum\n")
      file.write("\t\t\tvariable = tlp"+c_name+"_contact"+dir+"\n")
      file.write("\t\t\tboundary = 'TLP"+c_name+"'\n")
      file.write("\t\t[]\n")
    # get the magnitude
    file.write("\t\t["+ pp_composed_name +"]\n")
    file.write("\t\t\ttype = ParsedPostprocessor\n")
    xx = pp_name+"_x"
    yy = pp_name+"_y"
    zz = pp_name+"_z"
    if ((pc_name in half_force_set)):
      file.write("\t\t\texpression = '2* sqrt("+xx+"*"+xx+"+"+yy+"*"+yy+"+"+zz+"*"+zz +")'\n")
    else:
      file.write("\t\t\texpression = 'sqrt("+xx+"*"+xx+"+"+yy+"*"+yy+"+"+zz+"*"+zz +")'\n")
    file.write("\t\t\tpp_names = '"+xx + " "+yy + " "+zz+"' \n")
    file.write("\t\t[]\n\n")

  file.write("\t[]\n")
