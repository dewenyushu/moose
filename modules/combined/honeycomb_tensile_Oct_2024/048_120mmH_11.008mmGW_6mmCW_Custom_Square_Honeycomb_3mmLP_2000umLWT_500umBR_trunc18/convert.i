[Mesh]
  [file]
    type = FileMeshGenerator
    file = 048_120mmH_11.008mmGW_6mmCW_Custom_Square_Honeycomb_3mmLP_2000umLWT_500umBR_trunc18.unv
  []

  [traction_bottom]
    type = ParsedGenerateSideset
    combinatorial_geometry = 'abs(z+9.0)<10e-3' # 6mm from the center origin is the top/bottom face
    #normal = '0 0 -1'
    new_sideset_name = traction_bottom
    input = file
  []

  [traction_top]
    type = ParsedGenerateSideset
    combinatorial_geometry = 'abs(z-9.0)<10e-3'
    #normal = '0 0 1'
    new_sideset_name = traction_top
    input = traction_bottom
  []
[]
