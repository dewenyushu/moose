[Mesh]
  [file]
    type = FileMeshGenerator
    file = 032_1.1e+02mmH_11mmGW_6mmCW_Gyroid_3mmLP_1.05e+03umLWT_500umBR_trunc18.unv
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
