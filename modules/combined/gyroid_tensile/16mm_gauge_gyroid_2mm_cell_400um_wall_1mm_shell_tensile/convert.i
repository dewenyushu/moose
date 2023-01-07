[Mesh]
  [file]
    type = FileMeshGenerator
    file = 16mm_gauge_gyroid_2mm_cell_400um_wall_1mm_shell_tensile.unv
  []

  [separate_blocks]
    type = ParsedSubdomainMeshGenerator
    combinatorial_geometry = 'r:=sqrt(x*x+z*z); r<1.5' # 1.5mm is the inner radius of gyroid structure
    block_id = '1'
    input = file
  []

  [rename_blocks]
    type = RenameBlockGenerator
    old_block = '0        1'
    new_block = 'outside  gyroid'
    input = separate_blocks
  []

  [traction_bottom]
    type = ParsedGenerateSideset
    combinatorial_geometry = 'abs(y+24.9968)<2e-3' # 25mm from the center origin is the top/bottom face
    #normal = '0 0 -1'
    included_subdomains = outside
    new_sideset_name = traction_bottom
    input = rename_blocks
  []

  [traction_top]
    type = ParsedGenerateSideset
    combinatorial_geometry = 'abs(y-24.9968)<2e-3' # 25mm from the center origin is the top/bottom face
    #normal = '0 0 -1'
    included_subdomains = outside
    new_sideset_name = traction_top
    input = traction_bottom
  []
[]
