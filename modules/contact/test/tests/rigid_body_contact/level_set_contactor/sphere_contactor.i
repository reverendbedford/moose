# Probe a SphereContactor at every node of a small 2D grid and gold-file the
# signed distance and normal components against the analytical values for a
# circle of radius 0.5 centered at (1, 1).

[Mesh]
  [gen]
    type = GeneratedMeshGenerator
    dim = 2
    nx = 2
    ny = 2
    xmin = 0
    xmax = 2
    ymin = 0
    ymax = 2
  []
[]

[UserObjects]
  [sphere]
    type = SphereContactor
    center = '1 1 0'
    radius = 0.5
  []
[]

[AuxVariables]
  [sdf]
    family = LAGRANGE
    order = FIRST
  []
  [nx]
    family = LAGRANGE
    order = FIRST
  []
  [ny]
    family = LAGRANGE
    order = FIRST
  []
[]

[AuxKernels]
  [sdf]
    type = LevelSetContactorAux
    variable = sdf
    contactor = sphere
    quantity = signed_distance
    execute_on = 'INITIAL'
  []
  [nx]
    type = LevelSetContactorAux
    variable = nx
    contactor = sphere
    quantity = normal_x
    execute_on = 'INITIAL'
  []
  [ny]
    type = LevelSetContactorAux
    variable = ny
    contactor = sphere
    quantity = normal_y
    execute_on = 'INITIAL'
  []
[]

[VectorPostprocessors]
  [nodal]
    type = NodalValueSampler
    variable = 'sdf nx ny'
    sort_by = id
    execute_on = 'INITIAL'
  []
[]

[Problem]
  solve = false
[]

[Executioner]
  type = Steady
[]

[Outputs]
  [csv]
    type = CSV
    execute_on = 'INITIAL'
  []
[]
