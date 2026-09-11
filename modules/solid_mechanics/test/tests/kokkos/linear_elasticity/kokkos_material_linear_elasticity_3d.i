[Mesh]
  [generated]
    type = GeneratedMeshGenerator
    dim = 3
    nx = 2
    ny = 2
    nz = 2
  []
[]

[Variables]
  [disp_x]
  []
  [disp_y]
  []
  [disp_z]
  []
[]

[Kernels]
  [stress_x]
    type = KokkosStressDivergence
    variable = disp_x
    component = 0
    displacements = 'disp_x disp_y disp_z'
  []
  [stress_y]
    type = KokkosStressDivergence
    variable = disp_y
    component = 1
    displacements = 'disp_x disp_y disp_z'
  []
  [stress_z]
    type = KokkosStressDivergence
    variable = disp_z
    component = 2
    displacements = 'disp_x disp_y disp_z'
  []
[]

[Materials]
  [elasticity]
    type = KokkosComputeIsotropicElasticity
    displacements = 'disp_x disp_y disp_z'
    youngs_modulus = 10
    poissons_ratio = 0.25
  []
[]

[BCs]
  [disp_x_left]
    type = KokkosDirichletBC
    variable = disp_x
    boundary = left
    value = 0
  []
  [disp_x_right]
    type = KokkosDirichletBC
    variable = disp_x
    boundary = right
    value = 0
  []
  [disp_y]
    type = KokkosDirichletBC
    variable = disp_y
    boundary = bottom
    value = 0
  []
  [disp_z]
    type = KokkosDirichletBC
    variable = disp_z
    boundary = back
    value = 0
  []
[]

[Functions]
  [loading]
    type = ParsedFunction
    expression = t
  []
[]

[Controls]
  [loading]
    type = RealFunctionControl
    parameter = 'BCs/disp_x_right/value'
    function = loading
    execute_on = 'INITIAL TIMESTEP_BEGIN'
  []
[]

[Preconditioning]
  [smp]
    type = SMP
    full = true
  []
[]

[Executioner]
  type = Transient
  solve_type = NEWTON
  dt = 0.001
  num_steps = 5
  nl_abs_tol = 1e-12
[]

[Outputs]
  exodus = true
  file_base = kokkos_material_linear_elasticity_3d_out
[]
