[Mesh]
  [generated]
    type = GeneratedMeshGenerator
    dim = 2
    nx = 4
    ny = 4
  []
[]

[Variables]
  [disp_x]
  []
  [disp_y]
  []
[]

[NEML2]
  eager = true
  input = 'elasticity_neml2.i'
  [all]
    executor_name = neml2
    model = model
    input_kernels = neml2_strain
    output_backend = kokkos
    derivatives = 'neml2_stress neml2_strain Jacobian_mult'
  []
[]

[UserObjects]
  [assembly]
    type = TorchAssembly
  []
  [fe]
    type = TorchFEInterpolation
    assembly = assembly
  []
  [neml2_strain]
    type = TorchSmallStrain
    assembly = assembly
    fe = fe
    to_neml2 = neml2_strain
    displacements = 'disp_x disp_y'
  []
[]

[Kernels]
  [stress_x]
    type = KokkosStressDivergence
    variable = disp_x
    component = 0
    displacements = 'disp_x disp_y'
    stress = neml2_stress
  []
  [stress_y]
    type = KokkosStressDivergence
    variable = disp_y
    component = 1
    displacements = 'disp_x disp_y'
    stress = neml2_stress
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
    value = 0.1
  []
  [disp_y]
    type = KokkosDirichletBC
    variable = disp_y
    boundary = 'top bottom'
    value = 0
  []
[]

[Preconditioning]
  [smp]
    type = SMP
    full = true
  []
[]

[Executioner]
  type = Steady
  solve_type = NEWTON
  nl_abs_tol = 1e-12
[]

[Outputs]
  exodus = true
  file_base = neml2_action_kokkos_linear_elasticity_out
[]
