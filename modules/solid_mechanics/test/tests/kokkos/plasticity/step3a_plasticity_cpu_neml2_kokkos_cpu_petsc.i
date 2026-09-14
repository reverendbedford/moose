N = 16

# Step 3a: CPU NEML2, GPU Kokkos assembly, and CPU PETSc.
# NEML2 inputs move from the GPU to the CPU, and stress/tangent outputs move back to the GPU.
[Mesh]
  [generated]
    type = GeneratedMeshGenerator
    dim = 3
    nx = ${N}
    ny = ${N}
    nz = ${N}
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

[NEML2]
  eager = true
  input = '../../neml2/plasticity/perfect_neml2.i'
  [all]
    executor_name = neml2
    model = model
    device = cpu
    input_kernels = neml2_strain
    auto_output = false
    manage_state_advance = true
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
    displacements = 'disp_x disp_y disp_z'
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
  [stress]
    type = NEML2ToKokkosRankTwoMaterialProperty
    neml2_executor = neml2
    from_neml2 = neml2_stress
    to_moose = stress
  []
  [tangent]
    type = NEML2ToKokkosRankFourMaterialProperty
    neml2_executor = neml2
    from_neml2 = neml2_stress
    neml2_input_derivative = neml2_strain
    to_moose = Jacobian_mult
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
  petsc_options_iname = '-pc_type -ksp_type'
  petsc_options_value = 'gamg gmres'
  dt = 1e-3
  dtmin = 1e-3
  num_steps = 5
  nl_abs_tol = 1e-10
  residual_and_jacobian_together = true
[]

[Outputs]
  exodus = false
[]
