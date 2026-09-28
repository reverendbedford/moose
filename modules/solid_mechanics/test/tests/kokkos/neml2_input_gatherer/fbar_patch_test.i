# Homogeneous patch test for the total-mode F-bar Kokkos path.
#
# The applied displacement field is exactly reconstructible by linear hex shape functions:
#   u_x = alpha * x, u_y = beta * y, u_z = gamma * z
# with (alpha, beta, gamma) = (0.05, -0.02, 0.01). F_ust = diag(1.05, 0.98, 1.01) is CONSTANT
# across every quadrature point, so F_avg = F_ust and gamma_bar = cbrt(det F_avg / det F_ust) = 1
# identically. Under total-mode F-bar this makes F_stab = F_ust point-wise. The stabilized and
# unstabilized paths must therefore produce bit-for-bit identical stresses on this patch.
#
# The comparison is done externally by the test harness (`total_lagrangian_fbar_patch` +
# `total_lagrangian_fbar_patch_unstab` share the same *_out.csv gold file). Any drift indicates
# either the F-bar chain is being applied where gamma_bar should be 1 exactly, or that the
# rank-4 Jacobian chain leaks a nonzero non-local perturbation on constant F.
[Mesh]
  [gmg]
    type = GeneratedMeshGenerator
    dim = 3
    nx = 2
    ny = 2
    nz = 2
  []
[]

[Functions]
  [disp_x]
    type = ParsedFunction
    expression = '0.05*x'
  []
  [disp_y]
    type = ParsedFunction
    expression = '-0.02*y'
  []
  [disp_z]
    type = ParsedFunction
    expression = '0.01*z'
  []
[]

[Variables]
  [disp_x]
    [InitialCondition]
      type = FunctionIC
      function = disp_x
    []
  []
  [disp_y]
    [InitialCondition]
      type = FunctionIC
      function = disp_y
    []
  []
  [disp_z]
    [InitialCondition]
      type = FunctionIC
      function = disp_z
    []
  []
[]

[Kernels]
  [disp_x]
    type = KokkosTotalLagrangianStressDivergence
    variable = disp_x
    component = 0
    displacements = 'disp_x disp_y disp_z'
    stabilize_strain = ${STABILIZE}
  []
  [disp_y]
    type = KokkosTotalLagrangianStressDivergence
    variable = disp_y
    component = 1
    displacements = 'disp_x disp_y disp_z'
    stabilize_strain = ${STABILIZE}
  []
  [disp_z]
    type = KokkosTotalLagrangianStressDivergence
    variable = disp_z
    component = 2
    displacements = 'disp_x disp_y disp_z'
    stabilize_strain = ${STABILIZE}
  []
[]

[BCs]
  [disp_x]
    type = FunctionDirichletBC
    variable = disp_x
    boundary = 'left right'
    function = disp_x
    preset = true
  []
  [disp_y]
    type = FunctionDirichletBC
    variable = disp_y
    boundary = 'bottom top'
    function = disp_y
    preset = true
  []
  [disp_z]
    type = FunctionDirichletBC
    variable = disp_z
    boundary = 'back front'
    function = disp_z
    preset = true
  []
[]

[NEML2]
  eager = true
  input = 'full_tensor_bridge_neml2.i'
  [all]
    executor_name = neml2
    model = model
    device = cpu
    input_kernels = deformation_gradient
    auto_output = false
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
  [deformation_gradient]
    type = TorchDeformationGradient
    assembly = assembly
    fe = fe
    to_neml2 = deformation_gradient
    displacements = 'disp_x disp_y disp_z'
    stabilize_strain = ${STABILIZE}
  []
[]

[Materials]
  [f_bar_average]
    type = KokkosComputeFbarAverage
    displacements = 'disp_x disp_y disp_z'
  []
  [stress]
    type = NEML2ToKokkosFullRankTwoMaterialProperty
    neml2_executor = neml2
    from_neml2 = stress
    to_moose = stress
  []
  [tangent]
    type = NEML2ToKokkosFullRankFourMaterialProperty
    neml2_executor = neml2
    from_neml2 = stress
    neml2_input_derivative = deformation_gradient
    to_moose = tangent
  []
  [pk1]
    type = KokkosComputeLagrangianStressCustomPK2
    displacements = 'disp_x disp_y disp_z'
    custom_pk2_stress = stress
    custom_pk2_jacobian = tangent
    stabilize_strain = ${STABILIZE}
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
  nl_abs_tol = 1e-12
  nl_rel_tol = 1e-12
  solve_type = NEWTON
  petsc_options_iname = '-pc_type'
  petsc_options_value = 'lu'
  line_search = none
[]

# Patch-test points: pick the interior mid-cube node (0.5, 0.5, 0.5), which lies inside the
# mesh (not on any Dirichlet boundary). If the discrete solution matches the applied linear
# field on the interior, it satisfies the patch test.
#   Expected: u_x(0.5, 0.5, 0.5) = 0.05 * 0.5   =  0.025
#             u_y(0.5, 0.5, 0.5) = -0.02 * 0.5  = -0.010
#             u_z(0.5, 0.5, 0.5) =  0.01 * 0.5  =  0.005
[Postprocessors]
  [ux_center]
    type = PointValue
    variable = disp_x
    point = '0.5 0.5 0.5'
    execute_on = TIMESTEP_END
  []
  [uy_center]
    type = PointValue
    variable = disp_y
    point = '0.5 0.5 0.5'
    execute_on = TIMESTEP_END
  []
  [uz_center]
    type = PointValue
    variable = disp_z
    point = '0.5 0.5 0.5'
    execute_on = TIMESTEP_END
  []
[]

[Outputs]
  csv = true
[]
