# PetscJacobianTester input for the Kokkos F-bar path against the stateful crystal-plasticity
# NEML2 model (unlike total_lagrangian_jacobian_stabilized.i which uses the linear-elastic
# `full_tensor_bridge_neml2.i` bridge model). Steady solve on a small nonuniform-F configuration.
#
# Uses `[Variables]` with FunctionIC seeded from a mixed linear+quadratic displacement so that:
#   1. grad(u) varies across quadrature points => nontrivial element average
#   2. det F stays positive at every qp (avoids the sign-cbrt fallback)
#   3. Every internal node has a distinct displacement, so the cross-QP F-bar Jacobian coupling
#      is exercised.

N = 2

[Mesh]
  [gmg]
    type = GeneratedMeshGenerator
    dim = 3
    nx = ${N}
    ny = ${N}
    nz = ${N}
  []
[]

[Functions]
  [disp_x]
    type = ParsedFunction
    expression = '0.02*x + 0.005*y*y'
  []
  [disp_y]
    type = ParsedFunction
    expression = '-0.01*y + 0.004*x*z'
  []
  [disp_z]
    type = ParsedFunction
    expression = '0.003*z + 0.002*x*y'
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
    stabilize_strain = true
  []
  [disp_y]
    type = KokkosTotalLagrangianStressDivergence
    variable = disp_y
    component = 1
    displacements = 'disp_x disp_y disp_z'
    stabilize_strain = true
  []
  [disp_z]
    type = KokkosTotalLagrangianStressDivergence
    variable = disp_z
    component = 2
    displacements = 'disp_x disp_y disp_z'
    stabilize_strain = true
  []
[]

[BCs]
  [disp_x]
    type = FunctionDirichletBC
    variable = disp_x
    boundary = 'left right bottom top back front'
    function = disp_x
    preset = true
  []
  [disp_y]
    type = FunctionDirichletBC
    variable = disp_y
    boundary = 'left right bottom top back front'
    function = disp_y
    preset = true
  []
  [disp_z]
    type = FunctionDirichletBC
    variable = disp_z
    boundary = 'left right bottom top back front'
    function = disp_z
    preset = true
  []
[]

[NEML2]
  eager = true
  input = '../../neml2/crystal_plasticity/exact_kinematics_neml2.i'
  [all]
    executor_name = neml2
    model = model
    device = cpu
    input_kernels = deformation_gradient
    auto_output = false
    manage_state_advance = true
    initialize_outputs = 'plastic_deformation_gradient'
    initialize_output_values = 'initial_plastic_defgrad'
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
    stabilize_strain = true
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
    from_neml2 = neml2_stress
    to_moose = neml2_stress
  []
  [tangent]
    type = NEML2ToKokkosFullRankFourMaterialProperty
    neml2_executor = neml2
    from_neml2 = neml2_stress
    neml2_input_derivative = deformation_gradient
    to_moose = dneml2_stress_ddeformation_gradient
  []
  [pk1]
    type = KokkosComputeLagrangianStressCustomPK2
    displacements = 'disp_x disp_y disp_z'
    custom_pk2_stress = neml2_stress
    custom_pk2_jacobian = dneml2_stress_ddeformation_gradient
    stabilize_strain = true
  []
  [initial_orientation]
    type = GenericConstantRealVectorValue
    vector_name = orientation
    vector_values = '-0.54412095 -0.34931944 0.12600655'
  []
  [initial_plastic_defgrad]
    type = GenericConstantRankTwoTensor
    tensor_name = initial_plastic_defgrad
    tensor_values = '1 1 1'
  []
[]

[Preconditioning]
  [smp]
    type = SMP
    full = true
  []
[]

[Executioner]
  # Transient (not Steady) because the CP NEML2 model uses ScalarBackwardEulerTimeIntegration
  # and R2BackwardEulerTimeIntegration, which read the material `t` property that only the
  # Transient interface publishes. A single small step is enough for PetscJacobianTester to
  # exercise Jacobian assembly with a nonzero state increment.
  #
  # residual_and_jacobian_together is DISABLED here because it interferes with SNES's Jacobian
  # test: when residual+Jacobian are computed in one call, `-snes_test_jacobian`'s FD comparison
  # sees the analytic Jacobian assembled with the same shape it uses, producing the wrong verdict.
  type = Transient
  dt = 5e-3
  num_steps = 1
  solve_type = NEWTON
  line_search = none
  petsc_options_iname = '-pc_type'
  petsc_options_value = 'lu'
  automatic_scaling = true
[]
