N = 16

# Step 3 (LD counterpart of plasticity/step3): CPU NEML2, GPU Kokkos assembly,
# and CPU PETSc.
#
# Differences vs plasticity/step3_plasticity_cpu_neml2_kokkos_cpu_petsc.i:
#   * TorchDeformationGradient replaces TorchSmallStrain so NEML2 receives the
#     full non-symmetric F = I + grad(u).
#   * NEML2ToKokkosFullRankTwoMaterialProperty and NEML2ToKokkosFullRankFourMaterialProperty
#     transfer the full PK2 stress and full dS/dF without Mandel/permutation.
#   * KokkosComputeLagrangianStressCustomPK2 converts (PK2, dS/dF) to
#     (PK1, dP/d(grad u)) needed by the total-Lagrangian kernel.
#   * KokkosTotalLagrangianStressDivergence replaces KokkosStressDivergence.
#   * The crystal-plasticity exact-kinematics NEML2 model with an
#     initialize_outputs/initialize_output_values pair seeds
#     plastic_deformation_gradient before the first Newton evaluation.

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
    # Match host step1/step2 which use Physics/SolidMechanics volumetric_locking_correction=true
    # (default total F-bar mode). See ../neml2_input_gatherer/fbar_kokkos_vs_host.i for the
    # verification against the host reference solution.
    stabilize_strain = true
  []
[]

[Kernels]
  [stress_x]
    type = KokkosTotalLagrangianStressDivergence
    variable = disp_x
    component = 0
    displacements = 'disp_x disp_y disp_z'
    stabilize_strain = true
  []
  [stress_y]
    type = KokkosTotalLagrangianStressDivergence
    variable = disp_y
    component = 1
    displacements = 'disp_x disp_y disp_z'
    stabilize_strain = true
  []
  [stress_z]
    type = KokkosTotalLagrangianStressDivergence
    variable = disp_z
    component = 2
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
    preset = false
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

# hypre boomeramg + gmres, line_search=none: see step1_LD for the rationale.
[Executioner]
  type = Transient
  solve_type = NEWTON
  line_search = none
  petsc_options_iname = '-pc_type -pc_hypre_type -ksp_type'
  petsc_options_value = 'hypre boomeramg    gmres'
  dt = 1e-3
  dtmin = 1e-3
  num_steps = 5
  nl_abs_tol = 1e-10
  residual_and_jacobian_together = true
  automatic_scaling = true
[]

# Diagnostic sample points mirror plasticity/step3. ux_center probes the
# interior displacement; ux_right verifies that the Kokkos DirichletBC picks up
# the time-varying loading applied through RealFunctionControl.
[Postprocessors]
  [ux_right]
    type = PointValue
    variable = disp_x
    point = '1 0.5 0.5'
    execute_on = TIMESTEP_END
  []
  [ux_center]
    type = PointValue
    variable = disp_x
    point = '0.5 0.5 0.5'
    execute_on = TIMESTEP_END
  []
[]

[Outputs]
  exodus = false
  csv = true
[]
