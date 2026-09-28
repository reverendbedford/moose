# Nonuniform-deformation F-bar comparison: this input drives the Kokkos total-Lagrangian NEML2
# stack with `stabilize_strain = true` under a Dirichlet-pull setup where the interior displacement
# is genuinely solved (not just fixed by BCs). Its companion `fbar_kokkos_vs_host_reference.i`
# solves the same problem via the host `Physics/SolidMechanics/QuasiStatic` action with
# `formulation = TOTAL`, `strain = FINITE`, `volumetric_locking_correction = true` (which selects
# the same total-mode F-bar). The tests spec's CSVDiff compares interior point displacements.
#
# The problem is intentionally nonuniform: pulling the right face while fixing three orthogonal
# faces (left, bottom, back) produces an inhomogeneous deformation field, so F_ust varies across
# quadrature points inside each element and the F-bar element average is nontrivial. This is the
# case where the Kokkos non-local Jacobian term (Phase 3) is required for Newton to be consistent.

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

[Variables]
  [disp_x]
  []
  [disp_y]
  []
  [disp_z]
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
  [xfix]
    type = DirichletBC
    variable = disp_x
    boundary = left
    value = 0
    preset = true
  []
  [yfix]
    type = DirichletBC
    variable = disp_y
    boundary = bottom
    value = 0
    preset = true
  []
  [zfix]
    type = DirichletBC
    variable = disp_z
    boundary = back
    value = 0
    preset = true
  []
  [xpull]
    type = FunctionDirichletBC
    variable = disp_x
    boundary = right
    function = t
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
  type = Transient
  solve_type = NEWTON
  line_search = none
  petsc_options_iname = '-pc_type'
  petsc_options_value = 'lu'
  automatic_scaling = true
  dt = 5e-3
  dtmin = 1e-3
  num_steps = 3
  residual_and_jacobian_together = true
[]

# Compare against the host reference at nine interior points spanning the element interiors.
# These are the QP-adjacent Gauss-Legendre integration coords for a 2x2x2 mesh, so they capture
# the discretized displacement where F-bar's element average has the largest effect.
[Postprocessors]
  [u_c0]
    type = PointValue
    variable = disp_x
    point = '0.25 0.25 0.25'
  []
  [u_c1]
    type = PointValue
    variable = disp_x
    point = '0.75 0.25 0.25'
  []
  [u_c2]
    type = PointValue
    variable = disp_x
    point = '0.25 0.75 0.25'
  []
  [u_c3]
    type = PointValue
    variable = disp_x
    point = '0.25 0.25 0.75'
  []
  [u_c4]
    type = PointValue
    variable = disp_x
    point = '0.75 0.75 0.75'
  []
  [uy_c0]
    type = PointValue
    variable = disp_y
    point = '0.5 0.5 0.5'
  []
  [uz_c0]
    type = PointValue
    variable = disp_z
    point = '0.5 0.5 0.5'
  []
[]

[Outputs]
  csv = true
[]
