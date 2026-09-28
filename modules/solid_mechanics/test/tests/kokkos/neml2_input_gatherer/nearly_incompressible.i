# Nearly-incompressible linear-hex F-bar test. Runs the Kokkos LD stack with
# `stabilize_strain = true` by default. The tests spec compares this configuration against the
# host `nearly_incompressible_host.i` (Physics/SolidMechanics with
# `volumetric_locking_correction = true`, default total F-bar mode) via CSVDiff.
# Add `Kernels/*/stabilize_strain=false Materials/pk1/stabilize_strain=false
# UserObjects/deformation_gradient/stabilize_strain=false` on the CLI to compare against the
# unstabilized (locking-prone) control.
[Mesh]
  [gmg]
    type = GeneratedMeshGenerator
    dim = 3
    nx = 4
    ny = 4
    nz = 4
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

[BCs]
  [xfix]
    type = DirichletBC
    variable = disp_x
    boundary = left
    value = 0
  []
  [yfix]
    type = DirichletBC
    variable = disp_y
    boundary = left
    value = 0
  []
  [zfix]
    type = DirichletBC
    variable = disp_z
    boundary = left
    value = 0
  []
  # Small uniform shear along y on the right face
  [ypull]
    type = FunctionDirichletBC
    variable = disp_y
    boundary = right
    function = t
    preset = false
  []
[]

[NEML2]
  eager = true
  input = 'nearly_incompressible_neml2.i'
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
    stabilize_strain = true
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
  num_steps = 2
  residual_and_jacobian_together = true
[]

[Postprocessors]
  # Interior displacement at the mid-front-face node (2, 4, 4)/4 = (0.5, 1, 1) in physical coords.
  # For a shear-loaded cube, the tip displacement is the standard benchmark quantity.
  [uy_tip]
    type = PointValue
    variable = disp_y
    point = '1 0.5 0.5'
    execute_on = TIMESTEP_END
  []
  [ux_tip]
    type = PointValue
    variable = disp_x
    point = '1 0.5 0.5'
    execute_on = TIMESTEP_END
  []
[]

[Outputs]
  csv = true
[]
