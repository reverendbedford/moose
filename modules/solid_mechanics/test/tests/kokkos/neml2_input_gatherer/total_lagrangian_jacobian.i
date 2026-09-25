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
    expression = '0.1*x + 0.2*y + 0.3*z'
  []
  [disp_y]
    type = ParsedFunction
    expression = '0.4*x + 0.5*y + 0.6*z'
  []
  [disp_z]
    type = ParsedFunction
    expression = '0.7*x + 0.8*y + 0.9*z'
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
  []
  [disp_y]
    type = KokkosTotalLagrangianStressDivergence
    variable = disp_y
    component = 1
    displacements = 'disp_x disp_y disp_z'
  []
  [disp_z]
    type = KokkosTotalLagrangianStressDivergence
    variable = disp_z
    component = 2
    displacements = 'disp_x disp_y disp_z'
  []
[]

[BCs]
  [disp_x]
    type = FunctionDirichletBC
    variable = disp_x
    boundary = 'left right bottom top back front'
    function = disp_x
  []
  [disp_y]
    type = FunctionDirichletBC
    variable = disp_y
    boundary = 'left right bottom top back front'
    function = disp_y
  []
  [disp_z]
    type = FunctionDirichletBC
    variable = disp_z
    boundary = 'left right bottom top back front'
    function = disp_z
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
  []
[]

[Materials]
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
[]
