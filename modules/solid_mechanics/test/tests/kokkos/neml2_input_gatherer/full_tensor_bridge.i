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
    type = Diffusion
    variable = disp_x
  []
  [disp_y]
    type = Diffusion
    variable = disp_y
  []
  [disp_z]
    type = Diffusion
    variable = disp_z
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
  [check]
    type = KokkosFullTensorTestMaterial
    rank_two = stress
    rank_four = tangent
    pk1 = pk1_stress
    A = '1 2 3 5 7 11 13 17 19'
    B = '23 29 31 37 41 43 47 53 59'
    F = '1.1 0.2 0.3 0.4 1.5 0.6 0.7 0.8 1.9'
  []
[]

[AuxVariables]
  [rank_two_error]
    family = MONOMIAL
    order = CONSTANT
  []
  [rank_four_error]
    family = MONOMIAL
    order = CONSTANT
  []
  [pk1_error]
    family = MONOMIAL
    order = CONSTANT
  []
[]

[AuxKernels]
  [rank_two_error]
    type = KokkosMaterialRealAux
    variable = rank_two_error
    property = rank_two_error
  []
  [rank_four_error]
    type = KokkosMaterialRealAux
    variable = rank_four_error
    property = rank_four_error
  []
  [pk1_error]
    type = KokkosMaterialRealAux
    variable = pk1_error
    property = pk1_error
  []
[]

[Postprocessors]
  [rank_two_error]
    type = ElementAverageValue
    variable = rank_two_error
  []
  [rank_four_error]
    type = ElementAverageValue
    variable = rank_four_error
  []
  [pk1_error]
    type = ElementAverageValue
    variable = pk1_error
  []
[]

[Executioner]
  type = Steady
  nl_rel_tol = 1e-12
  nl_abs_tol = 1e-12
[]

[Outputs]
  csv = true
[]
