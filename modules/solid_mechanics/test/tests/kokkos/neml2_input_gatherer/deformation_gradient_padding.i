[Mesh]
  type = GeneratedMesh
  dim = 1
  nx = 1
[]

[Functions]
  [disp_x]
    type = ParsedFunction
    expression = '0.25*x'
  []
[]

[Variables]
  [disp_x]
    [InitialCondition]
      type = FunctionIC
      function = disp_x
    []
  []
[]

[Kernels]
  [disp_x]
    type = Diffusion
    variable = disp_x
  []
[]

[BCs]
  [disp_x]
    type = FunctionDirichletBC
    variable = disp_x
    boundary = 'left right'
    function = disp_x
  []
[]

[NEML2]
  eager = true
  input = 'deformation_gradient_neml2.i'
  [all]
    executor_name = neml2
    model = model
    device = cpu
    input_kernels = deformation_gradient
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
    displacements = disp_x
  []
[]

[AuxVariables]
  [F00]
    family = MONOMIAL
    order = CONSTANT
  []
  [F11]
    family = MONOMIAL
    order = CONSTANT
  []
  [F22]
    family = MONOMIAL
    order = CONSTANT
  []
  [F12]
    family = MONOMIAL
    order = CONSTANT
  []
[]

[AuxKernels]
  [F00]
    type = RankTwoAux
    variable = F00
    rank_two_tensor = copied_deformation_gradient
    index_i = 0
    index_j = 0
  []
  [F11]
    type = RankTwoAux
    variable = F11
    rank_two_tensor = copied_deformation_gradient
    index_i = 1
    index_j = 1
  []
  [F22]
    type = RankTwoAux
    variable = F22
    rank_two_tensor = copied_deformation_gradient
    index_i = 2
    index_j = 2
  []
  [F12]
    type = RankTwoAux
    variable = F12
    rank_two_tensor = copied_deformation_gradient
    index_i = 1
    index_j = 2
  []
[]

[Postprocessors]
  [F00]
    type = ElementAverageValue
    variable = F00
  []
  [F11]
    type = ElementAverageValue
    variable = F11
  []
  [F22]
    type = ElementAverageValue
    variable = F22
  []
  [F12]
    type = ElementAverageValue
    variable = F12
  []
[]

[Executioner]
  type = Steady
[]

[Outputs]
  csv = true
[]
