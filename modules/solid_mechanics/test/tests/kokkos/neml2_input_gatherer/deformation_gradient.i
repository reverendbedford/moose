[Mesh]
  [gmg]
    type = GeneratedMeshGenerator
    dim = 3
    nx = 1
    ny = 1
    nz = 1
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
    displacements = 'disp_x disp_y disp_z'
  []
[]

[AuxVariables]
  [F00]
    family = MONOMIAL
    order = CONSTANT
  []
  [F01]
    family = MONOMIAL
    order = CONSTANT
  []
  [F02]
    family = MONOMIAL
    order = CONSTANT
  []
  [F10]
    family = MONOMIAL
    order = CONSTANT
  []
  [F11]
    family = MONOMIAL
    order = CONSTANT
  []
  [F12]
    family = MONOMIAL
    order = CONSTANT
  []
  [F20]
    family = MONOMIAL
    order = CONSTANT
  []
  [F21]
    family = MONOMIAL
    order = CONSTANT
  []
  [F22]
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
  [F01]
    type = RankTwoAux
    variable = F01
    rank_two_tensor = copied_deformation_gradient
    index_i = 0
    index_j = 1
  []
  [F02]
    type = RankTwoAux
    variable = F02
    rank_two_tensor = copied_deformation_gradient
    index_i = 0
    index_j = 2
  []
  [F10]
    type = RankTwoAux
    variable = F10
    rank_two_tensor = copied_deformation_gradient
    index_i = 1
    index_j = 0
  []
  [F11]
    type = RankTwoAux
    variable = F11
    rank_two_tensor = copied_deformation_gradient
    index_i = 1
    index_j = 1
  []
  [F12]
    type = RankTwoAux
    variable = F12
    rank_two_tensor = copied_deformation_gradient
    index_i = 1
    index_j = 2
  []
  [F20]
    type = RankTwoAux
    variable = F20
    rank_two_tensor = copied_deformation_gradient
    index_i = 2
    index_j = 0
  []
  [F21]
    type = RankTwoAux
    variable = F21
    rank_two_tensor = copied_deformation_gradient
    index_i = 2
    index_j = 1
  []
  [F22]
    type = RankTwoAux
    variable = F22
    rank_two_tensor = copied_deformation_gradient
    index_i = 2
    index_j = 2
  []
[]

[Postprocessors]
  [F00]
    type = ElementAverageValue
    variable = F00
  []
  [F01]
    type = ElementAverageValue
    variable = F01
  []
  [F02]
    type = ElementAverageValue
    variable = F02
  []
  [F10]
    type = ElementAverageValue
    variable = F10
  []
  [F11]
    type = ElementAverageValue
    variable = F11
  []
  [F12]
    type = ElementAverageValue
    variable = F12
  []
  [F20]
    type = ElementAverageValue
    variable = F20
  []
  [F21]
    type = ElementAverageValue
    variable = F21
  []
  [F22]
    type = ElementAverageValue
    variable = F22
  []
[]

[Executioner]
  type = Steady
[]

[Outputs]
  csv = true
[]
