# Simple 1D plane strain test

[GlobalParams]
  displacements = 'disp_x disp_y disp_z'
  large_kinematics = true
[]

[Mesh]
  [msh]
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

[Kernels]
  [sdx]
    type = TotalLagrangianStressDivergence
    variable = disp_x
    component = 0
  []
  [sdy]
    type = TotalLagrangianStressDivergence
    variable = disp_y
    component = 1
  []
  [sdz]
    type = TotalLagrangianStressDivergence
    variable = disp_z
    component = 2
  []
[]

[Functions]
  [pullx]
    type = ParsedFunction
    expression = '4000 * t'
  []
  [pully]
    type = ParsedFunction
    expression = '-2000 * t'
  []
  [pullz]
    type = ParsedFunction
    expression = '3000 * t'
  []
[]

[BCs]
  [leftx]
    type = DirichletBC
    preset = true
    boundary = left
    variable = disp_x
    value = 0.0
  []
  [lefty]
    type = DirichletBC
    preset = true
    boundary = left
    variable = disp_y
    value = 0.0
  []
  [leftz]
    type = DirichletBC
    preset = true
    boundary = left
    variable = disp_z
    value = 0.0
  []
  [pull_x]
    type = FunctionNeumannBC
    boundary = right
    variable = disp_x
    function = pullx
  []
  [pull_y]
    type = FunctionNeumannBC
    boundary = top
    variable = disp_y
    function = pully
  []
  [pull_z]
    type = FunctionNeumannBC
    boundary = right
    variable = disp_z
    function = pullz
  []
[]

[Materials]
  [elastic_tensor]
    type = ComputeIsotropicElasticityTensor
    lambda = 5000
    shear_modulus = 5000
  []
  [compute_stress]
    type = ComputeLagrangianLinearElasticStress
    objective_rate = jaumann
  []
  [compute_strain]
    type = ComputeLagrangianStrain
    increment_approximation = linear
    alpha = 1.0
  []
[]

[Preconditioning]
  [precon]
    type = SMP
    full = true
  []
[]

[Postprocessors]
  [nonlin_it]
    type = NumNonlinearIterations
  []
  [cumulative_nonlin_it]
    type = CumulativeValuePostprocessor
    postprocessor = nonlin_it
  []
[]

[Executioner]
  type = Transient

  solve_type = 'newton'
  line_search = none

  petsc_options_iname = '-pc_type'
  petsc_options_value = 'lu'

  l_max_its = 2
  l_tol = 1e-14
  nl_max_its = 15
  nl_rel_tol = 1e-8
  nl_abs_tol = 1e-10

  start_time = 0.0
  dt = 0.1
  dtmin = 0.1
  end_time = 0.5
[]

[Outputs]
  exodus = false
  [out]
    type = CSV
    file_base = linear_1_jaumann
  []
[]
