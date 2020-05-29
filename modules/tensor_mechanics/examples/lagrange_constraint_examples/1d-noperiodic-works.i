# Simple 1D plane strain test

[GlobalParams]
  displacements = 'disp_x'
[]

[Mesh]
  [./msh]
    type = GeneratedMeshGenerator
    dim = 1
    nx = 4
    xmin = -1
    xmax = 3
  []

  [cnode]
    type = ExtraNodesetGenerator
    coord = '1'
    new_boundary = 1001
    input = msh
  []
[]

[Modules/TensorMechanics/Master]
  [./all]
    strain = SMALL
    add_variables = true
    generate_output = 'stress_xx'
  [../]
[../]

[Variables]
  [./lambda_11]
    family = SCALAR
    order = FIRST
  [../]
[]

[Kernels]
  [./lagrange_11]
    type = HomogenizationConstraint
    variable = disp_x
    lambda = lambda_11
    component = 0
    quantity = cauchy
    index_i = 0
    index_j = 0
    target = func
  [../]
[]

[ScalarKernels]
  [./force_11]
    type = DummyLagrange
    variable = lambda_11
  [../]
[]

[Functions]
  [./func]
    type = ParsedFunction
    value = '100*t'
  [../]
[]

[BCs]
#  [./Periodic]
#    [./all]
#      variable = disp_x
#      auto_direction = 'x'
#    [../]
#  [../]

  [./centerfix_x]
    type = DirichletBC
    boundary = 1001
    variable = disp_x
    value = 0
  [../]
[]

#[Constraints]
#  [./lr]
#    type = EqualValueBoundaryConstraint
#    variable = disp_x
#    master = '0'
#    slave = 'right'
#    penalty = 10e6
#  [../]
#[]

[Materials]
  [./elasticity_tensor]
    type = ComputeIsotropicElasticityTensor
    youngs_modulus = 1e5
    poissons_ratio = 0.3
  [../]
  [./stress]
    type = ComputeLinearElasticStress
  [../]
[]

[Preconditioning]
  [./smp]
    type = SMP
    full = true
  [../]
[]

[Postprocessors]
  [./sxx]
    type = ElementAverageValue
    variable = stress_xx
    execute_on = 'initial timestep_end'
  [../]
[]

[Executioner]
  type = Transient

  solve_type = 'newton'
  line_search = none

  petsc_options_iname = '-pc_type -ksp_view_rhs -ksp_view_mat'
  petsc_options_value = 'lu ascii ascii'

  l_max_its = 2
  l_tol = 1e-14
  nl_max_its = 15
  nl_rel_tol = 1e-8
  nl_abs_tol = 1e-10

  start_time = 0.0
  dt = 0.2
  dtmin = 0.2
  end_time = 1.0
[]

[Outputs]
  exodus = false
[]
