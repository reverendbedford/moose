Nelem = 3
Ngrain = 3

[Mesh]
  [gmg]
    type = GeneratedMeshGenerator
    dim = 3
    nx = ${fparse Nelem * Ngrain}
    ny = ${fparse Nelem * Ngrain}
    nz = ${fparse Nelem * Ngrain}
    xmin = -1
    ymin = -1
    zmin = -1
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
  input = 'perfect_neml2.i'
  [all]
    executor_name = neml2
    model = model
    #device = xpu
    device = cpu
    input_kernels = neml2_strain
    auto_output = false
    manage_state_advance = true
    parameters = 'sy'
    parameter_types = 'MATERIAL'
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
  [neml2_strain]
    type = TorchSmallStrain
    assembly = assembly
    fe = fe
    to_neml2 = neml2_strain
    displacements = 'disp_x disp_y disp_z'
  []
[]

[Kernels]
  [stress_x]
    type = KokkosStressDivergence
    variable = disp_x
    component = 0
    displacements = 'disp_x disp_y disp_z'
  []
  [stress_y]
    type = KokkosStressDivergence
    variable = disp_y
    component = 1
    displacements = 'disp_x disp_y disp_z'
  []
  [stress_z]
    type = KokkosStressDivergence
    variable = disp_z
    component = 2
    displacements = 'disp_x disp_y disp_z'
  []
[]

[Materials]
  [yield_stress]
    type = GenericConstantMaterial
    prop_names = 'sy'
    prop_values = '5.0'
  []

  [stress]
    type = NEML2ToKokkosRankTwoMaterialProperty
    neml2_executor = neml2
    from_neml2 = neml2_stress
    to_moose = stress
  []
  [tangent]
    type = NEML2ToKokkosRankFourMaterialProperty
    neml2_executor = neml2
    from_neml2 = neml2_stress
    neml2_input_derivative = neml2_strain
    to_moose = Jacobian_mult
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
    value = 0
  []
  [disp_y_bottom]
    type = KokkosDirichletBC
    variable = disp_y
    boundary = bottom
    value = 0
  []
  [disp_y_top]
    type = KokkosDirichletBC
    variable = disp_y
    boundary = top
    value = 0
    preset = false
  []
  [disp_z_back]
    type = KokkosDirichletBC
    variable = disp_z
    boundary = back
    value = 0
  []
  [disp_z_front]
    type = KokkosDirichletBC
    variable = disp_z
    boundary = front
    value = 0
    preset = false
  []
[]

[Functions]
  [loading_pos]
    type = ParsedFunction
    expression = t
  []
  [loading_neg]
    type = ParsedFunction
    expression = t
  []
[]

[Controls]
  [loading_top]
    type = RealFunctionControl
    parameter = 'BCs/disp_y_top/value'
    function = loading_pos
    execute_on = 'INITIAL TIMESTEP_BEGIN'
  []
#  [loading_front]
#    type = RealFunctionControl
#    parameter = 'BCs/disp_z_front/value'
#    function = loading_neg
#    execute_on = 'INITIAL TIMESTEP_BEGIN'
#  []
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
  petsc_options_iname = '-pc_type -ksp_type'
  petsc_options_value = 'gamg gmres'
  dt = 1e-3
  dtmin = 1e-3
  num_steps = 5
  nl_rel_tol = 1e-8
  nl_abs_tol = 1e-10

  automatic_scaling = false

  residual_and_jacobian_together = true

  l_tol = 1e-3
[]

[Outputs]
  file_base = 'results'
  exodus = true
  csv = true
  [pgraph]
    type = PerfGraphOutput
  []
[]
