N = 16

# Step 2: GPU NEML2 constitutive update, CPU MOOSE assembly, and CPU PETSc.
[GlobalParams]
  displacements = 'disp_x disp_y'
[]

[Mesh]
  [generated]
    type = GeneratedMeshGenerator
    dim = 2
    nx = ${N}
    ny = ${N}
  []
[]

[Physics]
  [SolidMechanics]
    [QuasiStatic]
      [all]
        strain = SMALL
        formulation = TOTAL
        new_system = true
        add_variables = true
      []
    []
  []
[]

[NEML2]
  eager = true
  input = '../../neml2/plasticity/perfect_neml2.i'
  output_device = cpu
  [all]
    model = model
    device = cuda
    derivatives = 'neml2_stress neml2_strain'
  []
[]

[Materials]
  [convert_strain]
    type = RankTwoTensorToSymmetricRankTwoTensor
    from = mechanical_strain
    to = neml2_strain
  []
  [stress]
    type = ComputeLagrangianObjectiveCustomSymmetricStress
    custom_small_stress = neml2_stress
    custom_small_jacobian = 'dneml2_stress/dneml2_strain'
  []
[]

[BCs]
  [disp_x_left]
    type = DirichletBC
    variable = disp_x
    boundary = left
    value = 0
  []
  [disp_x_right]
    type = FunctionDirichletBC
    variable = disp_x
    boundary = right
    function = t
    preset = false
  []
  [disp_y]
    type = DirichletBC
    variable = disp_y
    boundary = bottom
    value = 0
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
  automatic_scaling = true
  petsc_options_iname = '-pc_type'
  petsc_options_value = 'lu'
  dt = 1e-3
  dtmin = 1e-3
  num_steps = 5
  nl_abs_tol = 1e-10
  residual_and_jacobian_together = true
[]

[Postprocessors]
  [ux_right]
    type = PointValue
    variable = disp_x
    point = '1 0.5 0'
    execute_on = TIMESTEP_END
  []
  [ux_center]
    type = PointValue
    variable = disp_x
    point = '0.5 0.5 0'
    execute_on = TIMESTEP_END
  []
[]

[Outputs]
  exodus = false
  csv = true
[]
