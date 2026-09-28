# Host reference for `nearly_incompressible.i` -- same problem via Physics/SolidMechanics with
# volumetric_locking_correction = true (default mode = total).
[GlobalParams]
  displacements = 'disp_x disp_y disp_z'
[]

[Mesh]
  [gmg]
    type = GeneratedMeshGenerator
    dim = 3
    nx = 4
    ny = 4
    nz = 4
  []
[]

[Physics]
  [SolidMechanics]
    [QuasiStatic]
      [all]
        strain = FINITE
        formulation = TOTAL
        new_system = true
        add_variables = true
        volumetric_locking_correction = true
      []
    []
  []
[]

[NEML2]
  eager = true
  input = 'nearly_incompressible_neml2.i'
  output_device = cpu
  [all]
    model = model
    device = cpu
    derivatives = 'neml2_stress deformation_gradient'
  []
[]

[Materials]
  [stress]
    type = ComputeLagrangianStressCustomPK2
    custom_pk2_stress = neml2_stress
    custom_pk2_jacobian = 'dneml2_stress/ddeformation_gradient'
    large_kinematics = true
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
  [ypull]
    type = FunctionDirichletBC
    variable = disp_y
    boundary = right
    function = t
    preset = false
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
  file_base = nearly_incompressible_out
[]
