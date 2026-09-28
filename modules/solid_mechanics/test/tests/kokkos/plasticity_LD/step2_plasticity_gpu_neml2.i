N = 16

# Step 2 (LD counterpart of plasticity/step2): GPU NEML2 constitutive update,
# CPU MOOSE assembly, and CPU PETSc.
#
# Differences vs plasticity/step2_plasticity_gpu_neml2.i mirror the ones for
# step1: FINITE + TOTAL formulation, CP exact-kinematics NEML2 model,
# ComputeLagrangianStressCustomPK2, stateful initial conditions,
# volumetric locking correction. The only change vs step1_LD is
# `device = cuda` (so NEML2 runs on the GPU; output_device = cpu forces the
# results back to the host for the CPU MOOSE materials).

[GlobalParams]
  displacements = 'disp_x disp_y disp_z'
[]

[Mesh]
  [generated]
    type = GeneratedMeshGenerator
    dim = 3
    nx = ${N}
    ny = ${N}
    nz = ${N}
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
  input = '../../neml2/crystal_plasticity/exact_kinematics_neml2.i'
  output_device = cpu
  [all]
    model = model
    device = cuda
    derivatives = 'neml2_stress deformation_gradient'
    initialize_outputs = 'plastic_deformation_gradient'
    initialize_output_values = 'initial_plastic_defgrad'
  []
[]

[Materials]
  [stress]
    type = ComputeLagrangianStressCustomPK2
    custom_pk2_stress = neml2_stress
    custom_pk2_jacobian = 'dneml2_stress/ddeformation_gradient'
    large_kinematics = true
  []
  [initial_orientation]
    type = GenericConstantRealVectorValue
    vector_name = orientation
    vector_values = '-0.54412095 -0.34931944 0.12600655'
  []
  [initial_plastic_defgrad]
    type = GenericConstantRankTwoTensor
    tensor_name = initial_plastic_defgrad
    tensor_values = '1 1 1'
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
  [disp_z]
    type = DirichletBC
    variable = disp_z
    boundary = back
    value = 0
  []
[]

[Preconditioning]
  [smp]
    type = SMP
    full = true
  []
[]

# hypre boomeramg + gmres, line_search=none: see step1_LD for the rationale.
[Executioner]
  type = Transient
  solve_type = NEWTON
  line_search = none
  petsc_options_iname = '-pc_type -pc_hypre_type -ksp_type'
  petsc_options_value = 'hypre boomeramg    gmres'
  dt = 1e-3
  dtmin = 1e-3
  num_steps = 5
  nl_abs_tol = 1e-10
  residual_and_jacobian_together = true
  automatic_scaling = true
[]

[Postprocessors]
  [ux_right]
    type = PointValue
    variable = disp_x
    point = '1 0.5 0.5'
    execute_on = TIMESTEP_END
  []
  [ux_center]
    type = PointValue
    variable = disp_x
    point = '0.5 0.5 0.5'
    execute_on = TIMESTEP_END
  []
[]

[Outputs]
  exodus = false
  csv = true
[]
