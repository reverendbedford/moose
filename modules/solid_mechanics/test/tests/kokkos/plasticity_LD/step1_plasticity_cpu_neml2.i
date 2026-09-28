N = 16

# Step 1 (LD counterpart of plasticity/step1): CPU NEML2 constitutive update,
# CPU MOOSE assembly, and CPU PETSc.
#
# Differences vs plasticity/step1_plasticity_cpu_neml2.i:
#   * strain = FINITE / formulation = TOTAL exercises the total-Lagrangian path.
#   * The NEML2 model is the crystal-plasticity exact-kinematics model that
#     takes deformation_gradient and returns PK2 stress; ../plasticity/ uses
#     the small-strain J2 perfect_neml2 model. Material parameters do not
#     transfer between the two constitutive models.
#   * ComputeLagrangianStressCustomPK2 converts PK2/dS_dF into PK1 for the host
#     stress-divergence kernel; ../plasticity/ uses
#     ComputeLagrangianObjectiveCustomSymmetricStress for the small-strain stress.
#   * initial_orientation and initial_plastic_defgrad seed the stateful CP model.
#   * volumetric_locking_correction = true matches the host exact_kinematics.i
#     reference in ../../neml2/crystal_plasticity/.

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
    device = cpu
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

# Preconditioner: hypre boomeramg + gmres. SD plasticity/step1 uses gamg+gmres
# for scalable iterative solves; gamg stalls Newton on the non-symmetric CP LD
# tangent, but hypre boomeramg converges in the same Newton count as LU and is
# ~30% faster at N=16. See README.md for the pc comparison table. LU (the LD
# reference in this module) can be recovered with
#   -pc_type lu Executioner/petsc_options_iname=-pc_type Executioner/petsc_options_value=lu
# when an exact reference match is required.
# line_search=none matches the CP host reference; the CP return map already
# handles the material nonlinearity so backtracking is not needed.
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

# Sample points for direct comparison against step3/4/5 (same postprocessors,
# same locations). Used to verify the LD Kokkos assembly matches the host
# reference at machine precision on the smoke case and to a solver tolerance
# on the benchmark case.
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
