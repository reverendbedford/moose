# Host reference for `fbar_kokkos_vs_host.i`. Same mesh, BCs, NEML2 model, and solver settings,
# but the constitutive stack is the host `Physics/SolidMechanics/QuasiStatic` action with
# `formulation = TOTAL`, `strain = FINITE`, `volumetric_locking_correction = true`. Together with
# the default `volumetric_locking_correction_mode = total`, this selects the same total-mode F-bar
# formulation implemented by the Kokkos path.
#
# On the same problem, converged interior displacements should agree with the Kokkos stabilized
# output to solver tolerance. The tests spec's CSVDiff on the postprocessor columns is the
# acceptance criterion.

N = 2

[GlobalParams]
  displacements = 'disp_x disp_y disp_z'
[]

[Mesh]
  [gmg]
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
  [xfix]
    type = DirichletBC
    variable = disp_x
    boundary = left
    value = 0
    preset = true
  []
  [yfix]
    type = DirichletBC
    variable = disp_y
    boundary = bottom
    value = 0
    preset = true
  []
  [zfix]
    type = DirichletBC
    variable = disp_z
    boundary = back
    value = 0
    preset = true
  []
  [xpull]
    type = FunctionDirichletBC
    variable = disp_x
    boundary = right
    function = t
    preset = true
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
  num_steps = 3
  residual_and_jacobian_together = true
[]

[Postprocessors]
  [u_c0]
    type = PointValue
    variable = disp_x
    point = '0.25 0.25 0.25'
  []
  [u_c1]
    type = PointValue
    variable = disp_x
    point = '0.75 0.25 0.25'
  []
  [u_c2]
    type = PointValue
    variable = disp_x
    point = '0.25 0.75 0.25'
  []
  [u_c3]
    type = PointValue
    variable = disp_x
    point = '0.25 0.25 0.75'
  []
  [u_c4]
    type = PointValue
    variable = disp_x
    point = '0.75 0.75 0.75'
  []
  [uy_c0]
    type = PointValue
    variable = disp_y
    point = '0.5 0.5 0.5'
  []
  [uz_c0]
    type = PointValue
    variable = disp_z
    point = '0.5 0.5 0.5'
  []
[]

[Outputs]
  csv = true
  file_base = fbar_kokkos_vs_host_out
[]
