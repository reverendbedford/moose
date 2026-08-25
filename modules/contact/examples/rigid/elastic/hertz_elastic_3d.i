# Example: rigid sphere pressed into an elastic body, 3D quarter-symmetry,
# small strain.  Analytic level-set contact stack via `[RigidContact]`:
#
#   * SphereContactor supplies g_LS(x) = |x - c| - R and its normal.
#   * The action expands into: LowerDBlockFromSidesetGenerator,
#     RigidBodyContactSparsity, normal_lm variable + bounds,
#     RigidBodyNodalNCPKernel, RigidBodyNormalMechanicalContact (per
#     component), problem coverage flags, SMP preconditioning.
#   * PETSc SNESVINEWTONSSLS + ConstantBounds enforces lambda >= 0.
#
# Geometry (quarter of a sphere-on-sphere Hertz setup, symmetry planes
# at x = 0 and z = 0):
#   subdomain 1 = deformable quarter-sphere, radius 2, curved bottom on
#                 sideset 100.  Mesh's rigid indenter (subdomain 1000) is
#                 stripped; analytic sphere replaces it.
#   sideset 2   = top surface of deformable body, driven by function BC.
#   sidesets 1,3 = symmetry planes.
#
# Analytical Hertz (E = 1.40625e7, nu = 0.25, so E* = 1.5e7, R_eff = 1):
#   depth d = 0.01, contact radius a = sqrt(R d) = 0.1,
#   peak pressure p0 = 2 E* a / (pi R) = 9.55e5.

[GlobalParams]
  displacements = 'disp_x disp_y disp_z'
  large_kinematics = false
[]

[Mesh]
  [file]
    type = FileMeshGenerator
    file = ../../../test/tests/hertz_spherical/hertz_contact.e
  []
  [drop_rigid_indenter]
    type = BlockDeletionGenerator
    input = file
    block = 1000
  []
  allow_renumbering = false
[]

[Variables]
  # Disp variables span the lower-d block too, so RigidContact's NCP
  # kernel and BCs find them via block subset.  SolidMechanics adds
  # kernels + strain material on block 1 only.
  [disp_x]
    block = '1 contact_lower'
  []
  [disp_y]
    block = '1 contact_lower'
  []
  [disp_z]
    block = '1 contact_lower'
  []
[]

[Physics/SolidMechanics/QuasiStatic]
  [all]
    strain = SMALL
    add_variables = false
    new_system = true
    formulation = TOTAL
    block = 1
  []
[]

[UserObjects]
  [sphere]
    type = SphereContactor
    center = '0 -4 0'
    radius = 2.0
  []
[]

[RigidContact]
  [top]
    contactor = sphere
    boundary  = 100
    displacements = 'disp_x disp_y disp_z'
    lm_variable_name   = normal_lm
    lower_d_block_name = contact_lower
  []
[]

[AuxVariables]
  [stress_xx]
    order = CONSTANT
    family = MONOMIAL
    block = 1
  []
  [stress_yy]
    order = CONSTANT
    family = MONOMIAL
    block = 1
  []
  [stress_zz]
    order = CONSTANT
    family = MONOMIAL
    block = 1
  []
  [stress_xy]
    order = CONSTANT
    family = MONOMIAL
    block = 1
  []
  [stress_xz]
    order = CONSTANT
    family = MONOMIAL
    block = 1
  []
  [stress_yz]
    order = CONSTANT
    family = MONOMIAL
    block = 1
  []
[]

[AuxKernels]
  [stress_xx]
    type = RankTwoAux
    rank_two_tensor = cauchy_stress
    variable = stress_xx
    index_i = 0
    index_j = 0
    execute_on = 'TIMESTEP_END'
    block = 1
  []
  [stress_yy]
    type = RankTwoAux
    rank_two_tensor = cauchy_stress
    variable = stress_yy
    index_i = 1
    index_j = 1
    execute_on = 'TIMESTEP_END'
    block = 1
  []
  [stress_zz]
    type = RankTwoAux
    rank_two_tensor = cauchy_stress
    variable = stress_zz
    index_i = 2
    index_j = 2
    execute_on = 'TIMESTEP_END'
    block = 1
  []
  [stress_xy]
    type = RankTwoAux
    rank_two_tensor = cauchy_stress
    variable = stress_xy
    index_i = 0
    index_j = 1
    execute_on = 'TIMESTEP_END'
    block = 1
  []
  [stress_xz]
    type = RankTwoAux
    rank_two_tensor = cauchy_stress
    variable = stress_xz
    index_i = 0
    index_j = 2
    execute_on = 'TIMESTEP_END'
    block = 1
  []
  [stress_yz]
    type = RankTwoAux
    rank_two_tensor = cauchy_stress
    variable = stress_yz
    index_i = 1
    index_j = 2
    execute_on = 'TIMESTEP_END'
    block = 1
  []
[]

[Materials]
  [tensor]
    type = ComputeIsotropicElasticityTensor
    youngs_modulus = 1.40625e7
    poissons_ratio = 0.25
    block = 1
  []
  [stress]
    type = ComputeLagrangianLinearElasticStress
    block = 1
  []
[]

[Functions]
  [top_disp_y]
    type = PiecewiseLinear
    x = '0  1'
    y = '0 -0.01'
  []
[]

[BCs]
  [symm_x]
    type = DirichletBC
    variable = disp_x
    boundary = 1
    value = 0.0
  []
  [symm_z]
    type = DirichletBC
    variable = disp_z
    boundary = 3
    value = 0.0
  []
  [top_deform]
    type = FunctionDirichletBC
    variable = disp_y
    boundary = 2
    function = top_disp_y
  []
[]

[Executioner]
  type = Transient
  solve_type = NEWTON
  automatic_scaling = true

  petsc_options_iname = '-snes_type -pc_type -pc_factor_shift_type -pc_factor_shift_amount'
  petsc_options_value = 'vinewtonssls lu    NONZERO               1e-12'
  line_search = semismooth

  nl_rel_tol = 1e-9
  nl_abs_tol = 1e-8
  nl_max_its = 40
  l_max_its = 200

  start_time = 0.0
  end_time   = 1.0
  dt         = 0.1
[]

[Postprocessors]
  [max_lm]
    type = NodalExtremeValue
    variable = normal_lm
    block = contact_lower
    value_type = max
  []
  [num_nl]
    type = NumNonlinearIterations
  []
  [cumulative_nl]
    type = CumulativeValuePostprocessor
    postprocessor = num_nl
  []
[]

[Outputs]
  exodus = true
  csv = true
[]
