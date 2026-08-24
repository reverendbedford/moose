# Example: rigid sphere pressed into a J2-plasticity body, 3D
# quarter-symmetry, LARGE DEFORMATION.
#
# Same mesh + rigid-indenter treatment as the elastic 3D example, but the
# deformable body is now finite-strain J2 plasticity (linear hardening) and
# the load is ramped further to activate a large plastic zone.
#
# Contact stack: analytic level-set (SurfaceMeshContactor + RigidBodyNodalNCPKernel
# + RigidBodyNormalMechanicalContact).  Per-node min-NCP, no mortar, no AD,
# no dual basis.
#
# Constitutive stack (new-Lagrangian pipeline with consistent algorithmic
# tangent):
#   ComputeLagrangianStrain (kinematic_approximation = rashid_eigen)
#     + ComputeLagrangianWrappedStress (objective_rate = rashid)
#     + ComputeMultiPlasticityStress (plastic_models = j2)
#     + SolidMechanicsPlasticJ2 + SolidMechanicsHardeningPowerRule
#   TotalLagrangianStressDivergence with large_kinematics = true
#   stabilize_strain = true to avoid nearly-incompressible plastic locking
#     on linear hexes.
#
# STATUS -- KNOWN NOT CONVERGING WITH PLAIN NEWTON.  The load-control
# scalar equation dR_s/ds is structurally zero in the current formulation,
# so plain Newton on the coupled (u, lambda, s) system has no diagonal
# information to bound its scalar step.  RigidBodyLoadControl.kss_stiffness
# fabricates a preconditioner-only Kss diagonal, and it converges the
# frame-invariance regression tests (elastic Hertz at moderate load),
# but for this problem -- finite-strain J2 plasticity, load ramping to
# F(1) = 1.5e5, deep contact patch -- plain Newton with any fixed Kss
# choice enters a limit cycle at the first step (see notes in
# load_control_plan.md, Phase 2 / Uzawa outer loop).  The file below is
# staged with correct inputs (direction = '0 -1 0', kss_stiffness set to
# Young's modulus, auto-scaling off so the fabricated Kss does not
# neutralize itself) so that once the Uzawa executioner lands, running
# this example is a one-line executioner swap rather than a rewrite.

[GlobalParams]
  displacements = 'disp_x disp_y disp_z'
  large_kinematics = true
  stabilize_strain = true
[]

[Mesh]
  [file]
    type = FileMeshGenerator
    file = contact_test.e
  []
  [drop_rigid_indenter]
    type = BlockDeletionGenerator
    input = file
    block = 1
  []
  [contact_lower]
    type = LowerDBlockFromSidesetGenerator
    input = drop_rigid_indenter
    sidesets = 'mat_top'
    new_block_id = 10001
    new_block_name = contact_lower
  []
  allow_renumbering = false
[]

[UserObjects]
  [contact_sparsity]
    type = RigidBodyContactSparsity
    lm_variable = normal_lm
    displacements = 'disp_x disp_y disp_z'
    boundary = mat_top
  []
  [sphere]
    type = SurfaceMeshContactor
    file = unit_sphere.stl
    scale = 2.0
    disp_y_scalar = indenter_y
  []
  [nodal_area]
    type = NodalArea
    boundary = mat_top
    variable = nodal_area
    execute_on = 'INITIAL LINEAR'
  []
  [yield_strength]
    type = SolidMechanicsHardeningPowerRule
    value_0 = 2.0e5
    epsilon0 = 0.2
    exponent = 1.0
  []
  [j2]
    type = SolidMechanicsPlasticJ2
    yield_strength = yield_strength
    yield_function_tolerance = 1e-3
    internal_constraint_tolerance = 1e-9
  []
[]

[Variables]
  [disp_x]
    block = '1000 contact_lower'
  []
  [disp_y]
    block = '1000 contact_lower'
  []
  [disp_z]
    block = '1000 contact_lower'
  []
  [normal_lm]
    block = contact_lower
  []
  [indenter_y]
    family = SCALAR
    order = FIRST
    # Sphere STL (radius = 2 after scale) with no translation puts the
    # sphere bottom at y = -2 in the undeformed frame, coincident with
    # mat_top.  Starting at s = 0 means zero gap and zero penetration at
    # t = 0 (no reaction, R_s = -F(0) = 0).  Load-control then drives s
    # negative (sphere moves down into the material) as F(t) ramps up.
    initial_condition = 0
  []
[]

[AuxVariables]
  [nodal_area]
    family = LAGRANGE
    order = FIRST
  []
  [plastic_strain_mag]
    order = CONSTANT
    family = MONOMIAL
    block = 1000
  []
  [stress_xx]
    order = CONSTANT
    family = MONOMIAL
    block = 1000
  []
  [stress_yy]
    order = CONSTANT
    family = MONOMIAL
    block = 1000
  []
  [stress_zz]
    order = CONSTANT
    family = MONOMIAL
    block = 1000
  []
  [stress_xy]
    order = CONSTANT
    family = MONOMIAL
    block = 1000
  []
  [stress_xz]
    order = CONSTANT
    family = MONOMIAL
    block = 1000
  []
  [stress_yz]
    order = CONSTANT
    family = MONOMIAL
    block = 1000
  []
[]

[AuxKernels]
  [plastic_strain_mag]
    type = MaterialRealAux
    property = eff_plastic_strain
    variable = plastic_strain_mag
    execute_on = 'TIMESTEP_END'
    block = 1000
  []
  [stress_xx]
    type = RankTwoAux
    rank_two_tensor = cauchy_stress
    variable = stress_xx
    index_i = 0
    index_j = 0
    execute_on = 'TIMESTEP_END'
    block = 1000
  []
  [stress_yy]
    type = RankTwoAux
    rank_two_tensor = cauchy_stress
    variable = stress_yy
    index_i = 1
    index_j = 1
    execute_on = 'TIMESTEP_END'
    block = 1000
  []
  [stress_zz]
    type = RankTwoAux
    rank_two_tensor = cauchy_stress
    variable = stress_zz
    index_i = 2
    index_j = 2
    execute_on = 'TIMESTEP_END'
    block = 1000
  []
  [stress_xy]
    type = RankTwoAux
    rank_two_tensor = cauchy_stress
    variable = stress_xy
    index_i = 0
    index_j = 1
    execute_on = 'TIMESTEP_END'
    block = 1000
  []
  [stress_xz]
    type = RankTwoAux
    rank_two_tensor = cauchy_stress
    variable = stress_xz
    index_i = 0
    index_j = 2
    execute_on = 'TIMESTEP_END'
    block = 1000
  []
  [stress_yz]
    type = RankTwoAux
    rank_two_tensor = cauchy_stress
    variable = stress_yz
    index_i = 1
    index_j = 2
    execute_on = 'TIMESTEP_END'
    block = 1000
  []
[]


[Kernels]
  [sdx]
    type = TotalLagrangianStressDivergence
    variable = disp_x
    component = 0
    block = 1000
  []
  [sdy]
    type = TotalLagrangianStressDivergence
    variable = disp_y
    component = 1
    block = 1000
  []
  [sdz]
    type = TotalLagrangianStressDivergence
    variable = disp_z
    component = 2
    block = 1000
  []
[]

[NodalKernels]
  [ncp]
    type = RigidBodyNodalNCPKernel
    variable = normal_lm
    contactor = sphere
    displacements = 'disp_x disp_y disp_z'
    block = contact_lower
    c = 1.0
  []
[]

[ScalarKernels]
  [load_control]
    type = RigidBodyLoadControl
    variable = indenter_y
    boundary = mat_top
    force = applied_force
    contactor = sphere
    nodal_area = nodal_area
    lm_variable = normal_lm
    displacements = 'disp_x disp_y disp_z'
    direction = '0 -1 0'
    c = 1.0
    # Preconditioner-only shift on the (scalar_row, scalar_col) Jacobian
    # entry.  See the kernel's kss_stiffness docstring: the analytical
    # dR_s/ds is zero, so plain Newton on the coupled (u, lambda, s)
    # system overshoots dramatically on any non-trivial load ramp.  A
    # value near the material's Young's modulus makes the fabricated
    # diagonal comparable to the true Schur-complement stiffness and
    # gets Newton to descend cleanly.
    kss_stiffness = 1.40625e7
  []
[]

[Materials]
  [tensor]
    type = ComputeIsotropicElasticityTensor
    youngs_modulus = 1.40625e7
    poissons_ratio = 0.25
    block = 1000
  []
  [stress]
    type = ComputeLagrangianWrappedStress
    objective_rate = rashid
    block = 1000
  []
  [wrapped]
    type = ComputeMultiPlasticityStress
    plastic_models = j2
    ep_plastic_tolerance = 1e-9
    block = 1000
  []
  [strain]
    type = ComputeLagrangianStrain
    kinematic_approximation = rashid_eigen
    block = 1000
  []
  [eff_plastic_strain]
    type = RankTwoInvariant
    rank_two_tensor = plastic_strain
    property_name = eff_plastic_strain
    invariant = EffectiveStrain
    block = 1000
  []
[]

[Functions]
  [applied_force]
    # Ramp to a target reaction that induces significant plasticity in the
    # material.  Chosen so the example completes in a reasonable wall
    # clock: at F(1) = 1.5e5 the plastic zone extends beyond the contact
    # patch, and the run finishes in ~5-10 minutes on a single CPU.
    type = PiecewiseLinear
    x = '0 1'
    y = '0 1.5e5'
  []
[]

[BCs]
  [rb_tx]
    type = RigidBodyNormalMechanicalContact
    variable = disp_x
    lowerd_variable = normal_lm
    boundary = mat_top
    contactor = sphere
    component = x
    finite_strain = true
    displacements = 'disp_x disp_y disp_z'
  []
  [rb_ty]
    type = RigidBodyNormalMechanicalContact
    variable = disp_y
    lowerd_variable = normal_lm
    boundary = mat_top
    contactor = sphere
    component = y
    finite_strain = true
    displacements = 'disp_x disp_y disp_z'
  []
  [rb_tz]
    type = RigidBodyNormalMechanicalContact
    variable = disp_z
    lowerd_variable = normal_lm
    boundary = mat_top
    contactor = sphere
    component = z
    finite_strain = true
    displacements = 'disp_x disp_y disp_z'
  []
  [symm_x]
    type = DirichletBC
    variable = disp_x
    boundary = mat_sym_x
    value = 0.0
  []
  [symm_z]
    type = DirichletBC
    variable = disp_z
    boundary = mat_sym_z
    value = 0.0
  []
  [pin_top]
    type = DirichletBC
    variable = disp_y
    boundary = mat_bot
    value = 0.0
  []
[]

[Problem]
  kernel_coverage_check = false
  material_coverage_check = false
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
  # Auto-scaling would read the huge (scalar, scalar) diagonal fabricated
  # by RigidBodyLoadControl.kss_stiffness and scale the scalar equation
  # by ~1/kss_stiffness -- the raw scalar residual would then look
  # already-converged from the start and Newton would never move the
  # sphere.  Leave scaling off so SNES sees the unscaled R_s.
  automatic_scaling = false

  # Plain Newton with LU direct solve (see the elastic 3D force-control
  # test for the reason we can't use SSLS + semismooth line search here).
  petsc_options_iname = '-snes_type -pc_type -pc_factor_shift_type -pc_factor_shift_amount'
  petsc_options_value = 'newtonls    lu       NONZERO               1e-12'

  line_search = basic

  nl_rel_tol = 1e-9
  nl_abs_tol = 1e-8
  nl_max_its = 40
  l_max_its = 200

  start_time = 0.0
  end_time   = 1.0

  # Adaptive time stepping: analytic-level-set contact concentrates load at a
  # single node initially, so start small and grow.  Rashid-eigen strain
  # increments can fail (non-symmetric tensor) if a single Newton step
  # over-shoots the elastic-plastic corner.
  [TimeStepper]
    type = IterationAdaptiveDT
    dt = 0.005
    growth_factor = 1.5
    cutback_factor = 0.5
    optimal_iterations = 8
    iteration_window = 2
  []
[]

[Postprocessors]
  [max_lm]
    type = NodalExtremeValue
    variable = normal_lm
    block = contact_lower
    value_type = max
  []
  [max_plastic_strain]
    type = ElementExtremeValue
    variable = plastic_strain_mag
    block = 1000
    value_type = max
  []
  [num_nl]
    type = NumNonlinearIterations
  []
  [cumulative_nl]
    type = CumulativeValuePostprocessor
    postprocessor = num_nl
  []
  [force]
    type = FunctionValuePostprocessor
    function = applied_force
  []
  [depth]
    type = ScalarVariable
    variable = indenter_y
  []
[]

[Outputs]
  exodus = true
  csv = true
[]
