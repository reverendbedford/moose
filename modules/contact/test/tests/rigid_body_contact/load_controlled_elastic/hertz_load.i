# Load-controlled rigid-sphere-on-elastic-half-space Hertz test, 2D axisym.
#
# Physics is identical to hertz_sphere_elastic/hertz_elastic.i, but the
# contactor is driven by an applied force F(t) instead of by imposing a
# displacement on the top of the material.  Formulation:
#
#   * A Scalar variable `rigid_offset` (s, one DoF) is the sphere's
#     translation along `load_direction = (0, 1, 0)` (upward).
#   * SphereContactor is told about `offset_variable = rigid_offset`, so
#     every query point x is transformed to x - s * load_direction before
#     the SDF/normal is evaluated (rigid body moves up as s grows).
#   * NodalArea computes tributary weights w_i on the contact sideset.
#   * RigidBodyLoadControl enforces
#         F(t) = Sum_i w_i * lambda_i * (n_i . load_direction)
#     with F(t) ramped linearly from 0 to 1e4 over t in [0, 1].
#
# The material's top surface is now free in the load direction (no
# top_deform BC) — the top-face symmetry constraint on r = 0 is still
# there, and the "top" boundary (2) is left free.  The rigid body's
# equilibrium provides the mechanism that develops contact pressure.

[GlobalParams]
  displacements = 'disp_x disp_y'
  large_kinematics = false
[]

[Mesh]
  [file]
    type = FileMeshGenerator
    file = ../../hertz_spherical/hertz_contact_rz.e
  []
  [drop_rigid_indenter]
    type = BlockDeletionGenerator
    input = file
    block = 1000
  []
  [contact_lower]
    type = LowerDBlockFromSidesetGenerator
    input = drop_rigid_indenter
    sidesets = '100'
    new_block_id = 10001
    new_block_name = contact_lower
  []
  coord_type = RZ
  allow_renumbering = false
[]

[UserObjects]
  [sphere]
    type = SphereContactor
    center = '0 -4 0'                   # sphere top at y = -2 (tangent at t = 0)
    radius = 2.0
    offset_variable = rigid_offset
    load_direction = '0 1 0'            # applied force pushes sphere upward into material
  []
  [nodal_area]
    type = NodalArea
    boundary = 100
    variable = nodal_area
    execute_on = 'INITIAL LINEAR'
  []
[]

[Functions]
  [applied_force]
    type = PiecewiseLinear
    x = '0 1'
    y = '0 1.0e4'
  []
[]

[Variables]
  [disp_x]
    block = '1 contact_lower'
  []
  [disp_y]
    block = '1 contact_lower'
  []
  [normal_lm]
    block = contact_lower
  []
  [rigid_offset]
    family = SCALAR
    order = FIRST
    # Initial offset seeds active contact on multiple nodes at t = 0.
    # In 2D-axisym the r = 0 tip node has zero tributary area (NodalArea
    # includes the 2*pi*r factor), so it is decoupled from the load-control
    # equation; only nodes at r > 0 contribute.  s0 = 0.05 activates the
    # first several off-axis nodes, giving Newton enough coupling between
    # `s` and the LM DoFs at step 0.
    initial_condition = 0.05
  []
[]

[AuxVariables]
  [bounds_dummy]
    family = LAGRANGE
    order = FIRST
    block = contact_lower
  []
  [nodal_area]
    family = LAGRANGE
    order = FIRST
  []
[]

[Bounds]
  [lm_lo]
    type = ConstantBounds
    variable = bounds_dummy
    bounded_variable = normal_lm
    bound_type = lower
    bound_value = 0.0
  []
  [lm_hi]
    type = ConstantBounds
    variable = bounds_dummy
    bounded_variable = normal_lm
    bound_type = upper
    bound_value = 1e12
  []
[]

[Kernels]
  [sdx]
    type = TotalLagrangianStressDivergenceAxisymmetricCylindrical
    variable = disp_x
    component = 0
    block = 1
  []
  [sdy]
    type = TotalLagrangianStressDivergenceAxisymmetricCylindrical
    variable = disp_y
    component = 1
    block = 1
  []
[]

[NodalKernels]
  [ncp]
    type = RigidBodyNodalNCPKernel
    variable = normal_lm
    contactor = sphere
    displacements = 'disp_x disp_y'
    block = contact_lower
  []
[]

[ScalarKernels]
  [load_control]
    type = RigidBodyLoadControl
    variable = rigid_offset
    boundary = 100
    force = applied_force
    contactor = sphere
    nodal_area = nodal_area
    lm_variable = normal_lm
    displacements = 'disp_x disp_y'
    c = 1.0
    # Physical spring reacting the rigid body's offset (F_spring = -K*s
    # along load_direction).  Regularizes the otherwise-singular scalar
    # column of the Jacobian and keeps Newton from escaping the contact
    # region during the first few load steps.  Chosen to be large enough
    # that Newton cannot walk `s` past the material (K * s_escape >> F_peak)
    # but small enough that the fixed-point solution barely shifts: at
    # F = 1e4 and equilibrium s ≈ 5e-2, spring contributes 5e4 to R_s,
    # meaning ~10 % of the applied force ends up carried by the spring
    # rather than by contact.  In production this would be tuned down as
    # the equilibrium is approached, or replaced with a dynamic-relaxation
    # scheme; document that clearly.
    spring_stiffness = 1.0e6
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
  [strain]
    type = ComputeLagrangianStrainAxisymmetricCylindrical
    block = 1
  []
[]

[BCs]
  [rb_tx]
    type = RigidBodyNormalMechanicalContact
    variable = disp_x
    lowerd_variable = normal_lm
    boundary = 100
    contactor = sphere
    component = x
    displacements = 'disp_x disp_y'
  []
  [rb_ty]
    type = RigidBodyNormalMechanicalContact
    variable = disp_y
    lowerd_variable = normal_lm
    boundary = 100
    contactor = sphere
    component = y
    displacements = 'disp_x disp_y'
  []
  [symm_x]
    type = DirichletBC
    variable = disp_x
    boundary = 1
    value = 0.0
  []
  [pin_top]
    type = DirichletBC
    variable = disp_y
    boundary = 2                        # holds the material in place; rigid body
    value = 0.0                         # moves upward under the applied force
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

  # Adaptive dt: load-controlled contact off cold-start has to build up
  # both the rigid-body offset and the LM field from zero, which needs a
  # gentle first step so Newton can find equilibrium at F(t).
  [TimeStepper]
    type = IterationAdaptiveDT
    dt = 0.01
    growth_factor = 2.0
    cutback_factor = 0.5
    optimal_iterations = 10
    iteration_window = 4
  []
[]

[Postprocessors]
  [applied_force]
    type = FunctionValuePostprocessor
    function = applied_force
    execute_on = 'INITIAL TIMESTEP_END'
  []
  [offset]
    type = ScalarVariable
    variable = rigid_offset
    execute_on = 'INITIAL TIMESTEP_END'
  []
  [max_lm]
    type = NodalExtremeValue
    variable = normal_lm
    block = contact_lower
    value_type = max
    execute_on = 'INITIAL TIMESTEP_END'
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
  [csv]
    type = CSV
    execute_on = 'TIMESTEP_END'
  []
[]
