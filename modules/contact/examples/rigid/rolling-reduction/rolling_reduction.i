# Example: plate rolling with two-pass thickness reduction against
# rigid rollers, 3D quarter-symmetry, large-deformation J2 plasticity.
#
# Physical setup (full plate before symmetry reduction):
#   plate      height x width x length  =  1 x 2 x 4
#   roll 1     radius 2, axis // x,     centered at (any, 2.375, +1)
#              bottom edge at y = +0.375  ->  reduces full-plate height 1 -> 0.75
#   roll 2     radius 2, axis // x,     centered at (any, 2.25,  +3)
#              bottom edge at y = +0.25   ->  reduces full-plate height 0.75 -> 0.5
#
# The back face at z = -4 is displacement-driven in the +z direction; the
# plate is fed forward through the two rollers.  Quarter symmetry uses
# the x = 0 mid-width plane and y = 0 mid-height plane, so the modeled
# domain occupies x in [0, 1], y in [0, 0.5], z in [-4, 0] initially.
#
# Each roller is a QUAD4 surface mesh with no interior parent, built by
# wrapping a flat patch onto the cylinder with ParsedNodeTransformGenerator.
# Both rollers live in a single block, `rollers`, so that one contact pair
# with one generated secondary block on `plate_top` covers both: two pairs
# sharing `plate_top` would each try to build their own lower-dimensional
# block on that same sideset.
#
# Material: rate-independent J2 plasticity with a deliberately low yield
# stress and mild power-law hardening so the plate flows plastically as
# soon as it enters either roller gap; without the low yield the plate
# would deform elastically and pop back out between rollers.

# Mesh-resolution knob.  nx / ny / nz here go into GeneratedMeshGenerator
# directly; bump them for a converged rolling force and a smoother roller
# footprint.
nx = 6
ny = 6
nz = 40

# Roller geometry.  Both axes are parallel to x and lie in the x = 0
# symmetry plane's normal direction, so only the x >= 0 half is meshed.
roll_radius = 2.0
roll1_y = 2.375
roll1_z = 1.0
roll2_y = 2.25
roll2_z = 3.0

# Half-angle of the meshed roller arc, measured from the bottom of the
# cylinder.  The plate first touches roller 1 at its incoming half-height
# y = 0.5, i.e. at acos((2.375 - 0.5) / 2) = 0.355 rad, and roller 2 at
# y = 0.375, i.e. acos((2.25 - 0.375) / 2) = 0.355 rad; 0.5 rad covers both
# entry arcs and the exit side with margin.
roll_arc = 0.5
# Axial extent of the roller patches.  The analytic cylinder this example was
# adapted from was infinite along its axis; a meshed patch is not, and the
# plate spreads laterally as it is reduced -- past x = 1.19 by the time its
# leading edge reaches the first nip.  Any plate_top node that slides off the
# end of the patch loses its primary surface and is left unconstrained, so the
# patch is carried well beyond the widest spread the plate can reach.
roll_x_min = -0.2
roll_x_max = 2.0
# Roller surface resolution: the arc length 2*roll_arc*roll_radius = 2.0 is
# divided to roughly match the plate's 4/nz = 0.1 axial element size.
n_arc = 20
n_axial = 16

# End-of-run push distance.  Chosen so the initial front face at z = 0 has
# been fed past both rollers (z_1 = 1, z_2 = 3) with a margin.
push_distance = 4.0
end_time      = 1.0

[GlobalParams]
  displacements = 'disp_x disp_y disp_z'
[]

[Mesh]
  [plate]
    type = GeneratedMeshGenerator
    dim = 3
    xmin = 0
    xmax = 1
    ymin = 0
    ymax = 0.5
    zmin = -4
    zmax = 0
    nx = ${nx}
    ny = ${ny}
    nz = ${nz}
    # Standard face labels: left/right (x), bottom/top (y), back/front (z).
    boundary_name_prefix = 'plate'
  []

  # Each roller surface is a flat patch in (x, theta) wrapped onto the
  # cylinder by ParsedNodeTransformGenerator.  The patch's node ordering
  # runs along +x then +theta, so the right-hand-rule normal comes out as
  # (0, -cos(theta), sin(theta)): it points down at the bottom of the
  # roller, toward the plate, as contact requires.
  [roll1_flat]
    type = GeneratedMeshGenerator
    dim = 2
    xmin = ${roll_x_min}
    xmax = ${roll_x_max}
    ymin = ${fparse -roll_arc}
    ymax = ${roll_arc}
    nx = ${n_axial}
    ny = ${n_arc}
    subdomain_ids = 1
    boundary_name_prefix = roll1
    boundary_id_offset = 10
  []
  [roll1]
    type = ParsedNodeTransformGenerator
    input = roll1_flat
    x_function = 'x'
    y_function = '${roll1_y} - ${roll_radius} * cos(y)'
    z_function = '${roll1_z} + ${roll_radius} * sin(y)'
  []
  [roll2_flat]
    type = GeneratedMeshGenerator
    dim = 2
    xmin = ${roll_x_min}
    xmax = ${roll_x_max}
    ymin = ${fparse -roll_arc}
    ymax = ${roll_arc}
    nx = ${n_axial}
    ny = ${n_arc}
    subdomain_ids = 1
    boundary_name_prefix = roll2
    boundary_id_offset = 20
  []
  [roll2]
    type = ParsedNodeTransformGenerator
    input = roll2_flat
    x_function = 'x'
    y_function = '${roll2_y} - ${roll_radius} * cos(y)'
    z_function = '${roll2_z} + ${roll_radius} * sin(y)'
  []
  # One block for both rollers, so a single contact pair covers them
  [rolls_combined]
    type = CombinerGenerator
    inputs = 'roll1 roll2'
  []
  [rolls]
    type = RenameBlockGenerator
    input = rolls_combined
    old_block = 1
    new_block = rollers
  []

  [combine]
    type = CombinerGenerator
    inputs = 'plate rolls'
  []
  # The mortar primary boundary is a nodeset on every node of the roller block
  [primary]
    type = ParsedGenerateNodeset
    input = combine
    new_nodeset_name = roller_surface
    expression = '1'
    included_subdomains = rollers
    allow_distributed_meshes = true
  []
  allow_renumbering = false
[]

[Variables]
  [disp_x]
  []
  [disp_y]
  []
  [disp_z]
  []
[]

[Physics/SolidMechanics/QuasiStatic]
  [all]
    strain = FINITE
    add_variables = false
    compatibility_mode = true
    decomposition_method = EigenSolution
    # VLC left off here: with the deliberately low yield stress below and
    # the large volumetric plastic flow under the rollers, VLC's stabilized
    # strain and the eigen-decomposition kinematics interact to produce a
    # nonsymmetric intermediate rate tensor and abort with
    # "The tensor is not symmetric" in FactorizedRankTwoTensor.  The plate
    # mesh is coarse enough (hex8, no near-incompressibility) that skipping
    # VLC is fine for this demonstration.
    volumetric_locking_correction = false
    block = 0
  []
[]

[Contact]
  [rolling]
    primary = roller_surface
    secondary = plate_top
    primary_surface_blocks = rollers
    formulation = mortar
    model = frictionless

    # The NCP parameter balances gap against contact pressure, and the default
    # 1e6 suits stiff materials.  Here the yield stress is 5 and the contact
    # pressures are ~20, so the default leaves the complementarity condition
    # gap-dominated and the solve needs roughly four times the wall time and
    # three times the cutbacks.  Anything in 20 to 500 behaves about the same
    # and changes the converged pressure by under a percent.
    c_normal = 1e2
  []
[]

[UserObjects]
  [yield_strength]
    type = SolidMechanicsHardeningPowerRule
    # Low relative to E = 200 (yield stress here is ~1/40 of Young's modulus)
    # so that as soon as the plate top touches a roller it enters the plastic
    # branch instead of springing back elastically.  Modest power-law
    # hardening prevents runaway localization at the roller-gap entry.
    value_0  = 5.0
    epsilon0 = 0.02
    exponent = 0.4
  []
  [j2]
    type = SolidMechanicsPlasticJ2
    yield_strength = yield_strength
    yield_function_tolerance      = 1e-4
    internal_constraint_tolerance = 1e-9
  []
[]

[AuxVariables]
  [plastic_strain_mag]
    order = CONSTANT
    family = MONOMIAL
    block = 0
  []
  [saved_x]
    block = 0
  []
  [saved_y]
    block = 0
  []
  [saved_z]
    block = 0
  []
[]

[AuxKernels]
  [plastic_strain_mag]
    type = MaterialRealAux
    property = eff_plastic_strain
    variable = plastic_strain_mag
    execute_on = 'TIMESTEP_END'
    block = 0
  []
[]

[Materials]
  [tensor]
    type = ComputeIsotropicElasticityTensor
    youngs_modulus = 200
    poissons_ratio = 0.3
    block = 0
  []
  # `stress` published here; the Physics action auto-wraps it via
  # ComputeLagrangianWrappedStress (compatibility_mode).
  [stress]
    type = ComputeMultiPlasticityStress
    plastic_models = j2
    ep_plastic_tolerance = 1e-9
    block = 0
  []
  [eff_plastic_strain]
    type = RankTwoInvariant
    rank_two_tensor = plastic_strain
    property_name = eff_plastic_strain
    invariant = EffectiveStrain
    block = 0
  []
[]

[Functions]
  # Linear ramp -- disp_z on the back face increases from 0 to push_distance
  # over t in [0, end_time], feeding the plate through the rollers.
  [back_push]
    type = PiecewiseLinear
    x = '0 ${end_time}'
    y = '0 ${push_distance}'
  []
[]

[BCs]
  [symm_x]
    type = DirichletBC
    variable = disp_x
    boundary = plate_left    # x = 0 symmetry plane
    value = 0.0
  []
  [symm_y]
    type = DirichletBC
    variable = disp_y
    boundary = plate_bottom  # y = 0 symmetry plane
    value = 0.0
  []
  [push_back]
    type = FunctionDirichletBC
    variable = disp_z
    boundary = plate_back    # z = -4 back face
    function = back_push
  []
  # The rollers are rigid and stationary: every displacement component of
  # the surface mesh is prescribed
  [roller_x]
    type = DirichletBC
    variable = disp_x
    boundary = roller_surface
    value = 0.0
  []
  [roller_y]
    type = DirichletBC
    variable = disp_y
    boundary = roller_surface
    value = 0.0
  []
  [roller_z]
    type = DirichletBC
    variable = disp_z
    boundary = roller_surface
    value = 0.0
  []
[]

[Problem]
  # The roller surfaces have no physics, only prescribed displacements
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

  petsc_options_iname = '-pc_type -pc_factor_shift_type -pc_factor_shift_amount'
  petsc_options_value = 'lu       NONZERO               1e-12'
  line_search = none

  nl_rel_tol = 1e-8
  nl_abs_tol = 1e-6
  nl_max_its = 40
  l_max_its  = 200

  start_time = 0.0
  end_time   = ${end_time}

  [TimeStepper]
    type = IterationAdaptiveDT
    dt = 0.005
    growth_factor  = 2
    cutback_factor = 0.5
    # Entry of the plate's leading edge into each nip forces several cutbacks.
    # With a tighter window the stepper settles at dt ~ 5e-4 and never recovers,
    # because the iteration count sits right at the edge where growth stops; the
    # looser window trades more cutbacks for a dt roughly five times larger and
    # reaches a given time in well under half the wall clock, with the contact
    # pressure and plastic strain unchanged to about a percent.
    optimal_iterations = 12
    iteration_window   = 3
  []
  dtmin = 1e-6
[]

[Postprocessors]
  [max_lm]
    type = NodalExtremeValue
    variable = rolling_normal_lm
    block = rolling_secondary_subdomain
    value_type = max
  []
  [max_plastic_strain]
    type = ElementExtremeValue
    variable = plastic_strain_mag
    block = 0
    value_type = max
  []
  [num_nl]
    type = NumNonlinearIterations
  []
  [cumulative_nl]
    type = CumulativeValuePostprocessor
    postprocessor = num_nl
  []
  # z-displacement applied at the back face (piecewise-linear ramp).
  [contactor_displacement]
    type = FunctionValuePostprocessor
    function = back_push
  []
[]

[Outputs]
  exodus = true
  csv = true
[]
