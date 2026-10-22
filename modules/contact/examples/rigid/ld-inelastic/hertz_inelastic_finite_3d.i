# Example: rigid sphere pressed into a J2-plasticity body, 3D
# quarter-symmetry, LARGE DEFORMATION.
#
# A deformable sphere of radius 2 centered at the origin is pushed down onto a
# fixed rigid sphere of radius 2 centered at (0, -4, 0), so the two touch at
# (0, -2, 0) and the effective Hertz radius is 1.  Only the quarter x > 0,
# z > 0 of the lower hemisphere is modeled, with symmetry conditions on the
# cut planes.
#
# Contact stack: the rigid indenter is a QUAD4 *surface mesh* with no interior
# parent, used directly as the mortar primary through
# `Contact/.../primary_surface_blocks`.  The surface carries the displacement
# variables and is held fixed by Dirichlet conditions on every component.
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

# Top-level mesh-resolution knobs; passed to `nr` on SphereMeshGenerator.
sphere_refinement = 3
indenter_refinement = 3

# The indenter surface is trimmed to the polar cap above this height on the
# unit-center sphere, i.e. a 60 degree cap of planar radius 2*sin(60) = 1.73.
# The Hertz contact radius here is sqrt(R_eff * delta) = sqrt(1 * 0.1) ~ 0.32
# elastically and grows with the plastic zone, so the cap covers it severalfold.
indenter_cap_height = 1.0

[GlobalParams]
  displacements = 'disp_x disp_y disp_z'
[]

[Mesh]
  # --- Deformable body: quarter of the lower hemisphere ---
  [sphere]
    type = SphereMeshGenerator
    radius = 2.0
    nr = ${sphere_refinement}
  []
  [cut_upper_half]
    type = PlaneDeletionGenerator
    input = sphere
    point = '0 0 0'
    normal = '0 1 0'
  []
  [cut_neg_x]
    type = PlaneDeletionGenerator
    input = cut_upper_half
    point = '0 0 0'
    normal = '-1 0 0'
  []
  [cut_neg_z]
    type = PlaneDeletionGenerator
    input = cut_neg_x
    point = '0 0 0'
    normal = '0 0 -1'
  []
  [assign_block]
    type = RenameBlockGenerator
    input = cut_neg_z
    old_block = 0
    new_block = 1
  []
  [rename_curved_ss]
    type = RenameBoundaryGenerator
    input = assign_block
    old_boundary = 0
    new_boundary = 100
  []
  [flat_sidesets]
    type = SideSetsFromNormalsGenerator
    input = rename_curved_ss
    normals = '-1  0  0
                0  1  0
                0  0 -1'
    new_boundary = '1 2 3'
    normal_tol = 1e-6
    fixed_normal = true
  []

  # --- Rigid indenter: a surface mesh with no interior parent ---
  # The outer sides of a ball mesh are converted to a lower-dimensional block
  # and the ball itself is deleted, which leaves QUAD4 surface elements whose
  # normals, inherited from the sides they came from, point away from the ball
  # center and therefore toward the deformable body.
  [indenter_ball]
    type = SphereMeshGenerator
    radius = 2.0
    nr = ${indenter_refinement}
  []
  [indenter_shell]
    type = LowerDBlockFromSidesetGenerator
    input = indenter_ball
    sidesets = '0'
    new_block_id = 2
    new_block_name = indenter
  []
  [indenter_surface]
    type = BlockDeletionGenerator
    input = indenter_shell
    block = 0
    # Keep the surface elements and null their interior parents instead of
    # deleting them along with the ball
    delete_exteriors = false
  []
  [indenter_cap]
    type = PlaneDeletionGenerator
    input = indenter_surface
    point = '0 ${indenter_cap_height} 0'
    normal = '0 -1 0'
  []
  [indenter_cut_neg_x]
    type = PlaneDeletionGenerator
    input = indenter_cap
    point = '0 0 0'
    normal = '-1 0 0'
  []
  [indenter_cut_neg_z]
    type = PlaneDeletionGenerator
    input = indenter_cut_neg_x
    point = '0 0 0'
    normal = '0 0 -1'
  []
  # Place the indenter so that its pole touches the bottom of the body
  [indenter_place]
    type = TransformGenerator
    input = indenter_cut_neg_z
    transform = TRANSLATE
    vector_value = '0 -4 0'
  []

  [combine]
    type = CombinerGenerator
    inputs = 'flat_sidesets indenter_place'
  []
  # The mortar primary boundary is a nodeset on every node of the surface block
  [primary]
    type = ParsedGenerateNodeset
    input = combine
    new_nodeset_name = rigid_surface
    expression = '1'
    included_subdomains = indenter
    allow_distributed_meshes = true
  []
  allow_renumbering = false
[]

[Variables]
  # Displacements live on the body, the indenter surface and the generated
  # mortar block alike; SolidMechanics/QuasiStatic below adds kernels, the
  # strain material and the wrapping stress on block 1 only.
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
    compatibility_mode = true          # auto-wraps user's `stress` via ComputeLagrangianWrappedStress
    decomposition_method = EigenSolution  # -> kinematic_approximation = rashid_eigen on ComputeLagrangianStrain
    volumetric_locking_correction = true  # -> stabilize_strain = true
    block = 1
    # Assembled residual per component saved to aux vars for the
    # `contactor_force` reaction-sum postprocessor below.
    save_in = 'saved_x saved_y saved_z'
  []
[]

[Contact]
  [rigid]
    primary = rigid_surface
    secondary = 100
    primary_surface_blocks = indenter
    formulation = mortar
    model = frictionless
  []
[]

[UserObjects]
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

[AuxVariables]
  # Residual-save targets for the `contactor_force` reaction sum below.
  [saved_x]
    block = 1
  []
  [saved_y]
    block = 1
  []
  [saved_z]
    block = 1
  []
  [plastic_strain_mag]
    order = CONSTANT
    family = MONOMIAL
    block = 1
  []
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
  [plastic_strain_mag]
    type = MaterialRealAux
    property = eff_plastic_strain
    variable = plastic_strain_mag
    execute_on = 'TIMESTEP_END'
    block = 1
  []
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
  # `stress` is provided here; the Physics action auto-wraps it via
  # ComputeLagrangianWrappedStress (compatibility_mode) and creates
  # the ComputeLagrangianStrain material itself.
  [stress]
    type = ComputeMultiPlasticityStress
    plastic_models = j2
    ep_plastic_tolerance = 1e-9
    block = 1
  []
  [eff_plastic_strain]
    type = RankTwoInvariant
    rank_two_tensor = plastic_strain
    property_name = eff_plastic_strain
    invariant = EffectiveStrain
    block = 1
  []
[]

[Functions]
  [top_disp_y]
    type = PiecewiseLinear
    x = '0  1'
    y = '0 -0.1'
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
  # The indenter is rigid and stationary: every displacement component of the
  # surface mesh is prescribed
  [indenter_x]
    type = DirichletBC
    variable = disp_x
    boundary = rigid_surface
    value = 0.0
  []
  [indenter_y]
    type = DirichletBC
    variable = disp_y
    boundary = rigid_surface
    value = 0.0
  []
  [indenter_z]
    type = DirichletBC
    variable = disp_z
    boundary = rigid_surface
    value = 0.0
  []
[]

[Problem]
  # The indenter surface has no physics, only prescribed displacements
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

  nl_rel_tol = 1e-9
  nl_abs_tol = 1e-8
  nl_max_its = 20
  l_max_its = 200

  start_time = 0.0
  end_time   = 1.0

#  [TimeStepper]
#    type = IterationAdaptiveDT
#    dt = 0.005
#    growth_factor = 1.5
#    cutback_factor = 0.5
#    optimal_iterations = 10
#    iteration_window = 2
#  []
  dtmin = 1e-6

  dt = 0.025
[]

[Postprocessors]
  [max_lm]
    type = NodalExtremeValue
    variable = rigid_normal_lm
    block = rigid_secondary_subdomain
    value_type = max
  []
  [max_plastic_strain]
    type = ElementExtremeValue
    variable = plastic_strain_mag
    block = 1
    value_type = max
  []
  [num_nl]
    type = NumNonlinearIterations
  []
  [cumulative_nl]
    type = CumulativeValuePostprocessor
    postprocessor = num_nl
  []
  # Total y-force applied by the indenter at the contact face, obtained
  # by summing the assembled y-residual on the top BC and equilibrium.
  # `save_in` stores R_i = int (gradphi_i : sigma - phi_i*b), which is negative on a
  # top-pushed-down node under compression (sigma_yy < 0), so the raw sum is
  # flipped to give a positive compressive force.
  [contactor_force_raw]
    type = NodalSum
    variable = saved_y
    boundary = 2
    outputs = 'none'
  []
  [contactor_force]
    type = ScalePostprocessor
    value = contactor_force_raw
    scaling_factor = -1.0
  []
  # Indenter is fixed; the "contactor displacement" reported here is
  # the signed applied top displacement (negative going down).  Same
  # convention as the analytic-level-set version of this example.
  [contactor_displacement]
    type = FunctionValuePostprocessor
    function = top_disp_y
  []
[]

[Outputs]
  exodus = true
  csv = true
[]
