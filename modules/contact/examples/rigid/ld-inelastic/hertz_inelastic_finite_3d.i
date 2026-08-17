# Example: rigid sphere pressed into a J2-plasticity body, 3D
# quarter-symmetry, LARGE DEFORMATION.
#
# Same mesh + rigid-indenter treatment as the elastic 3D example, but the
# deformable body is now finite-strain J2 plasticity (linear hardening) and
# the load is ramped further to activate a large plastic zone.
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
# The mortar constraints run with use_displaced_mesh = true so the
# mortar-segment mesh sees the current-config positions of both bodies.

[GlobalParams]
  displacements = 'disp_x disp_y disp_z'
  large_kinematics = true
  stabilize_strain = true
[]

[Mesh]
  [file]
    type = FileMeshGenerator
    file = ../../../test/tests/hertz_spherical/hertz_contact.e
  []
  [secondary_lower]
    type = LowerDBlockFromSidesetGenerator
    input = file
    sidesets = '100'
    new_block_id = 10001
    new_block_name = secondary_lower
  []
  [primary_lower]
    type = LowerDBlockFromSidesetGenerator
    input = secondary_lower
    sidesets = '1000'
    new_block_id = 10000
    new_block_name = primary_lower
  []
  [refine_primary]
    type = RefineBlockGenerator
    input = primary_lower
    block = '1000 primary_lower'
    refinement = '2 2'
  []
  [rigid_all_nodes]
    type = ParsedGenerateNodeset
    input = refine_primary
    expression = '1'
    included_subdomains = '1000'
    new_nodeset_name = rigid_all_nodes
  []
  allow_renumbering = false
[]

[Variables]
  [disp_x]
    block = '1 1000'
  []
  [disp_y]
    block = '1 1000'
  []
  [disp_z]
    block = '1 1000'
  []
  [normal_lm]
    block = secondary_lower
    use_dual = true
  []
[]

[AuxVariables]
  [bounds_dummy]
    family = LAGRANGE
    order = FIRST
    block = secondary_lower
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
  [strain_xx]
    order = CONSTANT
    family = MONOMIAL
    block = 1
  []
  [strain_yy]
    order = CONSTANT
    family = MONOMIAL
    block = 1
  []
  [strain_zz]
    order = CONSTANT
    family = MONOMIAL
    block = 1
  []
  [strain_xy]
    order = CONSTANT
    family = MONOMIAL
    block = 1
  []
  [strain_xz]
    order = CONSTANT
    family = MONOMIAL
    block = 1
  []
  [strain_yz]
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
  [strain_xx]
    type = RankTwoAux
    rank_two_tensor = rotated_mechanical_strain
    variable = strain_xx
    index_i = 0
    index_j = 0
    execute_on = 'TIMESTEP_END'
    block = 1
  []
  [strain_yy]
    type = RankTwoAux
    rank_two_tensor = rotated_mechanical_strain
    variable = strain_yy
    index_i = 1
    index_j = 1
    execute_on = 'TIMESTEP_END'
    block = 1
  []
  [strain_zz]
    type = RankTwoAux
    rank_two_tensor = rotated_mechanical_strain
    variable = strain_zz
    index_i = 2
    index_j = 2
    execute_on = 'TIMESTEP_END'
    block = 1
  []
  [strain_xy]
    type = RankTwoAux
    rank_two_tensor = rotated_mechanical_strain
    variable = strain_xy
    index_i = 0
    index_j = 1
    execute_on = 'TIMESTEP_END'
    block = 1
  []
  [strain_xz]
    type = RankTwoAux
    rank_two_tensor = rotated_mechanical_strain
    variable = strain_xz
    index_i = 0
    index_j = 2
    execute_on = 'TIMESTEP_END'
    block = 1
  []
  [strain_yz]
    type = RankTwoAux
    rank_two_tensor = rotated_mechanical_strain
    variable = strain_yz
    index_i = 1
    index_j = 2
    execute_on = 'TIMESTEP_END'
    block = 1
  []
[]

[Bounds]
  [normal_lm_lower]
    type = ConstantBounds
    variable = bounds_dummy
    bounded_variable = normal_lm
    bound_type = lower
    bound_value = 0.0
  []
  [normal_lm_upper]
    type = ConstantBounds
    variable = bounds_dummy
    bounded_variable = normal_lm
    bound_type = upper
    bound_value = 1e12
  []
[]

[Kernels]
  [sdx_deform]
    type = TotalLagrangianStressDivergence
    variable = disp_x
    component = 0
    block = 1
  []
  [sdy_deform]
    type = TotalLagrangianStressDivergence
    variable = disp_y
    component = 1
    block = 1
  []
  [sdz_deform]
    type = TotalLagrangianStressDivergence
    variable = disp_z
    component = 2
    block = 1
  []
[]

[Materials]
  [elastic_deform]
    type = ComputeIsotropicElasticityTensor
    youngs_modulus = 1.40625e7
    poissons_ratio = 0.25
    block = 1
  []
  [stress_deform]
    type = ComputeLagrangianWrappedStress
    objective_rate = rashid
    block = 1
  []
  [wrapped_deform]
    type = ComputeMultiPlasticityStress
    plastic_models = j2
    ep_plastic_tolerance = 1e-9
    block = 1
  []
  [strain_deform]
    type = ComputeLagrangianStrain
    kinematic_approximation = rashid_eigen
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

[UserObjects]
  [yield_strength]
    type = SolidMechanicsHardeningPowerRule
    value_0 = 2.0e5                  # initial yield stress
    epsilon0 = 0.2                   # value = 2e5 * (1 + p/0.2)^1 -> linear hardening slope 1e6
    exponent = 1.0
  []
  [j2]
    type = SolidMechanicsPlasticJ2
    yield_strength = yield_strength
    yield_function_tolerance = 1e-3
    internal_constraint_tolerance = 1e-9
  []
  [weighted_gap_uo]
    type = LMWeightedGapUserObject
    primary_boundary = 1000
    secondary_boundary = 100
    primary_subdomain = primary_lower
    secondary_subdomain = secondary_lower
    lm_variable = normal_lm
    disp_x = disp_x
    disp_y = disp_y
    disp_z = disp_z
  []
[]

[Constraints]
  [weighted_gap_lm]
    type = ComputeWeightedGapLMMechanicalContact
    primary_boundary = 1000
    secondary_boundary = 100
    primary_subdomain = primary_lower
    secondary_subdomain = secondary_lower
    variable = normal_lm
    disp_x = disp_x
    disp_y = disp_y
    disp_z = disp_z
    c = 1
    use_displaced_mesh = true
    weighted_gap_uo = weighted_gap_uo
  []
  [normal_x]
    type = NormalMortarMechanicalContact
    primary_boundary = 1000
    secondary_boundary = 100
    primary_subdomain = primary_lower
    secondary_subdomain = secondary_lower
    variable = normal_lm
    secondary_variable = disp_x
    component = x
    compute_lm_residuals = false
    use_displaced_mesh = true
    weighted_gap_uo = weighted_gap_uo
  []
  [normal_y]
    type = NormalMortarMechanicalContact
    primary_boundary = 1000
    secondary_boundary = 100
    primary_subdomain = primary_lower
    secondary_subdomain = secondary_lower
    variable = normal_lm
    secondary_variable = disp_y
    component = y
    compute_lm_residuals = false
    use_displaced_mesh = true
    weighted_gap_uo = weighted_gap_uo
  []
  [normal_z]
    type = NormalMortarMechanicalContact
    primary_boundary = 1000
    secondary_boundary = 100
    primary_subdomain = primary_lower
    secondary_subdomain = secondary_lower
    variable = normal_lm
    secondary_variable = disp_z
    component = z
    compute_lm_residuals = false
    use_displaced_mesh = true
    weighted_gap_uo = weighted_gap_uo
  []
[]

[Functions]
  [top_disp_y]
    type = PiecewiseLinear
    x = '0  1'
    y = '0 -0.1'                   # 10x deeper indentation than the elastic example -> deep plastic zone
  []
[]

[BCs]
  [symm_x_deform]
    type = DirichletBC
    variable = disp_x
    boundary = 1
    value = 0.0
  []
  [symm_z_deform]
    type = DirichletBC
    variable = disp_z
    boundary = 3
    value = 0.0
  []
  [top_deform_dispy]
    type = FunctionDirichletBC
    variable = disp_y
    boundary = 2
    function = top_disp_y
  []
  [rigid_x]
    type = DirichletBC
    variable = disp_x
    boundary = rigid_all_nodes
    value = 0.0
    preset = true
  []
  [rigid_y]
    type = DirichletBC
    variable = disp_y
    boundary = rigid_all_nodes
    value = 0.0
    preset = true
  []
  [rigid_z]
    type = DirichletBC
    variable = disp_z
    boundary = rigid_all_nodes
    value = 0.0
    preset = true
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

  nl_rel_tol = 1e-9
  nl_abs_tol = 1e-8
  nl_max_its = 40
  l_max_its = 200

  # The semismooth (Fischer-Burmeister) line search absorbs the active-set
  # churn that stalls a plain Newton solve at this indentation depth, so we
  # can take much larger time steps than the earlier `line_search = 'none'`
  # + `dt = 6.25e-3` configuration.
  line_search = semismooth

  start_time = 0.0
  end_time   = 1.0
  dt         = 0.025
[]

[Postprocessors]
  [max_lm]
    type = NodalExtremeValue
    variable = normal_lm
    block = secondary_lower
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
[]

[Outputs]
  exodus = true
  csv = true
[]
