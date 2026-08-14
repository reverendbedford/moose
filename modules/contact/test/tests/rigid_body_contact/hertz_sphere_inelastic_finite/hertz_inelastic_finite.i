# Rigid sphere pressed into an inelastic half-space (J2 radial-return
# plasticity), 2D axisymmetric, LARGE DEFORMATION.
#
# Same mesh, mortar contact setup, and constitutive path as the small-strain
# inelastic test (rigid_body_contact/hertz_sphere_inelastic/), but:
#   * large_kinematics = true (kernel now uses the deformation gradient F
#     properly),
#   * use_displaced_mesh = true on the three mortar constraints so the
#     mortar-segment mesh sees the current-config positions of both bodies,
#   * ComputeLagrangianWrappedStress inherits ComputeLagrangianObjectiveStress
#     and applies its objective rate (Truesdell, the default) to advance the
#     Cauchy stress consistently with F.
#
# Uses TotalLagrangianStressDivergenceAxisymmetricCylindrical because the
# UpdatedLagrangian family in solid_mechanics has no axisymmetric variant; TL
# with large_kinematics = true is the same "new-Lagrangian" family and is
# large-deformation-correct.

[GlobalParams]
  displacements = 'disp_x disp_y'
  large_kinematics = true
[]

[Mesh]
  [file]
    type = FileMeshGenerator
    file = ../../hertz_spherical/hertz_contact_rz.e
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
  [refine]
    type = RefineBlockGenerator
    input = primary_lower
    block = '1 1000 secondary_lower primary_lower'
    refinement = '1 3 1 3'
  []
  [rigid_all_nodes]
    type = ParsedGenerateNodeset
    input = refine
    expression = '1'
    included_subdomains = '1000'
    new_nodeset_name = rigid_all_nodes
  []
  coord_type = RZ
  allow_renumbering = false
[]

[Variables]
  [disp_x]
    block = '1 1000'
  []
  [disp_y]
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
[]

[AuxKernels]
  [plastic_strain_mag]
    type = MaterialRealAux
    property = effective_plastic_strain
    variable = plastic_strain_mag
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
    type = TotalLagrangianStressDivergenceAxisymmetricCylindrical
    variable = disp_x
    component = 0
    block = 1
  []
  [sdy_deform]
    type = TotalLagrangianStressDivergenceAxisymmetricCylindrical
    variable = disp_y
    component = 1
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
    objective_rate = truesdell
    block = 1
  []
  [wrapped_deform]
    type = ComputeMultipleInelasticStress
    inelastic_models = 'j2'
    tangent_operator = elastic
    block = 1
  []
  [j2]
    type = IsotropicPlasticityStressUpdate
    yield_stress = 2.0e5
    hardening_constant = 1.0e6
    block = 1
  []
  [strain_deform]
    type = ComputeLagrangianStrainAxisymmetricCylindrical
    block = 1
  []
[]

[UserObjects]
  [weighted_gap_uo]
    type = LMWeightedGapUserObject
    primary_boundary = 1000
    secondary_boundary = 100
    primary_subdomain = primary_lower
    secondary_subdomain = secondary_lower
    lm_variable = normal_lm
    disp_x = disp_x
    disp_y = disp_y
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
[]

[Functions]
  [top_disp_y]
    type = PiecewiseLinear
    x = '0  1'
    y = '0 -0.05'   # 5x deeper indentation than small-strain tests -> genuinely large-def
  []
[]

[BCs]
  [symm_x_deform]
    type = DirichletBC
    variable = disp_x
    boundary = 1
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

  start_time = 0.0
  end_time   = 1.0
  dt         = 0.1
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

[VectorPostprocessors]
  [contact_lm]
    type = NodalValueSampler
    variable = normal_lm
    block = secondary_lower
    sort_by = x
    execute_on = 'TIMESTEP_END'
  []
[]

[Outputs]
  [csv]
    type = CSV
    execute_on = 'TIMESTEP_END'
  []
[]
