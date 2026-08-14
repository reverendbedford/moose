# Hertz-style rigid sphere pressed into an elastic half-space, 2D
# axisymmetric, small strain.
#
# Reuses the mesh from modules/contact/test/tests/hertz_spherical/
# hertz_contact_rz.e (subdomain 1 = deformable body, subdomain 1000 = rigid
# indenter with a very high stiffness).  Lower-d subdomains are added on
# each side of the interface at parse time.
#
# LM contact is enforced through the existing mortar path:
#   * ComputeWeightedGapLMMechanicalContact assembles the LM row on the
#     secondary lower-d block using MOOSE's mortar-segment machinery.
#   * NormalMortarMechanicalContact applies the -lambda*n traction on the
#     coupled displacement equations.
#   * PETSc SNESVINEWTONSSLS drives the complementarity through a lower
#     bound of 0 on the LM variable.
#
# Analytical Hertz for sphere-on-sphere with R1 = R2 = 2, both bodies:
#   E* = E / (2(1-nu^2))    (per body)
# With E = 1.40625e7 (deformable), nu = 0.25:
#   E* = 7.5e6, R = 1
# For depth of indentation d = 0.01:
#   a  = sqrt(R * d) = 0.1
#   p0 = 2 E* a / (pi R) = 4.775e5
#   P  = (4/3) E* R^(1/2) d^(3/2) = 1e4

[GlobalParams]
  displacements = 'disp_x disp_y'
  large_kinematics = false
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
    refinement = '1 3 1 3'     # more refinement on the coarse rigid indenter
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
  [sdx_rigid]
    type = TotalLagrangianStressDivergenceAxisymmetricCylindrical
    variable = disp_x
    component = 0
    block = 1000
  []
  [sdy_rigid]
    type = TotalLagrangianStressDivergenceAxisymmetricCylindrical
    variable = disp_y
    component = 1
    block = 1000
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
    type = ComputeLagrangianLinearElasticStress
    block = 1
  []
  [strain_deform]
    type = ComputeLagrangianStrainAxisymmetricCylindrical
    block = 1
  []
  [elastic_rigid]
    type = ComputeIsotropicElasticityTensor
    youngs_modulus = 1.40625e10   # 1000x stiffer than the deformable body -> effectively rigid
    poissons_ratio = 0.25
    block = 1000
  []
  [stress_rigid]
    type = ComputeLagrangianLinearElasticStress
    block = 1000
  []
  [strain_rigid]
    type = ComputeLagrangianStrainAxisymmetricCylindrical
    block = 1000
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
    weighted_gap_uo = weighted_gap_uo
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
  [symm_x_deform]
    type = DirichletBC
    variable = disp_x
    boundary = 1        # r = 0 symmetry axis (deformable body)
    value = 0.0
  []
  [top_deform_dispy]
    type = FunctionDirichletBC
    variable = disp_y
    boundary = 2        # top of deformable body: pushed down
    function = top_disp_y
  []
  [rigid_fix_x]
    type = DirichletBC
    variable = disp_x
    boundary = 1000     # bottom of rigid indenter: fixed
    value = 0.0
  []
  [rigid_fix_y]
    type = DirichletBC
    variable = disp_y
    boundary = 1000
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
  [num_nl]
    type = NumNonlinearIterations
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
