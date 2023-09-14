# Total number of blocks = nb**3
nb = 1

# Elements in nb = 1 problem = 2*ne**3
ne = 2

# So the total number of elements is nb**3 * 2 * ne**3
# The idea is to keep the number of MPI ranks the same as the
# number of blocks...

# Block 1 properties
E1 = 100000.0
nu1 = 0.3 
Sy1 = 150.0
H1 = 5000.0

# Block 2 properties
E2 = 150000.0
nu2 = 0.25
Sy2 = 200.0
H2 = 7500.0

[Mesh]
  [msh]
    type = CartesianMeshGenerator
    dim = 3
    
    dx = '1'
    dy = '1'
    dz = '0.5 0.5'

    ix = '${ne}'
    iy = '${ne}'
    iz = '${ne} ${ne}'

    subdomain_id = '1 2'
  []
  [tile]
    type = TiledMeshGenerator
    input = msh

    x_tiles = ${nb}
    y_tiles = ${nb}
    z_tiles = ${nb}
  []
[]

[GlobalParams]
  displacements = 'disp_x disp_y disp_z'
  large_kinematics = true
[]

[Modules]
  [TensorMechanics]
    [Master]
      [all]
        strain = FINITE
        add_variables = true
        new_system = true
        formulation = TOTAL
        volumetric_locking_correction = true
      []
    []
  []
[]

[Functions]
  [pfn]
    type = PiecewiseLinear
    x = '0    1'
    y = '0.0 0.1'
  []
[]

[BCs]
  [left]
    type = DirichletBC
    preset = true
    variable = disp_x
    boundary = left
    value = 0.0
  []

  [bottom]
    type = DirichletBC
    preset = true
    variable = disp_y
    boundary = bottom
    value = 0.0
  []

  [back]
    type = DirichletBC
    preset = true
    variable = disp_z
    boundary = back
    value = 0.0
  []

  [front]
    type = FunctionDirichletBC
    variable = disp_z
    boundary = front
    function = pfn
    preset = false
  []
[]

[Materials]
  [elastic_tensor1]
    type = ComputeIsotropicElasticityTensor
    youngs_modulus = ${E1}
    poissons_ratio = ${nu1}
    block = 1
  []
  [flow_stress1]
    type = DerivativeParsedMaterial
    property_name = flow_stress
    expression = '${Sy1}+${H1}*effective_plastic_strain'
    material_property_names = 'effective_plastic_strain'
    additional_derivative_symbols = 'effective_plastic_strain'
    derivative_order = 2
    compute = false
    block = 1
  []

  [elastic_tensor2]
    type = ComputeIsotropicElasticityTensor
    youngs_modulus = ${E2}
    poissons_ratio = ${nu2}
    block = 2
  []
  [flow_stress2]
    type = DerivativeParsedMaterial
    property_name = flow_stress
    expression = '${Sy2}+${H2}*effective_plastic_strain'
    material_property_names = 'effective_plastic_strain'
    additional_derivative_symbols = 'effective_plastic_strain'
    derivative_order = 2
    compute = false
    block = 2
  []

  [compute_stress1]
    type = ComputeSimoHughesJ2PlasticityStress
    flow_stress_material = flow_stress1
    block = 1
  []

  [compute_stress2]
    type = ComputeSimoHughesJ2PlasticityStress
    flow_stress_material = flow_stress2
    block = 2
  []
[]

[Preconditioning]
  [smp]
    type = SMP
    full = true
  []
[]

[Postprocessors]
  [nl_step]
    type = NumNonlinearIterations
  []
  [l_step]
    type = NumLinearIterations
  []
  [nl_total]
    type = CumulativeValuePostprocessor
    postprocessor = nl_step
  []
  [l_total]
    type = CumulativeValuePostprocessor
    postprocessor = l_step
  []
  [failed_steps]
    type = NumFailedTimeSteps
  []
  [elapsed]
    type = PerfGraphData
    section_name = "Root"
    data_type = total
  []
[]

[Executioner]
  type = Transient

  solve_type = 'newton'
  line_search = 'none'

  petsc_options = '-snes_converged_reason -ksp_converged_reason'
  petsc_options_iname = '-pc_type -pc_factor_mat_solver_package -ksp_gmres_restart -pc_hypre_boomeramg_strong_threshold -pc_hypre_boomeramg_interp_type -pc_hypre_boomeramg_coarsen_type -pc_hypre_boomeramg_agg_nl -pc_hypre_boomeramg_agg_num_paths -pc_hypre_boomeramg_truncfactor'
  petsc_options_value = 'hypre boomeramg 100 0.7 ext+i PMIS 4 2 0.4'

  l_max_its = 500
  nl_max_its = 25
  l_tol = 1e-4

  nl_abs_tol = 1e-10
  nl_rel_tol = 1e-8

  end_time = 1.0
  dtmin = 0.05
  dt = 0.05
[]

[Outputs]
  exodus = false
  csv = true
[]
