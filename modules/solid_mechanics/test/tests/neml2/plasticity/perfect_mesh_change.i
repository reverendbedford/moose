!include perfect.i

[Mesh]
  [initial_blocks]
    type = SubdomainBoundingBoxGenerator
    input = gmg
    block_id = 1
    bottom_left = '0 0.5 0'
    top_right = '1 1 1'
  []
[]

[NEML2]
  [all]
    block = 0
  []
  [moved]
    model = model
    block = 1
    device = cpu
    derivatives = 'neml2_stress neml2_strain'
  []
[]

[AuxVariables]
  [block_indicator]
    family = MONOMIAL
    order = CONSTANT
  []
[]

[AuxKernels]
  [block_indicator]
    type = ParsedAux
    variable = block_indicator
    expression = 'if(t < 0.0015, y, x)'
    use_xyzt = true
    execute_on = 'INITIAL TIMESTEP_BEGIN'
  []
[]

[MeshModifiers]
  [change_block]
    type = CoupledVarThresholdElementSubdomainModifier
    coupled_var = block_indicator
    criterion_type = ABOVE
    threshold = 0.5
    subdomain_id = 1
    complement_subdomain_id = 0
    execute_on = TIMESTEP_BEGIN
  []
[]

[Outputs]
  file_base = perfect_out
[]
