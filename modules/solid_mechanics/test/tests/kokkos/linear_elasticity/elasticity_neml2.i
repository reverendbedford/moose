[Models]
  [model]
    type = LinearIsotropicElasticity
    coefficients = '10 0.25'
    coefficient_types = 'YOUNGS_MODULUS POISSONS_RATIO'
    strain = 'neml2_strain'
    stress = 'neml2_stress'
  []
[]
