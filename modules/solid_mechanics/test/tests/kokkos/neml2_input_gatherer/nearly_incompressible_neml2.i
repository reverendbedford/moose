# Nearly-incompressible St. Venant-Kirchhoff elastic model for the F-bar locking-correction
# verification. Poisson's ratio 0.4999 makes the mesh acutely sensitive to volumetric locking with
# linear hex elements. Same PK2/deformation-gradient interface as
# ../../neml2/crystal_plasticity/exact_kinematics_neml2.i so it drops into the same MOOSE stack.
[Models]
  [gl_strain]
    type = GreenLagrangeStrain
    deformation_gradient = 'deformation_gradient'
    strain = 'elastic_strain'
  []
  [svk]
    type = LinearIsotropicElasticity
    coefficients = '1e5 0.49'
    coefficient_types = 'YOUNGS_MODULUS POISSONS_RATIO'
    strain = 'elastic_strain'
    stress = 'sr2_stress'
  []
  [full_stress]
    type = SR2ToR2
    input = 'sr2_stress'
    output = 'neml2_stress'
  []
  [model]
    type = ComposedModel
    models = 'gl_strain svk full_stress'
  []
[]
