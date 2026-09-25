[Models]
  [A]
    type = R2ConstantParameter
    value = A_tensor
    parameter = A_value
  []
  [B]
    type = R2ConstantParameter
    value = B_tensor
    parameter = B_value
  []
  [left]
    type = R2Multiplication
    A = A_value
    B = deformation_gradient
    to = left_product
  []
  [right]
    type = R2Multiplication
    A = left_product
    B = B_value
    to = stress
  []
  [model]
    type = ComposedModel
    models = 'A B left right'
  []
[]

[Tensors]
  [A_tensor]
    type = Python
    expr = 'R2(torch.tensor([[1.0, 2.0, 3.0], [5.0, 7.0, 11.0], [13.0, 17.0, 19.0]], dtype=torch.float64))'
  []
  [B_tensor]
    type = Python
    expr = 'R2(torch.tensor([[23.0, 29.0, 31.0], [37.0, 41.0, 43.0], [47.0, 53.0, 59.0]], dtype=torch.float64))'
  []
[]
