[StochasticTools]
[]

[Samplers]
  [sample]
    type = CartesianProduct
    linear_space_items = '0 1 10
                          0 1 10
                          0 1 10'
  []
  [test]
    type = CartesianProduct
    linear_space_items = '0.25 1 10
                          0.25 1 10
                          0.25 1 10'
  []
[]

[VectorPostprocessors]
  [values]
    type = GFunction
    sampler = sample
    q_vector = '0 0 0'
    execute_on = INITIAL
    outputs = none
  []
  [results]
    type = EvaluateSurrogate
    model = surrogate
    sampler = test
    execute_on = final
  []
[]

[Trainers]
  [train]
    type = NearestPointTrainer
    sampler = sample
    predictors = 'sampler/col_0 sampler/col_1 sampler/col_2'
    response = values/g_values
  []
[]

[Surrogates]
  [surrogate]
    type = NearestPointSurrogate
    trainer = train
  []
[]

[Outputs]
  csv = true
  execute_on = final
[]
