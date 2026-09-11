!include perfect.i

[Postprocessors]
  [dt]
    type = TimestepSize
    execute_on = 'INITIAL NONLINEAR'
  []
[]

[UserObjects]
  [reject_large_trial]
    type = Terminator
    expression = 'dt > 0.001'
    fail_mode = SOFT
    execute_on = NONLINEAR
  []
[]

[Executioner]
  [TimeStepper]
    type = ConstantDT
    dt = 0.002
    cutback_factor_at_failure = 0.5
    growth_factor = 1
  []
[]

[Outputs]
  file_base = perfect_out
  hide = dt
[]
