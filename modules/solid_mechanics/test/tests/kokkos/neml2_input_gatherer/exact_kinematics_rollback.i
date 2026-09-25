!include exact_kinematics.i

[Postprocessors]
  [dt]
    type = TimestepSize
    execute_on = 'INITIAL NONLINEAR'
  []
[]

[UserObjects]
  [reject_large_trial]
    type = Terminator
    expression = 'dt > 0.005'
    fail_mode = SOFT
    execute_on = NONLINEAR
  []
[]

[Executioner]
  [TimeStepper]
    type = ConstantDT
    dt = 0.01
    cutback_factor_at_failure = 0.5
    growth_factor = 1
  []
[]

[Outputs]
  file_base = exact_kinematics_out
  hide = dt
[]
