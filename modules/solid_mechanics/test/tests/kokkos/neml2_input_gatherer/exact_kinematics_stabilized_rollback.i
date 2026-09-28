# Time-step cutback test for the stabilized (total-mode F-bar) Kokkos LD path.
# Mirrors `exact_kinematics_rollback.i` but on `exact_kinematics_stabilized.i`. Sets a large
# time step (dt = 0.01) with a Terminator that fires whenever dt > 0.005 during any nonlinear
# iteration. The Terminator marks the step as failed, ConstantDT cuts back by cutback_factor
# 0.5, and MOOSE retries with dt = 0.005. If the Kokkos F-bar path correctly discards trial
# NEML2 state on the failed attempt (thanks to `manage_state_advance = true`), the retried step
# reproduces the exact_kinematics_stabilized trajectory bit-for-bit.
!include exact_kinematics_stabilized.i

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
  file_base = exact_kinematics_stabilized_out
  hide = dt
[]
