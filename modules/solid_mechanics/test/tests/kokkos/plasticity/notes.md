# Why did the CPU NEML2 test take more iterations than the GPU NEML2 test?

**Observation:**
In the plasticity benchmark results (`comparison.csv`), the `step1_plasticity_cpu_neml2` test took 72 Newton iterations and 15,702 KSP iterations, while the `step2_plasticity_gpu_neml2` test (and all subsequent GPU tests) took only 21 Newton iterations and 805 KSP iterations.

**Root Cause:**
This was not caused by a difference between CPU/GPU memory, device execution, or the NEML2 library itself. The CPU and GPU tests are mathematically identical. 

The discrepancy was traced to a single flag in the MOOSE input files: the `preset` parameter on the right boundary condition.

`step1_plasticity_cpu_neml2.i`:
```ini
  [disp_x_right]
    type = FunctionDirichletBC
    variable = disp_x
    boundary = right
    function = t
    preset = true   # <-- The culprit
  []
```

`step2_plasticity_gpu_neml2.i`:
```ini
  [disp_x_right]
    type = FunctionDirichletBC
    variable = disp_x
    boundary = right
    function = t
    preset = false  # <-- The correct configuration
  []
```

**Why `preset = true` causes an iteration spike:**
When `preset = true` is active, the boundary displacement (`u = t`) is aggressively enforced on the initial solution *before* the first residual evaluation. In the first timestep, the boundary nodes are instantly displaced while all internal nodes remain at `u = 0`. 

This creates a massive, artificial strain localized entirely in the single layer of elements adjacent to the boundary. For a nonlinear elastoplastic material model, this huge localized strain immediately triggers aggressive plastic yielding. This degrades the local stiffness matrix and puts the Newton solver in a very stiff, difficult state to recover from.

When `preset = false` is used, the initial solution is uniformly zero. The Newton solver gradually distributes the required boundary displacement across the entire domain over the course of its non-linear iterations. This smooths out the strain and plastic yielding, allowing the solver to converge in a fraction of the iterations (21 Newton / 805 KSP).

**Resolution:**
The `step1_plasticity_cpu_neml2.i` file was updated to use `preset = false` to match the GPU tests, ensuring an apples-to-apples performance comparison.
