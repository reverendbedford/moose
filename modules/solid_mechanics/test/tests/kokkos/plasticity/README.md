# NEML2/Kokkos/PETSc plasticity GPU benchmark

This directory compares the same three-dimensional perfect-plasticity problem
as progressively more work is moved to one GPU:

| Step | NEML2 | Assembly | PETSc |
| --- | --- | --- | --- |
| 1 | CPU | CPU | CPU |
| 2 | GPU | CPU | CPU |
| 3a | CPU | GPU (Kokkos) | CPU |
| 3b | GPU | GPU (Kokkos) | CPU |
| 4 | GPU | GPU (Kokkos) | GPU (AIJKokkos) |

Step 2 intentionally copies the NEML2 stress and tangent back to host memory.
Step 3a evaluates NEML2 on the CPU while Kokkos assembly remains on the GPU, so
the strain input moves device-to-host and the stress and tangent outputs move
host-to-device. Step 3b keeps NEML2 and assembly on the GPU. Both Step 3 cases
leave PETSc on standard CPU vectors and AIJ matrices. Step 4 also moves PETSc
vectors and matrices to the GPU.

The runner passes PETSc options on the MOOSE command line because MOOSE rebuilds
the PETSc options database before each solve.

All steps use one MPI rank. This makes the experiment usable with one GPU but
means that MPI time and CUDA-aware MPI performance are not measured. A later
multi-GPU experiment should compare `-use_gpu_aware_mpi 1` and `0` with at
least two MPI ranks and one rank per GPU.

## Run

First activate the CUDA-enabled MOOSE/NEML2 environment. For the stack in this
branch:

```bash
source kokkos-cuda-stack/scripts/activate.sh
```

Then run from this directory. The wrapper scripts preserve separate result
directories and encode the recommended progression:

```bash
./run_n32_smoke.sh
./run_n32_benchmarks.sh
./run_n64_smoke.sh
./run_n64_benchmarks.sh
```

The N=32 smoke run uses one unprofiled repetition. The formal N=32 run uses
five timing repetitions and one separate Nsight Systems profile per step. The
N=64 smoke run checks memory, runtime, and convergence with one repetition.
Only after it succeeds, the formal N=64 run performs three timing repetitions
without profiling.

## Coarse Exodus correctness check

Use the coarse-output runner before interpreting benchmark speedups:

```bash
./run_coarse_exodus.sh
```

It runs all five configurations once with `N=8`, forces Exodus and CSV output,
and writes separate `.e`, `.csv`, and `.log` files under
`coarse_exodus_n8/`. This run includes output overhead and is for correctness,
not timing.

For Steps 3a, 3b, and 4, the script also prints the final `ux_right` and
`ux_center` point values already defined in those inputs. At `t=0.005`, the
expected displacement is approximately `ux_right=0.005` and
`ux_center=0.0025`. Open the Exodus files in ParaView and compare `disp_x`
between all five steps.

Override the coarse mesh or output location when needed:

```bash
MESH_N=4 OUTPUT_DIR=/tmp/plasticity_exodus ./run_coarse_exodus.sh
```

Useful benchmark overrides:

```bash
EXE=/path/to/solid_mechanics-opt MESH_N=16 REPEATS=3 PROFILE=0 ./run_benchmarks.sh
```

`PROFILE=0` skips Nsight Systems. The normal repetitions measure wall time
without profiler overhead. When profiling is enabled, the runner separately
creates one `.nsys-rep` for each step.

The runner writes raw logs, PerfGraph JSON, Nsight reports, `comparison.csv`,
and `comparison.png` under `results/`. Run the analysis again with:

```bash
python3 analyze_results.py results
```

The analysis reports:

- total wall time;
- `NEML2::solve` self time from the MOOSE PerfGraph JSON;
- total CUDA kernel time from Nsight Systems;
- MPI time when available (not meaningful for the current one-rank runs);
- H2D and D2H copy time;
- total Newton and KSP iterations, used to check that timing differences are
  not caused by different convergence histories.

Start with `MESH_N=16` for a smoke run, then use 32 or 64 for timing. Compare
iteration totals and solution correctness before interpreting speedups.
