# NEML2/Kokkos/PETSc plasticity GPU benchmark

This directory compares the same three-dimensional perfect-plasticity problem
as progressively more work is moved to one GPU:

| Step | NEML2 | Assembly | PETSc |
| --- | --- | --- | --- |
| 1 | CPU | CPU | CPU |
| 2 | GPU | CPU | CPU |
| 3 | GPU | GPU (Kokkos) | CPU |
| 4 | GPU | GPU (Kokkos) | GPU (AIJKokkos) |

Step 2 intentionally copies the NEML2 stress and tangent back to host memory.
Step 3 leaves PETSc on standard CPU vectors and AIJ matrices. Step 4 selects
Kokkos vectors and AIJKokkos through `PETSC_OPTIONS` in the runner because
these unprefixed PETSc options cannot reliably be placed in a MOOSE input.

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

Then run from this directory:

```bash
MESH_N=32 REPEATS=5 ./run_benchmarks.sh
```

Useful overrides:

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
