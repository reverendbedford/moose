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
