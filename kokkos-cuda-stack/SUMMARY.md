# From-scratch MOOSE-CUDA — one-page summary

**Session:** 2026-08-26 00:14 – 04:35 (~4h 20m), fully autonomous.
**Directive:** build the whole PETSc → libmesh → MOOSE stack from source with no conda, so the `benchmark_cpu_and_gpu.py` GPU leg would actually work.
**Result:** it works.

## Deliverables

| Artefact | Where |
|---|---|
| Working CUDA-enabled `solid_mechanics-opt` | `moose-kokkos/modules/solid_mechanics/solid_mechanics-opt` — `--show-capabilities` reports `kokkos = 4.7.4`, `cuda = 12.4.0` |
| From-scratch stack | `moose-kokkos/kokkos-cuda-stack/prefix/` (~1.1 GB): PETSc 3.25.4 + libmesh + WASP; scripts to reproduce in `moose-kokkos/kokkos-cuda-stack/scripts/` |
| Full CPU + GPU benchmark data | `moose-kokkos/modules/solid_mechanics/test/tests/kokkos/linear_elasticity/cpu_and_gpu.{csv,png}` (mesh 256, 5 reps) + `moose-kokkos/kokkos-cuda-stack/logs/cpu_and_gpu_512.{csv,png}` (mesh 512, 3 reps) |
| Attempt-by-attempt journal | `moose-kokkos/kokkos-cuda-stack/ATTEMPTS.md` — 13 attempts with root cause + fix for each |
| Layout + rationale | `moose-kokkos/kokkos-cuda-stack/README.md` |
| Memory pointers for future sessions | `~/.claude/projects/.../memory/project_moose_cuda_from_scratch.md`, updated `project_petsc_kokkos_cpu_only.md` |

## Numbers at mesh 256 × 256

| Config | Wall time (s), mean ± s.d. |
|---|---|
| CPU (2 MPI × 1 thread) | 6.16 ± 0.29 |
| CPU (2 MPI × 4 threads) | 6.24 ± 0.30 |
| CPU (2 MPI × 16 threads) | 6.27 ± 0.34 |
| **GPU (1 MPI × 1 thread)** | **7.85 ± 0.04** |

GPU is ~1.25× SLOWER at this size. Two reasons:
1. **KSP solve stays on CPU** — the input doesn't set `-mat_type aijkokkos -vec_type kokkos`; only assembly ran on GPU (nsys confirmed `Moose::Kokkos::Dispatcher<…>` kernels). And trying those flags aborts with errorcode 76 because system OpenMPI is not CUDA-aware (`ompi_info` → `opal_built_with_cuda_support = false`).
2. **Mesh is too small** — 65 k elements can't amortise per-Newton kernel-launch and host↔device sync overhead against the CPU's ample DDR bandwidth. Assembly kernels are microseconds each.

Neither is a build defect. See ATTEMPTS.md attempts 11–13.

## Key gotchas from the 13 attempts (short list)

1. `--download-openblas` fails to link its .so because gcc-11 doesn't know about `/usr/lib/gcc/x86_64-linux-gnu/9/libgfortran.so`. **Use system BLAS/LAPACK** (`--with-blas-lib=/usr/lib/x86_64-linux-gnu/libblas.so`, `--with-lapack-lib=…liblapack.so`) and export `LDFLAGS=-L/usr/lib/gcc/x86_64-linux-gnu/9`.
2. Don't override `LIBMESH_BUILD_DIR`. MOOSE's two libmesh scripts disagree on where the build dir lives (`update_and_rebuild_libmesh.sh` respects the env var; `configure_libmesh.sh` hardcodes `${SRC_DIR}/build`).
3. Conda leaves build hints (`build_alias=x86_64-conda-linux-gnu`, `CMAKE_ARGS`, `CPP=x86_64-conda-linux-gnu-cpp`, `CXX_FOR_BUILD`, …) beyond `CONDA_*` that autoconf reads. **Whitelist-based env purge**, not unset-by-name.
4. Add WASP to the from-scratch stack (`scripts/build_wasp.sh`, `WASP_DIR=$PREFIX`). MOOSE's `moose.mk` requires it; conda supplied it via `moose-wasp` in the old env.
5. `nvidia-smi` polls too slowly (~1 s) to detect a MOOSE run whose total GPU kernel time is microseconds. **Use `nsys profile --stats=true` and grep for `Moose::Kokkos::`** as the definitive GPU-usage check.
6. System OpenMPI on this box is not CUDA-aware. GPU-KSP requires CUDA-aware MPI; not built in this session (would require rebuilding OpenMPI + PETSc + libmesh + MOOSE, ~3-4 h).

## Rerun

```bash
/home/chenghau.yang/packages/moose-kokkos/kokkos-cuda-stack/scripts/all.sh   # full chain, ~2 h
```

or later, in any shell:

```bash
. /home/chenghau.yang/packages/moose-kokkos/kokkos-cuda-stack/scripts/env.sh
export PETSC_DIR=$PREFIX PETSC_ARCH="" LIBMESH_DIR=$PREFIX WASP_DIR=$PREFIX
$MOOSE_DIR/modules/solid_mechanics/solid_mechanics-opt -i <input.i> --compute-device=cuda
```

To go back to the conda `moose` env: start a fresh shell (don't source `env.sh`), `conda activate moose`, then `cd framework && make clean && make -j` to rebuild against conda's PETSc/libmesh. The `moose` env itself was never touched during this session.

## Not done — sensible next steps (in priority order)

1. **Reproducibility check.** Run `kokkos-cuda-stack/scripts/all.sh` from scratch on a machine where the current prefix is deleted, to prove the scripts work end-to-end without human intervention. I did not do this on the current box because the artefacts are already there and the user wanted GPU numbers over 8h.
2. **CUDA-aware MPI + GPU KSP.** Rebuild OpenMPI with `--with-cuda=/usr/local/cuda` (or drop in a CUDA-aware MPICH), then rerun the from-scratch build. Then move the KSP solve to GPU (`-mat_type aijkokkos -vec_type kokkos …`). Only then can you get GPU winning wall time on this input.
3. **Bigger mesh sweep.** Add a `benchmark_mesh_scaling.py` (parametric over mesh_size) that plots wall-time vs problem size for CPU and GPU, so the crossover point is visible.
4. **Fix TestHarness deps.** `pip install pandas packaging pyyaml jinja2 numpy` (under system python) makes `./run_tests` work with the from-scratch env. Not needed for direct executable runs.
