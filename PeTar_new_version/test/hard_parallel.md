# Optional HARD_PARALLEL optimizations

This selectively ports [PR #3](https://github.com/joojie1212/PeTar_source/pull/3)
(`f74d968`, based on `5291b7b`). Nothing is enabled
by default. Define `HARD_PARALLEL` consistently for every C++ translation unit
to enable the selected implementation. Do not mix objects built with different
settings: some class layouts change.

## Build

Add `-DHARD_PARALLEL` to the `CXXFLAGS` passed to your existing configure command,
keeping your MPI, SIMD, BSE and other configure options. For example, replace
`CXXFLAGS='-g'` with `CXXFLAGS='-g -DHARD_PARALLEL'`, then clean and rebuild.

To revert, remove `-DHARD_PARALLEL` from the configure flags, reconfigure with
the same feature options, clean and rebuild. The switch is tested with
`#ifdef`: `-DHARD_PARALLEL=0` still enables it. Leave it undefined to disable it.
Changing `OMP_NUM_THREADS` does not disable these algorithms.

## Included behavior

- FDPS merge sort uses one parallel region per merge level; tree walks use
  `schedule(dynamic,1)` instead of chunks of four.
- Hard-cluster profile counters and energy contributions are accumulated in
  per-cluster slots and summed in cluster order. Interrupted integrators are
  registered in cluster order; removing completed entries preserves that order.
- H4 force work and real-particle force corrections use tasks within the existing
  OpenMP team. Neighbor step-update flags are snapshotted before H4 force work;
  both the single and group helpers use those snapshots. The single-neighbor
  snapshot bug and the missed `ARC_n_groups` accumulation in PR #3 are fixed.
- Isolated-particle drift runs in parallel, followed by serial stellar evolution
  in particle order and a search-radius refresh. This serializes the BSE calls
  in this drift routine, not every BSE call throughout the application.
- AR perturbation forces use four-target SIMD when compiling for x86-64 with
  AVX2 and FMA, double precision, and without `INTEGRATED_CUTOFF_FUNCTION`.
  `PETAR_AR_SIMD4` is internal and is derived from these conditions. Other
  configurations retain the scalar AR force path.

Without OpenMP, the new task helpers use a serial fallback and require no
OpenMP runtime. With `HARD_PARALLEL` undefined, the old algorithms, arithmetic,
flag propagation and interrupt-list ordering are retained, including their
existing limitations.

The group compute/commit split and compact-object stellar timestep replacement
from PR #3 are deliberately excluded. Group integration and its immediate
publication of mass/step changes keep the original implementation.

## Validation

Run from the repository root after building `PeTar_new_version/bse-interface/libbse.a`:

```sh
bash PeTar_new_version/test/bugfix/hard_parallel.sh
HARD_PARALLEL_SANITIZE=1 bash PeTar_new_version/test/bugfix/hard_parallel.sh
```

The runner requires GCC, OpenMP, AVX2/FMA hardware and the BSE library. It checks:

- sorting with duplicate keys, empty/small ranges, nonzero offsets and a
  non-power-of-two thread count;
- single-neighbor snapshot behavior when the saved and live flags disagree;
- isolated drift with stellar evolution disabled and enabled, including actual
  mass loss and stellar-type transitions, against the original serial ordering;
- H4 force tasks against direct pair forces, and cluster force corrections;
- scalar-reference AR forces with target counts that exercise partial SIMD
  blocks, plus the existing member-resolved group perturbation tests;
- 1/2/3/4-thread equality of drift/H4/correction results to the macro-off
  single-thread reference, and OpenMP-free linking of the task helpers.

AR SIMD is checked with numerical tolerances, not bitwise equivalence. Explicit
FMA and distance-sum ordering can change rounding. Neighbor flag propagation
with the macro enabled uses the start-of-evaluation snapshot. Neither the short
regressions nor successful compilation establish long-duration energy accuracy,
merger statistics or a speedup for a production workload. Validate those with
representative restarts before adopting the enabled build for production.
