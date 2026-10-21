# PeTar source bundle

This repository is a self-contained source snapshot of PeTar and the dependency
sources used with it. It includes the current local source changes present when
the snapshot was created, while excluding nested Git metadata, compiler output,
runtime results, caches, and machine-local editor settings.

## Layout

- `PeTar_new_version/` — PeTar
- `FDPS/` — Framework for Developing Particle Simulators
- `SDAR/` — slowdown algorithmic regularization library
- `mcluster/` — initial-condition generator

PeTar's default configure paths expect `FDPS` and `SDAR` beside the PeTar
directory, matching this layout.

The repository intentionally does not contain any directory named `sample` or
`tidletest`. Compiler products, simulation output, nested Git metadata, the
machine-local CUDA toolkit, and cached third-party binary packages are also
excluded. CUDA, MPI, Python/NumPy, Galpy, and system compiler runtimes must be
provided by the target machine.


## Added features in this snapshot

### Live collisionless dark matter

Configure with `--enable-darkmatter` to add an independently numbered live
dark-matter component. A run accepts the normal stellar snapshot plus a second
snapshot through `--dm-input`; synchronized outputs use the configurable
`--dm-output-prefix` (default `dmdata`). Stellar and dark-matter IDs occupy
separate namespaces, and paired restart snapshots must have the same snapshot
ID and time.

Dark matter uses a component-wide Plummer softening length and a tree-step
leapfrog. It participates in domain decomposition and MPI particle exchange,
but never enters neighbour/Hermite/SDAR, collision, or SSE/BSE paths. Two FDPS
distributed monopole trees evaluate DM self-gravity, DM-to-star forces, and
star-to-DM forces without a per-step all-gather. x86 builds can use the existing
Phantom-GRAPE SIMD kernels; `--enable-simd-64` is recommended when accurate
self-potential subtraction matters. CUDA kernels are not used for the DM force.

The implementation includes:

- independent stellar and DM snapshots with ASCII/binary restart support;
- complete star-star, DM-DM, and star-DM potential-energy accounting;
- consistent Galpy forces, moving-potential reactions, reference-frame shifts,
  and restart state shared by both particle components;
- MPI/OpenMP distributed-tree regression tests against direct summation and an
  independent Galpy trajectory integration.

### Black-hole inspiral, merger, and remnant placement

Configure with `--with-interrupt=bse --enable-bhmerger`. The build embeds
Python/NumPy and the bundled `bhmerger/precession-master` remnant-fit code.
Python is initialized once around the PeTar integration loop, rather than once
per merger, and merger calls and output are protected for OpenMP/MPI use.

For a bound BH-BH leaf pair, BSE is the sole owner of gravitational-wave
dissipation and evolves the semi-major axis and eccentricity. PeTar does not
apply a second Peters evolution. After a BSE update, PeTar uses the current
orbit in the constant-orbit Peters approximation to predict a merger deadline;
the prediction is refreshed after later perturbations or BSE updates.
Hyperbolic encounters retain the separate encounter prescription. Only an
actually due merger serializes remnant fitting and particle replacement.

The BH merger position/velocity handling differs materially from the upstream
path:

- the remnant is placed at the pre-merger binary centre-of-mass position;
- its velocity is the original centre-of-mass velocity plus the fitted recoil
  kick, and the fitted kick replaces rather than adds to a BSE kick;
- GW mass loss and kick momentum are accounted for separately as
  `-deltaM * Vcm + Mf * Vk`;
- remnant mass, Kerr-spin magnitude/direction, and recoil are validated before
  being applied, with the spin retained for hierarchical mergers;
- `petar.init --bhmerger` initializes the three Kerr-spin components, while
  ASCII snapshots write and restore them on restart (`--blackhole` enables the
  BHMERGER input format automatically);
- each accepted event records both progenitors, their pre-merger positions,
  velocities and spins, the binary orbit, and the final remnant in a shared
  `*.bhmerger` output protected by process and thread locks.

This is an operator-split, orbit-averaged prescription, not phase-resolved PN
integration or a waveform model. The callback timing, timestep convergence,
and long-duration production accuracy remain modeling limitations documented
under `PeTar_new_version/test/`.
10.12 修复了一个bug这个bug会导致并和时间判定的黑洞并和被无限延后


## Local stability fixes

### SDAR nearest-neighbor initialization

`SDAR/src/Hermite/neighbor.h` now initializes and resets `r_min_mass` in
the constructor and all neighbor reset/clear paths. Previously, the first
mass-weighted nearest-neighbor comparison could read this value before it
was assigned. Both single-particle and group searches now accept the first
positive-mass candidate explicitly, then retain the existing distance-squared
over mass comparison. Zero-mass candidates are excluded from nearest-neighbor
selection; neighbor-list membership is unchanged. This fixes an uninitialized
read, but has not been established as the cause of the `data.46` reproducibility
issue. Rebuild PeTar to include this header change.

The focused regression covers construction, every reset/clear path, both
candidate types, zero-mass candidates, and mass-weighted selection:

```sh
g++ -std=c++17 -O2 -Wall -Wextra -fopenmp -I SDAR/src SDAR/tests/neighbor_initialization.cpp -o /tmp/sdar-neighbor-test
/tmp/sdar-neighbor-test
```

### SDAR group-break mask bounds

`HermiteIntegrator::checkBreak()` now allocates a zero-initialized merger mask
over `groups.getSize()` slots instead of the active-group count. Active group
indices can have holes after groups are removed; for example, active slots
`[0, 2]` previously indexed past a two-element mask. The change fixes both the
merger-mask write and the break-check read, adds slot-bound assertions, and
replaces the variable-length stack array with `std::vector<unsigned char>`.
The group-break criteria are unchanged. Rebuild PeTar to apply the fix;
its effect on `data.46` reproducibility has not yet been validated.

The regression directly calls `checkBreak()` with active slots `[0, 2]`,
both with and without a merger at slot 2, and checks the resulting break list:

```sh
g++ -std=c++17 -O1 -g -fopenmp -fsanitize=address,undefined -I SDAR/src SDAR/tests/group_break_mask.cpp -o /tmp/sdar-group-break-test
/tmp/sdar-group-break-test
```

### SDAR initial group-member neighbor state

`initialSystemSingle()` now clears the neighbor record of every particle and
sets `r_neighbor_crit_sq` from its current `getRNeighbor()` before groups are
formed. Previously, particles placed into initial groups were removed from the
single-particle initialization list, leaving their member-level search radii at
the default value even though resolved-group force evaluation uses these records.
The later local initialization in `initialIntegration()` is unchanged, preserving
unrelated particles' neighbor lists and step-update flags during group changes.
This corrects neighbor initialization, not the force law; its effect on full-run
trajectories and `data.46` reproducibility has not yet been measured. Rebuild
PeTar to apply it. The `group_break_mask.cpp` regression above also checks this
initialization path and neighbor selection for an initial group member.

### Empty-list temporary arrays

Short-restart UBSan checks exposed zero-length variable-length arrays in
SDAR's `integrateToTimeList()` and `integrateGroupsOneStep()`, and PeTar's AR
perturber calculation. Empty integration requests now return immediately;
the selection, interruption-index, and single/group perturber-pointer arrays
use `std::vector`, which supports empty lists. This includes the case where
perturbers exist but all are singles or all are groups. Existing iteration
counts, force formulas, and integration criteria are unchanged. The focused
regression also covers an empty integration request and a request with no
group slots. Rebuild PeTar to apply these header changes.

Validation: the focused ASan/UBSan regression passed. A freshly compiled
`-O0` ASan/UBSan OpenMP/BSE build (BHMERGER disabled to match the 34-column
snapshot) completed a one-thread, 16-step restart from `data.46`, from
`t=4.598876953125` to `4.599365234375`, retaining 982 particles, with
`UBSAN_OPTIONS=halt_on_error=1` and leak detection disabled. No ASan/UBSan
diagnostics were reported in that run. This is not a multithreaded or long-term
accuracy validation; the existing BSE library was not sanitizer-instrumented.

### PeTar state consistency

The current PeTar source includes fixes for three state-consistency and
concurrency bugs that could propagate invalid particle data or crash a run:

- Hyperbolic gravitational-wave evolution now computes the new eccentricity
  from the newly updated semi-major axis instead of mixing the old and new
  orbital states.
- BSE orbital updates now rebuild particles from the new masses, semi-major
  axis, and eccentricity, then recalculate all cached binary quantities from
  that rebuilt state. Invalid or non-finite post-BSE orbits are detected and
  dumped immediately instead of being propagated into later integration.
- `Dynamic_merge` diagnostics and event records are serialized in an OpenMP
  critical section, preventing concurrent threads from corrupting the shared
  output stream.
- For ordinary stars with BSE types `kw=1-13`, predicted `Contact` or
  `Coalescence` events are deferred until the integrated instantaneous
  separation satisfies `r <= R1 + R2`. Detached stars remain in the
  integration and can respond to perturbations that prevent the predicted
  future collision, recovering objects that were previously removed by
  false-positive merger predictions based only on the osculating orbit or
  pericentre. The same physical-contact requirement is used for bound and
  unbound encounters, and stale delayed-collision flags are cleared.
- 10.18 修改了escaper.hpp中的能量判据，目前会给动能乘上1/2，

### Merger radius selection

By default, PeTar uses the physical surface-overlap distance `R1 + R2` for
instantaneous and delayed AR merger checks. For controlled collision-threshold
experiments, `--debug_lessmerger` replaces this pair distance with `1e-6` solar
radii. A positive custom distance can be supplied in solar radii, for example
`--debug_lessmerger 0.01` or `--debug_lessmerger=0.01`. Omitting the option
restores the physical `R1 + R2` criterion. The override also applies to the
contact gate for ordinary-star BSE `Contact` and `Coalescence` predictions, but
does not modify the stellar radii passed to BSE or used by dynamical tides.

### Single capture update and BSE tides (zhujie)

The additional dynamical-tide energy-loss prescription now acts only on
unbound encounters (`a < 0`). After a successful capture, further updates from
that prescription stop, including any remaining slowdown repetitions. Bound
binaries use BSE tidal evolution when `bse-tflag > 0`. This replaces repeated
extra tidal updates of captured binaries with the capture update followed by
BSE evolution, limiting overlap between the two prescriptions. The regression
checks bound-orbit exclusion, successful capture and post-capture exclusion.

### Optional frozen-binary optimization

Build with `--enable-frozen-binary` to enable the experimental, reversible
frozen state; it is disabled by default. Eligible hard, detached, weakly
perturbed binaries with `ecc < min(frozen-ecc-limit, 0.1)` can reduce frequent
BSE checks. Isolated frozen pairs advance their internal phase analytically;
frozen inner pairs in hierarchical groups use perturbation-safe SDAR slowdown.
Scheduled stellar updates, unsafe encounters and membership changes restore
normal integration. See [the PeTar documentation](PeTar_new_version/README.md#optional-frozen-binary-state)
for parameters and safety checks.

**Current local tests have not shown a significant overall performance
improvement from frozen mode.** No general speedup is claimed; the feature
remains optional and requires workload-specific validation.

## Opt-in hard-cluster OpenMP parallelism

The hard-cluster force loops can use the idle threads of the outer
hard-cluster team without opening a second parallel level for the FDPS tree.
The feature is opt-in: with no environment variables set the binary keeps the
serial hard integrator path (the outer cluster loop is still the same OpenMP
loop as before).

| variable | default | meaning |
|---|---|---|
| `PETAR_INNER_THREADS` | 1 | threads for the level-2 loops; 1 disables the feature |
| `PETAR_OUTER_THREADS` | `OMP_NUM_THREADS` | team size of the per-cluster loop |
| `PETAR_PAR_MODE=task` | team | reuse the current team with `taskloop` instead of nesting |
| `PETAR_PAR_H4` / `PETAR_PAR_AR` | 1 | enable the H4 force loop / AR group loop |
| `PETAR_AR_LPT` | 1 | largest-group-first scheduling for the AR loop |
| `PETAR_CLUSTER_LPT` | 0 | largest-cluster-first dispatch |
| `PETAR_PARSTAT` / `PETAR_PHASE_TIMER` | 1 / 0 | diagnostics at exit |

A typical task-mode run on 32 cores:

```sh
OMP_NUM_THREADS=32 PETAR_PAR_MODE=task \
PETAR_OUTER_THREADS=32 PETAR_INNER_THREADS=32 \
PETAR_PAR_H4=1 PETAR_PAR_AR=1 ./petar ...
```

The patch keeps every floating-point reduction in its original order: H4
force items are independent, AR groups are computed on private copies and
committed in ascending group order, and energy/profile counters are accumulated
per cluster and reduced serially. On the `data.46 + 16 steps` window,
1/2/4/8/16/32-thread runs are bit-identical in every snapshot and in the
`Physic:` energy lines. BSE/SSE Fortran calls are serialized with
`critical(petar_bse_fortran)` because their RNG and scratch state live in
COMMON blocks; long integrations with active BSE events still cannot be
promised bit-reproducible across thread counts.

The Hermite force loop now uses OpenMP, so the SDAR regression binaries above
must be linked with `-fopenmp`.

## Upstream starting points

The working-tree snapshot was based on these commits before including local
changes:

- PeTar: `e43cbdbcf3036be8597faab7f22e28e55c72bc9c`
- FDPS: `bce909e2bc7d89fe92853ba70758f84aebd14c8a`
- SDAR: `ac3e5d8aa2c630e69467b4f119b1442a04a48b22`
- mcluster: `8a4dacfe13b1ad88295737e269e47baff7cd3789`

The component directories retain their original license files and attribution.

## Difference from the GitHub PeTar source

This comparison was refreshed on 2026-10-06 against
`git@github.com:joojie1212/PeTar.git` branch `master`, commit
`45d1dce6d9bc6ed0e135d5ab7d0be385d9f14278` (`added merger part`). The local
PeTar working tree started from `b2ce49f8308348dc44a6603d92caa0e0b3c9d846`
and includes later uncommitted development. It is therefore a maintained
research snapshot, not a byte-for-byte mirror or a release tag.

Relative to that GitHub commit, this bundle adds or changes:

- the live, distributed dark-matter component (`darkmatter.hpp`,
  `darkmatter_force.hpp`, configure/build plumbing, I/O, restart, energy,
  Galpy, and regression tests);
- BSE-owned BH-BH orbital dissipation with approximate merger-time prediction,
  CM-preserving remnant placement, recoil/mass-loss momentum accounting,
  persistent remnant spin, and structured merger output;
- safer embedded-Python lifecycle and stronger validation of remnant-fit return
  values, including a direct/hyperbolic merger path;
- physical instantaneous-contact gating for ordinary-star BSE mergers and the
  configurable `--debug_lessmerger[=Rsun]` collision radius;
- consistent rebuilding of post-BSE orbital state, hyperbolic-GW eccentricity
  fixes, and serialized merger diagnostics;
- the optional reversible frozen-binary optimization and its safety checks;
- additional focused tests and validation notes under
  `PeTar_new_version/test/`.

The snapshot also omits generated Doxygen HTML, all `sample` and `tidletest`
trees, build products, and run data. FDPS, SDAR, and mcluster are committed as
ordinary source directories rather than Git submodules so a clone contains the
required dependency source immediately.

## Building PeTar

### lingsheng SIMD dependency (Arm/NEON)

The Arm SIMD compatibility backend requires an AArch64 CPU/compiler with
NEON and the external, header-only **SIMDe** library. The existing backend
notes use SIMDe v0.8.2; its headers are not bundled in this repository.
Point the include path at the directory containing `simde/`:

```sh
cd PeTar_new_version
./configure --with-simd=neon \
  CXXFLAGS="-O3 -march=armv8-a+simd -I/path/to/simde" \
  --with-fdps-prefix=../FDPS --with-sdar-prefix=../SDAR
make
```

Retain the MPI/BSE and installation options appropriate to the application.
Omit `--with-arch` for this backend (do not use `--with-arch=none`). It maps
the existing AVX-width kernel through SIMDe to 128-bit NEON operations;
SVE is not required. Native x86 AVX builds do not require SIMDe. This is
separate from McLuster's compiler/OpenMP SIMD optimization below.
See [Arm SIMD setup and numerical checks](PeTar_new_version/ARM_SIMD.md).
The saved backend validation covers portable and native x86 execution;
performance and numerical behavior still need validation on the target Arm CPU.

### McLuster pair-sum correction

The two original CPU potential loops (Plummer and the branch for profiles
`1`, `3`, and `-1`) and the accelerated Plummer loop now use `j < i`.
The old `j < i-1` range omitted `(0,1), (1,2), ..., (N-2,N-1)`: N-1
interactions. Each distinct pair must contribute once to
`U = -G sum(i>j) m_i m_j / r_ij`. The omission made the binding energy
too small in magnitude and affected subsequent model normalization.

This fix applies to the default CPU build as well as `mcluster_fast`;
the GPU implementation is unchanged. Existing initial-condition files are
not modified. Regenerate them when comparing corrected models; earlier
benchmark models/results predate this fix. A two-particle analytic check
guards against the old zero-potential result.

### McLuster acceleration algorithm

```sh
make -C mcluster mcluster_fast
OMP_NUM_THREADS=32 ./mcluster/mcluster_fast [normal McLuster options]
```

The optional CPU Plummer implementation copies masses and positions into
contiguous arrays, processes 256-row blocks for cache reuse, and uses OpenMP
parallelism and SIMD reductions. It still computes direct O(N²) gravity;
there is no tree or analytical-potential approximation. Extra coordinate
storage is 32 bytes per particle (96 MB for three million particles).
Random sampling and model parameters are unchanged, but floating-point
reduction order differs, so results need not be bitwise identical.

The build uses `-O3 -march=native -fno-math-errno` with OpenMP, without
`-ffast-math`. Rebuild on a different CPU rather than copying this binary.
Progress reports completed row blocks, not a linear estimate of elapsed
time. Before the correctness fix, a local 500,000-particle kernel test
measured about 1.48x speedup; this is not a universal performance guarantee.
See [algorithm details and numerical test command](mcluster/FAST_POTENTIAL.md).

McLuster's existing output routines have 50-byte filename buffers. Until
those routines are hardened, run inside the destination directory and pass
a short `-o` basename rather than a long absolute path.

### Main simulation build

From `PeTar_new_version`, configure and build using the options appropriate for
the target machine. The adjacent dependency directories are found by the
project's default configuration paths.

For a CPU build with BSE, BH mergers, and live dark matter, a representative
configuration is:

```sh
cd PeTar_new_version
./configure --prefix="$HOME" --with-interrupt=bse \
  --enable-bhmerger --enable-darkmatter --with-mpi=yes
make -j
make install
```

Python development headers and NumPy are required for `--enable-bhmerger`;
MPI and a Fortran compiler are required by this example and BSE respectively.
