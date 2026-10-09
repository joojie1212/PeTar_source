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
- each accepted event records both progenitors, their pre-merger positions,
  velocities and spins, the binary orbit, and the final remnant in a shared
  `*.bhmerger` output protected by process and thread locks.

This is an operator-split, orbit-averaged prescription, not phase-resolved PN
integration or a waveform model. The callback timing, timestep convergence,
and long-duration production accuracy remain modeling limitations documented
under `PeTar_new_version/test/`.

## Local stability fixes

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
