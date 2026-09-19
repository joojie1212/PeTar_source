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

## Building PeTar

From `PeTar_new_version`, configure and build using the options appropriate for
the target machine. The adjacent dependency directories are found by the
project's default configuration paths.
