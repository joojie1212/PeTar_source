# Contact orbit dispatch regression

With PeTar configured for BSE, run from the project directory:

```sh
make build/petar.tide_collision.test
OMP_NUM_THREADS=1 build/petar.tide_collision.test
```

The test calls `ARInteraction::modifyAndInterruptIter` with two main-sequence
stars on elliptic and hyperbolic encounters. Predicted pericentre overlap must
not merge stars while their instantaneous positions remain detached. Cases
placed physically inside the summed stellar radius at peri-centre must merge.
It checks merger status and finite output particle states as well as survival.

Nominal tangency can round to either side after converting the orbit to
particles. The tangency case uses the reconstructed pericentre; the cases
1e-12 to either side separately exercise both outcomes.

The tide model is skipped when its predicted pericentre is inside the stellar
surfaces, but that prediction alone no longer causes an irreversible merger.
These are controlled local encounters, not a replay of the original MPI job.

The eccentricity-switch regression verifies that the additional tide is
inactive at `e <= 0.9`, removes energy from a detached bound orbit at
`e > 0.9`, and still captures an initially unbound encounter. It checks
orbital-angular-momentum conservation, energy-loss consistency, and the
energy cap at `e = 0.9`; no further loss is applied after the handoff.

The BSE integration checks compare identical binaries with the per-call
cutoff, unrestricted tides, and `tflag=0`. Above the threshold the cutoff
agrees with tides disabled while the stellar clock advances. At and below
the threshold it agrees with unrestricted BSE. The configured `tflag`
remains unchanged, including when the user explicitly disabled BSE tides.

The production switch is enabled with `--stellar-evolution 2`. It changes
only tidal terms, preserving the other BSE evolution processes. Compact
pairs retain their existing gravitational-wave dispatch. The original
Fortran `evolv2` entry retains its unrestricted default via a wrapper.

For collision-threshold experiments, `--debug_lessmerger` replaces the
two-star surface-overlap criterion with a single pair separation of
`1e-6` solar radii. Specify a different threshold with
`--debug_lessmerger 0.01` (or `--debug_lessmerger=0.01`). Without this
option, PeTar uses the sum of the two stellar radii as before. The override
applies to PeTar's AR instantaneous and delayed merger checks, including
the contact gate for ordinary-star BSE merger predictions. It does not
change stellar radii used by BSE evolution or dynamical-tide calculations.
The collision test includes an overlap that must survive in this debug mode.
