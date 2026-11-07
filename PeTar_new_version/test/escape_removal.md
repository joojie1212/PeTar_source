# Atomic binary escape and synchronized deletion energy

`PeTar::removeParticles()` now gathers radial escape candidates and their named
companions, including across MPI ranks. A pair is protected only if its IDs are
reciprocal and its current relative two-body energy is negative. Selection uses
the pair CM radius and velocity relative to the stellar frame, and its external
potential (mass-weighted member potentials with mutual binding removed).
Both members are retained or removed together, including when one lies inside
the radius. Unbound/dangling pair IDs revert to individual selection; hierarchical
multiples are not automatically treated as one indivisible escape object.
Binaries without reciprocal partner IDs are not inferred by an all-pairs search.

Every candidate's stellar potential is recalculated directly from the current
real-particle state on all ranks. The method exchanges candidates, not the whole
cluster; cost is O(N_candidates*N_local) on steps with radial candidates.
Explicit positive-mass removals are included before zeroing masses; if one member
of a protected pair is forced out, both are removed. Already-booked zero-mass
remnants are simply cleaned up without a second energy adjustment.

For the full deletion set S, energy changes by

```
delta E = -(K_S + U_S,remaining + U_S,S + U_external,S)
        = -sum_S(K_i + m_i phi_i) + U_S,S.
```

The last term removes duplicate counting of mutual potential, including between
distinct escaping binaries and across ranks. All ranks apply the same global
energy-reference correction. K is evaluated in the same inertial frame used for
status energy. New escape output carries the recomputed pot_tot. The separate
pot_soft column remains the previous soft-force cache and is not used here.

The current isolated-star build is tested. Galpy uses a current-time potential
query with zero reaction mass. The optional DM path includes DM for selection
and excludes star-DM cross energy from stellar-reference corrections, consistent
with the separate stellar/DM status convention; combined DM energy auditing and
these optional configurations are not covered by the tests below.

## Tests

```
make build/petar.escape_removal.test
OMP_NUM_THREADS=4 mpiexec -n 1 ./build/petar.escape_removal.test
OMP_NUM_THREADS=4 mpiexec -n 2 ./build/petar.escape_removal.test
OMP_NUM_THREADS=4 mpiexec --oversubscribe -n 4 ./build/petar.escape_removal.test
```

The test invokes the actual `PeTar::removeParticles()` method with intentionally
wrong cached potentials. It checks the exact remaining IDs, global removal
counts, reference energy on every rank and independent before/after pairwise
energy. Cases cover the old half-binary failure, nonescaping bound CM, radius-only
selection, a boosted cluster frame, single escape, binary members on opposite
sides of the radius, unbound or missing partners, explicit removal, softening,
multiple simultaneous binaries, and zero-mass cleanup. MPI ownership interleaves
members and also exercises empty ranks. The normalized energy mismatch must be
less than 2e-12 in each case.

`diagnostics/escape_removal_fix/run_validation.py` additionally restarts the N500
AiywPI snapshots 159, 201 and 495 with both the preceding installed executable and
the rebuilt one, in isolated output directories. All targeted pairs must appear
in the fixed escape log at exactly the same timestamp. Expected maximum physical
energy residual limits, relative to each restart's initial energy, are 5e-4,
3e-5 and 3e-5 respectively. These short ASCII restarts do not reproduce the full
original trajectory or preserve every internal prediction cache. Full N500 runs
have not been repeated by this validation.

Build with `-D ESCAPE_REMOVAL_DEBUG` to emit an `ESCAPE_REMOVAL` line on rank 0 for
each positive-mass deletion batch, containing time, count, energy correction and
IDs with their recognized companions. The regression target enables this flag;
normal builds do not add these verbose records.
