# Significant merger-prediction orbit changes

With BSE_BASE and BHMERGER enabled, the normal output initialization creates
`<snapshot-prefix>.change` (normally `data.change`). It follows the existing
output/append settings. Rank 0 creates the header before the MPI barrier;
rows use the merger writer's OpenMP critical section and POSIX file lock.
Rows from different ranks are not guaranteed to be time ordered.

A row is written when the significant-orbit-change predicate causes an actual
prediction refresh for the same two members with a valid cached prediction.
Initial predictions, companion exchanges, and refreshes caused only by BSE
updates do not write rows. The existing prediction criterion is unchanged.

Columns, in order:

```
time_old time_new id1 id2 m1 m2 a_old e_old a_new e_new
tmerge_old tmerge_new rx ry rz vx vy vz
```

All dimensional quantities use internal simulation units, as in `.bhmerger`;
eccentricities are dimensionless. `m1/m2` are current masses. Relative vectors
are current `particle1 - particle2`. The new orbit is the one recalculated
for the prediction. Old orbit and time refer to the previous prediction
refresh, not necessarily the instant immediately before a physical encounter.
`tmerge_*` are absolute predicted merger times, not remaining lifetimes.
If the new orbit is invalid/unbound, `tmerge_new=-1` and the cache is cleared.
These rows describe orbit changes, not an attribution to external forces.

The prediction epoch travels with the particle through DataCopy and binary
snapshots. Adding this field changes the raw ParticleBase binary layout:
binary snapshots/hard dumps from older builds are not compatible. ASCII
snapshot columns are unchanged; ASCII restarts start fresh prediction caches,
so no change row bridges the restart until a new baseline has been established.

Build and run with a BHMERGER-enabled configuration:

```
make build/petar.orbit_change.test
./build/petar.orbit_change.test
```

The regression checks initial/exchanged-pair suppression, cached epoch copying,
old/new values and relative vectors, append behavior, and invalidation of an
unbound orbit. It does not run a full cluster simulation or a multi-rank run.
