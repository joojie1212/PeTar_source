# Dynamics-only merger-prediction changes

With BSE_BASE and BHMERGER enabled, `<snapshot-prefix>.change` records significant
changes measured on the gravitationally integrated orbit, before BSE modifies it.
Pure BSE evolution (including GW radiation), rejected BSE terminal-event proposals,
and hyperbolic GW capture do not themselves emit change rows. These are dynamical
orbit diagnostics, not proof of a particular third body's causal contribution:
integration error and osculating orbital variations are not force-decomposed.

## Clock and diagnostic state

`time_merger` is an absolute deadline. For an unchanged compact pair undergoing
pure BSE GW evolution it is preserved, even if already due. Internal updates
synchronize the separate orbital checkpoint on every accepted callback, without
an absolute eccentricity cutoff. Pair formation/exchange, mass/type changes and
GW capture establish fresh predictions silently. Disruption clears the state.
Rejected BH Roche/coalescence proposals retain the checked pre-BSE orbit/clock.

At each dynamical check, accumulate

```
dynamics_dlog_tgw += log[a_new^4 F(e_new)] - log[a_checkpoint^4 F(e_checkpoint)]
```

using the same full Peters eccentricity factor as the predictor. This signed net
accumulator excludes the BSE jump and survives intervening BSE callbacks. A net
ratio below 0.9 or above 1.1 triggers a change row and resets the accumulator.
Dynamics can advance OR delay the absolute deadline; it is not clamped with min.
Subthreshold increments of opposite signs cancel. The log evaluation avoids the
old absolute `abs(delta_e)>1e-3` problem at high eccentricity.

## Output schema (19 columns)

```
time_old time_new id1 id2 m1 m2 a_old e_old a_new e_new
tmerge_old tmerge_new rx ry rz vx vy vz dynamics_dlog_tgw
```

The first 18 column positions are retained; the final column is new. Use a new
output prefix/directory when switching from the old 18-column schema; do not
append new records to a historical change file. MPI rows may be time unordered.
Rank 0 initializes the header, and the writer uses OpenMP and POSIX file locks.

All dimensional quantities use internal simulation units. Relative vectors are
current `particle1 - particle2`. Old/new orbit and epoch now describe the **last
dynamical checkpoint and current integrated state**, with no BSE update between
them. The final column includes earlier subthreshold dynamical increments, so
it need not equal the log ratio of the two orbital endpoints in that row.
For a dynamical disruption the final column is NaN, `tmerge_new=-1`, and the
prediction is invalidated. Initial/exchanged-pair predictions are not changes.

`tmerge_*` remain absolute clock values, not isolated lifetimes of the adjacent
orbital columns. Since BSE preserves the clock while moving the checkpoint,
`tmerge_old-time_old` MUST NOT be used to calibrate the Peters coefficient from
`a_old,e_old` in new logs. Compute isolated lifetimes/contours using the run's
G, c, masses and full Peters factor. To measure an advance of the actual clock,
use `(tmerge_old-tmerge_new)/(tmerge_old-time_new)` only with a positive denominator.
A dynamical trigger and the direction of the absolute clock update are distinct
quantities if BSE's evolved orbit and the fixed Peters clock have drifted apart.

The new accumulator is copied with ParticleBase and reset with the checkpoint.
It changes raw binary particle/hard-dump layout: old binary snapshots and dumps
cannot be used with this build. ASCII columns are unchanged; ASCII restarts
establish a fresh clock/checkpoint (they do not preserve the old countdown).

## Regression

```
make build/petar.orbit_change.test
./build/petar.orbit_change.test
```

Run from the PeTar directory so embedded Python can find the merger module.
Checks cover output fields and particle-state copying, initial/exchange suppression, disruption,
2,000 internal updates without a clock reset or false row, high-e sensitivity,
real dynamical advances/delays, and weak changes accumulated across internal
updates. Actual Fortran BSE plus AR callback tests use e=0, 0.9, 0.9999 with 10
and 100 initial callbacks; they check zero GW-only rows, one row for an applied
velocity impulse, and merger at the fixed deadline with repeated BSE callbacks.
These dispatch tests do not substitute for a full star-cluster rerun or a
resolved three-body encounter convergence study.
