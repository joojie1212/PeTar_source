# PeTar BH–BH GW integration (experimental)

This implementation changes PeTar/BSE dispatch, not SDAR. Existing unrelated
working-tree edits are retained. It has not been installed into ~/bin or used
to restart the production simulation.

## Prescription

For a bound BH–BH leaf pair, PeTar bypasses BSE's immediate binary merger and
applies orbit-averaged Peters evolution to the current osculating a and e.
Conservative dynamics and external perturbations remain with SDAR. The
dissipative map retains eccentric anomaly and orbital plane, preserving the
original centre-of-mass position and velocity. This is operator splitting,
not phase-resolved PN integration.

Each callback recomputes the remaining inspiral time; it is not a permanent
deadline. Widening can postpone it and an unbound orbit cancels it. New or
legacy pairs initialize the dissipative clock without replaying old losses.
For BH–BH pairs only, negative time_interrupt = -(1 + last_GW_time) marks this
clock; particle layout is unchanged. Old binaries should not be used to
continue snapshots carrying this new marker.

The terminal condition is a(1+e) <= 10 G (m1+m2)/c^2, so the entire ellipse
is compact. This is an explicit numerical/modeling cutoff, not an exact GR
coalescence criterion. The merger is processed at the next available callback,
with step-level timing error recorded in BH_GW_terminal. There is no claim of
exact next-pericentre scheduling or waveform/phase accuracy.

The remnant is placed at the original binary CM. Its velocity is Vcm + Vk;
GW mass loss and kick are accounted for separately. The particle momentum
change is -deltaM * Vcm + Mf * Vk, not zero. BH_GW_budget logs both terms.
The remnant-fit kick replaces, rather than adds to, any BSE kick.

## Scope and limitations

Hyperbolic encounters retain the existing GW encounter prescription, not
the bound Peters equations. A request to merge outside the terminal scale
is rejected with dump_bh_gw_unresolved and an abort; strong unbound plunges
still require a phase-resolved PN treatment. Thus this is not a complete
replacement for all relativistic encounters. The old harmless SDAR diagnostic
NaN has not been repaired by changing SDAR.

Full production accuracy, timestep convergence, performance, and long-term
energy accounting still require validation. Short dump replay is not a full
cluster simulation validation.

## Build and tests

From the configured PeTar_new_version directory:

```sh
make -j2 build/petar.bh_peters.test build/petar.bh_gw_dispatch.test
build/petar.bh_peters.test
OMP_NUM_THREADS=1 PYTHONPATH=bhmerger/precession-master build/petar.bh_gw_dispatch.test
OMP_NUM_THREADS=1 build/petar.tide_collision.test
make -j2 build/petar.hard.debug build/petar.mpi.omp.avx2.bse
```

The numeric test checks the circular solution, eccentric Peters invariant,
convergence, split steps and rejection of unbound input. The dispatch test
checks legacy clock initialization, postponement, dissipation, CM preservation,
same-time idempotence, cancellation on unbinding, terminal merger and recoil
momentum accounting. These tests and both application builds passed locally.

The already-applied incremental patch artifacts are bh_gw_integration.patch,
bh_gw_clock.patch, bh_gw_build.patch and bh_gw_build_order.patch, in that order.
They accompany the new src/bh_peters.hpp, src/bh_gw_evolution.hpp and test
sources; the integration patch alone is not a standalone distribution.
Do not reapply these patches to this working tree.
