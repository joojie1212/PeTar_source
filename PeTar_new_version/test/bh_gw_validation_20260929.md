# BH GW validation — 2026-09-29

No production run was restarted and no SDAR source was changed. Tests used
the configured build with assertions enabled and OMP_NUM_THREADS=1.

## Passed

- `petar.bh_peters.test`: circular analytic solution, eccentric invariant,
  numerical convergence, split interval, widening and unbound rejection.
- `petar.bh_gw_dispatch.test`: legacy clock, delayed deadline, dissipation,
  CM preservation, repeated-time idempotence, cancellation after unbinding,
  terminal merger, remnant position and mass-loss/recoil momentum budget.
- `petar.tide_collision.test`: existing stellar tide/contact regressions.
- Real event 3 dump `dump_merger.0.2.2.1790589929`: replay completed with
  exit 0 within the 180-second limit, reaching its stored end time
  1.90735e-6. 598 particles; 144080 AR substeps, 411519 H4 steps.
  Final reported dE=2.27569 versus |Ekin+Epot| about 6.09316e8
  (ratio about 3.7e-9). This is an energy-bookkeeping diagnostic, not
  independent proof of the radiation model's accuracy.
- Real event 5 dump `dump_merger.0.3.4.1790592312`: exit 0 for all four
  AR step scaling parameters below. No merger during this short interval;
  therefore this replay does not validate a real terminal kick.
- Event 3 default replay logs and event 5 default/half/quarter stderr
  scans found no NaN/Inf tokens or wide-merger rejection.

## Step sensitivity: not yet accepted for production

All event 5 runs have the same stored end time, 1.90735e-6.

| AR ds scale | Last GW change callback | Reported dE_change | Reported dE |
|---|---:|---:|---:|
| 1 | 1.73573e-6 | -17106.3 | -2.98023e-8 |
| 0.5 | 1.73573e-6 | -17106.3 | 2.08616e-7 |
| 0.25 | 1.82176e-6 | -17954.4 | -2.68221e-7 |
| 0.125 | 1.86473e-6 | -18378.0 | -3.57628e-7 |

The magnitude of the GW loss changes by approximately 7.43% between scales
1 and 0.125 despite tiny corrected energy residuals. This means the
coupled evolution has not demonstrated sufficient timestep independence.
Logs indicate an unprocessed tail between the last dissipative callback
and the stored integration end. SDAR's callback is guarded by !time_end_flag
in symplectic_integrator.h; this is consistent with the observed endpoint
sensitivity. Whether the next production interval catches up correctly,
including pair regrouping and checkpoint/restart, still needs testing.
Do not infer that the entire tail is permanently lost in a continued run.

This verification did not implement a fix. The required next checks are
PeTar-side endpoint/clock handling without changing SDAR, cross-interval and
restart tests, and actual terminal-merger convergence with perturbations.
Strong unbound plunges remain outside the implemented Peters model.

## Reproduction and artifacts

Run the build/tests in bh_gw.md. Dump replays were isolated in
`/tmp/petar-gw-replay.Hyj5zV`, with PYTHONPATH set to the repository's
bhmerger/precession-master and OMP_NUM_THREADS=1. Commands used
`build/petar.hard.debug -s 1000000 <dump>` and added `-c 0.5`, `-c 0.25`,
or `-c 0.125` for the event 5 comparisons.

Logs: event3_verify.log/.err, event5_verify.log/.err, event5_half.log/.err,
event5_quarter.log/.err, event5_eighth.log/.err. Event 3 verification alone
produced about 2.0 GiB stdout and 60 MiB stderr; these artifacts were retained.
