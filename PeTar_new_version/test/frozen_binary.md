# Frozen-binary regression

Configure PeTar with `--enable-frozen-binary` to activate the optional
`FROZEN_BINARY` code path. The feature is disabled in normal builds.

The standalone policy regression can be built regardless of the current
configuration:

```sh
make build/petar.frozen_binary.test
build/petar.frozen_binary.test
```

The test checks the hard-binary energy threshold, external perturbation,
periapsis clearance, BSE event state, scheduling interval, pair-ID storage,
reversible frozen flag, analytic phase advance, and conservation of two-body
orbital energy. It also checks the circular detached override for a zero BSE
interval and the imposed finite recheck time. An eligible isolated frozen pair bypasses SDAR internal
substeps and follows the analytic Kepler solution while its centre of mass
continues in the normal hard/soft integration. A frozen inner binary in a
hierarchical group uses SDAR's maximum perturbation-safe slowdown while the
outer hierarchy remains integrated. A failed safety condition restores the
ordinary integration path.

Tight low-eccentricity pairs (`ecc < min(frozen-ecc-limit, 0.1)`, default **0.1**) with periapsis strictly
less than `10 * (R1 + R2)` bypass the configurable `frozen-radius-factor`
clearance. Both instantaneous separation and periapsis must still exceed
`R1 + R2`; stellar-event, binding-energy and perturbation checks remain active.
Other pairs retain the ordinary clearance policy. The regression covers this
override, its strict upper boundary, contact rejection and safety checks.

All frozen candidates now require strictly `ecc < 0.1`, including ordinary
clearance candidates. Older parameter files cannot relax this hard upper bound.
Existing numeric event columns are preserved; reason fields are appended.
Start reasons are `tight_low_ecc` or `ordinary_clearance`, followed by `ecc`
and `peri_over_rsum`. End reasons include `pair_mismatch`, `stellar_type`,
`invalid_orbit`, `eccentricity`, `invalid_bse_time`, `stellar_event`,
`contact_distance`, `clearance`, `binding_energy`, `external_perturbation`,
`bse_due`, `bse_update`, `tidal_orbit_change`, and `merger`.
When multiple eligibility checks fail, the first failed check is logged.
Thaw logging checks the flags before clearing them, avoiding duplicate ends.
