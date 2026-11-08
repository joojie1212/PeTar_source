# Fixed-step multiple retry and change-log shards (2026-11-08)

Stable multiples entering AR with `fix_step=later` retain that policy. The existing
post-interruption, pre-gravity backup is reused. A gravitational trial with an
invalid time increment, invalid final TTL denominator, or nonfinite/over-tolerance
extended-Hamiltonian step error is rejected. Both buffered ds values are halved;
the next iteration restores the existing integration backup and recomputes the
binary CM. That iteration skips the BSE/merger callback. Adaptive policy remains
active through the current call; the original `later` policy is restored on return.
Sixty-four consecutive failed trials terminate explicitly. Two-body and explicit
`always` policies retain their existing behavior. Existing TTL update recovery
and time-stall safeguards are retained.

This is step rollback, NOT full-interval restart. BSE, stellar metadata, random
state and event records precede the backup and are not undone or replayed. No
additional complete-particle copy is added per step. Error monitoring still runs;
there is no promise of a speedup for strongly interacting multiples. Repeated
fallback can cost more than continuous adaptation. AR_ADAPTIVE_MULTIPLE output
has been removed, including its formatting and critical section.

Production change logs are now `PREFIX.change.RANK.THREAD`. Each OpenMP worker
owns a persistent buffered stream, with no per-row OpenMP or OS file lock and no
per-row open/stat/close. Streams open once at run initialization, append without
repeating headers on restart, and flush at status/snapshot output and normal
stream destruction. Copies of ARInteraction share stream ownership. Nested
OpenMP writers are rejected because their thread IDs would collide. The existing
shared bhmerger event log and legacy openChange API remain unchanged.

Use fresh output directories for new runs. On restarting, preserve all shards,
including those from thread counts used earlier. An abrupt kill may lose buffered
rows since the most recent flush; snapshots trigger a flush. Do not assume the
shards are globally time ordered. Old analysis scripts expecting one data.change
must first merge the closed/quiescent logs:

    python3 PeTar_new_version/tools/merge_change_logs.py /run/data.change --output /run/data.change.merged

For runs continued from the old logger, add --include-legacy. Source files are
retained. The destination must be new. The merger sorts by time_new, time_old and
a deterministic tie-break; incomplete rows are rejected, not silently discarded.

## Validation

- Production MPI/OpenMP/AVX2/TT3/quad main binary and hard-debug build passed.
- Fixed-step hierarchical triple: retained `later`; callback times strictly
  increase. A temporary-header fault at step 8 forces fallback; finite oversized
  error and NaN cases pass without duplicate callbacks. No fault hooks are added
  to production code. Reproduce with `python3 PeTar_new_version/test/fixed_retry_fault.py`.
- Two simulated rank namespaces x four real OpenMP threads: 200 rows per shard,
  one header across append reopen, all 19 fields parse; copied logger ownership
  checked. Existing partial-write rollback test for shared merger logs retained.
- Merge utility: four input shards produce four correctly ordered records.
- Original stall dump replay: because the old particle binary record lacks
  merger_dyn_log_tgw, a COPY inserts an 8-byte zero at offset 128 in each
  280-byte real and 328-byte artificial record. Original data untouched.
  Checked original size 3808 and converted size 3904. Replay uses current source
  with TT2 matching the dump (production remains TT3).

Same-source reference differs only by forcing multiples to `none` at entry;
it also omits the removed log. Results for the same converted dump:

| AR tolerance | fixed/retry steps | continuous-adaptive steps | fixed/retry relative energy residual |
|---|---:|---:|---:|
| 1e-8 | 2041 | 2041 | 6.8853e-8 |
| 1e-9 | 4878 | 4287 | 6.1942e-8 |

No time stalls, step-limit warnings or slowdown emergency recovery occurred.
At tighter tolerance this case takes MORE steps, so these results establish
correctness on this case, not a performance improvement. No whole-cluster
speedup has been measured. Logs: diagnostics/ar_fixed_retry_validation/.

The build artifacts have not been installed into the user's bin directory;
existing running processes continue using their loaded code.
