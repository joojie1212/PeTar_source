#!/usr/bin/env python3
"""Four concurrent exclusive Slurm steps, one 4-thread MPI rank per model."""
import argparse
import csv
import math
import os
from pathlib import Path
import shutil
import signal
import subprocess
import tempfile
import time

HERE = Path(__file__).resolve().parent
NS = (500, 1000, 2000, 4000, 10000)
RH_REFS = ('0.0013047672', '0.013047672', '0.13047672')
FIXED = dict(REFERENCE_N='10000', PARTICLE_MASS_MSUN='30',
             PLUMMER_RH_OVER_A='1.30476603', G_PC3_MSUN_MYR2='0.00449830997959438',
             COULOMB_GAMMA='0.11', TARGET_TRH='300', OUTPUT_EVERY_TRH='0.5',
             ESCAPE_RADIUS_RH='100', MPI_RANKS='1', OMP_NUM_THREADS='4',
             INPUT_SNAPSHOT='input', OUTPUT_PREFIX='data',
             RUN_LOG_PREFIX='hard_reopt_run', FIND_DT_LOG_PREFIX='hard_reopt_find_dt',
             HARD_REOPT_RATIO='0.9', HARD_REOPT_WINDOWS='3', HARD_REOPT_COOLDOWN='10')


def models():
    for n in NS:
        for ref in RH_REFS:
            rh = float(ref) * (n / 10000) ** (1 / 3)
            rho = 3 * n * 30 / (4 * math.pi * (rh / 1.30476603) ** 3)
            trh = .138 * math.sqrt(n) * rh ** 1.5 / (math.sqrt(.00449830997959438 * 30) * math.log(.11*n))
            yield dict(case=f'N{n}_rhRef{ref}', N=n, rh_ref_pc=ref,
                       rh_pc=rh, rho0_Msun_pc3=rho, trh_Myr=trh, end_Myr=300*trh)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument('--print-config', action='store_true')
    group.add_argument('--run', action='store_true')
    args = parser.parse_args()
    rows = list(models())
    seed = os.environ.get('MCLUSTER_SEED', '72222')
    if not seed.isdigit() or int(seed) <= 0:
        parser.error('MCLUSTER_SEED must be a positive integer')
    escape_radius = os.environ.get('ESCAPE_RADIUS_RH', '100')
    try:
        if not math.isfinite(float(escape_radius)) or float(escape_radius) <= 0:
            raise ValueError
    except ValueError:
        parser.error('ESCAPE_RADIUS_RH must be a finite positive number')
    case_fixed = dict(FIXED, ESCAPE_RADIUS_RH=escape_radius)
    print('15 cases; 4 concurrent; 1 MPI rank x 4 OpenMP threads each; seed=' + seed, flush=True)
    print('Escape radius / initial rh = ' + escape_radius, flush=True)
    print('30 Msun BHs; Plummer; no primordial binaries; 300 initial trh; output every 0.5 trh', flush=True)
    print('Binary snapshots; no post-processing; requested walltime: 5 days', flush=True)
    print('N\trh_ref_pc\trh_pc\trho0_Msun_pc3\ttrh_Myr\tend_Myr', flush=True)
    for r in rows:
        print(f"{r['N']}\t{r['rh_ref_pc']}\t{r['rh_pc']:.12g}\t{r['rho0_Msun_pc3']:.12g}\t{r['trh_Myr']:.12g}\t{r['end_Myr']:.12g}", flush=True)
    if args.print_config:
        return 0
    if not os.environ.get('SLURM_JOB_ID'):
        parser.error('--run requires a Slurm allocation')
    root = Path(os.environ.get('OUTPUT_ROOT', str(Path.home()/'scratch/petar_n_rh_results'))).resolve()
    root.mkdir(parents=True, exist_ok=True)
    batch = Path(tempfile.mkdtemp(prefix=f"batch_{os.environ['SLURM_JOB_ID']}_seed{seed}_", dir=root))
    print('Batch directory: ' + str(batch), flush=True)
    # Execute preserved helpers, so edits to the submitted bundle cannot change later cases.
    scripts = batch / 'scripts'
    scripts.mkdir()
    for name in ('sweep.py', 'run_n_rh_series.slurm', 'run_case.sh', 'supervisor.sh'):
        shutil.copy2(HERE/name, scripts/name)
    fields = list(rows[0]) + ['seed', 'directory']
    with (batch/'manifest.tsv').open('w') as f:
        writer = csv.DictWriter(f, fields, delimiter='\t')
        writer.writeheader()
        for r in rows:
            r['seed'] = seed
            r['directory'] = str(batch/r['case'])
            Path(r['directory']).mkdir()
            writer.writerow(r)
    running = {}
    pending = rows.copy()
    failed = False
    stop = False

    def terminate(signum, frame):
        nonlocal stop
        stop = True
        # SIGINT to srun terminates the remote step; give it time to reap processes.
        for proc in running:
            proc.send_signal(signal.SIGINT)

    signal.signal(signal.SIGTERM, terminate)
    signal.signal(signal.SIGINT, terminate)
    with (batch/'status.tsv').open('w', buffering=1) as status:
        status.write('case\tstatus\texit_code\telapsed_seconds\ttimestamp\n')

        def record(r, state, code='', elapsed=''):
            status.write(f"{r['case']}\t{state}\t{code}\t{elapsed}\t{time.strftime('%Y-%m-%dT%H:%M:%S%z')}\n")

        for r in rows:
            record(r, 'queued')
        try:
            while pending or running:
                while pending and len(running) < 4 and not stop:
                    r = pending.pop(0)
                    env = os.environ.copy()
                    env.update(case_fixed, MCLUSTER_SEED=seed, TARGET_N=str(r['N']),
                               REFERENCE_HALF_MASS_RADIUS_PC=r['rh_ref_pc'],
                               CASE_DIR=r['directory'], MODEL_PREFIX=f"cluster_n{r['N']}")
                    case = Path(r['directory'])
                    with (case/'config.txt').open('w') as conf:
                        subprocess.run(['bash', str(scripts/'run_case.sh'), '--print-config'],
                                       env=env, stdout=conf, check=True)
                    cmd = ['srun', '--exclusive', '--exact', '--nodes=1', '--ntasks=1',
                           '--cpus-per-task=4', '--mem=4G', '--cpu-bind=cores', '--mpi=none',
                           '--job-name='+r['case'], '/usr/bin/time', '-v', '-o', str(case/'resources.txt'),
                           'bash', str(scripts/'run_case.sh'), '--foreground']
                    log = (case/'workflow.log').open('w')
                    try:
                        proc = subprocess.Popen(cmd, env=env, stdout=log, stderr=subprocess.STDOUT)
                    except BaseException:
                        log.close()
                        raise
                    running[proc] = (r, log, time.monotonic())
                    record(r, 'launched')
                    print('Launched ' + r['case'], flush=True)
                for proc, (r, log, start) in list(running.items()):
                    code = proc.poll()
                    if code is None:
                        continue
                    log.close()
                    state = 'completed' if code == 0 else 'failed'
                    record(r, state, code, f'{time.monotonic()-start:.1f}')
                    failed |= code != 0
                    print(f"{r['case']}: {state} ({code})", flush=True)
                    del running[proc]
                if stop and not running:
                    for r in pending:
                        record(r, 'cancelled')
                    break
                if running:
                    time.sleep(1)
        finally:
            for proc in running:
                if proc.poll() is None:
                    proc.send_signal(signal.SIGINT)
            for proc, (_, log, _) in running.items():
                try:
                    proc.wait(timeout=10)
                except subprocess.TimeoutExpired:
                    proc.kill()
                    proc.wait()
                log.close()
    print('Results: ' + str(batch), flush=True)
    return 1 if failed or stop else 0


if __name__ == '__main__':
    raise SystemExit(main())
