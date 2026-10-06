#!/usr/bin/env python3
"""Regression: python3 test/darkmatter_galpy.py /path/to/no-BSE-DM-Galpy-petar.

Requires numpy, scipy and mpirun. All generated fixtures stay in /tmp.
"""
import argparse
import os
from pathlib import Path
import subprocess
import tempfile

import numpy as np
from scipy.integrate import solve_ivp

G = 0.00449830997959438
DT = 1 / 1024
END = 1 / 16
OFFSET = np.array([3., .5, .1, .2, -.1, .03])


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('binary', type=Path)
    args = parser.parse_args()
    binary = args.binary.resolve()
    root = Path(tempfile.mkdtemp(prefix='petar-dm-galpy-'))
    print(f'Test outputs: {root}', flush=True)
    rng = np.random.default_rng(42)
    x = rng.uniform(-.6, .6, (24, 3))
    x[-1] = [-3., -.5, .2]  # DM exactly on the external potential's z axis.
    v = rng.uniform(-.03, .03, (24, 3))
    mass = np.r_[np.full(8, .01), np.full(16, .001)]
    stars = np.zeros((8, 21))
    stars[:, 0] = mass[:8]
    stars[:, 1:4], stars[:, 4:7] = x[:8], v[:8]
    stars[:, 9] = np.arange(1, 9)
    dm = np.zeros((16, 13))
    dm[:, 0] = mass[8:]
    dm[:, 1:4], dm[:, 4:7] = x[8:], v[8:]
    dm[:, 7] = np.arange(1, 17)
    star_fmt = ['%.17g']*21
    for column in (9, 11, 20):
        star_fmt[column] = '%d'
    dm_fmt = ['%.17g']*13
    dm_fmt[7] = '%d'
    np.savetxt(root/'initial.star', stars, fmt=star_fmt,
               header='0 8 0 '+' '.join(map(str, OFFSET)), comments='')
    np.savetxt(root/'initial.dm', dm, fmt=dm_fmt, header='0 16 0 0.05', comments='')
    for name, mode, changes in [('moving', 2, 'Nchange 0'),
                                 ('evolving', 0, 'Nchange 1 Index 0\nChangeMode 1\nChangeRate -0.2')]:
        (root/f'{name}.conf').write_text(
            f'Time 0 Task add\nNset 1\nSet 0\nNtype 1 Mode {mode}\n'
            f'GM 5 Pos 0 0 0 Vel 0.1 0 0\nType 17\nArg 5 2\n{changes}\n')

    def run(name, ranks=1, extra=(), source=None, fmt=1):
        directory = root/name
        directory.mkdir()
        starfile, dmfile = source or (root/'initial.star', root/'initial.dm')
        command = ['mpirun', '--bind-to', 'none', '-np', str(ranks), str(binary),
                   '-u', '1', '-i', str(fmt), '-a', '0', '-b', '0', '-s', str(DT),
                   '-r', '.001', '--soft-eps', '.02', '-T', '0',
                   '-o', str(END/2), '-t', str(END), *map(str, extra),
                   '--dm-input', str(dmfile), str(starfile)]
        env = dict(os.environ, OMP_NUM_THREADS='2', OMP_STACKSIZE='128M',
                   OPENBLAS_NUM_THREADS='1')
        with (directory/'run.log').open('w') as log:
            subprocess.run(command, cwd=directory, env=env, stdout=log,
                           stderr=subprocess.STDOUT, check=True, timeout=120)
        return directory

    def read(directory, index=2):
        with (directory/f'data.{index}').open() as stream:
            header = np.array(list(map(float, stream.readline().split())))
        s = np.loadtxt(directory/f'data.{index}', skiprows=1, ndmin=2)
        d = np.loadtxt(directory/f'dmdata.{index}', skiprows=1, ndmin=2)
        s, d = s[np.argsort(s[:, 9])], d[np.argsort(d[:, 7])]
        a = np.r_[s[:, :7], d[:, :7]]
        a[:, 1:4] += header[3:6]
        a[:, 4:7] += header[6:9]
        return a, d

    def check_close(a, b, label, tol=2e-9):
        error = np.max(np.abs(a-b))
        assert error < tol, (label, error)
        print(f'PASS {label}: max error {error:.3g}', flush=True)

    static = ('--galpy-type-arg', '17:5,2')
    one = run('static-1', extra=static)
    two = run('static-2', ranks=2, extra=static)
    initial, initial_dm = read(one, 0)
    positions = initial[:, 1:4]
    delta = positions[8:, None, :] - positions[None, :, :]
    r2 = np.sum(delta**2, axis=2) + .05**2
    internal_acc = -G*np.sum(mass[None, :, None]*delta/r2[:, :, None]**1.5, axis=1)
    external_acc = -5*positions[8:]/(np.sum(positions[8:]**2, axis=1)+4)[:, None]**1.5
    check_close(initial_dm[:, 8:11], internal_acc+external_acc, 'DM acceleration vs direct + analytic Plummer')
    status = np.loadtxt(one/'dmdata.status')
    expected_u = np.sum(-5*mass[8:]/np.sqrt(np.sum(positions[8:]**2, axis=1)+4))
    check_close(status[0, 7], expected_u, 'full external potential energy', 1e-12)
    check_close(status[0, 3], np.sum(.5*mass[8:]*np.sum(initial[8:, 4:7]**2, axis=1)),
                'inertial DM kinetic energy', 1e-12)
    check_close(status[:, 6], status[:, 3:6].sum(axis=1)+status[:, 7], 'DM energy sum', 1e-11)
    check_close(read(one)[0], read(two)[0], '1 vs 2 MPI ranks')

    # Independent inertial-frame integration tests kicks and multiple
    # moving-frame recentering events, including DM's frame corrections.
    eps2 = np.full((24, 24), .05**2)
    eps2[:8, :8] = .02**2
    def rhs(t, y):
        pos, vel = y[:72].reshape(24, 3), y[72:].reshape(24, 3)
        dr = pos[:, None, :]-pos[None, :, :]
        rr = np.sum(dr**2, axis=2)+eps2
        acc = -G*np.sum(mass[None, :, None]*dr/rr[:, :, None]**1.5, axis=1)
        acc -= 5*pos/(np.sum(pos**2, axis=1)+4)[:, None]**1.5
        return np.r_[vel.ravel(), acc.ravel()]
    y0 = np.r_[initial[:, 1:4].ravel(), initial[:, 4:7].ravel()]
    ref = solve_ivp(rhs, (0, END), y0, method='DOP853', rtol=1e-12, atol=1e-14).y[:, -1]
    check_close(read(one)[0][:, 1:7], np.c_[ref[:72].reshape(24, 3), ref[72:].reshape(24, 3)],
                'orbit vs inertial high-accuracy reference', 2e-7)
    resumed = run('static-restart', ranks=2, source=(one/'data.1', one/'dmdata.1'))
    check_close(read(one)[0], read(resumed)[0], 'static restart without repeating potential CLI')

    binary_run = run('binary', extra=static, fmt=2)
    assert (binary_run/'dmdata.1').stat().st_size == 32+16*104
    binary_restart = run('binary-restart', ranks=2, fmt=3,
                         source=(binary_run/'data.1', binary_run/'dmdata.1'))
    check_close(read(one)[0], read(binary_restart)[0], 'binary paired restart')
    zero = run('zero-potential')
    assert np.all(np.loadtxt(zero/'dmdata.status')[:, 7] == 0)
    print('PASS no external potential', flush=True)
    status_only = run('status-only', extra=(*static, '-w', '3'))
    check_close(np.loadtxt(status_only/'dmdata.status'), status,
                'status-only DM energy output')

    # The stellar one-particle shortcut must not bypass live DM evolution.
    np.savetxt(root/'one.star', stars[:1], fmt=star_fmt,
               header='0 1 0 '+' '.join(map(str, OFFSET)), comments='')
    single = run('one-star', extra=static,
                 source=(root/'one.star', root/'initial.dm'))
    assert np.loadtxt(single/'dmdata.status')[-1, 0] == END
    assert np.max(np.abs(read(single)[0][1:, 1:4]-positions[8:])) > 1e-4
    print('PASS one star with live DM', flush=True)

    for name in ('moving', 'evolving'):
        extra = ('--galpy-conf-file', root/f'{name}.conf')
        full = run(name+'-1', extra=extra)
        parallel = run(name+'-2', ranks=2, extra=extra)
        resumed = run(name+'-restart', ranks=2, extra=extra,
                      source=(full/'data.1', full/'dmdata.1'))
        check_close(read(full)[0], read(parallel)[0], name+' MPI agreement')
        check_close(read(full)[0], read(resumed)[0], name+' restart')
        check_close(np.fromstring((full/'data.2.galpy.state').read_text(), sep=' '),
                    np.fromstring((parallel/'data.2.galpy.state').read_text(), sep=' '),
                    name+' potential state across ranks')
    print('All DARKMATTER + GALPY regressions passed.', flush=True)


if __name__ == '__main__':
    main()
