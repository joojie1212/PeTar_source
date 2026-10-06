# Live DM and Galpy regression

Build PeTar with `--enable-darkmatter --with-external=galpy`, MPI, OpenMP,
AVX2 double-precision SIMD, and `--with-interrupt=off`. This test uses real
Galpy Plummer potentials, not a mocked force interface. Run:

```sh
python3 test/darkmatter_galpy.py /absolute/path/to/petar
```

Requires NumPy, SciPy and OpenMPI (`mpirun --bind-to none`). All generated
fixtures, logs and snapshots remain in a new `/tmp/petar-dm-galpy-*` directory.
The installed PeTar executable and existing runs are not modified.

Checks include:

- DM acceleration against direct softened gravity plus analytic Plummer force,
  including a particle on the external field's symmetry axis.
- Full external potential energy and inertial kinetic energy.
- Trajectories against an independent inertial-frame DOP853 integration,
  crossing multiple stellar frame-recentring events.
- One versus two MPI ranks, with two OpenMP threads per rank.
- Static CLI-potential restart without resupplying its force parameters.
- Paired binary restart, preserving the legacy DM binary record size.
- Zero external potential.
- Status-only output and one stellar particle with live DM.
- Moving-potential reactions and linearly evolving potential parameters,
  including MPI agreement and mid-run restart.

The pair of `data.N` and `dmdata.N` shares the stellar header's reference frame.
`data.N.galpy.state` is the common Galpy state; preserve the legacy
`data.N.galpy` and the input configuration schedule for evolving models.
