# Distributed dark-matter regression

From the PeTar source directory (FDPS in `../FDPS`):

```sh
mpicxx -std=c++17 -O2 -fopenmp \
  -DPARTICLE_SIMULATOR_MPI_PARALLEL -DPARTICLE_SIMULATOR_THREAD_PARALLEL \
  -I../FDPS/src -Isrc test/darkmatter_fdps.cxx -o /tmp/petar-dm-fdps-test
OMP_NUM_THREADS=2 mpirun -np 2 /tmp/petar-dm-fdps-test
```

For AVX2, add `-march=core-avx2 -DUSE_SIMD -DINTRINSIC_X86`.
For double precision, additionally add `-DCALC_EP_64bit`.
Run with 1, 2, and 4 ranks. The test uses direct summation as an independent
reference, distributed source/target ownership, probe-only cells, EP-SP forces,
the full cross energy of a star/DM pair, and softened self-potential removal.
The tolerances are `2e-11` for scalar/double SIMD and `3e-5` for float SIMD.

For end-to-end validation use paired ASCII stellar/DM initial conditions, run
with `-T 0` on one and two ranks, and compare snapshots sorted by each species'
ID. Repeat by restarting both snapshots from an intermediate output into a new
directory with `-a 0`. Check total energy as stellar status `Ekin + Epot` plus
DM status `E_dm_and_cross`; stellar status alone excludes the cross interaction.
