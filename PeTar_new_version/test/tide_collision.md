# Contact orbit dispatch regression

With PeTar configured for BSE, run from the project directory:

```sh
make build/petar.tide_collision.test
OMP_NUM_THREADS=1 build/petar.tide_collision.test
```

The test calls `ARInteraction::modifyAndInterruptIter` with two main-sequence
stars on elliptic and hyperbolic encounters. Predicted pericentre overlap must
not merge stars while their instantaneous positions remain detached. Cases
placed physically inside the summed stellar radius at peri-centre must merge.
It checks merger status and finite output particle states as well as survival.

Nominal tangency can round to either side after converting the orbit to
particles. The tangency case uses the reconstructed pericentre; the cases
1e-12 to either side separately exercise both outcomes.

The tide model is skipped when its predicted pericentre is inside the stellar
surfaces, but that prediction alone no longer causes an irreversible merger.
These are controlled local encounters, not a replay of the original MPI job.
