#!/usr/bin/env bash
set -euo pipefail
repo=$(cd -- "$(dirname -- "$0")/../.." && pwd)
cd "$repo/mcluster"
build=${ASTRA_MCLUSTER_BUILD_DIR:-/tmp/astra-mcluster-tests}
mkdir -p "$build"
make mcluster_sse
gcc -O0 -g -fopenmp -DSSE -Dmain=mcluster_entry -c main.c -o "$build/main.o"
objects=(deltat.o evolv1.o hrdiag.o kick.o mlwind.o mrenv.o ran3.o star.o zcnsts.o zfuncs.o comenv.o corerd.o dgcore.o evolv2.o gntage.o instar.o mix.o rl.o)
gcc -O0 -g -fopenmp tests/astra_binary_units.c "$build/main.o" "${objects[@]}" \
    -Wl,--wrap=evolv2_ -lm -lgfortran -o "$build/binary_units"
"$build/binary_units"

python3 tests/astra_smoke.py
gcc -std=gnu2x -Werror=implicit-function-declaration -fopenmp -DSSE -DDARKMATTER -fsyntax-only main.c
