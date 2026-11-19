#!/usr/bin/env bash
set -euo pipefail
repo=$(cd -- "$(dirname -- "$0")/../../.." && pwd)
cd "$repo"
output=${HARD_PARALLEL_BUILD_DIR:-/tmp/petar-hard-parallel-regressions}
mkdir -p "$output"
inc=(-I SDAR/src -I FDPS/src -I PeTar_new_version/src -I PeTar_new_version/bse-interface -I PeTar_new_version/bhmerger)
flags=(-std=c++17 -O1 -g -fopenmp -march=core-avx2 -DDEV_CODE -DUSE_SIMD -DINTRINSIC_X86
       -DPARTICLE_SIMULATOR_THREAD_PARALLEL -DSTELLAR_EVOLUTION -DBSE_BASE -DBSEBBF
       -DBHMERGER -DAR_TTL -DAR_SLOWDOWN_TREE -DAR_SLOWDOWN_TIMESCALE
       -DHARD_CHECK_ENERGY -DRECORD_CM_IN_HEADER -DCLUSTER_VELOCITY -DPROFILE)
if [[ ${HARD_PARALLEL_SANITIZE:-0} == 1 ]]; then
    flags+=(-fsanitize=address,undefined -fno-omit-frame-pointer -no-pie)
    export ASAN_OPTIONS=detect_leaks=0
fi
for mode in off on; do
    extra=(); threads=(1)
    if [[ $mode == on ]]; then extra=(-DHARD_PARALLEL); threads=(1 2 3 4); fi
    g++ "${flags[@]}" "${extra[@]}" "${inc[@]}" PeTar_new_version/test/bugfix/hard_parallel.cxx \
        PeTar_new_version/bse-interface/libbse.a -lgfortran -o "$output/hard-$mode"
    for n in "${threads[@]}"; do
        "$output/hard-$mode" "$n" "$output/$mode-$n.digest" > "$output/$mode-$n.log" 2>&1
        # Only drift/H4/cluster correction are compared bitwise. AR SIMD is
        # checked against direct forces with a numerical tolerance in C++.
        cmp "$output/off-1.digest" "$output/$mode-$n.digest"
        echo "HARD_PARALLEL=$mode, threads=$n: passed"
    done
done
# The helpers must not introduce a dependency on libgomp in serial builds.
cat > "$output/serial.cxx" <<'CPP'
#include "Common/inner_parallel.hpp"
#include <cassert>
int main() {
    assert(!InnerPar::inTeam());
    assert(!InnerPar::useTask(128));
    assert(InnerPar::taskCount(128)==1);
}
CPP
g++ -std=c++17 -DHARD_PARALLEL -I SDAR/src "$output/serial.cxx" -o "$output/serial"
"$output/serial"
echo 'Serial OpenMP fallback: passed'
