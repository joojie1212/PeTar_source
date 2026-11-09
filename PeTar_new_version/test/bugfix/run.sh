#!/usr/bin/env bash
set -euo pipefail
repo=$(cd -- "$(dirname -- "$0")/../../.." && pwd)
cd "$repo"
output=${BUGFIX_BUILD_DIR:-/tmp/petar-bugfix-tests}
mkdir -p "$output"
inc=(-I SDAR/src -I FDPS/src -I PeTar_new_version/src -I PeTar_new_version/bse-interface -I PeTar_new_version/bhmerger)
san=(-g -O1 -no-pie -fsanitize=address,undefined -fno-omit-frame-pointer)
g++ -std=c++17 "${san[@]}" -DINTEGRATED_CUTOFF_FUNCTION -DHARD_CHECK_ENERGY "${inc[@]}" PeTar_new_version/test/bugfix/core.cxx -o "$output/core"
# LeakSanitizer cannot enumerate processes inside the managed sandbox.
# Set BUGFIX_DETECT_LEAKS=1 in an unrestricted environment for the leak check.
timeout 60s env ASAN_OPTIONS=detect_leaks=${BUGFIX_DETECT_LEAKS:-0} "$output/core"
# Exercise both full internal gravity and directional changeover tie handling.
for cutoff in off on; do
    cutoff_flags=()
    if [[ "$cutoff" == on ]]; then cutoff_flags=(-DAR_CHANGEOVER); fi
g++ -std=c++17 "${san[@]}" "${cutoff_flags[@]}" -fopenmp -march=core-avx2 -DDEV_CODE -DUSE_SIMD -DINTRINSIC_X86 -DPARTICLE_SIMULATOR_THREAD_PARALLEL -DSTELLAR_EVOLUTION -DBSE_BASE -DBSEBBF -DBHMERGER -DAR_TTL -DAR_SLOWDOWN_TREE -DAR_SLOWDOWN_TIMESCALE -DHARD_CHECK_ENERGY -DRECORD_CM_IN_HEADER -DCLUSTER_VELOCITY -DPROFILE "${inc[@]}" PeTar_new_version/test/bugfix/integration.cxx PeTar_new_version/bse-interface/libbse.a -lgfortran -o "$output/integration"
tls_hex=$(readelf -SW "$output/integration" | awk '$2==".tbss" {print $6}')
tls_bytes=$((16#${tls_hex:-0}))
if (( tls_bytes > 1048576 )); then
    echo "static TLS regression: .tbss is $tls_bytes bytes" >&2
    exit 1
fi
timeout 60s env ASAN_OPTIONS=detect_leaks=${BUGFIX_DETECT_LEAKS:-0} "$output/integration"
done
g++ -std=c++17 -g -O1 "-DPETAR_PRECESSION_SOURCE_DIR=\"$repo/PeTar_new_version/bhmerger/precession-master\"" -I/usr/include/python3.11 -I/usr/lib/python3/dist-packages/numpy/core/include PeTar_new_version/test/bugfix/python_bridge.cxx -lpython3.11 -o "$output/python_bridge"
timeout 60s "$output/python_bridge"
timeout 60s "$output/python_bridge" --host
(cd /tmp && timeout 60s "$output/python_bridge" --import-only)
timeout 60s python3 PeTar_new_version/test/bugfix/precession_sign.py

# Trap invalid arithmetic in the actual production Fortran routines.
gfortran -g -O0 -ffpe-trap=invalid,zero,overflow -fcheck=all \
    PeTar_new_version/test/bugfix/sn_kick.f90 \
    PeTar_new_version/bse-interface/bse/kick.f \
    PeTar_new_version/bse-interface/bse/hrdiag.f \
    PeTar_new_version/bse-interface/libbse.a -o "$output/sn_kick"
timeout 60s "$output/sn_kick"

g++ -std=c++17 "${san[@]}" -DBSEBBF -DBSE_BASE \
    -I PeTar_new_version/bse-interface \
    PeTar_new_version/test/bugfix/stellar_merge.cxx \
    PeTar_new_version/bse-interface/libbse.a -lgfortran -o "$output/stellar_merge"
timeout 60s env ASAN_OPTIONS=detect_leaks=${BUGFIX_DETECT_LEAKS:-0} "$output/stellar_merge"
g++ -std=c++17 "${san[@]}" -fopenmp -DAR_TTL -DAR_SLOWDOWN_ARRAY \
    -DAR_SLOWDOWN_TIMESCALE -DHARD_CHECK_ENERGY "${inc[@]}" \
    PeTar_new_version/test/bugfix/array_merger.cxx -o "$output/array_merger"
timeout 60s env ASAN_OPTIONS=detect_leaks=${BUGFIX_DETECT_LEAKS:-0} "$output/array_merger"

for variant in default KDKDK_4TH; do
    flags=()
    if [[ "$variant" != default ]]; then flags+=("-D$variant"); fi
    g++ -std=c++17 "${san[@]}" "${flags[@]}" "${inc[@]}" \
        PeTar_new_version/test/bugfix/kickdrift.cxx -o "$output/kickdrift-$variant"
    env ASAN_OPTIONS=detect_leaks=${BUGFIX_DETECT_LEAKS:-0} "$output/kickdrift-$variant"
done
g++ -std=c++17 "${san[@]}" -fopenmp -DBHMERGER "${inc[@]}" \
    PeTar_new_version/test/bugfix/merger_log.cxx -Wl,--wrap=write -o "$output/merger_log"
env ASAN_OPTIONS=detect_leaks=${BUGFIX_DETECT_LEAKS:-0} "$output/merger_log"
g++ -std=c++17 "${san[@]}" -DINTEGRATED_CUTOFF_FUNCTION -DHARD_CHECK_ENERGY \
    "${inc[@]}" PeTar_new_version/test/bugfix/cm_prediction.cxx -o "$output/cm_prediction"
env ASAN_OPTIONS=detect_leaks=${BUGFIX_DETECT_LEAKS:-0} "$output/cm_prediction"

# Fixed-step triples retain their policy and rejected trials do not replay callbacks.
g++ -std=c++17 "${san[@]}" -fopenmp -DAR_TTL -DAR_SLOWDOWN_TREE \
    -DAR_SLOWDOWN_TIMESCALE -DHARD_CHECK_ENERGY "${inc[@]}" \
    PeTar_new_version/test/bugfix/fixed_retry.cxx -o "$output/fixed_retry"
env ASAN_OPTIONS=detect_leaks=${BUGFIX_DETECT_LEAKS:-0} "$output/fixed_retry"
