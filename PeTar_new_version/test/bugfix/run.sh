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
g++ -std=c++17 "${san[@]}" -fopenmp -march=core-avx2 -DDEV_CODE -DUSE_SIMD -DINTRINSIC_X86 -DPARTICLE_SIMULATOR_THREAD_PARALLEL -DSTELLAR_EVOLUTION -DBSE_BASE -DBSEBBF -DBHMERGER -DAR_TTL -DAR_SLOWDOWN_TREE -DAR_SLOWDOWN_TIMESCALE -DHARD_CHECK_ENERGY -DRECORD_CM_IN_HEADER -DCLUSTER_VELOCITY -DPROFILE "${inc[@]}" PeTar_new_version/test/bugfix/integration.cxx PeTar_new_version/bse-interface/libbse.a -lgfortran -o "$output/integration"
tls_hex=$(readelf -SW "$output/integration" | awk '$2==".tbss" {print $6}')
tls_bytes=$((16#${tls_hex:-0}))
if (( tls_bytes > 1048576 )); then
    echo "static TLS regression: .tbss is $tls_bytes bytes" >&2
    exit 1
fi
timeout 60s env ASAN_OPTIONS=detect_leaks=${BUGFIX_DETECT_LEAKS:-0} "$output/integration"
g++ -std=c++17 -g -O1 -I/usr/include/python3.11 -I/usr/lib/python3/dist-packages/numpy/core/include PeTar_new_version/test/bugfix/python_bridge.cxx -lpython3.11 -o "$output/python_bridge"
timeout 60s "$output/python_bridge"
timeout 60s python3 PeTar_new_version/test/bugfix/precession_sign.py
