#!/usr/bin/env bash
set -Eeuo pipefail

# Supervisor example based on the 1013_new_4 BH-merger run.  Build PeTar with
# ./configure --enable-hard-reopt [the other options used by your build]

SCRIPT_DIR=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)
PETAR_ROOT=$(cd -- "$SCRIPT_DIR/.." && pwd)
WORK_DIR=${WORK_DIR:-$PWD}
cd "$WORK_DIR"

GCC_HOME=${GCC_HOME:-"$HOME/tools/gcc-11"}
export PATH="$GCC_HOME/bin:$PATH"
export LD_LIBRARY_PATH="$GCC_HOME/lib64:$GCC_HOME/lib64/gcc/x86_64-linux-gnu/11${LD_LIBRARY_PATH:+:$LD_LIBRARY_PATH}"
export PYTHONPATH="$PETAR_ROOT/bhmerger/precession-master${PYTHONPATH:+:$PYTHONPATH}"
export UCX_VFS_ENABLE=${UCX_VFS_ENABLE:-n}
export OMP_STACKSIZE=${OMP_STACKSIZE:-128M}
export OMP_NUM_THREADS=${OMP_NUM_THREADS:-4}
export OMP_PROC_BIND=${OMP_PROC_BIND:-spread}
export OMP_PLACES=${OMP_PLACES:-cores}

MPI_RANKS=${MPI_RANKS:-1}
INPUT_SNAPSHOT=${INPUT_SNAPSHOT:-input}
END_TIME=${END_TIME:-0.01}
OUTPUT_INTERVAL=${OUTPUT_INTERVAL:-0.0001}
OUTPUT_PREFIX=${OUTPUT_PREFIX:-data}
HARD_REOPT_RATIO=${HARD_REOPT_RATIO:-0.9}
HARD_REOPT_WINDOWS=${HARD_REOPT_WINDOWS:-3}
HARD_REOPT_COOLDOWN=${HARD_REOPT_COOLDOWN:-10}
HARD_REOPT_EXIT_CODE=42
RUN_LOG_PREFIX=${RUN_LOG_PREFIX:-hard_reopt_run}
FIND_DT_LOG_PREFIX=${FIND_DT_LOG_PREFIX:-hard_reopt_find_dt}

require_command() {
    command -v "$1" >/dev/null 2>&1 || {
        printf 'Error: required command not found: %s\n' "$1" >&2
        exit 2
    }
}

request_value() {
    local key=$1 file=$2
    awk -F= -v key="$key" '$1 == key {value=substr($0, index($0, "=")+1)} END {print value}' "$file"
}

require_command mpirun
require_command petar
require_command petar.find.dt

petar_bin=$(readlink -f "$(command -v petar)")
current_snapshot=$INPUT_SNAPSHOT
selected_dt=
cycle=0

[[ -s "$current_snapshot" ]] || {
    printf 'Error: input snapshot is missing or empty: %s\n' "$current_snapshot" >&2
    exit 2
}

while :; do
    cycle=$((cycle+1))
    request_file=${OUTPUT_PREFIX}.hard_reopt
    run_log=${RUN_LOG_PREFIX}.${cycle}.log
    rm -f -- "$request_file" "$request_file.tmp"

    common_args=(
        -f "$OUTPUT_PREFIX"
        -a 1
        -u 1
        --write-group-info 0
        --stellar-evolution 2
        --energy-err-hard 1e-3
        --hard-reopt-ratio "$HARD_REOPT_RATIO"
        --hard-reopt-windows "$HARD_REOPT_WINDOWS"
        --hard-reopt-cooldown "$HARD_REOPT_COOLDOWN"
        -t "$END_TIME"
        -o "$OUTPUT_INTERVAL"
    )
    restart_args=()
    if [[ -n "$selected_dt" ]]; then
        restart_args=(
            -p input.par
            -s "$selected_dt"
            -r 0
            --r-search-min 0
            --r-bin 0
        )
    fi

    printf 'Starting PeTar cycle %d from %s%s\n' \
        "$cycle" "$current_snapshot" "${selected_dt:+ with dt=$selected_dt}"
    set +e
    mpirun --mca btl ^openib --bind-to none -np "$MPI_RANKS" "$petar_bin" \
        "${restart_args[@]}" "${common_args[@]}" "$current_snapshot" \
        2>&1 | tee "$run_log"
    petar_status=${PIPESTATUS[0]}
    set -e

    if [[ "$petar_status" -eq 0 ]]; then
        echo 'PeTar completed normally.'
        break
    fi
    if [[ "$petar_status" -ne "$HARD_REOPT_EXIT_CODE" ]]; then
        printf 'Error: PeTar failed with status %d; not restarting.\n' "$petar_status" >&2
        exit "$petar_status"
    fi
    [[ -s "$request_file" ]] || {
        printf 'Error: PeTar returned %d but %s is missing or empty.\n' \
            "$HARD_REOPT_EXIT_CODE" "$request_file" >&2
        exit 2
    }

    request_exit=$(request_value exit_code "$request_file")
    next_snapshot=$(request_value snapshot "$request_file")
    [[ "$request_exit" == "$HARD_REOPT_EXIT_CODE" && -s "$next_snapshot" ]] || {
        printf 'Error: invalid hard-reoptimization request in %s.\n' "$request_file" >&2
        exit 2
    }
    current_snapshot=$next_snapshot

    [[ -s input.par ]] || {
        echo 'Error: input.par is required for restart optimization.' >&2
        exit 2
    }
    find_log=${FIND_DT_LOG_PREFIX}.${cycle}.log
    printf 'Optimizing tree timestep from %s\n' "$current_snapshot"
    # find_dt uses -w 0 internally.  Explicitly disable the trigger in its
    # short benchmark processes because no restart snapshot is written there.
    petar.find.dt \
        -p "$petar_bin" \
        -m "$MPI_RANKS" \
        -o "$OMP_NUM_THREADS" \
        -a "-p input.par -u 1 --stellar-evolution 2 --energy-err-hard 1e-3 -r 0 --r-search-min 0 --r-bin 0 --hard-reopt-ratio -1 --hard-reopt-windows -1 --hard-reopt-cooldown -1" \
        "$current_snapshot" 2>&1 | tee "$find_log"

    selected_dt=$(awk '/Best performance choice: tree step:/ {value=$NF} END {print value}' "$find_log")
    if [[ ! "$selected_dt" =~ ^[0-9]+([.][0-9]*)?([eE][-+]?[0-9]+)?$ ]]; then
        printf 'Error: could not parse optimized timestep from %s.\n' "$find_log" >&2
        exit 2
    fi
    printf 'Selected tree timestep %s; restarting from %s.\n' \
        "$selected_dt" "$current_snapshot"
done
