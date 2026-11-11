#!/usr/bin/env bash
set -Eeuo pipefail

SCRIPT_DIR=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)
SUPERVISOR="$SCRIPT_DIR/supervisor.sh"
# The batch launcher supplies an existing, independent output directory.
CASE_DIR=${CASE_DIR:-$SCRIPT_DIR}
CASE_DIR=$(cd -- "$CASE_DIR" && pwd)

# Reference model: sample/bhmerger/1013_new_4/test6.sh
REFERENCE_N=${REFERENCE_N:-10000}
REFERENCE_HALF_MASS_RADIUS_PC=${REFERENCE_HALF_MASS_RADIUS_PC:-0.0013047672}
TARGET_N=${TARGET_N:-500}
PARTICLE_MASS_MSUN=${PARTICLE_MASS_MSUN:-30}
PLUMMER_RH_OVER_A=${PLUMMER_RH_OVER_A:-1.30476603}
G_PC3_MSUN_MYR2=${G_PC3_MSUN_MYR2:-0.00449830997959438}
COULOMB_GAMMA=${COULOMB_GAMMA:-0.02}
TARGET_TRH=${TARGET_TRH:-100}
OUTPUT_EVERY_TRH=${OUTPUT_EVERY_TRH:-0.2}
ESCAPE_RADIUS_RH=${ESCAPE_RADIUS_RH:-100}
MCLUSTER_SEED=${MCLUSTER_SEED:-72226}

MPI_RANKS=${MPI_RANKS:-1}
export OMP_STACKSIZE=${OMP_STACKSIZE:-128M}
export OMP_NUM_THREADS=${OMP_NUM_THREADS:-4}
export OMP_PROC_BIND=${OMP_PROC_BIND:-spread}
export OMP_PLACES=${OMP_PLACES:-cores}
export UCX_VFS_ENABLE=${UCX_VFS_ENABLE:-n}


MODEL_PREFIX=${MODEL_PREFIX:-cluster_n500}
MODEL_FILE="$MODEL_PREFIX.dat.10"
INPUT_SNAPSHOT=${INPUT_SNAPSHOT:-input}
OUTPUT_PREFIX=${OUTPUT_PREFIX:-data}
FOREGROUND=${FOREGROUND:-0}

require_command() {
    command -v "$1" >/dev/null 2>&1 || {
        printf 'Error: required command not found: %s\n' "$1" >&2
        exit 2
    }
}

validate_configuration() {
    awk -v n="$TARGET_N" -v nref="$REFERENCE_N" \
        -v rhref="$REFERENCE_HALF_MASS_RADIUS_PC" \
        -v mass="$PARTICLE_MASS_MSUN" -v grav="$G_PC3_MSUN_MYR2" \
        -v gamma="$COULOMB_GAMMA" -v ntrh="$TARGET_TRH" \
        -v output_trh="$OUTPUT_EVERY_TRH" -v escape_rh="$ESCAPE_RADIUS_RH" \
        'BEGIN {
            valid = (n == int(n) && n > 0 && nref == int(nref) && nref > 0 &&
                     rhref > 0 && mass > 0 && grav > 0 && gamma*n > 1 &&
                     ntrh > 0 && output_trh > 0 &&
                     escape_rh ~ /^[0-9]+([.][0-9]*)?([eE][-+]?[0-9]+)?$/ &&
                     escape_rh+0 > 0)
            exit !valid
        }' || {
        echo 'Error: invalid model parameters (also require COULOMB_GAMMA*TARGET_N > 1).' >&2
        exit 2
    }
}

validate_configuration

# For homologous Plummer models with the same particle mass,
# rho_0 proportional to N/r_h^3.  Thus fixed rho_0 requires r_h proportional
# to N^(1/3).
HALF_MASS_RADIUS_PC=$(awk -v rh="$REFERENCE_HALF_MASS_RADIUS_PC" \
    -v n="$TARGET_N" -v nref="$REFERENCE_N" \
    'BEGIN {printf "%.15g", rh*(n/nref)^(1.0/3.0)}')
# Fixed multiple of the initial rh; positive --r-escape also checks energy.
ESCAPE_RADIUS_PC=$(awk -v rh="$HALF_MASS_RADIUS_PC" -v factor="$ESCAPE_RADIUS_RH" \
    'BEGIN {printf "%.15g", rh*factor}')
PLUMMER_SCALE_RADIUS_PC=$(awk -v rh="$HALF_MASS_RADIUS_PC" \
    -v ratio="$PLUMMER_RH_OVER_A" \
    'BEGIN {printf "%.15g", rh/ratio}')
CENTRAL_DENSITY_MSUN_PC3=$(awk -v n="$TARGET_N" \
    -v mass="$PARTICLE_MASS_MSUN" -v scale="$PLUMMER_SCALE_RADIUS_PC" \
    'BEGIN {pi=atan2(0,-1); printf "%.15g", 3*n*mass/(4*pi*scale^3)}')

# Same definition as PeTar_new_version/tools/analysis/functions.py::calcTrh:
# trh = 0.138 sqrt(N) rh^(3/2) / [sqrt(G m) ln(gamma N)].
INITIAL_TRH_MYR=$(awk -v n="$TARGET_N" -v rh="$HALF_MASS_RADIUS_PC" \
    -v mass="$PARTICLE_MASS_MSUN" -v grav="$G_PC3_MSUN_MYR2" \
    -v gamma="$COULOMB_GAMMA" \
    'BEGIN {printf "%.15g", 0.138*sqrt(n)*rh^1.5/(sqrt(grav*mass)*log(gamma*n))}')
END_TIME_MYR=$(awk -v trh="$INITIAL_TRH_MYR" -v count="$TARGET_TRH" \
    'BEGIN {printf "%.15g", trh*count}')
OUTPUT_INTERVAL_MYR=$(awk -v trh="$INITIAL_TRH_MYR" -v count="$OUTPUT_EVERY_TRH" \
    'BEGIN {printf "%.15g", trh*count}')

print_configuration() {
    printf '%s\n' \
        "Reference N                  = $REFERENCE_N" \
        "Target N                     = $TARGET_N" \
        "Particle mass                = $PARTICLE_MASS_MSUN Msun" \
        "Target half-mass radius      = $HALF_MASS_RADIUS_PC pc" \
        "Escape radius ($ESCAPE_RADIUS_RH initial rh) = $ESCAPE_RADIUS_PC pc (energy >= 0 required)" \
        "Plummer scale radius         = $PLUMMER_SCALE_RADIUS_PC pc" \
        "Analytic central density     = $CENTRAL_DENSITY_MSUN_PC3 Msun/pc^3" \
        "Initial half-mass trh         = $INITIAL_TRH_MYR Myr" \
        "End time ($TARGET_TRH trh)          = $END_TIME_MYR Myr" \
        "Output interval ($OUTPUT_EVERY_TRH trh)    = $OUTPUT_INTERVAL_MYR Myr" \
        "McLuster seed                = $MCLUSTER_SEED"
}

prepare_case() {
    require_command mcluster
    require_command petar.init
    cd "$CASE_DIR"

    if [[ ! -s "$MODEL_FILE" ]]; then
        printf 'Generating %s with N=%s and r_h=%s pc...\n' \
            "$MODEL_FILE" "$TARGET_N" "$HALF_MASS_RADIUS_PC"
        mcluster \
            -N "$TARGET_N" -P 0 -m "$PARTICLE_MASS_MSUN" \
            -R "$HALF_MASS_RADIUS_PC" -f 0 -u 1 -C 5 \
            -s "$MCLUSTER_SEED" -o "$MODEL_PREFIX" >mc.log
    else
        printf 'Reusing existing McLuster model: %s\n' "$MODEL_FILE"
    fi

    if [[ ! -s "$INPUT_SNAPSHOT" ]]; then
        printf 'Creating PeTar input snapshot: %s\n' "$INPUT_SNAPSHOT"
        petar.init -s bse --blackhole -v kms2pcmyr \
            -f "$INPUT_SNAPSHOT" "$MODEL_FILE"
    else
        printf 'Reusing existing PeTar input snapshot: %s\n' "$INPUT_SNAPSHOT"
    fi
}

run_case() {
    require_command petar
    require_command petar.find.dt
    [[ -r "$SUPERVISOR" ]] || {
        printf 'Error: supervisor is missing or unreadable: %s\n' "$SUPERVISOR" >&2
        exit 2
    }

    prepare_case
    ulimit -S -s 131072
    print_configuration

    WORK_DIR="$CASE_DIR" \
    MPI_RANKS="$MPI_RANKS" \
    OMP_NUM_THREADS="$OMP_NUM_THREADS" \
    INPUT_SNAPSHOT="$INPUT_SNAPSHOT" \
    END_TIME="$END_TIME_MYR" \
    OUTPUT_INTERVAL="$OUTPUT_INTERVAL_MYR" \
    OUTPUT_PREFIX="$OUTPUT_PREFIX" \
    ESCAPE_RADIUS_PC="$ESCAPE_RADIUS_PC" \
        bash "$SUPERVISOR"

    echo "Simulation finished; automatic post-processing is disabled."
}

usage() {
    cat <<'EOF'
Usage: ./run.sh [--prepare-only | --print-config | --foreground | --run]

  --prepare-only  Generate the N=500 PeTar input without starting integration.
  --print-config  Print the scaled radius, trh, end time, and density.
  --foreground    Prepare and run attached to this terminal.
  --run           Internal entry point used by the detached workflow.

With no option, the workflow is started in the background and its PID is
written to workflow.pid.  Follow it with: tail -f workflow.log
EOF
}

case "${1:-}" in
    --prepare-only)
        prepare_case
        print_configuration
        exit
        ;;
    --print-config)
        print_configuration
        exit
        ;;
    --foreground)
        run_case
        exit
        ;;
    --run)
        run_case
        exit
        ;;
    -h|--help)
        usage
        exit
        ;;
    "")
        ;;
    *)
        usage >&2
        exit 2
        ;;
esac

cd "$CASE_DIR"
if [[ -s workflow.pid ]] && kill -0 "$(<workflow.pid)" 2>/dev/null; then
    printf 'Error: workflow is already running (PID %s).\n' "$(<workflow.pid)" >&2
    exit 1
fi

nohup "$0" --run >workflow.log 2>&1 &
workflow_pid=$!
printf '%s\n' "$workflow_pid" >workflow.pid
printf 'Started workflow PID %s; follow it with: tail -f %q\n' \
    "$workflow_pid" "$CASE_DIR/workflow.log"
