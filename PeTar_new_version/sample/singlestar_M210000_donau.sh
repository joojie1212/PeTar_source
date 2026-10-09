#!/bin/bash
#DSUB --mpi hmpi
#DSUB -q q_hpcapp
#DSUB -n petar_single_M210000
#DSUB -N 128
#DSUB -nn 1
#DSUB -rpn 128
#DSUB -R "cpu=1;mem=512M"
#DSUB -oo petar_single.%A.out
#DSUB -eo petar_single.%A.err

# Load HPCKit before enabling nounset: setvars.sh references unset variables.
source /work_ssd/software/HPCKit/25.2.1/setvars.sh >/dev/null
set -euo pipefail

export PATH="$HOME/opt/petar-arm-bse/bin:$HOME/bin:$PATH"
export PYTHONPATH="$HOME/opt/petar-arm-bse/include${PYTHONPATH:+:$PYTHONPATH}"
export UCX_VFS_ENABLE=n
export OMP_STACKSIZE=128M
export OMP_NUM_THREADS=4
export OMP_DYNAMIC=false
export OMP_PROC_BIND=false
unset OMP_PLACES
ulimit -s unlimited

PETAR_BIN=${PETAR_BIN:-$HOME/opt/petar-arm-bse/bin/petar}
PETAR_INIT=${PETAR_INIT:-$HOME/opt/petar-arm-bse/bin/petar.init}
MCLUSTER_BIN=${MCLUSTER_BIN:-$HOME/bin/mcluster}
PETAR_GETHER=${PETAR_GETHER:-$HOME/opt/petar-arm-bse/bin/petar.data.gether}
PETAR_PROCESS=${PETAR_PROCESS:-$HOME/opt/petar-arm-bse/bin/petar.data.process}
EXPECTED_CPU_TARGET=${EXPECTED_CPU_TARGET:-tsv110}

# BSE links PeTar to the GNU Fortran runtime.  Locate the matching GCC runtime
# explicitly because HPCKit may not expose it on compute nodes.
if [ -z "${GFORTRAN_LIBDIR:-}" ]; then
    GFORTRAN_LIB=""
    for candidate_dir in "$HOME/opt/petar-arm-bse/lib" /usr/lib64; do
        if [ -r "$candidate_dir/libgfortran.so.5" ]; then
            GFORTRAN_LIBDIR=$candidate_dir
            break
        fi
    done
fi
if [ -z "${GFORTRAN_LIBDIR:-}" ]; then
    GFORTRAN_LIB=""
    for compiler_driver in gfortran gcc; do
        if command -v "$compiler_driver" >/dev/null 2>&1; then
            candidate=$($compiler_driver -print-file-name=libgfortran.so.5)
            if [ "$candidate" != "libgfortran.so.5" ] && [ -r "$candidate" ]; then
                GFORTRAN_LIB=$candidate
                break
            fi
        fi
    done
    if [ -z "$GFORTRAN_LIB" ] && command -v ldconfig >/dev/null 2>&1; then
        GFORTRAN_LIB=$(ldconfig -p 2>/dev/null | awk '/libgfortran\.so\.5/{print $NF; exit}' || true)
    fi
    if [ -n "$GFORTRAN_LIB" ] && [ -r "$GFORTRAN_LIB" ]; then
        GFORTRAN_LIBDIR=$(dirname "$(readlink -f "$GFORTRAN_LIB")")
    fi
fi
if [ -z "${GFORTRAN_LIBDIR:-}" ] || [ ! -r "$GFORTRAN_LIBDIR/libgfortran.so.5" ]; then
    echo "ERROR: libgfortran.so.5 was not found." >&2
    echo "Install it in $HOME/opt/petar-arm-bse/lib or set GFORTRAN_LIBDIR explicitly." >&2
    exit 2
fi
export LD_LIBRARY_PATH="$GFORTRAN_LIBDIR${LD_LIBRARY_PATH:+:$LD_LIBRARY_PATH}"

# Positive custom instantaneous merger separation, in solar radii.
# 1e-6 matches the old bare --debug_lessmerger test behavior.
MERGER_RADIUS_RSUN=${MERGER_RADIUS_RSUN:-1e-6}

# Binary snapshots greatly reduce storage.  t=5000, o=10 gives about 501
# snapshots (about 331 by t=3300).  Allow roughly 10 GB including analysis.
END_TIME=10000
OUTPUT_INTERVAL=5
MPI_RANKS=32
OMP_THREADS=4

RUN_ROOT=${RUN_ROOT:-$HOME/petar_runs/singlestar_M210000_job${CCS_JOB_ID:-manual}}

for program in "$PETAR_BIN" "$PETAR_INIT" "$MCLUSTER_BIN" "$PETAR_GETHER" "$PETAR_PROCESS"; do
    if [ ! -x "$program" ]; then
        echo "ERROR: missing executable: $program" >&2
        exit 2
    fi
done
if ldd "$PETAR_BIN" 2>&1 | grep -q 'not found'; then
    echo "ERROR: unresolved shared libraries for $PETAR_BIN:" >&2
    ldd "$PETAR_BIN" >&2
    exit 2
fi

# -mcpu is a compile-time setting.  If the binary was built with
# -frecord-gcc-switches, verify that this job is using the intended build.
RECORDED_GCC_FLAGS=""
if command -v readelf >/dev/null 2>&1; then
    RECORDED_GCC_FLAGS=$(readelf --string-dump=.GCC.command.line "$PETAR_BIN" 2>/dev/null || true)
fi
if [[ "$RECORDED_GCC_FLAGS" == *"-mcpu=$EXPECTED_CPU_TARGET"* ]]; then
    echo "Verified PeTar CPU target: -mcpu=$EXPECTED_CPU_TARGET"
else
    echo "WARNING: cannot verify -mcpu=$EXPECTED_CPU_TARGET in $PETAR_BIN." >&2
    echo "Rebuild with CXXFLAGS='-O3 -mcpu=$EXPECTED_CPU_TARGET -ftree-vectorize -frecord-gcc-switches' to record and verify it." >&2
fi
if [ -z "${CCS_ALLOC_FILE:-}" ] || [ ! -r "$CCS_ALLOC_FILE" ]; then
    echo "ERROR: submit this script with dsub -s; CCS_ALLOC_FILE is missing." >&2
    exit 2
fi
if [ -e "$RUN_ROOT" ]; then
    echo "ERROR: run directory already exists: $RUN_ROOT" >&2
    exit 2
fi

mkdir -p "$RUN_ROOT"
cp "$0" "$RUN_ROOT/submitted_script.sh"
cd "$RUN_ROOT"

echo "Run directory: $RUN_ROOT"
echo "Generating McLuster initial conditions at $(date --iso-8601=seconds)"
export OMP_NUM_THREADS=8
taskset -c 0-31 "$MCLUSTER_BIN" \
  -M 210000 -m 0.7 -P 0 -R 1.621849 -f 0 -u 1 -C 5 \
  -B 0 -s 20261001 -o cluster > mc.log 2>&1

taskset -c 0 "$PETAR_INIT" -s bse -v kms2pcmyr \
  -f input cluster.dat.10 > petar_init.log 2>&1

NSTAR=$(awk 'NR==1 {print $2}' input)
echo "Initial condition ready: N=$NSTAR, M=210000 Msun, m=0.7 Msun, Rh=1.04381 pc"

export OMP_NUM_THREADS=$OMP_THREADS
echo "Starting PeTar with system scheduling: $MPI_RANKS MPI x $OMP_THREADS OpenMP, t=$END_TIME, o=$OUTPUT_INTERVAL"
echo "Custom merger radius: $MERGER_RADIUS_RSUN Rsun"
echo "Start time: $(date --iso-8601=seconds)"

mpirun \
  -x PATH -x LD_LIBRARY_PATH -x UCX_TLS=sm,rc \
  -x OMP_NUM_THREADS -x OMP_STACKSIZE -x OMP_DYNAMIC -x OMP_PROC_BIND \
  -x UCX_VFS_ENABLE \
  --bind-to none \
  --mca grpcomm_direct_priority 100 \
  --mca coll_tuned_use_dynamic_rules true \
  -x UCX_RNDV_THRESH=512K \
  -x BINDPROCNIC_NOSHOW=1 \
  -x UCX_RC_VERBS_ROCE_LOCAL_SUBNET=y \
  -x UCX_UD_VERBS_ROCE_LOCAL_SUBNET=y \
  --mca plm_rsh_agent /usr/bin/ssh \
  -np "$MPI_RANKS" \
  "$PETAR_BIN" \
    --debug_lessmerger "$MERGER_RADIUS_RSUN" \
    --energy-err-hard 1e-3 \
    -a 0 -u 1 -b 0 --stellar-evolution 2 \
    -i 2 -t "$END_TIME" -o "$OUTPUT_INTERVAL" input \
    > output.log 2>&1

echo "PeTar completed at $(date --iso-8601=seconds)"

"$PETAR_GETHER" data \
  > data_gether.log 2>&1

# Input snapshots are binary (-s binary); the executable was compiled with
# BSE fields (-i bse).  Binary processed snapshots avoid another large ASCII copy.
taskset -c 0-31 "$PETAR_PROCESS" \
  -i bse -s binary -o binary -G 0.00449830997959438 -n 32 \
  data.snap.lst > data_process.log 2>&1

echo "Post-processing completed at $(date --iso-8601=seconds)"
echo "Results: $RUN_ROOT"
