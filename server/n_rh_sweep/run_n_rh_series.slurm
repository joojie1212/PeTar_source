#!/usr/bin/env bash
#SBATCH --job-name=petar_n_rh
#SBATCH --partition=tiny
#SBATCH --nodes=1
#SBATCH --ntasks=4
#SBATCH --cpus-per-task=4
#SBATCH --hint=nomultithread
#SBATCH --mem=16G
#SBATCH --time=5-00:00:00
#SBATCH --output=petar_n_rh_%j.log
#SBATCH --error=petar_n_rh_%j.err
set -Eeuo pipefail

# sbatch spools this file, so BASH_SOURCE is NOT the original bundle directory.
# Submit from the directory containing this file, or export SWEEP_SCRIPT_DIR.
SCRIPT_DIR=${SWEEP_SCRIPT_DIR:-${SLURM_SUBMIT_DIR:-$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)}}
if [[ ${1:-} == --print-config ]]; then
    exec python3 "$SCRIPT_DIR/sweep.py" --print-config
fi
[[ ${1:-} == --run && $# -eq 1 ]] || {
    echo 'Preview: bash run_n_rh_series.slurm --print-config'
    echo 'Run:     sbatch run_n_rh_series.slurm --run'
    exit 2
}
: "${SLURM_JOB_ID:?Use sbatch, not bash, to run simulations}"
[[ ${SLURM_JOB_NUM_NODES:-0} == 1 && ${SLURM_NTASKS:-0} == 4 && ${SLURM_CPUS_PER_TASK:-0} == 4 ]] || {
    echo 'This launcher requires one node, 4 tasks, 4 CPUs/task.' >&2; exit 2;
}
# Some site module/Conda hooks reference unset variables.
set +u
source "${PETAR_ENV_FILE:-$HOME/petar_env.sh}"
set -u
conda activate petar-runtime
export OMP_NUM_THREADS=4 OMP_STACKSIZE=128M OMP_PLACES=cores OMP_PROC_BIND=spread
export OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1
export OUTPUT_ROOT=${OUTPUT_ROOT:-$HOME/scratch/petar_n_rh_results}
for cmd in srun mpirun mcluster petar petar.init petar.find.dt petar.format.transfer python3; do
    command -v "$cmd" >/dev/null || { echo "Missing command: $cmd" >&2; exit 2; }
done
[[ -x /usr/bin/time ]] || { echo 'GNU /usr/bin/time is required.' >&2; exit 2; }
python3 - <<'PY'
import numpy, scipy, precession
print('precession:', precession.__file__)
PY
exec python3 "$SCRIPT_DIR/sweep.py" --run
