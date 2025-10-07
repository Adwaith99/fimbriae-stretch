#!/usr/bin/env bash
#SBATCH --job-name=pulls
#SBATCH --output=logs/pulls_%A_%a.out
#SBATCH --partition=compute_full_node
#SBATCH --nodes=1
#SBATCH --gpus-per-node=4
#SBATCH --cpus-per-task=96
#SBATCH --time=24:00:00

set -euo pipefail

# (Optional) load your standard modules here if needed for the runner
# module purge; module load arch/avx512 StdEnv/2023 gcc/12.3 openmpi/4.1.5 cuda/12.2 gromacs/2024.4

ROOT="$(cd "$(dirname "$0")/.."; pwd)"
echo "[pulls_array_submit] SLURM_ARRAY_TASK_ID=${SLURM_ARRAY_TASK_ID}"
echo "[pulls_array_submit] FIM_PULL_MANIFEST=${FIM_PULL_MANIFEST:-<unset>}"

/usr/bin/env bash "$ROOT/scripts/runner.sh"
