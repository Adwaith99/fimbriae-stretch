#!/usr/bin/env bash
set -euo pipefail
ROOT="$(cd "$(dirname "$0")/.."; pwd)"
mkdir -p "$ROOT/logs"
python3 "$ROOT/scripts/generate_manifest.py"

CAP=$(python3 - <<PY
import yaml; print(yaml.safe_load(open("$ROOT/config.yaml"))["globals"]["slurm"]["array_cap"])
PY
)
PARTITION=$(python3 - <<PY
import yaml; print(yaml.safe_load(open("$ROOT/config.yaml"))["globals"]["slurm"]["partition"])
PY
)
GMXMOD=$(python3 - <<PY
import yaml; print(yaml.safe_load(open("$ROOT/config.yaml"))["globals"]["slurm"]["gromacs_module"])
PY
)
TLIM=$(python3 - <<PY
import yaml; print(yaml.safe_load(open("$ROOT/config.yaml"))["globals"]["slurm"]["time_limit"])
PY
)

# Only submit pulls for systems with BUILD_DONE
awk -F, 'NR>1 && $6=="PENDING" {print NR-1}' "$ROOT/manifests/manifest.csv" > /tmp/pending.idx || true
if [[ ! -s /tmp/pending.idx ]]; then
  echo "No PENDING pull rows."
  exit 0
fi

FIRST=$(head -n1 /tmp/pending.idx)
LAST=$(tail -n1 /tmp/pending.idx)

sbatch <<EOT
#!/bin/bash
#SBATCH --job-name=fimA_pulls
#SBATCH --partition=${PARTITION}
#SBATCH --gpus-per-node=4
#SBATCH --cpus-per-task=96
#SBATCH --time=${TLIM}
#SBATCH --array=${FIRST}-${LAST}%${CAP}
#SBATCH --requeue
#SBATCH --output=$ROOT/logs/pulls_%A_%a.out
module purge
module load ${GMXMOD}
bash "$ROOT/scripts/runner.sh"
EOT
echo "Submitted pull array indices ${FIRST}-${LAST} (cap ${CAP})."
