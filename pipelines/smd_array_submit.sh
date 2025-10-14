#!/usr/bin/env bash
set -euo pipefail

MAN=manifests/smd_manifest.csv
CFG=config.yaml

if [[ ! -f "$MAN" ]]; then
  echo "[smd-submit] ERROR: $MAN not found. Run: make smd-manifest" >&2
  exit 2
fi

readarray -t CFG_LINES < <(python3 - <<'PY'
import yaml
cfg = yaml.safe_load(open("config.yaml"))
slurm = cfg["globals"]["slurm"]
print(slurm["partition"])
print(slurm["gpus_per_node"])
print(slurm["cpus_per_task"])
print(slurm["time_limit"])
print(cfg["globals"]["smd"]["axis"])
print(slurm["gromacs_module"])
PY
)
PART="${CFG_LINES[0]}"
GPUS="${CFG_LINES[1]}"
CPUS="${CFG_LINES[2]}"
TLIM="${CFG_LINES[3]}"
AXIS="${CFG_LINES[4]}"
GMX_MOD="${CFG_LINES[5]}"

mkdir -p logs manifests/tmp

# systems present in the manifest
mapfile -t systems < <(awk -F, 'NR>1{print $1}' "$MAN" | sort -u)

for sys in "${systems[@]}"; do
  tmp="manifests/tmp/smd_manifest_${sys}.csv"
  # header
  head -n1 "$MAN" > "$tmp"
  # rows for this system
  awk -F, -v S="$sys" 'NR>1 && $1==S{print $0}' "$MAN" >> "$tmp"

  N=$(( $(wc -l < "$tmp") - 1 ))
  if [[ $N -le 0 ]]; then
    echo "[smd-submit] WARN: no rows for $sys, skipping."
    continue
  fi

  # Column 15 = array_cap (per our FIELDNAMES)
  ACAP=$(awk -F, 'NR==2{print $15}' "$tmp")
  [[ -z "$ACAP" ]] && ACAP=1

  echo "[smd-submit] Submitting $sys with $N tasks (cap %$ACAP)"

  sbatch <<SB
#!/bin/bash
#SBATCH --job-name=${sys}_smd
#SBATCH --partition=${PART}
#SBATCH --gpus-per-node=${GPUS}
#SBATCH --cpus-per-task=${CPUS}
#SBATCH --time=${TLIM}
#SBATCH --array=1-${N}%${ACAP}
#SBATCH --output=logs/${sys}_smd_%A_%a.out
#SBATCH --error=logs/${sys}_smd_%A_%a.err

module purge
module load ${GMX_MOD}

ROW=\$((SLURM_ARRAY_TASK_ID+1))
LINE=\$(sed -n "\${ROW}p" ${tmp})

python3 scripts/smd_runner.sh "\$LINE"
SB

done
