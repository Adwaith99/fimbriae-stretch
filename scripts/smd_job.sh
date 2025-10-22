#!/usr/bin/env bash
set -euo pipefail

# --- Slurm header intentionally minimal; flags come from sbatch CLI in the submitter ---
# (We still allow logs to go where sbatch says)

MANIFEST="${1:-manifests/smd_manifest.csv}"

# Load modules per config.yaml
readarray -t CFG_LINES < <(python3 - <<'PY'
import yaml
cfg = yaml.safe_load(open("config.yaml"))
slurm = cfg["globals"]["slurm"]
print(slurm.get("gromacs_module",""))
PY
)
GMX_MOD="${CFG_LINES[0]}"
if [[ -n "${GMX_MOD}" ]]; then
  module purge
  module load ${GMX_MOD}
fi

# Skip header (line 1)
if [[ "${SLURM_ARRAY_TASK_ID}" -eq 1 ]]; then
  echo "[smd-job] Task ${SLURM_ARRAY_TASK_ID} is header; skipping."
  exit 0
fi

# Read the CSV line whose index equals SLURM_ARRAY_TASK_ID
LINE="$(sed -n "${SLURM_ARRAY_TASK_ID}p" "${MANIFEST}")"
if [[ -z "${LINE}" ]]; then
  echo "[smd-job] No line for SLURM_ARRAY_TASK_ID=${SLURM_ARRAY_TASK_ID}" >&2
  exit 2
fi

# Compute mdrun -maxh from THIS job's time limit (leave headroom)
# Prefer SLURM_TIMELIMIT (minutes) or SLURM_TIMELIMIT_STR (D-HH:MM:SS)
MAXH_HOURS=""
if [[ -n "${SLURM_TIMELIMIT:-}" ]]; then
  # minutes → hours
  MAXH_HOURS=$(python3 - <<PY
import os, math
mins=float(os.environ["SLURM_TIMELIMIT"])
# keep 5% headroom, min 0.25h headroom
hours = mins/60.0
maxh = max(0.0, hours*0.95)
print(f"{maxh:.2f}")
PY
)
elif [[ -n "${SLURM_TIMELIMIT_STR:-}" ]]; then
  MAXH_HOURS=$(python3 - <<PY
import os
s=os.environ["SLURM_TIMELIMIT_STR"]
# formats: "D-HH:MM:SS" or "HH:MM:SS"
parts=s.split('-')
if len(parts)==2:
  d=int(parts[0]); hh,mm,ss=map(int,parts[1].split(':'))
  hours=d*24+hh+mm/60+ss/3600
else:
  hh,mm,ss=map(int,parts[0].split(':'))
  hours=hh+mm/60+ss/3600
maxh=max(0.0, hours*0.95)
print(f"{maxh:.2f}")
PY
)
fi
export MAXH_HOURS

echo "[smd-job] Node=$(hostname)  Task=${SLURM_ARRAY_TASK_ID}  maxh=${MAXH_HOURS:-unset}"

# Run the per-row driver
bash scripts/smd_runner.sh "${LINE}"
