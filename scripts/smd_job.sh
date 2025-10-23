#!/usr/bin/env bash
set -euo pipefail

MANIFEST="${1:?Usage: smd_job.sh <ABS_PATH_TO_MANIFEST>}"

# We are started in repo root because submitter used --chdir=<repo_root>
ROOT="$(pwd)"
echo "[smd-job] CWD=${ROOT}"
echo "[smd-job] MANIFEST=${MANIFEST}"

# Load modules per config.yaml
readarray -t CFG_LINES < <(python3 - <<'PY'
import yaml
cfg = yaml.safe_load(open("config.yaml"))
print(cfg["globals"]["slurm"].get("gromacs_module",""))
PY
)
GMX_MOD="${CFG_LINES[0]}"
if [[ -n "${GMX_MOD}" ]]; then
  module purge || true
  # --- HPC arch + GPU-direct comms ---------------------------------
  # Load AVX-512 tuning module (safe to ignore if unavailable)
  module load arch/avx512 || true
  module load ${GMX_MOD}
fi

# Enable FPU↔GPU direct communications in GROMACS (CUDA-aware MPI path)
export GMX_ENABLE_DIRECT_GPU_COMM=ON

# Log what we actually have
echo "[smd-job] GMX_ENABLE_DIRECT_GPU_COMM=${GMX_ENABLE_DIRECT_GPU_COMM}"
module list 2>&1 | sed 's/^/[smd-job] module: /'


# Skip header line if array index == 1
if [[ "${SLURM_ARRAY_TASK_ID}" -eq 1 ]]; then
  echo "[smd-job] Task ${SLURM_ARRAY_TASK_ID} is header; skipping."
  exit 0
fi

# Read the exact CSV line == SLURM_ARRAY_TASK_ID
LINE="$(sed -n "${SLURM_ARRAY_TASK_ID}p" "${MANIFEST}")"
if [[ -z "${LINE}" ]]; then
  echo "[smd-job] No line for SLURM_ARRAY_TASK_ID=${SLURM_ARRAY_TASK_ID}" >&2
  exit 2
fi

# Prefer MAXH_HOURS exported by submitter; fall back to SLURM if present
if [[ -n "${MAXH_HOURS:-}" ]]; then
  echo "[smd-job] maxh=${MAXH_HOURS} (from submitter)"
elif [[ -n "${SLURM_TIMELIMIT:-}" ]]; then
  export MAXH_HOURS="$(python3 - <<PY
import os
mins=float(os.environ["SLURM_TIMELIMIT"])
print(f"{(mins/60.0)*0.95:.2f}")
PY
)"
  echo "[smd-job] maxh=${MAXH_HOURS} (from SLURM_TIMELIMIT)"
elif [[ -n "${SLURM_TIMELIMIT_STR:-}" ]]; then
  export MAXH_HOURS="$(python3 - <<PY
import os
s=os.environ["SLURM_TIMELIMIT_STR"]
parts=s.split('-')
if len(parts)==2:
  d=int(parts[0]); hh,mm,ss=map(int,parts[1].split(':'))
  hours=d*24+hh+mm/60+ss/3600
else:
  hh,mm,ss=map(int,parts[0].split(':'))
  hours=hh+mm/60+ss/3600
print(f"{hours*0.95:.2f}")
PY
)"
  echo "[smd-job] maxh=${MAXH_HOURS} (from SLURM_TIMELIMIT_STR)"
else
  echo "[smd-job] maxh=unset (no export and no SLURM vars)"
fi


# Launch the per-row runner
bash "scripts/smd_runner.sh" "${LINE}"
