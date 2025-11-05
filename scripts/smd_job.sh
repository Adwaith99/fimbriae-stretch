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

############################
# Auto-restart logic: check completion and resubmit if needed
############################

# Parse manifest line to get run directory path
IFS=',' read -r system variant replicate speed_nm_per_ns k_kj dt_ps target_extension_nm axis perf_ns_per_day start_time_ps final_tpr final_xtc start_id anchor_chain array_cap <<< "${LINE}"

# Strip quotes
stripq() { local v="${1}"; v="${v#\"}"; v="${v%\"}"; printf "%s" "${v}"; }
system="$(stripq "${system}")"
variant="$(stripq "${variant}")"
replicate="$(stripq "${replicate}")"
speed_nm_per_ns="$(stripq "${speed_nm_per_ns}")"
start_id="$(stripq "${start_id}")"

RUN_DIR="${ROOT}/smd/${system}/${variant}/v$(printf "%.3f" "${speed_nm_per_ns}")/rep${replicate}/${start_id}"

# Check if simulation completed
COMPLETE=$(python3 - <<'PY'
import re, os, sys

run_dir = os.getenv("RUN_DIR")
log_file = f"{run_dir}/pull.log"
expected_file = f"{run_dir}/expected_nsteps.txt"

if not os.path.exists(log_file) or not os.path.exists(expected_file):
    print("0")
    sys.exit(0)

expected_nsteps = int(open(expected_file).read().strip())

# Parse pull.log for max step reached
max_step = 0
with open(log_file) as f:
    for line in f:
        # Match "Step 12345" or "Statistics over 12345 steps"
        m = re.search(r'(?:^Step\s+|Statistics over\s+)(\d+)', line)
        if m:
            max_step = max(max_step, int(m.group(1)))

# Complete if we've done >= expected steps
complete = max_step >= expected_nsteps
print("1" if complete else "0")
PY
)

RESTART_COUNT="${RESTART_COUNT:-0}"
MAX_RESTARTS=50

if [[ "${COMPLETE}" == "1" ]]; then
  echo "[smd-job] Simulation COMPLETE (step count reached expected)"
elif [[ "${RESTART_COUNT}" -ge "${MAX_RESTARTS}" ]]; then
  echo "[smd-job] WARNING: Max restarts (${MAX_RESTARTS}) reached; not resubmitting" >&2
else
  echo "[smd-job] Simulation INCOMPLETE; resubmitting (restart ${RESTART_COUNT}/${MAX_RESTARTS})"
  sleep 10
  
  # Build sbatch command with same parameters but as single job (no array)
  NEXT_RESTART=$((RESTART_COUNT + 1))
  
  # Inherit current job's parameters
  SBATCH_ARGS=(
    --chdir="${ROOT}"
    --job-name="${SLURM_JOB_NAME:-smd-restart}"
    --partition="${SLURM_JOB_PARTITION:-default}"
    --time="${SLURM_JOB_TIME:-24:00:00}"
    --output="${ROOT}/logs/restart_${system}_${variant}_v${speed_nm_per_ns}_rep${replicate}_${start_id}_%j.log"
    --error="${ROOT}/logs/restart_${system}_${variant}_v${speed_nm_per_ns}_rep${replicate}_${start_id}_%j.log"
  )
  
  # Add CPU or GPU specific args
  if [[ -n "${GMX_CMD:-}" && "${GMX_CMD}" == *"gmx_mpi"* ]]; then
    # CPU mode
    SBATCH_ARGS+=(
      --nodes="${SLURM_JOB_NUM_NODES:-1}"
      --ntasks-per-node="${SLURM_NTASKS_PER_NODE:-1}"
      --cpus-per-task="${SLURM_CPUS_PER_TASK:-8}"
    )
  else
    # GPU mode
    [[ -n "${SLURM_JOB_GPUS:-}" ]] && SBATCH_ARGS+=(--gres="${SLURM_JOB_GPUS}")
    SBATCH_ARGS+=(--cpus-per-task="${SLURM_CPUS_PER_TASK:-8}")
  fi
  
  # Export necessary variables
  EXPORTS="ALL,RESTART_COUNT=${NEXT_RESTART},MAXH_HOURS=${MAXH_HOURS:-}"
  [[ -n "${GMX_CMD:-}" ]] && EXPORTS="${EXPORTS},GMX_CMD=${GMX_CMD}"
  [[ -n "${GMX_MODULE:-}" ]] && EXPORTS="${EXPORTS},GMX_MODULE=${GMX_MODULE}"
  
  # Create a temporary single-line manifest for this specific job
  TEMP_MANIFEST="${ROOT}/manifests/tmp/restart_${system}_${variant}_v${speed_nm_per_ns}_rep${replicate}_${start_id}.csv"
  mkdir -p "$(dirname "${TEMP_MANIFEST}")"
  # Write header + this line (as line 2 so SLURM_ARRAY_TASK_ID=2 works, but we'll use sed -n 2p)
  head -n1 "${MANIFEST}" > "${TEMP_MANIFEST}"
  echo "${LINE}" >> "${TEMP_MANIFEST}"
  
  # Submit as array with single task (index=2, the data line)
  JID=$(sbatch --parsable --array=2 "${SBATCH_ARGS[@]}" --export="${EXPORTS}" "${ROOT}/scripts/smd_job.sh" "${TEMP_MANIFEST}")
  echo "[smd-job] Resubmitted as job ${JID}"
fi
