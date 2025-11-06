#!/usr/bin/env bash
set -euo pipefail

# Submit SMD runs that are NOT completed on disk.
# Default behavior (USE_LEDGER=0): Only skips runs with "Finished mdrun" in pull.log
# - Will resubmit incomplete/failed runs (even if in ledger)
# - Will resubmit queued/running jobs (user must cancel manually first)
# - Ledger is updated but not used to block submissions
#
# Alternative (USE_LEDGER=1): Also skip runs in smd_submitted.csv
# - Prevents duplicate submissions if you trust the ledger
# - export USE_LEDGER=1 before running to enable
#
# Typical workflow:
# 1. scancel <jobid>  # Cancel unwanted GPU jobs
# 2. make smd-submit-new-cpu  # Resubmit incomplete runs to CPU

MANIFEST="${1:-manifests/smd_manifest.csv}"
LEDGER="manifests/smd_submitted.csv"      # Tracks what was submitted (informational only)
SAFETY="${SMD_TIME_SAFETY_FACTOR:-1.20}"  # walltime pad; override via env
JOB_SCRIPT="${SMD_JOB_SCRIPT:-scripts/smd_job.sh}"
USE_LEDGER="${USE_LEDGER:-0}"             # Set USE_LEDGER=1 to prevent resubmission of ledger entries

[[ -f "${MANIFEST}" ]] || { echo "[smd-submit-new] ERROR: manifest not found: ${MANIFEST}" >&2; exit 2; }
[[ -f "${JOB_SCRIPT}" ]] || { echo "[smd-submit-new] ERROR: job script not found: ${JOB_SCRIPT}" >&2; exit 2; }

# Read Slurm defaults from config (partition, gpus, cpus, nodes, ntasks_per_node, max_wall_hours)
readarray -t CFG < <(python3 - <<'PY'
import yaml
cfg=yaml.safe_load(open("config.yaml"))
g=cfg.get("globals",{})
s=g.get("slurm",{})
print(s.get("partition", "default"))
print(s.get("gpus_per_node", ""))
print(s.get("cpus_per_task", 8))
print(s.get("nodes", 1))
print(s.get("ntasks_per_node", 1))
print(s.get("gromacs_module", ""))
print(g.get("max_wall_hours", 168))
PY
)
PART="${CFG[0]}"
GRES_SPEC="${CFG[1]}"
CPUS="${CFG[2]}"
NODES="${CFG[3]}"
NTASKS_PER_NODE="${CFG[4]}"
GMX_MOD="${CFG[5]}"
MAX_WALL_HOURS="${CFG[6]}"

# CPU mode toggle
CPU_MODE="${CPU:-0}"

mkdir -p logs manifests/tmp
touch "${LEDGER}"

# Build a ledger set (only used if USE_LEDGER=1)
declare -A submitted
if [[ "${USE_LEDGER}" = "1" ]]; then
  while IFS=, read -r s v r sp sid || [[ -n "${s:-}" ]]; do
    [[ -z "${s:-}" ]] && continue
    submitted["$s,$v,$r,$sp,$sid"]=1
  done < "${LEDGER}"
  echo "[smd-submit-new] Ledger blocking enabled (USE_LEDGER=1): will skip ${#submitted[@]} previously submitted runs."
else
  echo "[smd-submit-new] Ledger blocking disabled (default): only completed runs will be skipped."
fi

# Parse manifest rows to TSV (TAB-separated!)
tmpdir=$(mktemp -d); trap 'rm -rf "${tmpdir}"' EXIT
rows="${tmpdir}/rows.tsv"
# lineno  system  variant  replicate  speed  target  start_id  array_cap
# Note: we no longer extract perf from manifest (column 9); we'll look it up from config per system
awk -F, -v OFS='\t' 'NR>1{
  gsub(/^"|"$/, "", $0);
  print NR, $1, $2, $3, $4, $7, $13, $15
}' "${MANIFEST}" > "${rows}"

############################
# Build a set of queued array indices from SLURM (avoid double-queuing)
############################
declare -A queued_indices
if command -v squeue >/dev/null 2>&1; then
  # Only check jobs for current user and this manifest
  squeue -u "$USER" --name smd: --format="%j %A %a" | awk '{print $1,$3}' | while read jname arrayidx; do
    # Extract array index from job name and array index
    # Job name format: smd:<SYS>:v<SPD>
    # Array index: <arrayidx>
    [[ "$jname" =~ ^smd: ]] || continue
    [[ -n "$arrayidx" ]] || continue
    queued_indices["$arrayidx"]=1
  done
fi
declare -A system_perf
readarray -t PERF_MAP < <(python3 - <<'PY'
import yaml
cfg = yaml.safe_load(open("config.yaml"))
global_perf = float(cfg.get("globals", {}).get("perf_ns_per_day", 10.0))
for sys in cfg.get("systems", []):
    name = sys.get("name")
    perf = float(sys.get("perf_ns_per_day", global_perf))
    print(f"{name}:{perf}")
PY
)
for entry in "${PERF_MAP[@]}"; do
  IFS=':' read -r sname sperf <<< "${entry}"
  system_perf["${sname}"]="${sperf}"
done

# Group NEW rows by (system, speed)
declare -A grp_idxs grp_time grp_cap
BAD=0
while IFS=$'\t' read -r L SYS VAR REP SPD TGT START CAP; do
  # Skip empty lines defensively
  [[ -z "${L}" ]] && continue

  KEY="$SYS,$VAR,$REP,$SPD,$START"

  # Skip if in ledger AND ledger blocking is enabled
  if [[ "${USE_LEDGER}" = "1" && -n "${submitted[$KEY]:-}" ]]; then
    continue
  fi

  # Lookup perf from config for this system
  PERF="${system_perf[${SYS}]:-10.0}"

  # Validate required numeric fields
  if [[ -z "${SPD}" || -z "${TGT}" ]]; then
    echo "[smd-submit-new] WARN: skipping line ${L} (missing numeric: speed='${SPD}', target='${TGT}')" >&2
    BAD=$((BAD+1))
    continue
  fi

  # Numeric sanity
  if ! python3 - <<PY >/dev/null 2>&1
spd=float("${SPD}"); tgt=float("${TGT}"); perf=float("${PERF}")
assert spd>0 and tgt>0 and perf>0
PY
  then
    echo "[smd-submit-new] WARN: skipping line ${L} (invalid numbers: speed=${SPD}, target=${TGT}, perf=${PERF})" >&2
    BAD=$((BAD+1))
    continue
  fi

  # Check if run is truly DONE on disk
  RUN="smd/${SYS}/${VAR}/v$(printf "%0.3f" "${SPD}")/rep${REP}/${START}"
  IS_DONE=0
  if [[ -f "${RUN}/pull.xtc" && -f "${RUN}/pull.log" ]]; then
    if grep -q "Finished mdrun" "${RUN}/pull.log" 2>/dev/null; then
      IS_DONE=1
    fi
  fi

  # Check if this array index is already queued
  if [[ -n "${queued_indices[${L}]:-}" ]]; then
    echo "[smd-submit-new] SKIP queued: ${RUN} (array index ${L})" >&2
    continue
  fi

  if [[ "${IS_DONE}" -eq 1 ]]; then
    echo "[smd-submit-new] SKIP completed: ${RUN}" >&2
    continue
  fi

  # Determine restart status
  RESTART_MSG="start from scratch"
  if [[ -f "${RUN}/pull.cpt" ]]; then
    # Try to extract last completed step from pull.log
    LAST_STEP=""
    if [[ -f "${RUN}/pull.log" ]]; then
      LAST_STEP=$(awk '/Step/ {s=$2} END{print s}' "${RUN}/pull.log")
    fi
    if [[ -n "$LAST_STEP" ]]; then
      RESTART_MSG="restart from checkpoint (step $LAST_STEP)"
    else
      RESTART_MSG="restart from checkpoint (step unknown)"
    fi
  fi

  # Not completed = submit (even if in ledger or has partial files)
  if [[ -f "${RUN}/pull.log" ]]; then
    echo "[smd-submit-new] RESUBMIT incomplete/failed: ${RUN} -- ${RESTART_MSG}" >&2
  elif [[ -n "${submitted[$KEY]:-}" ]]; then
    echo "[smd-submit-new] RESUBMIT (was in ledger but not done): ${RUN} -- ${RESTART_MSG}" >&2
  else
    echo "[smd-submit-new] SUBMIT new: ${RUN} -- ${RESTART_MSG}" >&2
  fi

  GRP="${SYS}|${SPD}"
  # Collect sparse array indices (CSV line numbers)
  if [[ -n "${grp_idxs[$GRP]:-}" ]]; then
    grp_idxs[$GRP]+=","${L}
  else
    grp_idxs[$GRP]="${L}"
  fi

  # walltime hours per row = (target/speed)/perf * 24 * SAFETY, capped at MAX_WALL_HOURS
  ht=$(python3 - <<PY
import math
spd=float("${SPD}"); tgt=float("${TGT}"); perf=float("${PERF}"); pad=float("${SAFETY}")
max_hrs=float("${MAX_WALL_HOURS}")
ns=tgt/spd; days=(ns/perf)
h=days*24*pad
h=max(1.0, min(h, max_hrs))
print("{:.3f}".format(h))
PY
)
  cur="${grp_time[$GRP]:-0.0}"
  grp_time[$GRP]=$(python3 - <<PY
a=float("${cur}"); b=float("${ht}")
print("{:.3f}".format(max(a,b)))
PY
)
  grp_cap[$GRP]="${CAP:-5}"
done < "${rows}"

if [[ "${BAD}" -gt 0 ]]; then
  echo "[smd-submit-new] NOTE: skipped ${BAD} invalid manifest row(s)." >&2
fi

any=0
for GRP in "${!grp_idxs[@]}"; do
  IFS='|' read -r SYS SPD <<< "${GRP}"
  IDXS="${grp_idxs[$GRP]}"
  CAP="${grp_cap[$GRP]:-5}"
  # Convert hours → HH:MM:SS (ceil to minutes)
  TSTR=$(python3 - <<PY
import math
h=float("${grp_time[$GRP]}")
m=math.ceil(h*60)
print(f"{m//60:02d}:{m%60:02d}:00")
PY
)
  # ---- Resolve absolute paths & chdir to repo root ----
  MANIFEST_ABS="$(readlink -f "${MANIFEST}")"
  ROOT="$(dirname "${MANIFEST_ABS}")/.."
  ROOT="$(readlink -f "${ROOT}")"
  JOB_SCRIPT_ABS="$(readlink -f "${JOB_SCRIPT}")"
  mkdir -p "${ROOT}/logs"

  # Compute maxh directly from the same hours used to build TSTR (keep 5% headroom)
  MAXH=$(python3 - <<PY
import math
h=float("${grp_time[$GRP]}")   # hours
print(f"{h*0.95:.2f}")
PY
)

  JNAME="smd:${SYS}:v$(printf "%0.3f" "${SPD}")"
  OUT="${ROOT}/logs/${SYS}_v$(printf "%0.3f" "${SPD}")_%A_%a.log"   # unified stdout+stderr

  # Build sbatch args based on CPU/GPU mode
  SBATCH_ARGS=(
    --chdir="${ROOT}"
    --job-name="${JNAME}"
    --partition="${PART}"
    --cpus-per-task="${CPUS}"
    --time="${TSTR}"
    --array="${IDXS}%${CAP}"
    --output="${OUT}"
    --error="${OUT}"
  )
  EXPORTS="ALL,MAXH_HOURS=${MAXH},GMX_MODULE=${GMX_MOD}"

  if [[ "${CPU_MODE}" = "1" ]]; then
    # CPU mode
    SBATCH_ARGS+=(--nodes="${NODES}" --ntasks-per-node="${NTASKS_PER_NODE}")
    EXPORTS="${EXPORTS},GMX_CMD=srun gmx_mpi mdrun"
    echo "[smd-submit-new] CPU submit: ${JNAME} rows=${IDXS} time=${TSTR} cap=${CAP} nodes=${NODES} ntasks/node=${NTASKS_PER_NODE}"
  else
    # GPU mode
    [[ -n "${GRES_SPEC}" ]] && SBATCH_ARGS+=(--gres="gpu:${GRES_SPEC}")
    echo "[smd-submit-new] GPU submit: ${JNAME} rows=${IDXS} time=${TSTR} cap=${CAP} gres=${GRES_SPEC}"
  fi

  if [[ -n "${SLURM_TEST_ONLY:-}" ]]; then
    echo "[smd-submit-new] DRY-RUN: sbatch --test-only (no submission)"
    any=1
    continue
  fi

  JID=$(
    sbatch --parsable \
      "${SBATCH_ARGS[@]}" \
      --export="${EXPORTS}" \
      "${JOB_SCRIPT_ABS}" "${MANIFEST_ABS}"
  )
  echo "[smd-submit-new] submitted: ${JID}"



  # Append submitted keys for this group to the ledger (real submissions only)
  pat="^($(echo "${IDXS}" | sed 's/,/|/g'))$"
  awk -F $'\t' -v pat="${pat}" ' $1 ~ pat { printf "%s,%s,%s,%s,%s\n",$2,$3,$4,$5,$8 } ' "${rows}" >> "${LEDGER}"

  any=1
done

if [[ "${any}" -eq 0 ]]; then
  echo "[smd-submit-new] Nothing new to submit (ledger up to date or all invalid/finished)."
fi
