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
# Build a set of queued array tasks from SLURM (avoid double-queuing)
# Key = "<JOB_NAME>|<ARRAY_INDEX>" (e.g., "smd:fimA_WT:v0.020|17")
############################
declare -A queued              # keyed by JNAME|index when full job name present
declare -A queued_index        # keyed by plain array index (also used for 'sq' PD parsing)

# Prefer site alias 'sq' if available: parse PD (pending) array ranges from JOBID column
if command -v sq >/dev/null 2>&1; then
  [[ -n "${SMD_DEBUG:-}" ]] && echo "[smd-submit-new][dbg] using 'sq' to detect pending array indices" >&2
  # Capture sq output and feed to Python to expand ranges
  sq_output=$(sq 2>&1)
  if [[ -n "${SMD_DEBUG:-}" ]]; then
    echo "[smd-submit-new][dbg] sq returned $(echo "$sq_output" | wc -l) lines" >&2
  fi
  while read -r idx; do
    [[ -z "$idx" ]] && continue
    queued_index["$idx"]=1
    [[ -n "${SMD_DEBUG:-}" ]] && echo "[smd-submit-new][dbg] sq pending index: $idx" >&2
  done < <(echo "$sq_output" | python3 - <<'PY'
import sys,re,os
lines=sys.stdin.read().splitlines()
debug=os.environ.get('SMD_DEBUG','')
if debug:
    print(f"[sq-parse] Python received {len(lines)} lines", file=sys.stderr)
for i,line in enumerate(lines):
    if i==0 and 'JOBID' in line:
        continue
    if not line.strip():
        continue
    parts=re.split(r'\s+', line.strip())
    # Expect columns: JOBID USER ACCOUNT NAME ST ...
    if len(parts) < 5:
        if debug:
            print(f"[sq-parse] skipping line {i} (only {len(parts)} parts): {line[:60]}", file=sys.stderr)
        continue
    jobid,name,state = parts[0], parts[3], parts[4]
    if debug:
        print(f"[sq-parse] line {i}: jobid={jobid} name={name} state={state}", file=sys.stderr)
    # Match any smd: job (name may be truncated)
    if not name.startswith('smd:'):
        continue
    if state != 'PD':
        continue
    m=re.search(r'_\[([^\]]+)', jobid)  # allow truncated closing bracket
    if not m:
        if debug:
            print(f"[sq-parse] no bracket range in jobid: {jobid}", file=sys.stderr)
        continue
    payload=m.group(1)
    if debug:
        print(f"[sq-parse] extracted payload from {jobid}: {payload}", file=sys.stderr)
    for tok in payload.split(','):
        core=tok.split('%',1)[0]
        if re.match(r'^\d+-\d+$', core):
            a,b=core.split('-')
            a=int(a); b=int(b)
            if a<=b:
                for k in range(a,b+1):
                    print(k)
        elif re.match(r'^\d+$', core):
            print(int(core))
PY
  )
fi

# Also query squeue for exact task indices and names
if command -v squeue >/dev/null 2>&1; then
  u="${USER:-$(id -un 2>/dev/null)}"
  # Try per-task first; some clusters don't show %i, so we fall back to parsing array ranges from %A.
  got_any=0
  while read -r jname arrayidx state; do
    [[ "$jname" =~ ^smd: ]] || continue
    [[ "$arrayidx" =~ ^[0-9]+$ ]] || continue
    queued["$jname|$arrayidx"]=1; got_any=1
    queued_index["$arrayidx"]=1
    [[ -n "${SMD_DEBUG:-}" ]] && echo "[smd-submit-new][dbg] queued task: $jname|$arrayidx (state=$state)" >&2
  done < <(squeue -h -u "$u" -o "%j %i %T" 2>/dev/null)

  if [[ "${got_any}" -eq 0 ]]; then
    # Fallback: parse job array ranges from the JOBID field like 436236_[17-21%5,30,33-35]
    while read -r jobid jname; do
      [[ "$jname" =~ ^smd: ]] || continue
      # Extract content inside brackets []
      rng=$(python3 - "$jobid" <<'PY'
import re,sys
jobid=sys.argv[1]
m=re.search(r'_\[(.+?)\]$', jobid)
print(m.group(1) if m else "")
PY
)
      [[ -z "$rng" ]] && continue
      # Split by comma and expand numeric ranges; strip %limits
      IFS=',' read -ra parts <<< "$rng"
      for p in "${parts[@]}"; do
        p_no=${p%%\%*}
        if [[ "$p_no" =~ ^[0-9]+-[0-9]+$ ]]; then
          IFS='-' read -r a b <<< "$p_no"
          for ((k=a; k<=b; k++)); do
            queued["$jname|$k"]=1
            queued_index["$k"]=1
            [[ -n "${SMD_DEBUG:-}" ]] && echo "[smd-submit-new][dbg] queued range: $jname|$k" >&2
          done
        elif [[ "$p_no" =~ ^[0-9]+$ ]]; then
          queued["$jname|$p_no"]=1
          queued_index["$p_no"]=1
          [[ -n "${SMD_DEBUG:-}" ]] && echo "[smd-submit-new][dbg] queued single: $jname|$p_no" >&2
        fi
      done
    done < <(squeue -h -u "$u" -o "%A %j" 2>/dev/null)
  fi
fi
if [[ -n "${SMD_DEBUG:-}" ]]; then
  # Print a compact sorted list of pending/queued indices detected
  if [[ "${!queued_index[@]}" ]]; then
    mapfile -t _qtmp < <(for k in "${!queued_index[@]}"; do echo "$k"; done | sort -n)
    echo "[smd-submit-new][dbg] final queued indices: ${_qtmp[*]}" >&2
  else
    echo "[smd-submit-new][dbg] no queued indices detected" >&2
  fi
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

  # Compute job name and run directory
  JNAME="smd:${SYS}:v$(printf "%0.3f" "${SPD}")"
  RUN="smd/${SYS}/${VAR}/v$(printf "%0.3f" "${SPD}")/rep${REP}/${START}"
  [[ -n "${SMD_DEBUG:-}" ]] && echo "[smd-submit-new][dbg] row=${L} JNAME=${JNAME} RUN=${RUN}" >&2
  # Skip if queued (prefer exact JNAME|idx, fall back to index-only)
  if [[ -n "${queued[${JNAME}|${L}]:-}" || -n "${queued_index[${L}]:-}" ]]; then
    echo "[smd-submit-new] SKIP queued: ${RUN} (job ${JNAME} idx ${L})" >&2
    continue
  fi
  
  # Step-based completion only
  IS_DONE=0
  LAST_STEP=""
  EXPECTED_STEPS=""
  if [[ -f "${RUN}/pull.log" ]]; then
    LAST_STEP=$(LOG_PATH="${RUN}/pull.log" python3 - <<'PY'
import os,re
log=os.environ.get('LOG_PATH','')
max_step=0
try:
  with open(log) as f:
    for line in f:
      m=re.search(r'(?:^Step\s+|Statistics over\s+)(\d+)', line)
      if m:
        v=int(m.group(1))
        if v>max_step:
          max_step=v
except Exception:
  pass
print(max_step if max_step>0 else "")
PY
)
  fi
  if [[ -f "${RUN}/expected_nsteps.txt" ]]; then
    EXPECTED_STEPS=$(cat "${RUN}/expected_nsteps.txt")
  fi
  # Fallback: parse nsteps from pull.mdp if expected file is missing
  if [[ -z "${EXPECTED_STEPS}" && -f "${RUN}/pull.mdp" ]]; then
    EXPECTED_STEPS=$(awk 'tolower($1)=="nsteps"{print $3}' "${RUN}/pull.mdp" | tail -n1)
  fi
  # Only do arithmetic if both values are integers
  if [[ "${LAST_STEP}" =~ ^[0-9]+$ && "${EXPECTED_STEPS}" =~ ^[0-9]+$ ]]; then
    if (( LAST_STEP >= EXPECTED_STEPS )); then
      IS_DONE=1
    fi
  fi
  # No longer rely on 'Finished mdrun' to mark completion (handles maxh terminations correctly)
  [[ -n "${SMD_DEBUG:-}" ]] && echo "[smd-submit-new][dbg] steps last=${LAST_STEP} expected=${EXPECTED_STEPS} done=${IS_DONE}" >&2

  if [[ "${IS_DONE}" -eq 1 ]]; then
    echo "[smd-submit-new] SKIP completed: ${RUN}" >&2
    continue
  fi

  # Determine restart status
  RESTART_MSG="start from scratch"
  if [[ -f "${RUN}/pull.cpt" ]]; then
    if [[ "${LAST_STEP}" =~ ^[0-9]+$ && "${EXPECTED_STEPS}" =~ ^[0-9]+$ ]]; then
      RESTART_MSG="restart from checkpoint (step ${LAST_STEP} of ${EXPECTED_STEPS})"
    elif [[ "${LAST_STEP}" =~ ^[0-9]+$ ]]; then
      RESTART_MSG="restart from checkpoint (step ${LAST_STEP})"
    else
      RESTART_MSG="restart from checkpoint (step unknown)"
    fi
  fi

  # Not completed = submit (even if in ledger or has partial files)
  echo "[smd-submit-new] SUBMIT: ${RUN} -- ${RESTART_MSG}" >&2

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
