#!/usr/bin/env bash
set -euo pipefail

MANIFEST="${1:-manifests/smd_manifest.csv}"
LEDGER="manifests/smd_submitted.csv"      # Keys of runs already submitted
SAFETY="${SMD_TIME_SAFETY_FACTOR:-1.20}"  # walltime pad; override via env
JOB_SCRIPT="${SMD_JOB_SCRIPT:-scripts/smd_job.sh}"

[[ -f "${MANIFEST}" ]] || { echo "[smd-submit-new] ERROR: manifest not found: ${MANIFEST}" >&2; exit 2; }
[[ -f "${JOB_SCRIPT}" ]] || { echo "[smd-submit-new] ERROR: job script not found: ${JOB_SCRIPT}" >&2; exit 2; }

# Read Slurm defaults from config (partition, gpus, cpus)
readarray -t CFG < <(python3 - <<'PY'
import yaml
cfg=yaml.safe_load(open("config.yaml"))
s=cfg["globals"]["slurm"]
print(s["partition"])
print(s["gpus_per_node"])
print(s["cpus_per_task"])
PY
)
PART="${CFG[0]}"; GPUS="${CFG[1]}"; CPUS="${CFG[2]}"

mkdir -p logs manifests/tmp
touch "${LEDGER}"

# Build a ledger set
declare -A submitted
while IFS=, read -r s v r sp sid || [[ -n "${s:-}" ]]; do
  [[ -z "${s:-}" ]] && continue
  submitted["$s,$v,$r,$sp,$sid"]=1
done < "${LEDGER}"

# Parse manifest rows to TSV (TAB-separated!)
tmpdir=$(mktemp -d); trap 'rm -rf "${tmpdir}"' EXIT
rows="${tmpdir}/rows.tsv"
# lineno  system  variant  replicate  speed  perf  target  start_id  array_cap
awk -F, -v OFS='\t' 'NR>1{
  gsub(/^"|"$/, "", $0);
  print NR, $1, $2, $3, $4, $9, $7, $13, $15
}' "${MANIFEST}" > "${rows}"

# Group NEW rows by (system, speed)
declare -A grp_idxs grp_time grp_cap
BAD=0
while IFS=$'\t' read -r L SYS VAR REP SPD PERF TGT START CAP; do
  # Skip empty lines defensively
  [[ -z "${L}" ]] && continue

  KEY="$SYS,$VAR,$REP,$SPD,$START"

  # Skip if already submitted
  if [[ -n "${submitted[$KEY]:-}" ]]; then
    continue
  fi

  # Validate required numeric fields
  if [[ -z "${SPD}" || -z "${TGT}" || -z "${PERF}" ]]; then
    echo "[smd-submit-new] WARN: skipping line ${L} (missing numeric: speed='${SPD}', target='${TGT}', perf='${PERF}')" >&2
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

  # Skip if already finished on disk
  RUN="smd/${SYS}/${VAR}/v$(printf "%0.3f" "${SPD}")/rep${REP}/${START}"
  if [[ -f "${RUN}/pull.xtc" || -f "${RUN}/pull.log" ]]; then
    echo "[smd-submit-new] SKIP finished: ${RUN}" >&2
    # do NOT mark in ledger on dry-run; the real submission loop handles ledger updates
    continue
  fi

  GRP="${SYS}|${SPD}"
  # Collect sparse array indices (CSV line numbers)
  if [[ -n "${grp_idxs[$GRP]:-}" ]]; then
    grp_idxs[$GRP]+=","${L}
  else
    grp_idxs[$GRP]="${L}"
  fi

  # walltime hours per row = (target/speed)/perf * 24 * SAFETY
  ht=$(python3 - <<PY
import math
spd=float("${SPD}"); tgt=float("${TGT}"); perf=float("${PERF}"); pad=float("${SAFETY}")
ns=tgt/spd; days=(ns/perf)
print("{:.3f}".format(days*24*pad))
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
  # ---- NEW: resolve absolute paths & chdir to repo root ----
  MANIFEST_ABS="$(readlink -f "${MANIFEST}")"
  ROOT="$(dirname "${MANIFEST_ABS}")/.."
  ROOT="$(readlink -f "${ROOT}")"
  JOB_SCRIPT_ABS="$(readlink -f "${JOB_SCRIPT}")"
  mkdir -p "${ROOT}/logs"

  OUT="${ROOT}/logs/${SYS}_v$(printf "%0.3f" "${SPD}")_%A_%a.log"   # unified stdout+stderr

  echo "[smd-submit-new] sbatch --chdir='${ROOT}' --job-name=${JNAME} --partition=${PART} --gpus-per-node=${GPUS} --cpus-per-task=${CPUS} --time=${TSTR} --array=${IDXS}%${CAP} --output='${OUT}' ${JOB_SCRIPT_ABS} ${MANIFEST_ABS}"

  if [[ -n "${SLURM_TEST_ONLY:-}" ]]; then
    echo "[smd-submit-new] DRY-RUN: sbatch --test-only (no submission)"
    any=1
    continue
  fi

  JID=$(
    sbatch --parsable \
      --chdir="${ROOT}" \
      --job-name="${JNAME}" \
      --partition="${PART}" \
      --gpus-per-node="${GPUS}" \
      --cpus-per-task="${CPUS}" \
      --time="${TSTR}" \
      --array="${IDXS}%${CAP}" \
      --output="${OUT}" \
      --error="${OUT}" \
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
