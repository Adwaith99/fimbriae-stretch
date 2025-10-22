#!/usr/bin/env bash
set -euo pipefail

MANIFEST="${1:-manifests/smd_manifest.csv}"
LEDGER="manifests/smd_submitted.csv"      # Keys of runs already submitted
SAFETY="${SMD_TIME_SAFETY_FACTOR:-1.20}"  # walltime pad; override via env
JOB_SCRIPT="${SMD_JOB_SCRIPT:-scripts/smd_job.sh}"

[[ -f "${MANIFEST}" ]] || { echo "[smd-submit-new] ERROR: manifest not found: ${MANIFEST}" >&2; exit 2; }
[[ -f "${JOB_SCRIPT}" ]] || { echo "[smd-submit-new] ERROR: job script not found: ${JOB_SCRIPT}" >&2; exit 2; }

# Read Slurm defaults from config (partition, gpus, cpus); module is loaded inside job script
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

# Parse manifest rows to TSV for grouping
tmpdir=$(mktemp -d); trap 'rm -rf "${tmpdir}"' EXIT
rows="${tmpdir}/rows.tsv"
# lineno  system  variant  replicate  speed  perf  target  start_id  array_cap
awk -F, 'NR>1{
  gsub(/^"|"$/, "", $0);
  print NR, $1, $2, $3, $4, $9, $7, $13, $15
}' "${MANIFEST}" > "${rows}"

# Group NEW rows by (system, speed)
declare -A grp_idxs grp_time grp_cap
while IFS=$'\t' read -r L SYS VAR REP SPD PERF TGT START CAP; do
  KEY="$SYS,$VAR,$REP,$SPD,$START"
  # Skip if already submitted
  [[ -n "${submitted[$KEY]:-}" ]] && continue
  # Skip if already finished
  RUN="smd/${SYS}/${VAR}/v$(printf "%0.3f" "${SPD}")/rep${REP}/${START}"
  if [[ -f "${RUN}/pull.xtc" || -f "${RUN}/pull.log" ]]; then
    echo "$KEY" >> "${LEDGER}"
    continue
  fi
  GRP="${SYS}|${SPD}"
  grp_idxs[$GRP]="${grp_idxs[$GRP]:-}${grp_idxs[$GRP]:+,}${L}"
  # walltime hours per row = (target/speed)/perf * 24 * SAFETY
  ht=$(python3 - <<PY
import math
spd=float("${SPD}"); tgt=float("${TGT}"); perf=float("${PERF}"); pad=float("${SAFETY}")
ns=tgt/spd; days=(ns/perf) if perf>0 else 1e6
print("{:.3f}".format(days*24*pad))
PY
)
  cur="${grp_time[$GRP]:-0.0}"
  grp_time[$GRP]=$(python3 - <<PY
a=float("${cur}"); b=float("${ht}")
print("{:.3f}".format(max(a,b)))
PY
)
  grp_cap[$GRP]="${CAP}"
done < "${rows}"

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
  JNAME="smd:${SYS}:v$(printf "%0.3f" "${SPD}")"
  OUT="logs/${SYS}_v$(printf "%0.3f" "${SPD}")_%A_%a.out"
  ERR="logs/${SYS}_v$(printf "%0.3f" "${SPD}")_%A_%a.err"
  echo "[smd-submit-new] sbatch --job-name=${JNAME} --partition=${PART} --gpus-per-node=${GPUS} --cpus-per-task=${CPUS} --time=${TSTR} --array=${IDXS}%${CAP}"

  JID=$(
    sbatch --parsable \
      --job-name="${JNAME}" \
      --partition="${PART}" \
      --gpus-per-node="${GPUS}" \
      --cpus-per-task="${CPUS}" \
      --time="${TSTR}" \
      --array="${IDXS}%${CAP}" \
      --output="${OUT}" \
      --error="${ERR}" \
      "${JOB_SCRIPT}" "${MANIFEST}"
  )
  echo "[smd-submit-new] submitted: ${JID}"

  # Append submitted keys for this group to the ledger
  pat="^($(echo "${IDXS}" | sed 's/,/|/g'))$"
  awk -F $'\t' -v pat="${pat}" ' $1 ~ pat { printf "%s,%s,%s,%s,%s\n",$2,$3,$4,$5,$8 } ' "${rows}" >> "${LEDGER}"

  any=1
done

if [[ "${any}" -eq 0 ]]; then
  echo "[smd-submit-new] Nothing new to submit (ledger up to date)."
fi
