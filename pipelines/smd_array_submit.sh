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
g = cfg.get("globals", {})
slurm = g.get("slurm", {})
print(g.get("target_extension_nm", 5.0))           # 0
print(g.get("walltime_safety", 1.20))              # 1
print(g.get("max_wall_hours", 168))                # 2
print(slurm.get("partition", "default"))           # 3
print(slurm.get("array_cap", 5))                   # 4
print(slurm.get("cpus_per_task", 8))               # 5
print(slurm.get("gromacs_module", ""))             # 6
print(slurm.get("gpus_per_node", ""))              # 7
print(slurm.get("nodes", 1))                       # 8
print(slurm.get("ntasks_per_node", 1))             # 9
PY
)
TARGET_NM="${CFG_LINES[0]}"
SAFETY="${CFG_LINES[1]}"
MAX_HRS="${CFG_LINES[2]}"
PART="${CFG_LINES[3]}"
CAP="${CFG_LINES[4]}"
CPUS="${CFG_LINES[5]}"
GMX_MOD="${CFG_LINES[6]}"
GRES_SPEC="${CFG_LINES[7]}"
NODES="${CFG_LINES[8]}"
NTASKS_PER_NODE="${CFG_LINES[9]}"

# CPU mode toggle (export CPU=1 to use CPU)
CPU_MODE="${CPU:-0}"

# CPU mode toggle (export CPU=1 to use CPU)
CPU_MODE="${CPU:-0}"

# Normalize Windows line endings
sed -i 's/\r$//' "$MAN" 2>/dev/null || true

compute_wall() {
  local sys="$1" minspeed="$2"
  python3 - "$CFG" "$sys" "$minspeed" "$TARGET_NM" "$SAFETY" "$MAX_HRS" <<'PY'
import sys, yaml, math
cfg, sysname, mins, target_nm, safety, max_hrs = sys.argv[1], sys.argv[2], float(sys.argv[3]), float(sys.argv[4]), float(sys.argv[5]), int(float(sys.argv[6]))
data = yaml.safe_load(open(cfg))
globals_perf = float(data.get("globals", {}).get("perf_ns_per_day", 10.0))
perf = globals_perf
for s in data.get("systems", []):
    if s.get("name") == sysname:
        perf = float(s.get("perf_ns_per_day", globals_perf))
        break
needed_ns = target_nm / mins if mins > 0 else 0.0
days = needed_ns / perf if perf > 0 else 1e9
hours = math.ceil(days * 24 * float(safety))
hours = max(1, min(hours, max_hrs))
print(f"{hours}:00:00")
PY
}

mkdir -p logs manifests/tmp

# Collect systems with PENDING rows
mapfile -t SYSTEMS < <(awk -F, 'NR>1 && tolower($7)=="pending"{print $2}' "$MAN" | sort -u)
[ "${#SYSTEMS[@]}" -eq 0 ] && { echo "[smd-submit] No PENDING rows."; exit 0; }

echo "[smd-submit] Systems to submit: ${#SYSTEMS[@]}"
for sys in "${SYSTEMS[@]}"; do
  RUNMAN="manifests/tmp/smd_${sys}.csv"
  awk -F, -v S="$sys" 'NR==1||($2==S&&tolower($7)=="pending"){print $0}' "$MAN" > "$RUNMAN"
  N=$(awk -F, 'NR>1{c++} END{print c+0}' "$RUNMAN")
  [ "$N" -eq 0 ] && { echo "[smd-submit] Skip ${sys} (no rows)"; continue; }

  MINS=$(awk -F, 'NR>1{print $4}' "$RUNMAN" | sort -g | head -n1)
  WALL=$(compute_wall "$sys" "$MINS")

  SBATCH_ARGS=(--job-name="pulls_${sys}" --partition="${PART}" --array="1-${N}%${CAP}" --time="${WALL}" --output="logs/pulls_%A_%a.out")
  EXPORTS="ALL,FIM_PULL_MANIFEST=${RUNMAN},SYSTEM=${sys},GMX_MODULE=${GMX_MOD}"

  if [ "$CPU_MODE" = "1" ]; then
    # CPU mode
    SBATCH_ARGS+=(--nodes="${NODES}" --ntasks-per-node="${NTASKS_PER_NODE}" --cpus-per-task="${CPUS}")
    EXPORTS="${EXPORTS},GMX_CMD=srun gmx_mpi mdrun -v -deffnm pull"
    echo "[smd-submit] CPU submit sys=${sys} rows=${N} wall=${WALL} part=${PART} nodes=${NODES} ntasks/node=${NTASKS_PER_NODE} cpus/task=${CPUS}"
  else
    # GPU mode
    [ -n "$GRES_SPEC" ] && SBATCH_ARGS+=(--gres="gpu:${GRES_SPEC}")
    SBATCH_ARGS+=(--cpus-per-task="${CPUS}")
    echo "[smd-submit] GPU submit sys=${sys} rows=${N} wall=${WALL} part=${PART} gres=${GRES_SPEC} cpus/task=${CPUS}"
  fi

  sbatch "${SBATCH_ARGS[@]}" --export="${EXPORTS}" scripts/smd_runner.sh
done

echo "[smd-submit] Submission complete."
