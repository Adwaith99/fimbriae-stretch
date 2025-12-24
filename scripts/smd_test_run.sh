#!/usr/bin/env bash
set -euo pipefail

# Single-line SMD test wrapper: submit a short Slurm job for GPU or CPU
# - Prepares the run via smd_runner.sh
# - Runs mdrun for ~10 min using -maxh 0.1 (request 15 min walltime)
# - Intended to gauge performance before full SMD arrays
#
# Usage:
#   scripts/smd_test_run.sh gpu --line 2 [--clean]
#   scripts/smd_test_run.sh cpu --line 2 [--clean]
#
# Notes:
# - GPU/CPU resources and modules are taken from config.yaml (globals.slurm)
#   * For GPU tests: set a GPU partition and uncomment the CUDA-enabled gromacs_module.
#   * For CPU tests: use a CPU partition and keep the non-CUDA gromacs_module.
# - This wrapper does NOT alter runner behavior; it only packages a short sbatch.
# - CLEAN (optional): when set, removes the smd/<sys>/<variant>/v<speed>/rep<rep>/<start> dir after the test.

ROOT_DIR="$(pwd)"
MAN="manifests/smd_manifest.csv"
CFG="config.yaml"

if [[ $# -lt 2 ]]; then
  echo "Usage: $0 <gpu|cpu> --line <N> [--clean]" >&2
  exit 2
fi

MODE="$1"; shift
CLEAN=0
LINE_NUM=""

while [[ $# -gt 0 ]]; do
  case "$1" in
    --line)
      LINE_NUM="$2"; shift 2;;
    --clean)
      CLEAN=1; shift;;
    *)
      echo "Unknown arg: $1" >&2; exit 2;;
  esac
done

if [[ -z "$LINE_NUM" ]]; then
  echo "ERROR: specify --line <N>" >&2
  exit 2
fi

# Normalize Windows line endings
sed -i 's/\r$//' "$MAN" 2>/dev/null || true

readarray -t CFG_LINES < <(python3 - <<'PY'
import yaml
cfg = yaml.safe_load(open("config.yaml"))
g = cfg.get("globals", {})
slurm = g.get("slurm", {})
print(slurm.get("partition", "default"))           # 0
print(slurm.get("array_cap", 1))                     # 1 (unused here)
print(slurm.get("cpus_per_task", 4))                # 2
print(slurm.get("gromacs_module", ""))             # 3
print(slurm.get("gpus_per_node", ""))              # 4
print(slurm.get("nodes", 1))                         # 5
print(slurm.get("ntasks_per_node", 1))               # 6
PY
)
PART="${CFG_LINES[0]}"
CPUS="${CFG_LINES[2]}"
GMX_MOD="${CFG_LINES[3]}"
GRES_SPEC="${CFG_LINES[4]}"
NODES="${CFG_LINES[5]}"
NTASKS_PER_NODE="${CFG_LINES[6]}"

echo "[smd-test] Testing manifest line $LINE_NUM (MODE=$MODE)" >&2

mkdir -p logs tmp

# Extract CSV line using sed (1-based line number)
LINE_TXT=$(sed -n "${LINE_NUM}p" "$MAN")
if [[ -z "$LINE_TXT" ]]; then
  echo "ERROR: manifest line $LINE_NUM is empty or missing" >&2
  exit 1
fi

# Skip header line
if [[ "$LINE_TXT" =~ ^system,variant, ]]; then
  echo "ERROR: line $LINE_NUM is the CSV header; use line 2+ for data rows" >&2
  exit 1
fi

echo "[smd-test] Processing manifest line $LINE_NUM: $(echo "$LINE_TXT" | cut -d, -f1-3)" >&2

# Extract run directory path for optional cleanup using Python CSV parser
echo "[smd-test] Extracting run_dir from Python CSV parser..." >&2
run_dir=$(python3 - "$MAN" "$LINE_NUM" "$ROOT_DIR" <<'PYEOF'
import sys, csv
manifest, line_num, root = sys.argv[1], int(sys.argv[2]), sys.argv[3]
with open(manifest, 'r') as f:
    reader = csv.DictReader(f)
    for i, row in enumerate(reader, start=2):  # start=2 because line 1 is header
        if i == line_num:
            speed = float(row['speed_nm_per_ns'])
            vtag = f"v{speed:.3f}"
            path = f"{root}/smd/{row['system']}/{row['variant']}/{vtag}/rep{row['replicate']}/{row['start_id']}"
            print(path)
            break
PYEOF
)

if [[ -z "$run_dir" ]]; then
  echo "ERROR: failed to parse run directory from line $LINE_NUM" >&2
  exit 1
fi

echo "[smd-test] run_dir='$run_dir'" >&2

# Build sbatch script
job="tmp/smd_test_${MODE}_L${LINE_NUM}.sbatch"
cat > "$job" <<SB
#!/usr/bin/env bash
#SBATCH --job-name=smdt_L${LINE_NUM}
#SBATCH --partition=${PART}
#SBATCH --time=00:15:00
#SBATCH --output=logs/smdt_%j.out
SB
if [[ "$MODE" == "cpu" ]]; then
  cat >> "$job" <<SB
#SBATCH --nodes=${NODES}
#SBATCH --ntasks-per-node=${NTASKS_PER_NODE}
#SBATCH --cpus-per-task=${CPUS}
SB
else
  if [[ -n "$GRES_SPEC" ]]; then
    cat >> "$job" <<SB
#SBATCH --gres=gpu:${GRES_SPEC}
SB
  fi
  cat >> "$job" <<SB
#SBATCH --cpus-per-task=${CPUS}
SB
fi

cat >> "$job" <<'SB'
set -euo pipefail
module purge || true
SB
if [[ -n "$GMX_MOD" ]]; then
  echo "module load ${GMX_MOD}" >> "$job"
fi
cat >> "$job" <<SB
export GMX_MODULE="${GMX_MOD}"
SB
if [[ "$MODE" == "cpu" ]]; then
  cat >> "$job" <<'SB'
export GMX_CMD="srun gmx_mpi mdrun -maxh 0.1"
SB
else
  cat >> "$job" <<'SB'
export GMX_CMD="gmx mdrun -maxh 0.1"
SB
fi

# Invoke runner with the manifest line
printf '%s\n' "bash scripts/smd_runner.sh \"$LINE_TXT\"" >> "$job"

# Optional cleanup
if [[ $CLEAN -eq 1 ]]; then
  printf '%s\n' "rm -rf \"$run_dir\"" >> "$job"
fi

# Submit job
echo "[smd-test] Job script: $job" >&2
echo "[smd-test] Submitting sbatch..." >&2

sbatch_output=$(sbatch "$job" 2>&1) || sbatch_status=$?
sbatch_status=${sbatch_status:-$?}

# Print sbatch output
echo "========== sbatch output ==========" >&2
echo "$sbatch_output" >&2
echo "==================================" >&2

if [[ $sbatch_status -ne 0 ]]; then
  echo "[smd-test] ERROR: sbatch failed with status $sbatch_status" >&2
  exit 1
fi

jid=$(echo "$sbatch_output" | awk '{print $NF}')
echo "[smd-test] SUCCESS: Submitted ${MODE} test for line ${LINE_NUM} as job ${jid}" >&2
