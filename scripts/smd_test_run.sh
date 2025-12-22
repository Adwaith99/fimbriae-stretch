#!/usr/bin/env bash
set -euo pipefail

# Single-line SMD test wrapper: submit a short Slurm job for GPU or CPU
# - Prepares the run via smd_runner.sh
# - Runs mdrun for ~6 min using -maxh 0.1 (request 10 min walltime)
# - Intended to gauge performance before full SMD arrays
#
# Usage:
#   scripts/smd_test_run.sh gpu --lines 12[,34,...] [--clean]
#   scripts/smd_test_run.sh cpu --lines 12[,34,...] [--clean]
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
  echo "Usage: $0 <gpu|cpu> --lines <N[,N2,...]> [--clean]" >&2
  exit 2
fi

MODE="$1"; shift
CLEAN=0
LINE_SPEC=""

while [[ $# -gt 0 ]]; do
  case "$1" in
    --lines)
      LINE_SPEC="$2"; shift 2;;
    --clean)
      CLEAN=1; shift;;
    *)
      echo "Unknown arg: $1" >&2; exit 2;;
  esac
done

if [[ -z "$LINE_SPEC" ]]; then
  echo "ERROR: specify --lines <N[,N2,...]>" >&2
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

IFS=',' read -r -a LINES <<< "$LINE_SPEC"

mkdir -p logs tmp

for LN in "${LINES[@]}"; do
  # Extract CSV line (skip header)
  LINE_TXT=$(awk -F, -v N="$LN" 'NR==N{print $0}' "$MAN")
  if [[ -z "$LINE_TXT" ]]; then
    echo "WARN: no manifest row at line $LN; skipping" >&2
    continue
  fi

  # Parse fields to derive run dir for optional cleanup
  read -r system variant replicate speed_nm_per_ns _dt _k _axis _perf _start _tpr _xtc start_id _anch _cap <<<$(echo "$LINE_TXT" | awk -F, '{print $1,$2,$3,$4,$6,$5,$8,$9,$10,$11,$12,$13,$14,$15}')
  # Format speed tag
  vtag=$(python3 - <<PY
s=float("${speed_nm_per_ns}")
print(f"v{s:.3f}")
PY
)
  run_dir="$ROOT_DIR/smd/${system}/${variant}/${vtag}/rep${replicate}/${start_id}"

  # Build sbatch script
  job="tmp/smd_test_${MODE}_L${LN}.sbatch"
  cat > "$job" <<SB
#!/usr/bin/env bash
#SBATCH --job-name="smdt_${system}_${variant}_L${LN}"
#SBATCH --partition="${PART}"
#SBATCH --time="00:10:00"
#SBATCH --output="logs/smdt_%j.out"
SB
  if [[ "$MODE" == "cpu" ]]; then
    cat >> "$job" <<SB
#SBATCH --nodes="${NODES}"
#SBATCH --ntasks-per-node="${NTASKS_PER_NODE}"
#SBATCH --cpus-per-task="${CPUS}"
SB
  else
    if [[ -n "$GRES_SPEC" ]]; then
      cat >> "$job" <<SB
#SBATCH --gres="gpu:${GRES_SPEC}"
SB
    fi
    cat >> "$job" <<SB
#SBATCH --cpus-per-task="${CPUS}"
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

  jid=$(sbatch "$job" | awk '{print $NF}')
  echo "[smd-test] Submitted ${MODE} test for line ${LN} as job ${jid}"
done
