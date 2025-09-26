#!/usr/bin/env bash
# Submit staged build: EM → (NVT,NPT)x5 ↓POSRES → Final NPT (equilibrium MD @ dt from config)
# USAGE: fimA_build_submit_staged.sh <SYSTEM> <PDB> <BOX_X> <BOX_Y> <BOX_Z> <GMX_MODULE>

set -euo pipefail
SYS=${1:?}
PDB_IN=${2:?}
BOX_X=${3:?}
BOX_Y=${4:?}
BOX_Z=${5:?}
GMX_MOD=${6:?}

ROOT="$(cd "$(dirname "$0")/.."; pwd)"
WD="$ROOT/systems/${SYS}/00_build"
LOG="$ROOT/logs"
mkdir -p "$WD" "$LOG"

# === Pull scheduler & equilibrium settings from config ===
read -r PARTITION GPUS CPUS TIME_EM TIME_NVT TIME_NPT TIME_NPT_FINAL EQ_NS <<< "$(python3 - <<'PY'
import yaml
c=yaml.safe_load(open("config.yaml"))
g=c["globals"]; s=g["slurm"]; eq=g["equilibrium_md"]
print(
    s["partition"],
    s["gpus_per_node"],
    s["cpus_per_task"],
    eq["time_em"],
    eq["time_nvt"],
    eq["time_npt"],
    eq["time_npt_final"],
    eq["length_ns"],
)
PY
)"

# Restraint schedule (unchanged)
FC_LIST=(1000 500 200 100 50)

# Stage step counts (equilibration phases keep 1 fs as in equil_stage.sh)
NVT_STEPS=200000
NPT_STEPS=500000

# >>> Final NPT steps: use dt_ps from config (intended 0.002 ps)
FINAL_NPT_STEPS=$(python3 - <<'PY'
import yaml, math
c=yaml.safe_load(open("config.yaml"))
dt_ps=float(c["globals"]["dt_ps"])           # e.g., 0.002
length_ns=float(c["globals"]["equilibrium_md"]["length_ns"])
print(int(round(length_ns*1_000_000/dt_ps))) # steps = (ns*1e6 ps/ns)/dt_ps
PY
)

# 0) Prepare inputs (first time only)
[[ -f "$WD/input.pdb" ]] || cp -f "$ROOT/$PDB_IN" "$WD/input.pdb"

submit_stage () {
  local NAME=$1 MODE=$2 TEMP=$3 NSTEPS=$4 FC=$5 FINALFLAG=$6 TIME=$7 DEP=$8
  sbatch ${DEP:+--dependency=afterok:$DEP} <<EOT | awk '{print $4}'
#!/bin/bash
#SBATCH --job-name=${NAME}_${SYS}
#SBATCH --partition=${PARTITION}
#SBATCH --gpus-per-node=${GPUS}
#SBATCH --cpus-per-task=${CPUS}
#SBATCH --time=${TIME}
#SBATCH --output=$LOG/${NAME}_${SYS}.out
module purge
module load ${GMX_MOD}
bash "$ROOT/pipelines/equil_stage.sh" "$SYS" "$MODE" "$TEMP" "$NSTEPS" "$FC" "$FINALFLAG"
EOT
}

# EM
EM_ID=$(submit_stage "em" "em" 0 50000 0 0 "${TIME_EM}" "")
echo "EM: $EM_ID"
DEP="$EM_ID"

# (NVT,NPT) × 5 with decreasing restraints (1 fs stages)
for FC in "${FC_LIST[@]}"; do
  NVT_ID=$(submit_stage "nvt_fc${FC}" "nvt" 303.15 ${NVT_STEPS} ${FC} 0 "${TIME_NVT}" "$DEP")
  echo "NVT fc=${FC}: $NVT_ID"; DEP="$NVT_ID"
  NPT_ID=$(submit_stage "npt_fc${FC}" "npt" 303.15 ${NPT_STEPS} ${FC} 0 "${TIME_NPT}" "$DEP")
  echo "NPT fc=${FC}: $NPT_ID"; DEP="$NPT_ID"
done

# Final NPT equilibrium MD (dt from config, longer walltime)
FINAL_ID=$(submit_stage "npt_final" "npt" 303.15 ${FINAL_NPT_STEPS} 0 1 "${TIME_NPT_FINAL}" "$DEP")
echo "Final NPT: $FINAL_ID"
echo "Chain submitted. Final job id: $FINAL_ID"
