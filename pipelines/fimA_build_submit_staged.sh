#!/usr/bin/env bash
# Submit a staged build: EM → [ (NVT,NPT) x 5 with decreasing POSRES ] → final NPT (no POSRES)
# Each stage is a separate Slurm job with its own --time.
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

# ===== Stage parameters you can tune =====
# Restraint schedule (fcx=fcy=fcz)
FC_LIST=(1000 500 200 100 50)

# Steps per stage
NVT_STEPS=200000     # 200 ps at 1 fs
NPT_STEPS=500000     # 500 ps at 1 fs
FINAL_NPT_STEPS=1000000  # 1 ns no-restraint

# Walltimes per stage (HH:MM:SS)
# Keep these tight to avoid asking the scheduler for unnecessary time.
TIME_EM="00:30:00"
TIME_NVT="01:00:00"
TIME_NPT="02:00:00"
TIME_NPT_FINAL="06:00:00"   # longer for the final NPT

# ========================================

# 0) Prepare inputs (only first time)
if [[ ! -f "$WD/input.pdb" ]]; then
  cp -f "$ROOT/$PDB_IN" "$WD/input.pdb"
fi

# Helper: submit one stage; returns JOBID
# args: STAGE_NAME MODE TEMPK NSTEPS FC TIME DEPEND_JOBID
submit_stage () {
  local STAGE_NAME=$1 MODE=$2 TEMPK=$3 NSTEPS=$4 FC=$5 TIME=$6 DEP=$7
  sbatch ${DEP:+--dependency=afterok:$DEP} <<EOT | awk '{print $4}'
#!/bin/bash
#SBATCH --job-name=${STAGE_NAME}_${SYS}
#SBATCH --partition=$(python3 - <<PY
import yaml; c=yaml.safe_load(open("$ROOT/config.yaml"))
print(c["globals"]["slurm"]["partition"])
PY
)
#SBATCH --gpus-per-node=$(python3 - <<PY
import yaml; c=yaml.safe_load(open("$ROOT/config.yaml"))
print(c["globals"]["slurm"]["gpus_per_node"])
PY
)
#SBATCH --cpus-per-task=$(python3 - <<PY
import yaml; c=yaml.safe_load(open("$ROOT/config.yaml"))
print(c["globals"]["slurm"]["cpus_per_task"])
PY
)
#SBATCH --time=${TIME}
#SBATCH --output=$LOG/${STAGE_NAME}_${SYS}.out
module purge
module load ${GMX_MOD}
bash "$ROOT/pipelines/equil_stage.sh" "$SYS" "$MODE" "$TEMPK" "$NSTEPS" "$FC"
EOT
}

# ---- submit chain ----

# EM (short)
EM_ID=$(submit_stage "em" "em" 0 50000 0 "${TIME_EM}" "")
echo "Submitted EM: $EM_ID"

# Five pairs: (NVT, NPT) with decreasing restraints
DEP="$EM_ID"
for FC in "${FC_LIST[@]}"; do
  NVT_ID=$(submit_stage "nvt_fc${FC}" "nvt" 303.15 ${NVT_STEPS} ${FC} "${TIME_NVT}" "$DEP")
  echo "Submitted NVT fc=${FC}: $NVT_ID"
  DEP="$NVT_ID"

  NPT_ID=$(submit_stage "npt_fc${FC}" "npt" 303.15 ${NPT_STEPS} ${FC} "${TIME_NPT}" "$DEP")
  echo "Submitted NPT fc=${FC}: $NPT_ID"
  DEP="$NPT_ID"
done

# Final NPT without restraints (longer walltime)
FINAL_ID=$(submit_stage "npt_final" "npt" 303.15 ${FINAL_NPT_STEPS} 0 "${TIME_NPT_FINAL}" "$DEP")
echo "Submitted final NPT (no restraints): $FINAL_ID"
echo "Chain submitted. Final job id: $FINAL_ID"
