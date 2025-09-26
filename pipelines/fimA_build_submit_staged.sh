#!/usr/bin/env bash
# Submit a staged build: EM → [ (NVT,NPT) x 5 with decreasing POSRES ] → final NPT (no POSRES)
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

# Stage force constants (fcx=fcy=fcz) for POSRES stages
# e.g., strong→weak: 1000, 500, 200, 100, 50 (tune as you like)
FC_LIST=(1000 500 200 100 50)

# 0) Prepare inputs (only first time)
if [[ ! -f "$WD/input.pdb" ]]; then
  cp -f "$ROOT/$PDB_IN" "$WD/input.pdb"
fi

# Helper to submit a stage and return JOBID
submit_stage () {
  local STAGE_NAME=$1 MODE=$2 TEMPK=$3 NSTEPS=$4 FC=$5 DEP=$6
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
#SBATCH --time=04:00:00
#SBATCH --output=$LOG/${STAGE_NAME}_${SYS}.out
module purge
module load ${GMX_MOD}
bash "$ROOT/pipelines/equil_stage.sh" "$SYS" "$MODE" "$TEMPK" "$NSTEPS" "$FC"
EOT
}

# ---- submit chain ----

# EM (no T)
EM_ID=$(submit_stage "em" "em" 0 50000 0 "")
echo "Submitted EM: $EM_ID"

# Five pairs: (NVT, NPT) with decreasing restraints
DEP="$EM_ID"
# You can tweak nsteps if you want longer/shorter segments
NVT_STEPS=200000     # 200 ps at 1 fs
NPT_STEPS=500000     # 500 ps at 1 fs
for FC in "${FC_LIST[@]}"; do
  NVT_ID=$(submit_stage "nvt_fc${FC}" "nvt" 303.15 $NVT_STEPS $FC "$DEP")
  echo "Submitted NVT fc=${FC}: $NVT_ID"
  DEP="$NVT_ID"
  NPT_ID=$(submit_stage "npt_fc${FC}" "npt" 303.15 $NPT_STEPS $FC "$DEP")
  echo "Submitted NPT fc=${FC}: $NPT_ID"
  DEP="$NPT_ID"
done

# Final NPT without restraints (longer, e.g., 1 ns)
FINAL_NPT_STEPS=1000000
FINAL_ID=$(submit_stage "npt_final" "npt" 303.15 $FINAL_NPT_STEPS 0 "$DEP")
echo "Submitted final NPT (no restraints): $FINAL_ID"

echo "Chain submitted. Final job id: $FINAL_ID"