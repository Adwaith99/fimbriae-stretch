#!/usr/bin/env bash
# Submit staged build: EM → (NVT,NPT)x5 (paired temps per FC) → Final NPT
# USAGE: fimA_build_submit_staged.sh <SYSTEM> <PDB> <BOX_X> <BOX_Y> <BOX_Z> <GMX_MODULE> [--force]
# Env override: BUILD_FORCE=1  (re-run even if outputs exist)

set -euo pipefail
SYS=${1:?}
PDB_IN=${2:?}
BOX_X=${3:?}
BOX_Y=${4:?}
BOX_Z=${5:?}
GMX_MOD=${6:?}
BUILD_FORCE=${BUILD_FORCE:-0}
if [[ "${7:-}" == "--force" ]]; then BUILD_FORCE=1; fi

ROOT="$(cd "$(dirname "$0")/.."; pwd)"
WD="$ROOT/systems/${SYS}/00_build"
LOG="$ROOT/logs"
mkdir -p "$WD" "$LOG"

# ---- read slurm + eq settings from config ----
read -r PARTITION GPUS CPUS TIME_EM TIME_NVT TIME_NPT TIME_NPT_FINAL EQ_NS <<< "$(python3 - <<'PY'
import yaml
c=yaml.safe_load(open("config.yaml"))
g=c["globals"]; s=g["slurm"]; eq=g["equilibrium_md"]
print(s["partition"], s["gpus_per_node"], s["cpus_per_task"],
      eq["time_em"], eq["time_nvt"], eq["time_npt"], eq["time_npt_final"], eq["length_ns"])
PY
)"

# Restraint schedule (high → low)
FC_LIST=(1000 500 200 100 50)

# Step counts (restrained stages @ 1 fs)
NVT_STEPS=200000   # 200 ps
NPT_STEPS=500000   # 500 ps

# Final NPT steps @ globals.dt_ps
FINAL_NPT_STEPS=$(python3 - <<'PY'
import yaml
c=yaml.safe_load(open("config.yaml"))
dt=float(c["globals"]["dt_ps"])
length_ns=float(c["globals"]["equilibrium_md"]["length_ns"])
print(int(round(length_ns*1_000_000/dt)))
PY
)

# Temperature schedule for BOTH NVT and NPT per FC (defaults = your original)
temp_for_fc() {
python3 - "$1" <<'PY'
import sys, yaml
fc=sys.argv[1]
cfg=yaml.safe_load(open("config.yaml"))
default={"1000":100.0,"500":200.0,"200":303.15,"100":303.15,"50":303.15}
tmap=cfg.get("globals",{}).get("build",{}).get("temp_schedule_by_fc",default)
print(tmap.get(str(fc), default.get(str(fc), 303.15)))
PY
}

# Helper: check if a stage's expected output exists
is_done() {
  # $1 = filename to check (relative to $WD)
  [[ -f "$WD/$1" ]]
}

# Helper: submit a stage unless already done; handle dependencies
submit_stage () {
  local NAME=$1 MODE=$2 TEMP=$3 NSTEPS=$4 FC=$5 FINALFLAG=$6 TIME=$7 DEP=$8 OUTFILE=$9

  if [[ "$BUILD_FORCE" == "0" ]] && is_done "$OUTFILE"; then
    echo "SKIP  ${NAME}_${SYS} (found $OUTFILE)"
    # If we skip, clear dependency so the next stage doesn't wait on us
    echo ""
    return
  fi

  echo "SUBMIT ${NAME}_${SYS} (temp=${TEMP}K, fc=${FC}, steps=${NSTEPS}, final=${FINALFLAG})"
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
# Pass a stable root to the stage script for config.yaml resolution
export PIPE_ROOT="$ROOT"
bash "$ROOT/pipelines/equil_stage.sh" "$SYS" "$MODE" "$TEMP" "$NSTEPS" "$FC" "$FINALFLAG" "${BOX_X}" "${BOX_Y}" "${BOX_Z}"
EOT
}

# Ensure input present once
[[ -f "$WD/input.pdb" ]] || cp -f "$PDB_IN" "$WD/input.pdb"

# ==== Chain submission with resume ====
DEP=""

# 0) EM
EM_JOBID=$(submit_stage "em" "em" 0 50000 0 0 "${TIME_EM}" "$DEP" "em.gro")
if [[ -n "${EM_JOBID}" ]]; then DEP="$EM_JOBID"; fi

# 1..5) Paired NVT/NPT at same T for each FC
for FC in "${FC_LIST[@]}"; do
  T=$(temp_for_fc "$FC")

  # NVT (output: nvt_fcXX.gro)
  NVT_JOBID=$(submit_stage "nvt_fc${FC}" "nvt" "${T}" ${NVT_STEPS} ${FC} 0 "${TIME_NVT}" "$DEP" "nvt_fc${FC}.gro")
  if [[ -n "${NVT_JOBID}" ]]; then DEP="$NVT_JOBID"; fi

  # NPT (output: npt_fcXX.gro)
  NPT_JOBID=$(submit_stage "npt_fc${FC}" "npt" "${T}" ${NPT_STEPS} ${FC} 0 "${TIME_NPT}" "$DEP" "npt_fc${FC}.gro")
  if [[ -n "${NPT_JOBID}" ]]; then DEP="$NPT_JOBID"; fi
done

# Final unrestrained NPT (output: npt_final.gro)
FINAL_JOBID=$(submit_stage "npt_final" "npt" 303.15 ${FINAL_NPT_STEPS} 0 1 "${TIME_NPT_FINAL}" "$DEP" "npt_final.gro")
if [[ -n "${FINAL_JOBID}" ]]; then
  echo "Final NPT: $FINAL_JOBID"
else
  echo "Final NPT: SKIPPED (npt_final.gro present)"
fi
