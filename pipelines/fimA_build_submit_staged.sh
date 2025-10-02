#!/usr/bin/env bash
# Submit staged build: EM → (NVT,NPT)x5 (paired temps per FC) → Final NPT
# USAGE: fimA_build_submit_staged.sh <SYSTEM> <PDB> <BOX_X> <BOX_Y> <BOX_Z> <GMX_MODULE> [--force]
# Env: BUILD_FORCE=1 to force rerun even if outputs exist.

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

# Restraint sequence and per-FC temperature (NVT & NPT share the same T)
FC_LIST=(1000 500 200 100 50)
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

# Step counts
NVT_STEPS=200000   # 200 ps @ 1 fs
NPT_STEPS=500000   # 500 ps @ 1 fs
FINAL_NPT_STEPS=$(python3 - <<'PY'
import yaml
c=yaml.safe_load(open("config.yaml"))
dt=float(c["globals"]["dt_ps"])
length_ns=float(c["globals"]["equilibrium_md"]["length_ns"])
print(int(round(length_ns*1_000_000/dt)))
PY
)

# Small helpers
have() { [[ -f "$WD/$1" ]]; }
must_run() {  # return 0 if we should submit (either force, or outfile missing)
  [[ "$BUILD_FORCE" == "1" ]] && return 0
  ! have "$1"
}

submit_stage () {
  local NAME=$1 MODE=$2 TEMP=$3 NSTEPS=$4 FC=$5 FINALFLAG=$6 TIME=$7 DEP="$8" OUTFILE=$9

  if ! must_run "$OUTFILE"; then
    echo "SKIP  ${NAME}_${SYS} (found $OUTFILE)"
    printf "%s" ""  # print empty job id
    return 0
  fi

  echo "SUBMIT ${NAME}_${SYS} (T=${TEMP} K, FC=${FC}, steps=${NSTEPS}, final=${FINALFLAG})"
  # Build the sbatch command; use --parsable to return the job ID cleanly.
  local cmd=( sbatch --parsable --job-name="${NAME}_${SYS}" --partition="${PARTITION}"
              --gpus-per-node="${GPUS}" --cpus-per-task="${CPUS}" --time="${TIME}"
              --output="$LOG/${NAME}_${SYS}.out" )
  if [[ -n "$DEP" ]]; then cmd+=( --dependency="afterok:${DEP}" ); fi

  # Show a preview line (helpful for debugging)
  echo "[PREVIEW] ${cmd[*]}  <<< HEREDOC"

  # Feed the job script via heredoc. IMPORTANT: delimiter must be at start of line with no spaces.
  local JID
  JID="$(
    "${cmd[@]}" <<EOT
#!/bin/bash
set -euo pipefail
module purge
module load ${GMX_MOD}
# Export project root for equil_stage to resolve config.yaml reliably
export PIPE_ROOT="$ROOT"
bash "$ROOT/pipelines/equil_stage.sh" "$SYS" "$MODE" "$TEMP" "$NSTEPS" "$FC" "$FINALFLAG" "${BOX_X}" "${BOX_Y}" "${BOX_Z}"
EOT
  )"

  # sbatch --parsable returns just the job id (or ARRAYID_TASKS); print it
  echo "$JID"
}

# Seed input once
[[ -f "$WD/input.pdb" ]] || cp -f "$PDB_IN" "$WD/input.pdb"

DEP=""

# 0) EM
EM_JID="$(submit_stage "em" "em" 0 50000 0 0 "${TIME_EM}" "$DEP" "em.gro")"
[[ -n "$EM_JID" ]] && DEP="$EM_JID"

# 1..5) NVT/NPT pairs
for FC in "${FC_LIST[@]}"; do
  T="$(temp_for_fc "$FC")"
  NVT_JID="$(submit_stage "nvt_fc${FC}" "nvt" "${T}" ${NVT_STEPS} ${FC} 0 "${TIME_NVT}" "$DEP" "nvt_fc${FC}.gro")"
  [[ -n "$NVT_JID" ]] && DEP="$NVT_JID"

  NPT_JID="$(submit_stage "npt_fc${FC}" "npt" "${T}" ${NPT_STEPS} ${FC} 0 "${TIME_NPT}" "$DEP" "npt_fc${FC}.gro")"
  [[ -n "$NPT_JID" ]] && DEP="$NPT_JID"
done

# Final NPT (unrestrained)
FINAL_JID="$(submit_stage "npt_final" "npt" 303.15 ${FINAL_NPT_STEPS} 0 1 "${TIME_NPT_FINAL}" "$DEP" "npt_final.gro")"
if [[ -n "$FINAL_JID" ]]; then
  echo "Final NPT: ${FINAL_JID}"
else
  echo "Final NPT: SKIPPED (npt_final.gro present)"
fi
