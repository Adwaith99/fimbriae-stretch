#!/usr/bin/env bash
set -euo pipefail

# Preprocess SMD trajectories: extract protein-only, center, fix PBC, and stride
# Usage: preproc_traj.sh <smd_run_dir_relative_path>
# Example: preproc_traj.sh smd/fimA_WT/AtoD/v0.020/rep1/s075

# Find repo root
find_root() {
  local p="$PWD"
  while [[ "$p" != "/" ]]; do
    [[ -f "$p/config.yaml" ]] && { echo "$p"; return 0; }
    p="$(dirname "$p")"
  done
  return 1
}
ROOT="$(find_root)"
if [[ -z "$ROOT" ]]; then
  echo "[preproc-traj] ERROR: cannot locate repo root (config.yaml)" >&2
  exit 2
fi

# User settings (can be overridden by environment)
RAW_XTC_NAME="${RAW_XTC_NAME:-pull.xtc}"      # trajectory filename in run dir
PROT_GROUP="${PROT_GROUP:-1}"                 # Protein group index for gmx trjconv
GMX="${GMX:-gmx}"                              # gmx command
DX_NM="${DX_NM:-0.02}"                         # target extension step (nm) between frames

if [[ $# -ne 1 ]]; then
    echo "Usage: $0 <smd_run_dir_relative_path>" >&2
    echo "Example: $0 smd/fimA_WT/AtoD/v0.020/rep1/s075" >&2
    exit 1
fi

SDIR_REL="$1"
SDIR="${ROOT}/${SDIR_REL}"

if [[ ! -d "${SDIR}" ]]; then
    echo "[preproc-traj] WARN: Directory ${SDIR} does not exist, skipping" >&2
    exit 0
fi

# Extract speed, rep, start_id from path: smd/<system>/<variant>/v<speed>/rep<rep>/<start_id>
# Parse: strip leading "smd/" and extract components
path_after_smd="${SDIR_REL#smd/}"
IFS='/' read -r system variant vtag rep start_id <<< "$path_after_smd"

if [[ -z "${system}" || -z "${variant}" || -z "${vtag}" || -z "${rep}" || -z "${start_id}" ]]; then
  echo "[preproc-traj] ERROR: cannot parse path ${SDIR_REL} (expected smd/<system>/<variant>/v<speed>/rep<rep>/<start_id>)" >&2
  exit 2
fi

# Extract speed from vtag (e.g., "v0.020" -> "0.020")
speed="${vtag#v}"

# Look up dt_ps for this speed from config, or fall back to formula
DT_PS_THIS=$(python3 - <<PY
import yaml, math
cfg = yaml.safe_load(open(r"${ROOT}/config.yaml"))
speed_nm_per_ns = float("${speed}")
dx_nm = float("${DX_NM}")

# dt_ps = (dx_nm / speed_nm_per_ns) * 1000
dt_ps = (dx_nm / speed_nm_per_ns) * 1000.0
print(f"{dt_ps:.1f}")
PY
)

xtc_in="${SDIR}/${RAW_XTC_NAME}"
tpr_in="${SDIR}/pull.tpr"

if [[ ! -f "${xtc_in}" ]]; then
    echo "[preproc-traj] SKIP: Missing trajectory ${xtc_in}" >&2
    exit 0
fi
if [[ ! -f "${tpr_in}" ]]; then
    echo "[preproc-traj] SKIP: Missing tpr ${tpr_in}" >&2
    exit 0
fi

# Output directory: analysis/preproc_traj/<system>/<variant>/v<speed>/rep<rep>/<start_id>
OUTROOT="${ROOT}/analysis/preproc_traj/${system}/${variant}"
outdir="${OUTROOT}/${vtag}/${rep}/${start_id}"
mkdir -p "${outdir}"

xtc_out="${outdir}/traj_protein_dt${DT_PS_THIS}.xtc"
gro_out="${outdir}/start_protein.gro"

if [[ -f "${xtc_out}" && -f "${gro_out}" ]]; then
    echo "[preproc-traj] SKIP: Already processed ${SDIR_REL}" >&2
    exit 0
fi

echo "[preproc-traj] INFO: Processing ${SDIR_REL} (dt=${DT_PS_THIS} ps)" >&2
echo "[preproc-traj] INFO:   -> ${xtc_out}" >&2
echo "[preproc-traj] INFO:   -> ${gro_out}" >&2

# Optionally load GROMACS module from config if available
GMX_MODULE=$(python3 - <<PY
import yaml
cfg=yaml.safe_load(open(r"${ROOT}/config.yaml"))
m=cfg.get('globals',{}).get('slurm',{}).get('gromacs_module','')
print(m)
PY
)
if [[ -n "${GMX_MODULE}" ]] && type module >/dev/null 2>&1; then
  module purge >/dev/null 2>&1 || true
  module load ${GMX_MODULE} >/dev/null 2>&1 || true
  echo "[preproc-traj] INFO: loaded module: ${GMX_MODULE}" >&2
fi

# Be nice: single-thread for OpenMP
export OMP_NUM_THREADS=1

# Use index.ndx if present (to ensure group names exist)
NDX="${SDIR}/index.ndx"; NDX_FLAG=()
[[ -f "${NDX}" ]] && NDX_FLAG=( -n "${NDX}" )

# 1) Protein-only, PBC-fixed, centered, strided trajectory
printf "%s\n%s\n" "${PROT_GROUP}" "${PROT_GROUP}" | ${GMX} trjconv \
    -s "${tpr_in}" \
    -f "${xtc_in}" \
    -o "${xtc_out}" \
    "${NDX_FLAG[@]}" \
    -pbc nojump \
    -center \
    -dt "${DT_PS_THIS}" \
    >/dev/null 2>&1 || {
  echo "[preproc-traj] ERROR: gmx trjconv failed for ${SDIR_REL}" >&2
  exit 2
}

# 2) Start structure at t=0 (same PBC/centering)
printf "%s\n%s\n" "${PROT_GROUP}" "${PROT_GROUP}" | ${GMX} trjconv \
    -s "${tpr_in}" \
    -f "${xtc_in}" \
    -o "${gro_out}" \
    "${NDX_FLAG[@]}" \
    -pbc nojump \
    -center \
    -dump 0 \
    >/dev/null 2>&1 || {
  echo "[preproc-traj] ERROR: gmx trjconv (start frame) failed for ${SDIR_REL}" >&2
  exit 2
}

echo "[preproc-traj] INFO: Done ${SDIR_REL}" >&2
