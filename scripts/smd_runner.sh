#!/usr/bin/env bash
set -euo pipefail

# Input: one CSV line (no header)
LINE="$1"

# CSV → vars (ordered per FIELDNAMES in smd_build_manifest.py)
IFS=',' read -r system variant replicate speed_nm_per_ns k_kj dt_ps target_ext_nm axis perf_ns_per_day start_time_ps final_tpr final_xtc start_id anchor_chain array_cap <<< "$LINE"

# Resolve repo root by crawling up to find config.yaml
find_root() {
  local p="$PWD"
  while [[ "$p" != "/" ]]; do
    if [[ -f "$p/config.yaml" ]]; then echo "$p"; return 0; fi
    p="$(dirname "$p")"
  done
  return 1
}
ROOT="$(find_root)"
if [[ -z "$ROOT" ]]; then
  echo "[smd-runner] ERROR: cannot locate repo root (config.yaml)" >&2
  exit 2
fi

# Output directory
run_root="$ROOT/smd/${system}/${variant}/v$(printf "%.3f" "${speed_nm_per_ns}")/rep${replicate}/${start_id}"
mkdir -p "${run_root}"
cd "${run_root}"

echo "[smd-runner] System=${system} Variant=${variant} Rep=${replicate} Speed=${speed_nm_per_ns} nm/ns Start=${start_time_ps} ps"

# Anchor ITP path from template in config
anchor_itp=$(python3 - <<PY
import yaml, os
cfg = yaml.safe_load(open("${ROOT}/config.yaml"))
tpl = cfg["globals"]["smd"]["anchor_chain_template"]
print(os.path.join("${ROOT}", tpl.format(system="${system}", chain="${anchor_chain}")))
PY
)
if [[ ! -f "${anchor_itp}" ]]; then
  echo "[smd-runner] ERROR: anchor ITP not found: ${anchor_itp}" >&2
  exit 2
fi

# Resolve final_tpr/final_xtc relative to repo root if they are not absolute
resolve_path() {
  case "$1" in
    /*) printf "%s\n" "$1" ;;
    *)  printf "%s/%s\n" "${ROOT}" "$1" ;;
  esac
}
TPR_ABS="$(resolve_path "${final_tpr}")"
XTC_ABS="$(resolve_path "${final_xtc}")"

if [[ ! -f "${TPR_ABS}" ]]; then
  echo "[smd-runner] ERROR: final_tpr not found at ${TPR_ABS}" >&2
  exit 2
fi
if [[ ! -f "${XTC_ABS}" ]]; then
  echo "[smd-runner] ERROR: final_xtc not found at ${XTC_ABS}" >&2
  exit 2
fi

ln -sf "${TPR_ABS}" final_npt.tpr
ln -sf "${XTC_ABS}" final_npt.xtc


# Extract start structure
echo 0 | gmx trjconv -s final_npt.tpr -f final_npt.xtc -dump "${start_time_ps}" -o start.gro

# Ensure index with [ Anchor ] and [ Pulled ]
NDX_SRC="${ROOT}/systems/${system}/00_build/index.ndx"
if [[ -f "${NDX_SRC}" ]]; then ln -sf "${NDX_SRC}" index.ndx; fi
if ! grep -q "^\[ *Anchor *\]" index.ndx 2>/dev/null || ! grep -q "^\[ *Pulled *\]" index.ndx 2>/dev/null; then
  echo "[smd-runner] INFO: [Anchor]/[Pulled] missing; attempting to build indices..."
  if [[ -f "${ROOT}/scripts/build_indices_and_posres.py" ]]; then
  echo "[smd-runner] INFO: running build_indices_and_posres.py ${system} ${variant}"

  # Make sure prints flush immediately, and prevent silent indefinite hangs.
  # Requires coreutils `timeout` (present on most clusters).
  export PYTHONUNBUFFERED=1
  export GMX_MAXBACKUP=-1

  # Run with a 5-minute timeout; stream stdout/stderr line-by-line.
  set +e
  timeout 300s stdbuf -oL -eL python3 "${ROOT}/scripts/build_indices_and_posres.py" "${system}" "${variant}"
  rc=$?
  set -e

  if [[ $rc -eq 124 ]]; then
    echo "[smd-runner] ERROR: build_indices_and_posres.py timed out after 300s."
    echo "[smd-runner] INFO: Active gmx processes at timeout:"
    ps -o pid,etime,cmd -u "$USER" | grep -E "gmx( |$)" | grep -v grep || echo "[smd-runner] (none)"
    exit 2
  elif [[ $rc -ne 0 ]]; then
    echo "[smd-runner] ERROR: build_indices_and_posres.py exited with rc=${rc}" >&2
    exit $rc
  fi

  # try again
  NDX_SRC="${ROOT}/systems/${system}/00_build/index.ndx"
  if [[ -f "${NDX_SRC}" ]]; then
    ln -sf "${NDX_SRC}" index.ndx
  fi
else
  echo "[smd-runner] ERROR: ${ROOT}/scripts/build_indices_and_posres.py not found." >&2
fi


fi
if ! grep -q "^\[ *Anchor *\]" index.ndx || ! grep -q "^\[ *Pulled *\]" index.ndx; then
  echo "[smd-runner] ERROR: Could not find/create index with [Anchor] and [Pulled]." >&2
  exit 2
fi

# Compute COM(Anchor->Pulled) along x
# Write plain numeric output (no xvg headers/legends)
gmx select -s start.gro -n index.ndx -select 'com of group "Anchor"'  -os com_anchor.xvg  -xvg none
gmx select -s start.gro -n index.ndx -select 'com of group "Pulled"'  -os com_pulled.xvg  -xvg none


# Files now have numeric rows only; columns are: time  x  y  z
ax=$(awk 'NF>=2{last=$2} END{print last}' com_anchor.xvg 2>/dev/null)
px=$(awk 'NF>=2{last=$2} END{print last}' com_pulled.xvg 2>/dev/null)
# Fallback if empty
[[ -z "${ax}" ]] && ax="NaN"
[[ -z "${px}" ]] && px="NaN"

if [[ "${ax}" == "NaN" || "${px}" == "NaN" ]]; then
  echo "[smd-runner] ERROR: Failed to parse COM x-coordinates." >&2
  exit 2
fi

dx=$(python3 - <<PY
ax=float("${ax}"); px=float("${px}")
print(px-ax)
PY
)
absdx=$(python3 - <<PY
import math
print(abs(float("${dx}")))
PY
)

# nsteps from target extension and speed
speed_nm_per_ps=$(python3 - <<PY
print(float("${speed_nm_per_ns}")/1000.0)
PY
)
nsteps=$(python3 - <<PY
dt=float("${dt_ps}"); target=float("${target_ext_nm}"); rate=float("${speed_nm_per_ps}")
print(int(round( (target / rate) / dt )))
PY
)

# Choose vec sign so extension is positive with positive rate along x
sign=$(python3 - <<PY
dx=float("${dx}")
print("neg" if dx < 0 else "pos")
PY
)
if [[ "${sign}" == "neg" ]]; then
  vec="-1 0 0"
else
  vec="1 0 0"
fi


# Write pull.mdp
cat > pull.mdp <<MDP
integrator              = md
dt                      = ${dt_ps}
nsteps                  = ${nsteps}
nstxout-compressed      = 1000
nstenergy               = 1000
nstlog                  = 1000
continuation            = no
constraints             = h-bonds
constraint_algorithm    = lincs
lincs_iter              = 1
lincs_order             = 4
tcoupl                  = v-rescale
tc-grps                 = System
tau_t                   = 1.0
ref_t                   = 303
pcoupl                  = Parrinello-Rahman
pcoupltype              = isotropic
tau_p                   = 5.0
ref_p                   = 1.0
compressibility         = 4.5e-5

pull                    = yes
pull-ncoords            = 1
pull-ngroups            = 2
pull-group1-name        = Anchor
pull-group2-name        = Pulled
pull-coord1-groups      = 1 2
pull-coord1-type        = umbrella
pull-coord1-geometry    = direction-periodic
pull-coord1-vec         = ${vec}
pull-coord1-init        = ${absdx}
pull-coord1-k           = ${k_kj}
pull-coord1-rate        = ${speed_nm_per_ps}
pull-print-components   = yes
pull-nstxout            = 100
pull-nstfout            = 100
MDP

# Topology (prefer system master top)
TOP_SRC="${ROOT}/systems/${system}/00_build/topol.top"
if [[ -f "${TOP_SRC}" ]]; then
  cp "${TOP_SRC}" topol.top
  if ! grep -q "$(basename "${anchor_itp}")" topol.top; then
    echo "#include \"$(realpath --relative-to=. "${anchor_itp}")\"" >> topol.top
  fi
else
  echo "[smd-runner] WARN: ${TOP_SRC} not found — proceeding with local topology if present." >&2
fi

# Grompp + MDrun
# Tell the preprocessor where to find the forcefield folder referenced by topol.top
INC="-I ${ROOT}/systems/${system}/00_build"
if [[ ! -d "${ROOT}/systems/${system}/00_build" ]]; then
  echo "[smd-runner] WARN: include dir not found: ${ROOT}/systems/${system}/00_build"
fi

ffdir_guess=$(grep -oE '"[^"]+\.ff/forcefield\.itp"' topol.top | sed -E 's/^"([^"]+)\/forcefield\.itp"$/\1/')
echo "[smd-runner] topol includes FF: ${ffdir_guess:-<unknown>}"


gmx grompp ${INC} -f pull.mdp -c start.gro -p topol.top -n index.ndx -o pull.tpr


# >>> ADD THIS BLOCK (DRY RUN SWITCH) <<<
if [[ "${DRY_RUN:-}" == "1" ]]; then
  echo "[smd-runner] DRY RUN: grompp succeeded; skipping mdrun."
  exit 0
fi
# >>> END ADD <<<
gmx mdrun  -deffnm pull -ntmpi 1


# Log resolved parameters
python3 - <<PY
import json
d = {
  "system":"${system}","variant":"${variant}","replicate":int("${replicate}"),
  "speed_nm_per_ns":float("${speed_nm_per_ns}"),"dt_ps":float("${dt_ps}"),
  "k_kj_mol_nm2":float("${k_kj}"),"target_extension_nm":float("${target_ext_nm}"),
  "axis":"${axis}","start_time_ps":float("${start_time_ps}"),
  "anchor_itp":"${anchor_itp}","dx_x":float("${dx}"),
  "pull_vec":"${vec}","nsteps":int("${nsteps}")
}
open("run.json","w").write(json.dumps(d, indent=2))
PY
