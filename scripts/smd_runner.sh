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

# Link inputs
ln -sf "$(realpath "${final_tpr}")" final_npt.tpr
ln -sf "$(realpath "${final_xtc}")" final_npt.xtc

# Extract start structure
echo 0 | gmx trjconv -s final_npt.tpr -f final_npt.xtc -dump "${start_time_ps}" -o start.gro

# Ensure index with [ Anchor ] and [ Pulled ]
NDX_SRC="${ROOT}/systems/${system}/00_build/index.ndx"
if [[ -f "${NDX_SRC}" ]]; then ln -sf "${NDX_SRC}" index.ndx; fi
if ! grep -q "^\[ *Anchor *\]" index.ndx 2>/dev/null || ! grep -q "^\[ *Pulled *\]" index.ndx 2>/dev/null; then
  echo "[smd-runner] INFO: [Anchor]/[Pulled] missing; attempting to build indices..."
  if [[ -f "${ROOT}/scripts/build_indices_and_posres.py" ]]; then
    python3 "${ROOT}/scripts/build_indices_and_posres.py" --system "${system}" --variant "${variant}"
    [[ -f "${NDX_SRC}" ]] && ln -sf "${NDX_SRC}" index.ndx
  fi
fi
if ! grep -q "^\[ *Anchor *\]" index.ndx || ! grep -q "^\[ *Pulled *\]" index.ndx; then
  echo "[smd-runner] ERROR: Could not find/create index with [Anchor] and [Pulled]." >&2
  exit 2
fi

# Compute COM(Anchor->Pulled) along x
gmx select -s start.gro -n index.ndx -select 'com of group "Anchor"' -os com_anchor.xvg
gmx select -s start.gro -n index.ndx -select 'com of group "Pulled"' -os com_pulled.xvg

ax=$(awk 'NF==4{t=$1;x=$2;y=$3;z=$4} END{print x}' com_anchor.xvg 2>/dev/null || echo "NaN")
px=$(awk 'NF==4{t=$1;x=$2;y=$3;z=$4} END{print x}' com_pulled.xvg 2>/dev/null || echo "NaN")
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

# Flip vec sign to ensure extension proceeds positive along x with positive rate
if python3 - <<PY
print("neg" if float("${dx}")<0 else "pos")
PY
 | grep -q neg; then
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
gmx grompp -f pull.mdp -c start.gro -p topol.top -n index.ndx -o pull.tpr
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
