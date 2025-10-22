#!/usr/bin/env bash
set -euo pipefail

# Input: one CSV line (no header)
LINE="$1"

############################
# Parse manifest CSV (positional)
# Header order:
# system,variant,replicate,speed_nm_per_ns,k_kj_mol_nm2,dt_ps,target_extension_nm,axis,perf_ns_per_day,start_time_ps,final_tpr,final_xtc,start_id,anchor_chain,array_cap
############################
IFS=',' read -r \
  system \
  variant \
  replicate \
  speed_nm_per_ns \
  k_kj \
  dt_ps \
  target_extension_nm \
  axis \
  perf_ns_per_day \
  start_time_ps \
  final_tpr \
  final_xtc \
  start_id \
  anchor_chain \
  array_cap \
  <<< "${LINE}"

: "${system:?}"; : "${variant:?}"; : "${replicate:?}"
: "${speed_nm_per_ns:?}"; : "${k_kj:?}"; : "${dt_ps:?}"
: "${target_extension_nm:?}"; : "${axis:?}"; : "${start_time_ps:?}"
: "${final_tpr:?}"; : "${final_xtc:?}"; : "${start_id:?}"; : "${anchor_chain:?}"

# strip any accidental quotes
stripq(){ local v="${1}"; v="${v#\"}"; v="${v%\"}"; printf "%s" "${v}"; }
replicate="$(stripq "${replicate}")"
speed_nm_per_ns="$(stripq "${speed_nm_per_ns}")"
k_kj="$(stripq "${k_kj}")"
dt_ps="$(stripq "${dt_ps}")"
target_extension_nm="$(stripq "${target_extension_nm}")"
start_time_ps="$(stripq "${start_time_ps}")"
perf_ns_per_day="$(stripq "${perf_ns_per_day}")"
array_cap="$(stripq "${array_cap}")"

############################
# Locate repo root (config.yaml)
############################
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
  echo "[smd-runner] ERROR: cannot locate repo root (config.yaml)" >&2
  exit 2
fi

############################
# Paths and run directory
############################
# Run dir: smd/<sys>/<var>/v<speed>/rep<rep>/<start_id>
run_root="$ROOT/smd/${system}/${variant}/v$(printf "%.3f" "${speed_nm_per_ns}")/rep${replicate}/${start_id}"
mkdir -p "${run_root}"
cd "${run_root}"

echo "[smd-runner] System=${system} Variant=${variant} Rep=${replicate} Speed=${speed_nm_per_ns} nm/ns Start=${start_time_ps} ps"

# Resolve build dir once (FF/top live here)
BUILD_DIR="${ROOT}/systems/${system}/00_build"
if [[ ! -d "${BUILD_DIR}" ]]; then
  echo "[smd-runner] ERROR: build dir missing: ${BUILD_DIR}" >&2
  exit 2
fi

############################
# Anchor ITP and anchor residue range from config
############################
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

anchor_res_range=$(python3 - <<PY
import yaml
cfg=yaml.safe_load(open("${ROOT}/config.yaml"))
res=None
for s in cfg["systems"]:
    if s["name"]=="${system}":
        for v in s.get("variants", []):
            if v["id"]=="${variant}":
                res=v["anchor"]["res"]; break
print(res if res is not None else "")
PY
)
if [[ -z "${anchor_res_range}" ]]; then
  echo "[smd-runner] ERROR: anchor residue range not found in config for ${system}/${variant}" >&2
  exit 2
fi

############################
# Variant-specific posre (run-local)
############################
POSRE_ITP="${run_root}/posre_anchor_${system}_${variant}_chain_${anchor_chain}_${anchor_res_range}.itp"
python3 - <<PY
import re
itp_path = r"""${anchor_itp}"""
out_path = r"""${POSRE_ITP}"""
rng = r"""${anchor_res_range}"""
k = float(${k_kj})
m=re.match(r'^\s*(\d+)\s*-\s*(\d+)\s*$', rng)
if m:
    lo,hi=sorted(map(int,m.groups()))
else:
    x=int(rng.strip()); lo=hi=x
txt=open(itp_path).read()
in_atoms=False; rows=[]
for line in txt.splitlines():
    s=line.strip()
    if not s or s.startswith(";"): continue
    if s.startswith("["):
        in_atoms = "atoms" in s
        continue
    if in_atoms:
        parts=s.split()
        if len(parts)>=6 and parts[0].isdigit() and parts[2].isdigit():
            ai=int(parts[0]); resnr=int(parts[2])
            if lo<=resnr<=hi: rows.append(ai)
if not rows:
    raise SystemExit(f"[smd-runner] ERROR: no atoms in range {rng} in {itp_path}")
with open(out_path,"w") as f:
    f.write("[ position_restraints ]\n")
    f.write("; ai  funct  fc_x   fc_y   fc_z\n")
    for ai in rows:
        f.write(f"{ai:6d}   1   {k:.3f}  {k:.3f}  {k:.3f}\n")
print(f"[smd-runner] Wrote posre: {out_path} ({len(rows)} atoms; k={k:.3f})")
PY

############################
# Link final NPT files and dump start structure
############################
resolve_path() { case "$1" in /*) printf "%s\n" "$1" ;; *) printf "%s/%s\n" "${ROOT}" "$1" ;; esac; }
TPR_ABS="$(resolve_path "${final_tpr}")"
XTC_ABS="$(resolve_path "${final_xtc}")"
[[ -f "${TPR_ABS}" ]] || { echo "[smd-runner] ERROR: final_tpr not found: ${TPR_ABS}" >&2; exit 2; }
[[ -f "${XTC_ABS}" ]] || { echo "[smd-runner] ERROR: final_xtc not found: ${XTC_ABS}" >&2; exit 2; }
ln -sf "${TPR_ABS}" final_npt.tpr
ln -sf "${XTC_ABS}" final_npt.xtc

echo 0 | gmx trjconv -s final_npt.tpr -f final_npt.xtc -dump "${start_time_ps}" -o start.gro

############################
# Ensure index.ndx has [ Anchor ] and [ Pulled ] (build if needed)
############################
if [[ -f "${BUILD_DIR}/index.ndx" ]]; then
  ln -sf "${BUILD_DIR}/index.ndx" index.ndx
fi
if ! grep -q "^\[ *Anchor *\]" index.ndx 2>/dev/null || ! grep -q "^\[ *Pulled *\]" index.ndx 2>/dev/null; then
  echo "[smd-runner] INFO: [Anchor]/[Pulled] missing; building indices..."
  export PYTHONUNBUFFERED=1 GMX_MAXBACKUP=-1
  timeout 300s python3 "${ROOT}/scripts/build_indices_and_posres.py" "${system}" "${variant}" || {
    rc=$?
    if [[ $rc -eq 124 ]]; then
      echo "[smd-runner] ERROR: index builder timed out after 300s." >&2
    else
      echo "[smd-runner] ERROR: index builder failed rc=${rc}." >&2
    fi
    exit 2
  }
  ln -sf "${BUILD_DIR}/index.ndx" index.ndx
fi
if ! grep -q "^\[ *Anchor *\]" index.ndx || ! grep -q "^\[ *Pulled *\]" index.ndx; then
  echo "[smd-runner] ERROR: Could not find/create index with [ Anchor ] and [ Pulled ]." >&2
  exit 2
fi

############################
# Build run-local topol.top (absolute includes) + insert posre include after anchor chain include
############################
SRC_TOP="${BUILD_DIR}/topol.top"
RUN_TOP="${run_root}/topol.top"
SRC_TOP="${SRC_TOP}" BUILD_DIR="${BUILD_DIR}" RUN_TOP="${RUN_TOP}" ANCHOR_ITP="${anchor_itp}" POSRE_ITP="${POSRE_ITP}" python3 - <<'PY'
import os, re
src_top   = os.environ["SRC_TOP"]
build_dir = os.environ["BUILD_DIR"]
run_top   = os.environ["RUN_TOP"]
anchor_itp= os.environ["ANCHOR_ITP"]
posre_itp = os.environ["POSRE_ITP"]
anchor_base = os.path.basename(anchor_itp)

txt=open(src_top,"r").read()
lines=txt.splitlines()
rewritten=[]
inc_re = re.compile(r'^\s*#\s*include\s+"([^"]+)"\s*$')
for line in lines:
    m = inc_re.match(line)
    if not m:
        rewritten.append(line); continue
    inc = m.group(1)
    cand = os.path.join(build_dir, inc)
    rewritten.append(f'#include "{cand}"' if os.path.isfile(cand) else line)

out=[]; inserted=False
for line in rewritten:
    out.append(line)
    m = inc_re.match(line)
    if (not inserted) and m:
        inc_path = m.group(1)
        if os.path.basename(inc_path) == anchor_base:
            out.append(f'#include "{posre_itp}"'); inserted=True
if not inserted:
    out.append(f'#include "{posre_itp}"')
open(run_top,"w").write("\n".join(out)+"\n")
print(f"[smd-runner] Wrote run-local topol: {run_top}")
PY

############################
# Probe grompp to get SIGNED initial distance; derive signed rate
############################
vec="1 0 0"  # always +x; sign is carried by init and rate

# Base (positive) rate in nm/ps, fixed decimal
base_rate=$(python3 - <<PY
rate = float("${speed_nm_per_ns}") / 1000.0
print("{:.6f}".format(rate))
PY
)

# Minimal probe mdp
cat > pull_probe.mdp <<PROBE
integrator              = md
dt                      = ${dt_ps}
nsteps                  = 10
tcoupl                  = no
pcoupl                  = no
cutoff-scheme           = Verlet
coulombtype             = PME
rcoulomb                = 1.2
vdwtype                 = cutoff
vdw-modifier            = Force-switch
rlist                   = 1.2
rvdw                    = 1.2
rvdw-switch             = 1.0
pbc                     = xyz
constraints             = h-bonds
constraint-algorithm    = lincs

pull                    = yes
pull-ncoords            = 1
pull-ngroups            = 2
pull-group1-name        = Anchor
pull-group2-name        = Pulled
pull-coord1-type        = umbrella
pull-coord1-geometry    = direction-periodic
pull-coord1-vec         = ${vec}
pull-coord1-groups      = 1 2
pull-coord1-init        = 0.000000
pull-coord1-rate        = ${base_rate}
PROBE

# Run probe grompp from BUILD_DIR; capture log
pushd "${BUILD_DIR}" >/dev/null
gmx grompp \
  -f "${run_root}/pull_probe.mdp" \
  -c "${run_root}/start.gro" \
  -r "${run_root}/start.gro" \
  -p "${run_root}/topol.top" \
  -n "${run_root}/index.ndx" \
  -o "${run_root}/probe.tpr" \
  2>&1 | tee "${run_root}/probe.grompp.log"
popd >/dev/null

# Extract signed "distance at start" for group 2
signed_init=$(awk '/^ *2[[:space:]]/{for(i=1;i<=NF;i++){if($i=="nm"){print $(i-1); exit}}}' "${run_root}/probe.grompp.log")
if [[ -z "${signed_init}" ]]; then
  signed_init=$(awk '/distance at start/{for(i=1;i<=NF;i++){if($i=="nm"){print $(i-1); exit}}}' "${run_root}/probe.grompp.log")
fi
if [[ -z "${signed_init}" ]]; then
  echo "[smd-runner] ERROR: could not extract signed start distance from probe.grompp.log" >&2
  exit 2
fi

init_signed=$(python3 - <<PY
v=float("${signed_init}")
print("{:.6f}".format(v))
PY
)
rate_signed=$(python3 - <<PY
base=float("${base_rate}"); v=float("${signed_init}")
print("{:.6f}".format(-base if v<0.0 else base))
PY
)
echo "[smd-runner] grompp start distance (signed): ${init_signed} nm; pull rate: ${rate_signed} nm/ps"

############################
# Steps based on magnitude only (uses base_rate)
############################
total_time_ps=$(python3 - <<PY
speed=float("${base_rate}")     # nm/ps
t_ns = float("${target_extension_nm}")/speed
print("{:.3f}".format(t_ns*1000.0))
PY
)
nsteps=$(python3 - <<PY
import math
dt=float("${dt_ps}")
Tps=float("${total_time_ps}")
print(int(math.ceil(Tps/dt)))
PY
)
if ! [[ "${nsteps}" =~ ^[0-9]+$ ]] || [[ "${nsteps}" -le 0 ]]; then
  echo "[smd-runner] ERROR: bad nsteps='${nsteps}' (dt_ps=${dt_ps}, total_time_ps=${total_time_ps})" >&2
  exit 2
fi

############################
# Final pull.mdp (Option B: ALL atoms, 10 ps stride)
############################
cat > pull.mdp <<MDP
; ---------- pull.mdp (SMD) ----------
integrator              = md
dt                      = ${dt_ps}
nsteps                  = ${nsteps}

; Thermostat
tcoupl                  = V-rescale
tc-grps                 = Protein NonProtein
tau_t                   = 0.5   0.5
ref_t                   = 303.15 303.15

; Barostat (NVT for SMD)
pcoupl                  = no

; Non-bonded
cutoff-scheme           = Verlet
coulombtype             = PME
rcoulomb                = 1.2
vdwtype                 = cutoff
vdw-modifier            = Force-switch
rlist                   = 1.2
rvdw                    = 1.2
rvdw-switch             = 1.0
DispCorr                = no
pbc                     = xyz

; Constraints
constraints             = h-bonds
constraint-algorithm    = lincs
lincs_iter              = 1
lincs_order             = 4

; Pull code
pull                    = yes
pull-ncoords            = 1
pull-ngroups            = 2
pull-group1-name        = Anchor
pull-group2-name        = Pulled
pull-coord1-type        = umbrella
pull-coord1-geometry    = direction-periodic
pull-coord1-vec         = ${vec}          ; fixed axis (+x)
pull-coord1-groups      = 1 2
pull-coord1-init        = ${init_signed}  ; SIGNED from probe grompp
pull-coord1-rate        = ${rate_signed}  ; SIGNED, same sign as init
pull-coord1-k           = ${k_kj}
pull-print-components   = yes

; Output control (ALL atoms)
nstxout                 = 0
nstvout                 = 0
nstfout                 = 0
; 10 ps stride at dt=0.002 → every 5000 steps
nstxout-compressed      = 5000
nstenergy               = 5000
nstlog                  = 5000
nstcalcenergy           = 5000
; Pull traces: keep GROMACS defaults (~100 steps) by omitting pull-nstxout/pull-nstfout
MDP

############################
# Final grompp (from BUILD_DIR so FF includes resolve)
############################
pushd "${BUILD_DIR}" >/dev/null
gmx grompp \
  -f "${run_root}/pull.mdp" \
  -c "${run_root}/start.gro" \
  -r "${run_root}/start.gro" \
  -p "${run_root}/topol.top" \
  -n "${run_root}/index.ndx" \
  -o "${run_root}/pull.tpr"
popd >/dev/null

# Dry-run switch
if [[ "${DRY_RUN:-}" == "1" ]]; then
  echo "[smd-runner] DRY RUN: grompp succeeded; skipping mdrun."
  exit 0
fi

# Production mdrun (single-rank launch here; Slurm script will map resources)
gmx mdrun -deffnm pull -ntmpi 1

############################
# Minimal run.json (for auditing)
############################
python3 - <<PY
import json
d = {
  "system": "${system}",
  "variant": "${variant}",
  "replicate": int("${replicate}"),
  "start_time_ps": float("${start_time_ps}"),
  "axis": "${axis}",
  "speed_nm_per_ns": float("${speed_nm_per_ns}"),
  "rate_nm_per_ps_signed": float("${rate_signed}"),
  "dt_ps": float("${dt_ps}"),
  "k_kj_mol_nm2": float("${k_kj}"),
  "target_extension_nm": float("${target_extension_nm}"),
  "init_nm_signed": float("${init_signed}"),
  "nsteps": int("${nsteps}"),
  "topol": "topol.top",
  "posre_itp": "${POSRE_ITP}",
  "index": "index.ndx"
}
open("run.json","w").write(json.dumps(d, indent=2))
PY
