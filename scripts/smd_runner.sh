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
if [[ -n "${GMX_MODULE:-}" ]]; then
  echo "[smd-runner] GMX_MODULE=${GMX_MODULE}"
fi
if [[ -n "${SLURM_JOB_ID:-}" ]]; then
  echo "[smd-runner] SLURM job: ${SLURM_JOB_ID} task=${SLURM_ARRAY_TASK_ID:-} cpus=${SLURM_CPUS_PER_TASK:-} gres=${SLURM_JOB_GRES:-}"
fi

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
SRC_TOP="${SRC_TOP}" BUILD_DIR="${BUILD_DIR}" RUN_TOP="${RUN_TOP}" ANCHOR_ITP="${anchor_itp}" POSRE_ITP="${POSRE_ITP}" ROOT_DIR="${ROOT}" python3 - <<'PY'
import os, re
src_top   = os.environ["SRC_TOP"]
build_dir = os.environ["BUILD_DIR"]
run_top   = os.environ["RUN_TOP"]
anchor_itp= os.environ["ANCHOR_ITP"]
posre_itp = os.environ["POSRE_ITP"]
root_dir  = os.environ["ROOT_DIR"]
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
    if os.path.isfile(cand):
        # If the include exists in build_dir, make it absolute to avoid CWD sensitivity
        rewritten.append(f'#include "{os.path.abspath(cand)}"')
    else:
        # Try resolving relative include against repo root
        inc_clean = inc[2:] if inc.startswith("./") else inc
        root_cand = os.path.join(root_dir, inc_clean)
        if os.path.isfile(root_cand):
            rewritten.append(f'#include "{os.path.abspath(root_cand)}"')
        else:
            # Last resort: keep line but normalize leading ./ to avoid odd relative prefixes
            if inc.startswith("./"):
                rewritten.append(f'#include "{inc_clean}"')
            else:
                rewritten.append(line)

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

# No probe grompp: let GROMACS compute the initial distance at grompp by
# enabling pull-coord1-start in the mdp. This removes an expensive probe
# grompp and avoids parsing logs for a signed distance. We still compute
# nsteps based on the (positive) base_rate and the target extension magnitude.

# Base (positive) rate in nm/ps, fixed decimal
base_rate=$(python3 - <<PY
rate = float("${speed_nm_per_ns}") / 1000.0
print("{:.6f}".format(rate))
PY
)

# Pull vector (direction): keep +x as before; GROMACS will determine the
# signed initial distance when pull-coord1-start = yes is used in the mdp.
vec="1 0 0"

############################
# Steps based on magnitude only (uses base_rate)
############################
# --- Steps based on target extension and speed (CORRECTED UNITS) ---
# base_rate is nm/ps (we already divided by 1000 earlier)
total_time_ps=$(python3 - <<PY
speed_nm_per_ps=float("${base_rate}")   # nm/ps
target_nm=float("${target_extension_nm}")
t_ps = target_nm / speed_nm_per_ps      # (nm) / (nm/ps) = ps
print("{:.3f}".format(t_ps))
PY
)

nsteps=$(python3 - <<PY
import math
dt_ps=float("${dt_ps}")
t_ps=float("${total_time_ps}")
print(int(math.ceil(t_ps/dt_ps)))
PY
)

if ! [[ "${nsteps}" =~ ^[0-9]+$ ]] || [[ "${nsteps}" -le 0 ]]; then
  echo "[smd-runner] ERROR: bad nsteps='${nsteps}' (dt_ps=${dt_ps}, total_time_ps=${total_time_ps})" >&2
  exit 2
fi

# Save expected nsteps for completion checking
echo "${nsteps}" > expected_nsteps.txt
echo "[smd-runner] Expected nsteps: ${nsteps} (saved to expected_nsteps.txt)"


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
pull-coord1-start      = yes               ; let grompp compute initial distance
pull-coord1-rate        = ${base_rate}     ; positive base rate (nm/ps)
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
nstcalcenergy           = 100
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

# --- mdrun parameters from config + maxh from Slurm timelimit ---
readarray -t MDR < <(python3 - <<PY
import yaml
cfg=yaml.safe_load(open(r"${ROOT}/config.yaml"))
m=cfg.get("globals",{}).get("smd",{}).get("mdrun",{})
def g(k,default=""):
    v=m.get(k,default)
    return "" if v is None else str(v)
print(g("ntmpi",""))       # 0
print(g("ntomp",""))       # 1
print(g("nb",""))          # 2
print(g("bonded",""))      # 3
print(g("pme",""))         # 4
print(g("update",""))      # 5
print(g("npme",""))        # 6
print(g("ntomp_pme",""))   # 7
PY
)
NTMPI="${MDR[0]}"; NTOMP="${MDR[1]}"; NB="${MDR[2]}"; BONDED="${MDR[3]}"; PME="${MDR[4]}"; UPDATE="${MDR[5]}"; NPME="${MDR[6]}"; NTOMP_PME="${MDR[7]}"

## Build args list (runner supplies -v and -deffnm pull; submitters set only the base command)
# We'll decide CPU/GPU mode after determining GMX_CMD and then append flags accordingly.

# Decide mdrun command: CPU mode can export GMX_CMD="srun gmx_mpi mdrun"
# Default to standard gmx mdrun if not provided by submitter
GMX_CMD="${GMX_CMD:-gmx mdrun}"

# Determine mode from GMX_CMD
CPU_MODE=0
if [[ "${GMX_CMD}" == *"gmx_mpi mdrun"* ]]; then
  CPU_MODE=1
fi

# Base mdargs
mdargs=( -v -deffnm pull )

# Check for checkpoint (restart scenario)
if [[ -f "pull.cpt" ]]; then
  echo "[smd-runner] Found checkpoint pull.cpt; continuing from restart"
  # For SMD/pull code, also ensure we continue writing pull traces consistently
  mdargs+=( -cpi pull.cpt -px pull_pullx.xvg -pf pull_pullf.xvg )
  # Mark this as a restart (increment counter)
  RESTART_NUM="${RESTART_COUNT:-0}"
  touch ".restart_${RESTART_NUM}"
  echo "[smd-runner] Created restart marker: .restart_${RESTART_NUM}"
fi

if [[ ${CPU_MODE} -eq 1 ]]; then
  echo "[smd-runner] Detected CPU mode (GMX_CMD contains gmx_mpi). Omitting -ntmpi/-ntomp and GPU offload flags."
  # Optionally set OpenMP threads if not provided; align with Slurm cpus-per-task
  if [[ -z "${OMP_NUM_THREADS:-}" && -n "${SLURM_CPUS_PER_TASK:-}" ]]; then
    export OMP_NUM_THREADS="${SLURM_CPUS_PER_TASK}"
    echo "[smd-runner] OMP_NUM_THREADS=${OMP_NUM_THREADS} (from SLURM_CPUS_PER_TASK)"
  fi
  # Do NOT append -ntmpi/-ntomp or any GPU offload flags in CPU mode
else
  # GPU/default mode: honor config flags
  [[ -n "${NTMPI}"      ]] && mdargs+=( -ntmpi "${NTMPI}" )
  [[ -n "${NTOMP}"      ]] && mdargs+=( -ntomp "${NTOMP}" )
  [[ -n "${NB}"         ]] && mdargs+=( -nb "${NB}" )
  [[ -n "${BONDED}"     ]] && mdargs+=( -bonded "${BONDED}" )
  [[ -n "${PME}"        ]] && mdargs+=( -pme "${PME}" )
  [[ -n "${UPDATE}"     ]] && mdargs+=( -update "${UPDATE}" )
  [[ -n "${NPME}"       ]] && mdargs+=( -npme "${NPME}" )
  [[ -n "${NTOMP_PME}"  ]] && mdargs+=( -ntomp_pme "${NTOMP_PME}" )
fi

# -maxh from job script (computed or exported by submitter) with a buffer for resubmission logic
if [[ -n "${MAXH_HOURS:-}" ]]; then
  # Subtract 10 min (0.167 hrs) to leave time for completion check and resubmit
  MAXH_ADJUSTED=$(python3 - <<PY
h=float("${MAXH_HOURS}")
adjusted=max(0.5, h-0.167)
print(f"{adjusted:.2f}")
PY
)
  mdargs+=( -maxh "${MAXH_ADJUSTED}" )
  echo "[smd-runner] Using -maxh ${MAXH_ADJUSTED} (original: ${MAXH_HOURS}, buffer: 10 min)"
fi

# Dry-run switch
if [[ "${DRY_RUN:-}" == "1" ]]; then
  echo "[smd-runner] DRY RUN: would run: ${GMX_CMD} with args: ${mdargs[*]}"
  exit 0
fi

# Log command for diagnostics
echo "[smd-runner] Executing: ${GMX_CMD} ${mdargs[*]}"

# Build command string with proper quoting and execute
cmd="${GMX_CMD}"
for a in "${mdargs[@]}"; do
  cmd+=" "$(printf %q "$a")
done
eval "$cmd"


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
  "rate_nm_per_ps": float("${base_rate}"),
  "dt_ps": float("${dt_ps}"),
  "k_kj_mol_nm2": float("${k_kj}"),
  "target_extension_nm": float("${target_extension_nm}"),
  "nsteps": int("${nsteps}"),
  "topol": "topol.top",
  "posre_itp": "${POSRE_ITP}",
  "index": "index.ndx"
}
open("run.json","w").write(json.dumps(d, indent=2))
PY
