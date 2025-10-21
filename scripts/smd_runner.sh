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

# Get anchor residue range from config (e.g., "135-150" or "140")
anchor_res_range=$(python3 - <<PY
import yaml
cfg=yaml.safe_load(open("${ROOT}/config.yaml"))
res=None
for s in cfg["systems"]:
    if s["name"]=="${system}":
        for v in s.get("variants", []):
            if v["id"]=="${variant}":
                res=v["anchor"]["res"]
                break
print(res if res is not None else "")
PY
)
if [[ -z "${anchor_res_range}" ]]; then
  echo "[smd-runner] ERROR: anchor residue range not found in config for ${system}/${variant}" >&2
  exit 2
fi

# Create a variant-specific posre file in the RUN dir from the anchor chain ITP
POSRE_ITP="${run_root}/posre_anchor_${system}_${variant}_chain_${anchor_chain}_${anchor_res_range}.itp"
python3 - <<PY
import re, yaml
itp_path = r"""${anchor_itp}"""
out_path = r"""${POSRE_ITP}"""
rng = r"""${anchor_res_range}"""
k = float(${k_kj})
# parse range "a-b" or "x"
m=re.match(r'^\s*(\d+)\s*-\s*(\d+)\s*$', rng)
if m:
    lo,hi=sorted(map(int,m.groups()))
else:
    x=int(rng.strip()); lo=hi=x
txt=open(itp_path).read()
in_atoms=False; rows=[]
for line in txt.splitlines():
    s=line.strip()
    if not s or s.startswith(";"):
        continue
    if s.startswith("["):
        in_atoms = "atoms" in s
        continue
    if in_atoms:
        parts=s.split()
        if len(parts)>=6 and parts[0].isdigit() and parts[2].isdigit():
            ai=int(parts[0]); resnr=int(parts[2])
            if lo<=resnr<=hi:
                rows.append(ai)
if not rows:
    raise SystemExit(f"[smd-runner] ERROR: no atoms in range {rng} in {itp_path}")
with open(out_path,"w") as f:
    f.write("[ position_restraints ]\n")
    f.write("; ai  funct  fc_x   fc_y   fc_z\n")
    for ai in rows:
        f.write(f"{ai:6d}   1   {k:.3f}  {k:.3f}  {k:.3f}\n")
print(f"[smd-runner] Wrote posre: {out_path} ({len(rows)} atoms; k={k:.3f})")
PY




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
# Headerless numeric output (columns: time  x  y  z)
gmx select -s start.gro -n index.ndx -select 'com of group "Anchor"'  -os com_anchor.xvg  -xvg none
gmx select -s start.gro -n index.ndx -select 'com of group "Pulled"'  -os com_pulled.xvg  -xvg none

# For a single-frame input, these files have exactly one numeric row; grab x (col 2)
ax=$(awk 'NF>=4{print $2; exit}' com_anchor.xvg 2>/dev/null)
px=$(awk 'NF>=4{print $2; exit}' com_pulled.xvg 2>/dev/null)

if [[ -z "${ax}" || -z "${px}" ]]; then
  echo "[smd-runner] ERROR: Failed to read COM x from com_*.xvg (got ax='${ax}', px='${px}')" >&2
  exit 2
fi

# Compute dx and |dx| with fixed 6-decimal formatting
dx=$(awk -v ax="$ax" -v px="$px" 'BEGIN{printf "%.6f", (px-ax)}')
absdx=$(awk -v v="$dx" 'BEGIN{if (v<0) v=-v; printf "%.6f", v}')

# If nearly zero, warn (indices may be wrong)
awk -v v="$absdx" 'BEGIN{if (v+0.0 < 1e-6) { print "[smd-runner] WARNING: pull-coord1-init ~ 0; check Anchor/Pulled groups" > "/dev/stderr"; }}'


# nsteps from target extension and speed
speed_nm_per_ps=$(python3 - <<PY
print(float("${speed_nm_per_ns}")/1000.0)
PY
)
rate_fmt=$(python3 - <<PY
rate=float("${speed_nm_per_ns}")/1000.0
print(f"{rate:.6f}")
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
; ---------- pull.mdp (SMD) ----------
;define                  = -DPOSRES_ANCHOR   ; (unused since posre is included unconditionally in topol.top)

integrator              = md
dt                      = ${dt_ps}
nsteps                  = ${nsteps}

; ── Thermostat (two groups) ──
tcoupl                  = V-rescale
tc-grps                 = Protein NonProtein
tau_t                   = 0.5   0.5
ref_t                   = 303.15 303.15

; ── Barostat (NVT for SMD) ──
pcoupl                  = no
; If you ever enable barostat, prefer light isotropic:
;pcoupl                  = C-rescale
;pcoupltype              = isotropic
;tau_p                   = 2.0
;ref_p                   = 1.0
;compressibility         = 4.5e-5

; ── Non-bonded ──
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

; ── Constraints ──
constraints             = h-bonds
constraint-algorithm    = lincs
lincs_iter              = 1
lincs_order             = 4

; ─────────────── Pull code ───────────────
pull                    = yes
pull-ncoords            = 1
pull-ngroups            = 2

pull-group1-name        = Anchor
pull-group2-name        = Pulled

pull-coord1-type        = umbrella
pull-coord1-geometry    = direction-periodic
pull-coord1-vec         = ${vec}
pull-coord1-groups      = 1 2
pull-coord1-init        = ${absdx}
pull-coord1-rate        = ${rate_fmt}
pull-coord1-k           = ${k_kj}
pull-print-components   = yes

; ── Outputs (compact, analysis-friendly) ──
; No .trr (huge):
nstxout                 = 0
nstvout                 = 0
nstfout                 = 0
; Keep compressed coords (XTC) — every 2 ps at dt=0.002 & 1000:
nstxout-compressed      = 1000
; Energies/logs:
nstenergy               = 500
nstlog                  = 500
nstcalcenergy           = 100
; Pull outputs: use GROMACS defaults (typically 100) by omitting pull-nstxout/pull-nstfout
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


# Run grompp from the system's 00_build dir so #include paths (e.g., charmm36-*.ff) resolve correctly
BUILD_DIR="${ROOT}/systems/${system}/00_build"
if [[ ! -d "${BUILD_DIR}" ]]; then
  echo "[smd-runner] ERROR: build dir missing: ${BUILD_DIR}" >&2
  exit 2
fi

# Log what FF topol.top references (optional)
ffline=$(grep -oE '"[^"]+\.ff/forcefield\.itp"' "${BUILD_DIR}/topol.top" | head -n1 || true)
echo "[smd-runner] topol includes FF: ${ffline:-<unknown>}"

# Create run-local topol.top with absolute includes + posre include just after the anchor chain itp
BUILD_DIR="${ROOT}/systems/${system}/00_build"
SRC_TOP="${BUILD_DIR}/topol.top"
RUN_TOP="${run_root}/topol.top"

# (Optional) show what FF the source top references
ffdir_guess=$(grep -oE '"[^"]+\.ff/forcefield\.itp"' "${SRC_TOP}" | head -n1 || true)
echo "[smd-runner] topol includes FF: ${ffdir_guess:-<unknown>}"

SRC_TOP="${SRC_TOP}" BUILD_DIR="${BUILD_DIR}" RUN_TOP="${RUN_TOP}" ANCHOR_ITP="${anchor_itp}" POSRE_ITP="${POSRE_ITP}" python3 - <<'PY'
import os, re

src_top   = os.environ["SRC_TOP"]
build_dir = os.environ["BUILD_DIR"]
run_top   = os.environ["RUN_TOP"]
anchor_itp= os.environ["ANCHOR_ITP"]
posre_itp = os.environ["POSRE_ITP"]

anchor_base = os.path.basename(anchor_itp)

# 1) Rewrite #include "X" → absolute path if X exists under BUILD_DIR
txt   = open(src_top, "r").read()
lines = txt.splitlines()
rewritten = []
for line in lines:
    m = re.match(r'^\s*#\s*include\s+"([^"]+)"\s*$', line)
    if not m:
        rewritten.append(line)
        continue
    inc = m.group(1)
    candidate = os.path.join(build_dir, inc)
    if os.path.isfile(candidate):
        rewritten.append(f'#include "{candidate}"')
    else:
        # keep as-is if we can’t resolve under BUILD_DIR
        rewritten.append(line)

# 2) Insert posre include immediately after the *anchor chain* include.
#    Match by the BASENAME of the included file, so it works whether it’s relative or absolute.
out = []
inserted = False
inc_re = re.compile(r'^\s*#\s*include\s+"([^"]+)"\s*$')
for line in rewritten:
    out.append(line)
    m = inc_re.match(line)
    if (not inserted) and m:
        inc_path = m.group(1)
        if os.path.basename(inc_path) == anchor_base:
            out.append(f'#include "{posre_itp}"')
            inserted = True

if not inserted:
    # Fall back to appending (but it’s better if it goes after the anchor include)
    out.append(f'#include "{posre_itp}"')

open(run_top, "w").write("\n".join(out) + "\n")
print(f"[smd-runner] Wrote run-local topol: {run_top}")
PY



# Run grompp from BUILD_DIR (so forcefield folder resolves), but point to RUN files and RUN_TOP
pushd "${BUILD_DIR}" >/dev/null
gmx grompp \
  -f "${run_root}/pull.mdp" \
  -c "${run_root}/start.gro" \
  -r "${run_root}/start.gro" \
  -p "${run_root}/topol.top" \
  -n "${run_root}/index.ndx" \
  -o "${run_root}/pull.tpr"
popd >/dev/null




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
