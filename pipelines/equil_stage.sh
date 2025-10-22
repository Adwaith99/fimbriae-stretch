#!/usr/bin/env bash
# MODE: em|nvt|npt
# FC: restraint level (0 => none)
# FINALFLAG: 0 normal stage, 1 => final NPT (equilibrium MD) using dt_ps from config
# BOX_X/Y/Z: box dimensions in nm (passed from submitter)
set -euo pipefail

SYS=${1:?}; MODE=${2:?}; TEMPK=${3:?}; NSTEPS=${4:?}; FC=${5:?}; FINAL=${6:-0}
BOX_X=${7:-}; BOX_Y=${8:-}; BOX_Z=${9:-}

ROOT="$(cd "$(dirname "$0")/.."; pwd)"
export PIPE_ROOT="$ROOT"
WD="$ROOT/systems/${SYS}/00_build"
cd "$WD"

# ---------- helpers ----------
pick_prev_conf() {
  for c in \
    npt_final.gro \
    npt_fc1000.gro npt_fc500.gro npt_fc200.gro npt_fc100.gro npt_fc50.gro \
    nvt_fc1000.gro nvt_fc500.gro nvt_fc200.gro nvt_fc100.gro nvt_fc50.gro \
    em_all.gro em_posres.gro em.gro ions.gro solv.gro boxed.gro aligned.pdb clean.pdb
  do [[ -f "$c" ]] && { echo "$c"; return; }; done
  echo "clean.pdb"
}

# Build GPU flags deterministically; robust GPU/CPU detection
build_mdrun_flags() {
  if [[ -n "${GMX_MDRUN_FLAGS:-}" ]]; then echo "${GMX_MDRUN_FLAGS}"; return; fi
  _int_last(){ local n; n="$(echo "$1" | sed -E 's/.*[^0-9]([0-9]+)[^0-9]*$/\1/;t;d')"||true; [[ -z "$n" ]]&&echo 0||echo "$n"; }
  local RAW_GPUS="${SLURM_GPUS_PER_NODE:-${SLURM_GPUS_ON_NODE:-${SLURM_GPUS:-}}}"
  [[ -z "$RAW_GPUS" ]] && RAW_GPUS="$(python3 - <<'PY'
import yaml, os; print(yaml.safe_load(open(os.environ["PIPE_ROOT"] + "/config.yaml"))["globals"]["slurm"]["gpus_per_node"])
PY
)"
  local GPUS="$(_int_last "$RAW_GPUS")"; (( GPUS<1 )) && GPUS=1
  local CPUS="$(_int_last "${SLURM_CPUS_PER_TASK:-0}")"
  if (( CPUS<1 )) && command -v taskset >/dev/null 2>&1; then
    local aff; aff="$(taskset -pc $$ 2>/dev/null | awk -F': ' 'NR==1{print $2}')"||aff=""
    if [[ -n "$aff" ]]; then
      local cnt=0 part; IFS=',' read -ra parts <<< "$aff"
      for part in "${parts[@]}"; do
        if [[ "$part" =~ ^([0-9]+)-([0-9]+)$ ]]; then cnt=$((cnt+${BASH_REMATCH[2]}-${BASH_REMATCH[1]}+1))
        elif [[ "$part" =~ ^[0-9]+$ ]]; then cnt=$((cnt+1)); fi
      done; CPUS="$cnt"
    fi
  fi
  if (( CPUS<1 )) && command -v getconf >/dev/null 2>&1; then CPUS="$(getconf _NPROCESSORS_ONLN 2>/dev/null || echo 0)"; fi
  if (( CPUS<1 )); then
    CPUS="$(python3 - <<'PY'
import yaml, os; print(yaml.safe_load(open(os.environ["PIPE_ROOT"] + "/config.yaml"))["globals"]["slurm"]["cpus_per_task"])
PY
)"; CPUS="$(_int_last "$CPUS")"
  fi
  (( CPUS<1 )) && CPUS=1
  local NTMPI=$(( GPUS * 2 )); (( NTMPI<1 )) && NTMPI=1
  local NTOMP=$(( CPUS / NTMPI )); (( NTOMP<1 )) && NTOMP=1
  echo "[MDLAYOUT] gpus=${GPUS} cpus=${CPUS} -> ntmpi=${NTMPI} ntomp=${NTOMP}" 1>&2
  export OMP_NUM_THREADS=$NTOMP OMP_PLACES=cores OMP_PROC_BIND=close GMX_ENABLE_DIRECT_GPU_COMM=1
  echo "-ntmpi $NTMPI -ntomp $NTOMP -nb gpu -bonded gpu -pme gpu -update gpu -npme 1"
}

# ---------- preparation (only if missing) ----------
if [[ ! -f clean.pdb ]]; then
  FF="charmm36-jul2022"; [[ -d ${FF}.ff ]] || ln -s "../../../${FF}.ff" .
  [[ -f input.pdb ]] || { echo "ERROR: input.pdb missing in $WD" >&2; exit 2; }
  gmx pdb2gmx -f input.pdb -o clean.pdb -p topol.top -ff ${FF} -water tip3p
  echo "1" | gmx editconf -f clean.pdb -o aligned.pdb -princ -center 0 0 0
  if [[ -z "$BOX_X" || -z "$BOX_Y" || -z "$BOX_Z" ]]; then
    read -r BOX_X BOX_Y BOX_Z <<<"$(python3 - <<PY
import yaml, os
s = next((x for x in yaml.safe_load(open(os.environ["ROOT"] + "/config.yaml"))["systems"]
         if x["name"] == os.environ["SYS"]), None)
print(*s["box"] if s and "box" in s else (10.0,10.0,10.0))
PY
)"; fi
  gmx editconf -f aligned.pdb -o boxed.gro -c -box ${BOX_X} ${BOX_Y} ${BOX_Z} -bt triclinic
  gmx solvate -cp boxed.gro -cs spc216.gro -o solv.gro -p topol.top
  cat > ions.mdp <<'EOF'
integrator=steep
nsteps=2000
emtol=1000
cutoff-scheme=Verlet
rcoulomb=1.2
rvdw=1.2
EOF
  gmx grompp -f ions.mdp -c solv.gro -p topol.top -o ions.tpr

  # ---- NEW: read SOL_index from config.yaml (per-system), default = 13 ----
  SOLVENT_GROUP_INDEX="$(python3 - <<'PY'
import yaml, os, sys
cfg = yaml.safe_load(open(os.environ["ROOT"] + "/config.yaml"))
sys_name = os.environ.get("SYS")
sysrec = next((s for s in cfg.get("systems", []) if s.get("name")==sys_name), {})
# primary key: systems[].SOL_index  (integer)
val = sysrec.get("SOL_index", 13)
# allow optional nested location if you later move it under prep.genion.solvent_group_index
if not isinstance(val, int):
    val = ((sysrec.get("prep", {}) or {}).get("genion", {}) or {}).get("solvent_group_index", 13)
print(val if isinstance(val, int) else 13)
PY
)"

  # tiny guard: ensure it's a positive integer
  case "$SOLVENT_GROUP_INDEX" in
    ''|*[!0-9]*)
      echo "WARNING: SOL_index invalid ('$SOLVENT_GROUP_INDEX'); falling back to 13" >&2
      SOLVENT_GROUP_INDEX=13
      ;;
  esac

  printf "%s\n" "$SOLVENT_GROUP_INDEX" | gmx genion -s ions.tpr -o ions.gro -p topol.top \
    -pname NA -nname CL -conc 0.15 -neutral

  python3 "$ROOT/scripts/patch_posre_macros.py"
fi

prev_conf="$(pick_prev_conf)"

# ---------- timestep selection ----------
DT_PS=0.001
if [[ "$FINAL" == "1" && "$MODE" == "npt" && "$FC" == "0" ]]; then
  DT_PS=$(python3 - <<'PY'
import yaml, os
print(yaml.safe_load(open(os.environ["PIPE_ROOT"] + "/config.yaml"))["globals"]["dt_ps"])
PY
); fi


# ---------- write stage MDP ----------
gen_mdp () {
  local mode=$1 temp=$2 nsteps=$3 finalflag=$4
  if [[ "$mode" == "em" ]]; then
    cat > stage.mdp <<'EOF'
integrator=steep
nsteps=50000
emtol=100
emstep=0.001
constraints=none
cutoff-scheme=Verlet
coulombtype=PME
rcoulomb=1.2
vdwtype=Cut-off
vdw-modifier=Force-switch
rvdw_switch=1.0
rvdw=1.2
pbc=xyz
EOF
    return
  fi

  local DEF=""; [[ "$FC" != "0" ]] && DEF="-DPOSRES -DPOSRE_FC=${FC}"

  # Defaults for intermediate stages
  local nstcalc=200 nstxoutc=50000 nstenergy=10000 nstlog=10000
  # Final equilibrium NPT: compute strides with actual dt
  if [[ "$finalflag" == "1" && "$mode" == "npt" && "$FC" == "0" ]]; then
    export DT_PS
    read -r nstxoutc nstcalc nstenergy nstlog <<< "$(python3 - <<'PY'
import os,yaml
g=yaml.safe_load(open(os.environ["PIPE_ROOT"] + "/config.yaml"))["globals"]["equilibrium_md"]
dt=float(os.environ.get("DT_PS","0.001"))
def steps(ps): return int(round(ps/dt))
print(steps(g["xtc_write_ps"]), 100, steps(100.0), steps(100.0))
PY
)"
  fi

  {
    echo "define=${DEF}"
    echo "integrator=md"
    echo "dt=${DT_PS}"
    echo "nsteps=${nsteps}"
    echo "tcoupl=V-rescale"
    echo "tc-grps=Protein Non-Protein"
    echo "tau_t=1.0 1.0"
    echo "ref_t=${temp} ${temp}"
    echo "constraints=h-bonds"
    echo "cutoff-scheme=Verlet"
    echo "rlist=1.2"
    echo "rcoulomb=1.2"
    echo "vdwtype=Cut-off"
    echo "vdw-modifier=Force-switch"
    echo "rvdw_switch=1.0"
    echo "rvdw=1.2"
    echo "coulombtype=PME"
    echo "pbc=xyz"
    echo "nstxout-compressed = ${nstxoutc}"
    echo "nstenergy           = ${nstenergy}"
    echo "nstlog              = ${nstlog}"
    echo "nstcalcenergy       = ${nstcalc}"
    echo "nstcomm             = ${nstcalc}"
    echo "xtc-precision       = 1000"
  } > stage.mdp

  if [[ "$mode" == "nvt" ]]; then
    {
      echo "gen-vel=yes"
      echo "gen-temp=${temp}"
      echo "gen-seed=-1"
      echo "pcoupl=no"
    } >> stage.mdp
  else
    {
      echo "gen-vel=no"
      echo "pcoupl=C-rescale"
      echo "pcoupltype=isotropic"
      echo "tau_p=2.0"
      echo "ref_p=1.0"
      echo "compressibility=4.5e-5"
      if [[ "$FC" != "0" ]]; then echo "refcoord_scaling = com"; fi
    } >> stage.mdp
  fi
}

TOP_USE="topol.top"
gen_mdp "$MODE" "$TEMPK" "$NSTEPS" "$FINAL"

OUT="${MODE}$([[ "$FC" != "0" ]] && echo "_fc${FC}" || echo "")"

# grompp
gmx grompp -f stage.mdp -c "$prev_conf" -r "$prev_conf" -p "$TOP_USE" -o "${OUT}.tpr"

# mdrun
if [[ "$MODE" == "em" ]]; then
  echo "[INFO] EM: launching 'gmx mdrun -v -deffnm ${OUT}'"
  gmx mdrun -v -deffnm "${OUT}"
else
  FLAGS="$(build_mdrun_flags)"
  echo "[INFO] ${MODE}: launching 'gmx mdrun -v -deffnm ${OUT} ${FLAGS}'"
  gmx mdrun -v -deffnm "${OUT}" ${FLAGS}
fi

# finalize
if [[ "$MODE" == "npt" && "$FC" == "0" && "$FINAL" == "1" ]]; then
  cp -f "${OUT}.gro" npt_final.gro
  cp -f "${OUT}.tpr" npt_final.tpr
  [[ -f "${OUT}.xtc" ]] && cp -f "${OUT}.xtc" npt_final.xtc
  touch BUILD_DONE
fi
