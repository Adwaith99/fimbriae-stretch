#!/usr/bin/env bash
# MODE: em|nvt|npt
# FC: restraint level (0 => none)
# FINALFLAG: 0 normal stage, 1 => final NPT (equilibrium MD) using dt_ps from config
# BOX_X/Y/Z: box dimensions in nm (passed from submitter)
set -euo pipefail

SYS=${1:?}
MODE=${2:?}
TEMPK=${3:?}
NSTEPS=${4:?}
FC=${5:?}
FINAL=${6:-0}
BOX_X=${7:-}
BOX_Y=${8:-}
BOX_Z=${9:-}

ROOT="$(cd "$(dirname "$0")/.."; pwd)"
WD="$ROOT/systems/${SYS}/00_build"
cd "$WD"

# ---------- helpers ----------
pick_prev_conf() {
  for c in \
    npt_final.gro \
    npt_fc1000.gro npt_fc500.gro npt_fc200.gro npt_fc100.gro npt_fc50.gro \
    nvt_fc1000.gro nvt_fc500.gro nvt_fc200.gro nvt_fc100.gro nvt_fc50.gro \
    em_all.gro em_posres.gro em.gro ions.gro solv.gro boxed.gro aligned.pdb clean.pdb
  do
    [[ -f "$c" ]] && { echo "$c"; return; }
  done
  echo "clean.pdb"
}

build_mdrun_flags() {
  # If user forced flags, use them
  if [[ -n "${GMX_MDRUN_FLAGS:-}" ]]; then
    echo "${GMX_MDRUN_FLAGS}"
    return
  fi

  # helper: extract last integer from a string like "h100:4" -> 4
  _int_last () {
    local s="$1"
    local n
    n="$(echo "$s" | sed -E 's/.*[^0-9]([0-9]+)[^0-9]*$/\1/;t;d')" || true
    if [[ -z "$n" ]]; then echo 1; else echo "$n"; fi
  }

  # GPUs: SLURM may give "gpu:h100:4" or "h100:4" or just "4"
  local RAW_GPUS="${SLURM_GPUS_PER_NODE:-${SLURM_GPUS_ON_NODE:-${SLURM_GPUS:-}}}"
  if [[ -z "$RAW_GPUS" ]]; then
    RAW_GPUS="$(python3 - <<'PY'
import yaml
print(yaml.safe_load(open("config.yaml"))["globals"]["slurm"]["gpus_per_node"])
PY
)"
  fi
  local GPUS="$(_int_last "$RAW_GPUS")"
  (( GPUS < 1 )) && GPUS=1

  # CPUs:
  # Prefer per-task; fall back to per-gpu × gpus; else config.
  local RAW_CPUS="${SLURM_CPUS_PER_TASK:-}"
  if [[ -z "$RAW_CPUS" && -n "${SLURM_CPUS_PER_GPU:-}" ]]; then
    RAW_CPUS=$(( $(_int_last "${SLURM_CPUS_PER_GPU}") * GPUS ))
  fi
  if [[ -z "$RAW_CPUS" ]]; then
    RAW_CPUS="$(python3 - <<'PY'
import yaml
print(yaml.safe_load(open("config.yaml"))["globals"]["slurm"]["cpus_per_task"])
PY
)"
  fi
  local CPUS="$(_int_last "$RAW_CPUS")"
  (( CPUS < 1 )) && CPUS=1

  # Layout (your preference):
  #   ntmpi = gpus * 2
  #   ntomp = floor(cpus / ntmpi), min 1
  local NTMPI=$(( GPUS * 2 ))
  (( NTMPI < 1 )) && NTMPI=1
  local NTOMP=$(( CPUS / NTMPI ))
  (( NTOMP < 1 )) && NTOMP=1

  # Expose to log so we can verify
  echo "[MDLAYOUT] gpus=${GPUS} cpus=${CPUS} -> ntmpi=${NTMPI} ntomp=${NTOMP}" 1>&2

  export OMP_NUM_THREADS=$NTOMP
  export OMP_PLACES=cores
  export OMP_PROC_BIND=close
  export GMX_ENABLE_DIRECT_GPU_COMM=1

  echo "-ntmpi $NTMPI -ntomp $NTOMP -nb gpu -bonded gpu -pme gpu -update gpu -npme 1"
}



# ---------- preparation (only if missing) ----------
if [[ ! -f clean.pdb ]]; then
  FF="charmm36-jul2022"
  [[ -d ${FF}.ff ]] || ln -s "../../../${FF}.ff" .
  [[ -f input.pdb ]] || { echo "ERROR: input.pdb missing in $WD" >&2; exit 2; }

  gmx pdb2gmx -f input.pdb -o clean.pdb -p topol.top -ff ${FF} -water tip3p
  echo "1" | gmx editconf -f clean.pdb -o aligned.pdb -princ -center 0 0 0

  if [[ -z "$BOX_X" || -z "$BOX_Y" || -z "$BOX_Z" ]]; then
    read -r BOX_X BOX_Y BOX_Z <<<"$(python3 - <<PY
import yaml
c=yaml.safe_load(open("$ROOT/config.yaml"))
s=next(x for x in c["systems"] if x["name"]=="$SYS")
print(*s["box"])
PY
)"
  fi
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
  printf "13\n" | gmx genion -s ions.tpr -o ions.gro -p topol.top -pname NA -nname CL -conc 0.15 -neutral

  # Ensure per-chain posre files use POSRE_FC macro (run once)
  python3 "$ROOT/scripts/patch_posre_macros.py"
fi

prev_conf="$(pick_prev_conf)"

# ---------- timestep selection ----------
DT_PS=0.001
if [[ "$FINAL" == "1" && "$MODE" == "npt" && "$FC" == "0" ]]; then
  DT_PS=$(python3 - <<'PY'
import yaml
print(yaml.safe_load(open("config.yaml"))["globals"]["dt_ps"])
PY
)
fi

# ---------- write stage MDP ----------
gen_mdp () {
  local mode=$1; local temp=$2; local nsteps=$3; local finalflag=$4

  if [[ "$mode" == "em" ]]; then
    # EM: compile & run with defaults; include constraints=none (your preference)
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

  # MDP defines for restrained vs unrestrained
  local DEF=""
  if [[ "$FC" != "0" ]]; then
    DEF="-DPOSRES -DPOSRE_FC=${FC}"
  fi

  # Defaults for intermediate stages
  local nstcalc=200
  local nstxoutc=50000
  local nstenergy=10000
  local nstlog=10000

  # Final equilibrium NPT: compute strides with actual dt
  if [[ "$finalflag" == "1" && "$mode" == "npt" && "$FC" == "0" ]]; then
    export DT_PS
    read -r nstxoutc nstcalc nstenergy nstlog <<< "$(python3 - <<'PY'
import os,yaml
g=yaml.safe_load(open("config.yaml"))["globals"]["equilibrium_md"]
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
    echo "nstcomm             = ${nstcalc}"   # fix NOTE 1
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
    } >> stage.mdp
  fi
}

TOP_USE="topol.top"
gen_mdp "$MODE" "$TEMPK" "$NSTEPS" "$FINAL"

OUT="${MODE}$([[ "$FC" != "0" ]] && echo "_fc${FC}" || echo "")"

# ---------- grompp (always) ----------
gmx grompp -f stage.mdp -c "$prev_conf" -r "$prev_conf" -p "$TOP_USE" -o "${OUT}.tpr"

# ---------- mdrun ----------
if [[ "$MODE" == "em" ]]; then
  # EM: run with defaults (no GPU/tMPI/OpenMP overrides), but keep -v for progress
  echo "[INFO] EM: launching 'gmx mdrun -v -deffnm ${OUT}'"
  gmx mdrun -v -deffnm "${OUT}"
else
  # NVT/NPT: use GPU+tMPI/OpenMP flags (or user's GMX_MDRUN_FLAGS if provided)
  FLAGS="$(build_mdrun_flags)"
  echo "[INFO] ${MODE}: launching 'gmx mdrun -v -deffnm ${OUT} ${FLAGS}'"
  gmx mdrun -v -deffnm "${OUT}" ${FLAGS}
fi

# ---------- finalize ----------
if [[ "$MODE" == "npt" && "$FC" == "0" && "$FINAL" == "1" ]]; then
  cp -f "${OUT}.gro" npt_final.gro
  cp -f "${OUT}.tpr" npt_final.tpr
  [[ -f "${OUT}.xtc" ]] && cp -f "${OUT}.xtc" npt_final.xtc
  touch BUILD_DONE
fi
