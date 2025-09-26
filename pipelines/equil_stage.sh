#!/usr/bin/env bash
# MODE: em|nvt|npt
# FC: restraint level (0 => none)
# FINALFLAG: 0 normal stage, 1 => final NPT (equilibrium MD) with tuned outputs
set -euo pipefail

SYS=${1:?}
MODE=${2:?}
TEMPK=${3:?}
NSTEPS=${4:?}
FC=${5:?}
FINAL=${6:?0}

ROOT="$(cd "$(dirname "$0")/.."; pwd)"
WD="$ROOT/systems/${SYS}/00_build"
cd "$WD"

# Prepare initial topology/solvent if missing
if [[ ! -f clean.pdb ]]; then
  FF="charmm36-jul2022"
  [[ -d ${FF}.ff ]] || ln -s "../../../${FF}.ff" .
  gmx pdb2gmx -f input.pdb -o clean.pdb -p topol.top -ff ${FF} -water tip3p
  echo "1" | gmx editconf -f clean.pdb -o aligned.pdb -princ -center 0 0 0
  read -r BX BY BZ <<<"$(python3 - <<PY
import yaml; c=yaml.safe_load(open("$ROOT/config.yaml"))
sys=next(s for s in c["systems"] if s["name"]=="$SYS")
print(*sys["box"])
PY
)"
  gmx editconf -f aligned.pdb -o boxed.gro -c -box ${BX} ${BY} ${BZ} -bt triclinic
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
fi

# choose previous conf
prev_conf="npt_final.gro"
if [[ ! -f $prev_conf ]]; then
  for c in npt_fc*\.gro nvt_fc*\.gro npt.gro nvt.gro em.gro ions.gro solv.gro boxed.gro aligned.pdb clean.pdb; do
    [[ -f $c ]] && prev_conf="$c"
  done
fi

# build-time dt (keep 1 fs for robust equilibration)
DT_PS=0.001

# Generate stage MDP
gen_mdp () {
  local mode=$1; local temp=$2; local nsteps=$3; local finalflag=$4
  if [[ "$mode" == "em" ]]; then
    cat > stage.mdp <<'EOF'
integrator=steep
nsteps=50000
emtol=100
emstep=0.001
cutoff-scheme=Verlet
coulombtype=PME
rcoulomb=1.2
vdwtype=PME
rvdw=1.2
fourierspacing=0.12
pbc=xyz
DispCorr=no
EOF
    return
  fi

  # Common for nvt/npt
  {
    [[ "$FC" != "0" ]] && echo "define=-DPOSRES" || echo "define="
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

  # Output control
  if [[ "$finalflag" == "1" && "$mode" == "npt" && "$FC" == "0" ]]; then
    # Equilibrium MD outputs for sampling — read from config
    python3 - <<'PY' >> stage.mdp
import yaml
g=yaml.safe_load(open("config.yaml"))["globals"]["equilibrium_md"]
def to_steps(ps): 
    dt=0.001
    return int(round(ps/dt))
print("nstxout-compressed =", to_steps(g["xtc_write_ps"]))  # XTC stride
print("nstenergy           =", to_steps(100.0))              # energies every 100 ps
print("nstlog              =", to_steps(100.0))              # logs every 100 ps
print("nstcalcenergy       =", 100)
print("xtc-precision       =", 1000)
PY
  else
    # Tighter but modest output for intermediate stages
    cat >> stage.mdp <<'EOF'
nstxout-compressed = 50000
nstenergy           = 10000
nstlog              = 10000
nstcalcenergy       = 200
xtc-precision       = 1000
EOF
  fi
}

# Optional: stage-local posre scaling
TOP_USE="topol.top"
if [[ "$FC" != "0" ]]; then
  python3 "$ROOT/scripts/scale_posre.py" posre.itp "posre_stage_${FC}.itp" "$FC"
  awk '
    BEGIN{done=0}
    /^\#include \"posre\.itp\"/ && !done {print "#include \"posre_stage_'$FC'.itp\""; done=1; next}
    {print}
  ' topol.top > topol_stage.top
  TOP_USE="topol_stage.top"
fi

gen_mdp "$MODE" "$TEMPK" "$NSTEPS" "$FINAL"

OUT="${MODE}$([[ "$FC" != "0" ]] && echo "_fc${FC}" || echo "")"
gmx grompp -f stage.mdp -c "$prev_conf" -r "$prev_conf" -p "$TOP_USE" -o "${OUT}.tpr"
gmx mdrun -deffnm "${OUT}"

# If final equilibrium NPT: mark and leave artifacts for sampling
if [[ "$MODE" == "npt" && "$FC" == "0" && "$FINAL" == "1" ]]; then
  cp -f "${OUT}.gro" npt_final.gro
  cp -f "${OUT}.tpr" npt_final.tpr
  cp -f "${OUT}.xtc" npt_final.xtc
  touch BUILD_DONE
fi
