#!/usr/bin/env bash
# Run one equilibration stage for a system inside systems/<SYS>/00_build
# MODE: em|nvt|npt
# FC:   0 (no restraints) or integer (e.g. 1000, 500, ...)
set -euo pipefail

SYS=${1:?}
MODE=${2:?}
TEMPK=${3:?}
NSTEPS=${4:?}
FC=${5:?}

ROOT="$(cd "$(dirname "$0")/.."; pwd)"
WD="$ROOT/systems/${SYS}/00_build"
cd "$WD"

# Locate prior structure/checkpoint to continue from
prev_conf=""
if [[ -f npt_final.gro ]]; then
  prev_conf="npt_final.gro"
else
  # find most recent stage output as starting conf
  for c in npt_fc*\.gro nvt_fc*\.gro npt.gro nvt.gro em.gro ions.gro solv.gro boxed.gro aligned.pdb clean.pdb; do
    [[ -f $c ]] && prev_conf="$c"
  done
  [[ -z "$prev_conf" ]] && prev_conf="clean.pdb"
fi

# Prepare initial generation if first time
if [[ ! -f clean.pdb ]]; then
  # First-time build prep (pdb2gmx etc.)
  FF="charmm36-jul2022"
  [[ -d ${FF}.ff ]] || ln -s "../../../${FF}.ff" .
  GMX_FORCE_FIELD="$FF" gmx pdb2gmx -f input.pdb -o clean.pdb -p topol.top -ff ${FF} -water tip3p
  echo "1" | gmx editconf -f clean.pdb -o aligned.pdb -princ -center 0 0 0
  gmx editconf -f aligned.pdb -o boxed.gro -c -box \
    $(python3 - <<PY
import yaml; c=yaml.safe_load(open("$ROOT/config.yaml"))
sys=next(s for s in c["systems"] if s["name"]=="$SYS")
print(*sys["box"])
PY
) -bt triclinic
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
  prev_conf="ions.gro"
fi

# MDP generator
gen_mdp () {
  local mode=$1 temp=$2 nsteps=$3
  case "$mode" in
    em)
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
      ;;
    nvt)
      cat > stage.mdp <<EOF
define=$( [[ "$FC" != "0" ]] && echo "-DPOSRES" || echo "" )
integrator=md
dt=0.001
nsteps=${nsteps}
tcoupl=V-rescale
tc-grps=Protein Non-Protein
tau_t=1.0 1.0
ref_t=${temp} ${temp}
constraints=h-bonds
cutoff-scheme=Verlet
rlist=1.2
rcoulomb=1.2
vdwtype=Cut-off
vdw-modifier=Force-switch
rvdw_switch=1.0
rvdw=1.2
coulombtype=PME
pbc=xyz
gen-vel=yes
gen-temp=${temp}
gen-seed=-1
EOF
      ;;
    npt)
      cat > stage.mdp <<EOF
define=$( [[ "$FC" != "0" ]] && echo "-DPOSRES" || echo "" )
integrator=md
dt=0.001
nsteps=${nsteps}
tcoupl=V-rescale
tc-grps=Protein Non-Protein
tau_t=1.0 1.0
ref_t=${temp} ${temp}
pcoupl=C-rescale
pcoupltype=isotropic
tau_p=2.0
ref_p=1.0
compressibility=4.5e-5
constraints=h-bonds
cutoff-scheme=Verlet
rlist=1.2
rcoulomb=1.2
vdwtype=Cut-off
vdw-modifier=Force-switch
rvdw_switch=1.0
rvdw=1.2
coulombtype=PME
pbc=xyz
gen-vel=no
EOF
      ;;
  esac
}

# If restraints requested, generate a scaled posre and a temporary top
TOP_USE="topol.top"
if [[ "$FC" != "0" ]]; then
  python3 "$ROOT/scripts/scale_posre.py" posre.itp "posre_stage_${FC}.itp" "$FC"
  # make a temp top that includes the staged posre instead of default
  awk '
    BEGIN{done=0}
    /^\#include \"posre\.itp\"/ && !done {print "#include \"posre_stage_'$FC'.itp\""; done=1; next}
    {print}
  ' topol.top > topol_stage.top
  TOP_USE="topol_stage.top"
fi

# Generate MDP for this mode
gen_mdp "$MODE" "$TEMPK" "$NSTEPS"

# Grompp + mdrun
outbase="${MODE}$([[ "$FC" != "0" ]] && echo "_fc${FC}" || echo "")"
gmx grompp -f stage.mdp -c "$prev_conf" -r "$prev_conf" -p "$TOP_USE" -o "${outbase}.tpr"
gmx mdrun -deffnm "${outbase}"

# Mark final done if this was the final NPT (no FC)
if [[ "$MODE" == "npt" && "$FC" == "0" ]]; then
  cp -f "${outbase}.gro" npt_final.gro
  touch BUILD_DONE
fi
