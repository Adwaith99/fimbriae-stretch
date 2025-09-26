#!/usr/bin/env bash
# Build & equilibrate once per system. Minimal, based on your ladder.
# USAGE: fimA_build_slurm.sh <SYSTEM_NAME> <PDB> <BOX_X> <BOX_Y> <BOX_Z> <GMX_MODULE_SPEC>

set -euo pipefail
SYS=${1:?}
PDB_IN=${2:?}
BOX_X=${3:?}
BOX_Y=${4:?}
BOX_Z=${5:?}
GMX_MOD=${6:?}

WD="systems/${SYS}/00_build"
mkdir -p "$WD"
cp -f "$PDB_IN" "$WD/input.pdb"

module purge
# shellcheck disable=SC2086
module load ${GMX_MOD}

cd "$WD"

# Force field (adjust if needed)
FF="charmm36-jul2022"
[[ -d ${FF}.ff ]] || ln -s "../../../${FF}.ff" .

# pdb2gmx
GMX_FORCE_FIELD="$FF" gmx pdb2gmx -f input.pdb -o clean.pdb -p topol.top -ff ${FF} -water tip3p

# align & center
echo "1" | gmx editconf -f clean.pdb -o aligned.pdb -princ -center 0 0 0

# box/solvate/ions
gmx editconf -f aligned.pdb -o boxed.gro -c -box ${BOX_X} ${BOX_Y} ${BOX_Z} -bt triclinic
gmx solvate  -cp boxed.gro -cs spc216.gro -o solv.gro -p topol.top

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

# EM
cat > em.mdp <<'EOF'
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

gmx grompp -f em.mdp -c ions.gro -p topol.top -o em.tpr
gmx mdrun -v -deffnm em -ntmpi 4 -ntomp 24

# NVT then NPT at 303 K with mild posres on protein
cat > nvt.mdp <<'EOF'
define=-DPOSRES
integrator=md
dt=0.001
nsteps=200000
tcoupl=V-rescale
tc-grps=Protein Non-Protein
tau_t=1.0 1.0
ref_t=303.15 303.15
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
EOF

gmx grompp -f nvt.mdp -c em.gro -r em.gro -p topol.top -o nvt.tpr
gmx mdrun -v -deffnm nvt -ntmpi 4 -ntomp 24

cat > npt.mdp <<'EOF'
define=-DPOSRES
integrator=md
dt=0.001
nsteps=500000
tcoupl=V-rescale
tc-grps=Protein Non-Protein
tau_t=1.0 1.0
ref_t=303.15 303.15
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
EOF

gmx grompp -f npt.mdp -c nvt.gro -t nvt.cpt -r nvt.gro -p topol.top -o npt.tpr
gmx mdrun -v -deffnm npt -ntmpi 4 -ntomp 24

# Final relaxed structure for pulls
cp -f npt.gro npt_303K_0.gro
cp -f npt.tpr topol.tpr

touch BUILD_DONE
echo "Build for ${SYS} complete."
