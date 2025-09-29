#!/usr/bin/env bash
# Prep-only on login node: pdb2gmx → editconf → solvate → genion.
# Optionally compile EM (grompp only) to validate topology/MDP/includes.
#
# USAGE:
#   build_preflight.sh <SYSTEM> [--em-grompp]
#
# This DOES NOT run any mdrun. It only prepares files and, if requested,
# compiles em.tpr so you can see errors immediately without queueing.

set -euo pipefail

SYS=${1:? "Usage: build_preflight.sh <SYSTEM> [--em-grompp]"}
EM_CHECK=${2:-""}

ROOT="$(cd "$(dirname "$0")/.."; pwd)"
WD="$ROOT/systems/${SYS}/00_build"
LOG="$ROOT/logs"
mkdir -p "$WD" "$LOG"

REPORT="$LOG/prep_${SYS}.txt"
: > "$REPORT"

say() { echo -e "$*" | tee -a "$REPORT"; }

# ---- Pull config bits ----
read -r BOX_X BOX_Y BOX_Z PDB_REL <<<"$(python3 - <<PY
import yaml
c=yaml.safe_load(open("$ROOT/config.yaml"))
s=next(x for x in c["systems"] if x["name"]=="$SYS")
print(*s["box"], s["pdb"])
PY
)"

PDB_IN="$ROOT/$PDB_REL"
[[ -f "$PDB_IN" ]] || { say "ERROR: PDB not found: $PDB_IN"; exit 2; }

say "== PREP for system: $SYS =="
say "PDB: $PDB_IN"
say "Box: ${BOX_X} ${BOX_Y} ${BOX_Z} nm"
cd "$WD"

# Clean previous partial prep (but keep final NPT artifacts if they exist)
rm -f clean.pdb aligned.pdb boxed.gro solv.gro ions.gro ions.tpr topol.top ions.mdp

# Force field (adjust if needed)
FF="charmm36-jul2022"
[[ -d ${FF}.ff ]] || ln -s "../../../${FF}.ff" .

# pdb2gmx
say "\n-- pdb2gmx --"
gmx pdb2gmx -f "$PDB_IN" -o clean.pdb -p topol.top -ff ${FF} -water tip3p 2>&1 | tee -a "$REPORT"

# align & center
say "\n-- editconf: align/principal axes --"
echo "1" | gmx editconf -f clean.pdb -o aligned.pdb -princ -center 0 0 0 2>&1 | tee -a "$REPORT"

# box
say "\n-- editconf: define triclinic box --"
gmx editconf -f aligned.pdb -o boxed.gro -c -box ${BOX_X} ${BOX_Y} ${BOX_Z} -bt triclinic 2>&1 | tee -a "$REPORT"

# solvate
say "\n-- solvate --"
gmx solvate -cp boxed.gro -cs spc216.gro -o solv.gro -p topol.top 2>&1 | tee -a "$REPORT"

# ions pre-grompp
say "\n-- grompp: ions --"
cat > ions.mdp <<'EOF'
integrator=steep
nsteps=2000
emtol=1000
cutoff-scheme=Verlet
rcoulomb=1.2
rvdw=1.2
EOF
gmx grompp -f ions.mdp -c solv.gro -p topol.top -o ions.tpr 2>&1 | tee -a "$REPORT"

# genion
say "\n-- genion --"
printf "13\n" | gmx genion -s ions.tpr -o ions.gro -p topol.top -pname NA -nname CL -conc 0.15 -neutral 2>&1 | tee -a "$REPORT"

say "\nPrep complete: produced ions.gro (ready for EM)."

# Optional: EM grompp only (no run) to surface MDP/include issues quickly
if [[ "$EM_CHECK" == "--em-grompp" ]]; then
  say "\n-- OPTIONAL: EM grompp (compile only) --"
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
  gmx grompp -f em.mdp -c ions.gro -r ions.gro -p topol.top -o em.tpr 2>&1 | tee -a "$REPORT"
  say "EM grompp OK → em.tpr compiled."
fi

say "\nReport: $REPORT"
