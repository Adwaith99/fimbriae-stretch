#!/usr/bin/env bash
# One array task == one replicate (single pass, no self-requeue)
set -euo pipefail

ROOT="$(cd "$(dirname "$0")/.."; pwd)"
CFG="$ROOT/config.yaml"
MAN="${FIM_PULL_MANIFEST:-$ROOT/manifests/manifest.csv}"
IDX="${SLURM_ARRAY_TASK_ID:-1}"

say(){ echo "[$(date +'%F %T')] $*"; }

# ---------- 0) Row lookup ----------

echo "[runner] Using manifest: $MAN"
echo "[runner] Array index: $IDX"

ROW=$(awk -F, -v n="$IDX" '
BEGIN{IGNORECASE=1}
NR==1{next}
{
  # trim whitespace + CR from each field
  for(i=1;i<=NF;i++){ gsub(/^[ \t\r\n]+|[ \t\r\n]+$/, "", $i) }
  status=$7
  if (status=="" || tolower(status)=="pending") {
    pending[++cnt]=$0
  }
}
END{
  if (n in pending) print pending[n]
}' "$MAN")

[[ -n "${ROW:-}" ]] || { echo "[runner] ERROR: No pending row for array index $IDX in $MAN"; exit 3; }

IFS=, read -r uid system variant speed k rep status <<<"$ROW"
echo "[runner] ROW: uid=$uid system=$system variant=$variant speed=$speed k=$k rep=$rep status=$status"


# ---------- 1) Ensure build done ----------
BD="$ROOT/systems/$system/00_build"
[[ -f "$BD/npt_final.tpr" && -f "$BD/npt_final.gro" ]] || { say "Build artifacts missing in $BD"; exit 2; }

# ---------- 2) Ensure indices + local posres exist for this variant ----------
IDXDIR="$ROOT/indices/$system"
FINAL_NDX="$IDXDIR/${variant}_pull.ndx"
POSRE_ITP="$IDXDIR/${variant}_posre_anchor_local.itp"
if [[ ! -f "$FINAL_NDX" || ! -f "$POSRE_ITP" ]]; then
  say "Building indices/posres for $system/$variant …"
  python3 "$ROOT/scripts/build_indices_and_posres.py" "$system" "$variant"
fi
[[ -f "$FINAL_NDX" && -f "$POSRE_ITP" ]] || { say "Failed to produce $FINAL_NDX / $POSRE_ITP"; exit 3; }

# ---------- 3) Run directory + materials ----------
RUN_DIR="$ROOT/systems/$system/20_pulls/${variant}/v$(printf "%.2f" "$speed")/k$k/rep$rep"
mkdir -p "$RUN_DIR"; cd "$RUN_DIR"

STARTS_DIR="$ROOT/systems/$system/20_pulls/_starts"
START_GRO="$STARTS_DIR/rep$(printf "%03d" "$rep").gro"
BASE_GRO_FALLBACK="$BD/npt_final.gro"
TOP_TPR="$BD/topol.tpr"
TOP_TOP="$BD/topol.top"

if [[ -f "$START_GRO" ]]; then
  say "Using sampled start: $START_GRO"
  cp -f "$START_GRO" start.gro
else
  say "Sampled start not found; using fallback: $BASE_GRO_FALLBACK"
  cp -f "$BASE_GRO_FALLBACK" start.gro
fi

# ---------- 4) Compute pull init distance + correct vec sign on x-axis ----------
# We will use direction-periodic and a pure x-axis vec (+/- 1 0 0); we decide sign by COM_x(Pulled) - COM_x(Anchor).
read -r PULL_INIT PULL_SIGN <<EOF
$(gmx distance -s "$BD/npt_final.tpr" -f start.gro -n "$FINAL_NDX" -select 'com of group "Pulled" plus com of group "Anchor"' -oall /dev/null 2>/dev/null \
 | awk '
   /#/ {next}
   NF>=4 {dx=$2; dy=$3; dz=$4; d=sqrt(dx*dx+dy*dy+dz*dz); sign=(dx>=0?"+":"-");}
   END{ if(d>0){printf "%.6f %s\n", d, sign} }'
)
EOF
if [[ -z "${PULL_INIT:-}" ]]; then
  # Fallback: let grompp compute init distance; assume positive x
  PULL_INIT="0.0"
  VEC_X="+"
  say "Warning: could not compute COM distance; falling back to vec +1 0 0 with pull_coord1_init=0"
else
  VEC_X="${PULL_SIGN}"
  say "Init distance = ${PULL_INIT} nm ; vec sign on x = ${VEC_X}"
fi
[[ "$VEC_X" == "+" ]] && PULL_VEC="1 0 0" || PULL_VEC="-1 0 0"

# ---------- 5) Render pull.mdp from template with correct numbers ----------
DT_PS=$(python3 - <<PY
import yaml; print(yaml.safe_load(open("$CFG"))["globals"]["dt_ps"])
PY
)
KVAL="$k"
RATE="$speed"
# Choose a reasonable cap on nsteps from target_extension if present
TARGET_EXT=$(python3 - <<PY
import yaml; print(yaml.safe_load(open("$CFG"))["globals"]["target_extension_nm"])
PY
)
# Conservative steps: equal to target distance / (rate*dt) with 1.1 safety
NSTEPS=$(python3 - <<PY
dt=$DT_PS; v=$RATE; targ=$TARGET_EXT
steps=max(1,int((targ/(v*dt))*1.1))
print(steps)
PY
)

# simple Jinja-ish replacement
python3 - "$ROOT/templates/pull.tpl.mdp" > pull.mdp <<PY
import sys
tpl=open(sys.argv[1]).read()
def repl(s, **kw):
    for k,v in kw.items():
        s=s.replace('{{ '+k+' }}', str(v))
    return s
print(repl(tpl,
    dt_ps=${DT_PS},
    nsteps=${NSTEPS},
    pull_init_nm=${PULL_INIT},
    rate_nm_per_ns=${RATE},
    k_kj_mol_nm2=${KVAL}
))
PY

# Patch vec sign in-place
awk -v v1="$(echo "$PULL_VEC")" '
  /^pull_coord1_vec/ {print "pull_coord1_vec         = " v1; next}
  {print}
' pull.mdp > pull.mdp.tmp && mv pull.mdp.tmp pull.mdp

# ---------- 6) Ensure anchor local posre is included (define is already in template) ----------
# We include a one-liner include file to avoid editing topol.top
echo '#include "'"$POSRE_ITP"'"' > includes_posre_anchor.itp

# ---------- 7) Grompp + mdrun ----------
# GPU layout similar to equil
build_mdrun_flags() {
  # derive GPUs/CPUs like your equil_stage.sh (simplified but robust)
  local GPUS="${SLURM_GPUS_PER_NODE:-${SLURM_GPUS_ON_NODE:-${SLURM_GPUS:-4}}}"
  GPUS=$(echo "$GPUS" | sed -E 's/[^0-9]*([0-9]+).*/\1/'); [[ -z "$GPUS" ]] && GPUS=1
  local CPUS="${SLURM_CPUS_PER_TASK:-96}"
  CPUS=$(echo "$CPUS" | sed -E 's/[^0-9]*([0-9]+).*/\1/'); [[ -z "$CPUS" ]] && CPUS=16
  local NTMPI=$(( GPUS * 2 ))
  (( NTMPI<1 )) && NTMPI=1
  local NTOMP=$(( CPUS / NTMPI )); (( NTOMP<1 )) && NTOMP=1
  echo "-ntmpi $NTMPI -ntomp $NTOMP -nb gpu -bonded gpu -pme gpu -update gpu -npme 1"
}

module purge
module load $(python3 - <<PY
import yaml; print(yaml.safe_load(open("$CFG"))["globals"]["slurm"]["gromacs_module"])
PY
)

# Compose a minimal topol include file
cp -f "$TOP_TOP" topol.top
echo '#include "includes_posre_anchor.itp"' >> topol.top

# Use variant index that contains [ Anchor ], [ Pulled ], [ NonProtein ]
gmx grompp -f pull.mdp -c start.gro -p topol.top -n "$FINAL_NDX" -o pull.tpr -maxwarn 1

FLAGS="$(build_mdrun_flags)"
say "Launching mdrun with $FLAGS"
gmx mdrun -v -deffnm pull ${FLAGS}

# ---------- 8) Completion tag (single pass) ----------
CUR_EXT=0
if [[ -f pullx.xvg ]]; then
  CUR_EXT=$(awk '!/^[@#]/{x=$2} END{print (x?x:0)}' pullx.xvg || echo 0)
fi
TARGET=$(python3 - <<PY
import yaml; print(yaml.safe_load(open("$CFG"))["globals"]["target_extension_nm"])
PY
)
STATUS="PARTIAL"
awk -v U="$uid" -v cur="$CUR_EXT" -v tgt="$TARGET" -F, 'BEGIN{OFS=","} NR==1{print;next} {if($1==U){if(cur>=tgt){$7="DONE"} else {$7="PARTIAL"}} print}' "$MAN" > "$MAN.tmp" && mv "$MAN.tmp" "$MAN"
if awk -v e="$CUR_EXT" -v t="$TARGET" 'BEGIN{exit(e>=t?0:1)}'; then
  say "[$uid] DONE (extension ${CUR_EXT} ≥ ${TARGET})"
else
  say "[$uid] PARTIAL (extension ${CUR_EXT} < ${TARGET})"
fi
