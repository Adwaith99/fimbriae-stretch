#!/usr/bin/env bash
# One array task == one replicate (self-requeue; single-file -append)
set -euo pipefail

ROOT="$(cd "$(dirname "$0")/.."; pwd)"
MAN="$ROOT/manifests/manifest.csv"
CFG="$ROOT/config.yaml"
IDX="${SLURM_ARRAY_TASK_ID:-1}"
JOBID="${SLURM_JOB_ID:-0}"

# 1) row lookup
ROW=$(awk -F, -v i="$IDX" 'NR==1{next} NR-1==i{print; exit}' "$MAN")
if [[ -z "$ROW" ]]; then echo "No manifest row for index $IDX"; exit 1; fi
IFS=, read -r uid system variant speed k rep status <<<"$ROW"

# 2) ensure build done
if [[ ! -f "$ROOT/systems/$system/00_build/BUILD_DONE" ]]; then
  echo "Build not done for $system; exiting."
  exit 0
fi

RUN_DIR="$ROOT/systems/$system/20_pulls/${variant}/v$(printf "%.2f" "$speed")/k$k/rep$rep"
mkdir -p "$RUN_DIR"; cd "$RUN_DIR"

# 3) materials
STARTS_DIR="$ROOT/systems/$system/20_pulls/_starts"
START_GRO="$STARTS_DIR/rep$(printf "%03d" "$rep").gro"
BASE_GRO_FALLBACK="$ROOT/systems/$system/00_build/npt_final.gro"
TOP_TPR="$ROOT/systems/$system/00_build/topol.tpr"
TOP_TOP="$ROOT/systems/$system/00_build/topol.top"

if [[ -f "$START_GRO" ]]; then
  echo "Using sampled start: $START_GRO"
  cp -f "$START_GRO" start.gro
else
  echo "Sampled start not found; using fallback: $BASE_GRO_FALLBACK"
  cp -f "$BASE_GRO_FALLBACK" start.gro
fi

# 4) build indices & local posre (once per system+variant)
if [[ ! -s "$ROOT/indices/$system/${variant}_pull.ndx" ]] || [[ ! -s "$ROOT/indices/$system/${variant}_posre_anchor_local.itp" ]]; then
  python3 "$ROOT/scripts/build_indices_and_posres.py" "$system" "$variant"
fi

# Copy index & posre locally
cp -f "$ROOT/indices/$system/${variant}_pull.ndx" pull.ndx
cp -f "$ROOT/indices/$system/${variant}_posre_anchor_local.itp" posre_anchor_local.itp

# 5) COM distance for init
gmx distance -s start.gro -f start.gro -n pull.ndx \
  -select 'com of group "Anchor" plus com of group "Pulled"' -oall comdist.xvg >/dev/null 2>&1 || true
PINIT=$(awk '!/^[@#]/{last=$2} END{printf("%.3f\n",(last?last:0))}' comdist.xvg)

# 6) make MDP from template
python3 - "$ROOT" "$speed" "$PINIT" <<'PY'
import sys, yaml, pathlib
ROOT, rate, pinit = pathlib.Path(sys.argv[1]), float(sys.argv[2]), float(sys.argv[3])
cfg = yaml.safe_load(open(ROOT/"config.yaml"))
tpl = (ROOT/"templates"/"pull.tpl.mdp").read_text()
dt = float(cfg["globals"]["dt_ps"])
target = float(cfg["globals"]["target_extension_nm"])
k = float(cfg["globals"]["k_kj_mol_nm2"])
nsteps = int(round((target/rate)*1000.0/dt))
out = (tpl.replace("{{ dt_ps }}", f"{dt:.3f}")
          .replace("{{ nsteps }}", str(nsteps))
          .replace("{{ pull_init_nm }}", f"{pinit:.3f}")
          .replace("{{ rate_nm_per_ns }}", f"{rate:.4f}")
          .replace("{{ k_kj_mol_nm2 }}", f"{k:.1f}"))
open("pull.mdp","w").write(out)
print(f"Wrote pull.mdp (nsteps={nsteps})")
PY

# 7) grompp + mdrun (-append)
gmx grompp -f pull.mdp -c start.gro -p "$TOP_TOP" -n pull.ndx -o pull.tpr -D POSRES_ANCHOR_LOCAL
if [[ -f pull.cpt ]]; then
  gmx mdrun -deffnm pull -cpi pull.cpt -append -maxh 24
else
  gmx mdrun -deffnm pull -append -maxh 24
fi

# 8) progress & requeue
CUR_EXT=$(awk '!/^[@#]/{x=$2} END{print (x?x:0)}' pullx.xvg)
TARGET=$(python3 - <<PY
import yaml; print(yaml.safe_load(open("$CFG"))["globals"]["target_extension_nm"])
PY
)
if awk -v e="$CUR_EXT" -v t="$TARGET" 'BEGIN{exit(e>=t?0:1)}'; then
  echo "[$uid] DONE ($CUR_EXT >= $TARGET)"
  awk -v U="$uid" -F, 'BEGIN{OFS=","} NR==1{print;next} $1==U{$6="DONE"}{print}' "$MAN" > "$MAN.tmp" && mv "$MAN.tmp" "$MAN"
  exit 0
else
  echo "[$uid] continue ($CUR_EXT < $TARGET) → requeue"
  scontrol requeue "$JOBID"
  exit 0
fi
