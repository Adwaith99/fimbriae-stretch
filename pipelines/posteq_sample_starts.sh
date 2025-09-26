#!/usr/bin/env bash
# Extract evenly spaced frames from final equilibrium NPT (npt_final.tpr/xtc)
# and materialize per-replicate starting .gro files that runner.sh will use.
# USAGE: postequ_sample_starts.sh <SYSTEM>

set -euo pipefail
SYS=${1:?}
ROOT="$(cd "$(dirname "$0")/.."; pwd)"
BD="$ROOT/systems/${SYS}/00_build"
OUTDIR="$ROOT/systems/${SYS}/20_pulls/_starts"
MAN="$ROOT/manifests/manifest.csv"

[[ -f "$BD/npt_final.tpr" && -f "$BD/npt_final.xtc" ]] || { echo "Final NPT artifacts missing"; exit 1; }
mkdir -p "$OUTDIR"

# Pull sampling params and #reps per system (fallback to globals)
read -r SAMPLE_EVERY_PS WARMUP_PS NREPS <<< "$(python3 - <<'PY'
import yaml,sys,csv
c=yaml.safe_load(open("config.yaml"))
g=c["globals"]["equilibrium_md"]
sample=int(g["sample_every_ps"]); warm=int(g["warmup_ns"]*1000)
# pick reps for this system (max over its variants)
import pandas as pd
import pathlib
try:
    df=pd.read_csv(pathlib.Path("manifests/manifest.csv"))
    n=int(df[df["system"]==sys.argv[1]]["rep"].astype(int).max())
except Exception: n=int(c["globals"]["n_reps_default"])
print(sample, warm, n)
PY
"$SYS")"

# Use trjconv to write frames at SAMPLE_EVERY_PS after WARMUP
# Create a temporary .ndx selecting System
echo -e "q\n" | gmx make_ndx -f "$BD/npt_final.tpr" -o /tmp/sys.ndx >/dev/null
# Write a time range list for -dt (we’ll rely on -dt with starting time)
# We’ll set -b (begin) to WARMUP_PS and -dt to SAMPLE_EVERY_PS
gmx trjconv -s "$BD/npt_final.tpr" -f "$BD/npt_final.xtc" -o "$OUTDIR/frame_.gro" -sep -b "$WARMUP_PS" -dt "$SAMPLE_EVERY_PS" -n /tmp/sys.ndx <<EOF >/dev/null
System
EOF

# Enumerate produced frames and symlink first NREPS to rep###.gro
i=0
for f in "$OUTDIR"/frame_*.gro; do
  [[ -e "$f" ]] || break
  i=$((i+1))
  rep=$(printf "%03d" "$i")
  ln -sf "$(basename "$f")" "$OUTDIR/rep${rep}.gro"
done

echo "Wrote starts in $OUTDIR (frames: $i; linked up to repNNN)"
