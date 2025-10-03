#!/usr/bin/env bash
# Usage: pipelines/analyze_eq.sh <SYSTEM>
set -euo pipefail
SYS=${1:?}

ROOT="$(cd "$(dirname "$0")/.."; pwd)"
WD="$ROOT/systems/${SYS}/00_build"
LOG="$ROOT/logs"
mkdir -p "$LOG"

GMX_MOD="$(python3 - <<'PY'
import yaml
c=yaml.safe_load(open("config.yaml"))
print(c["globals"]["slurm"]["gromacs_module"])
PY
)"

module --force purge || true
module load ${GMX_MOD}

echo "[EQ] extracting → ${WD}/xvg/"
bash "$ROOT/scripts/extract_eq.sh" "$WD" "$WD/xvg"

echo "[EQ] plotting grid with your script"
cd "$WD"
python3 "$ROOT/scripts/plot_eq_grid.py" | tee "$LOG/eq_plot_${SYS}.log"
echo "[EQ] done → ${WD}/eq_grid.png , ${WD}/eq_grid.pdf"
