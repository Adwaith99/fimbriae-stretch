#!/usr/bin/env bash
# Usage: pipelines/analyze_eq.sh <SYSTEM>
set -euo pipefail
[[ "${DEBUG:-0}" == "1" ]] && set -x

SYS=${1:?}

ROOT="$(cd "$(dirname "$0")/.."; pwd)"
WD="$ROOT/systems/${SYS}/00_build"
LOG="$ROOT/logs"
mkdir -p "$LOG"

# Load GROMACS from config
GMX_MOD="$(python3 - <<'PY'
import yaml
c=yaml.safe_load(open("config.yaml"))
print(c["globals"]["slurm"]["gromacs_module"])
PY
)"

# Tolerate sticky base modules
module --force purge || true
module load ${GMX_MOD}

# Sanity: gmx present?
if ! command -v gmx >/dev/null 2>&1; then
  echo "[ERR] 'gmx' not found on PATH after loading '${GMX_MOD}'."
  echo "      Try: module load gromacs/2024.4 (or the site’s gromacs module)."
  exit 2
fi

echo "[ENV] which gmx: $(command -v gmx)"
echo "[ENV] gmx version:"
(gmx --version || true)

echo "[EQ] extracting energies for ${SYS} → ${WD}/xvg/"
bash "$ROOT/scripts/extract_eq.sh" "$WD" "$WD/xvg" | tee "$LOG/extract_${SYS}.log" || true

echo "[EQ] plotting grid with your script"
cd "$WD"
python3 "$ROOT/scripts/plot_eq_grid.py" | tee "$LOG/eq_plot_${SYS}.log"

echo "[EQ] done → ${WD}/eq_grid.png , ${WD}/eq_grid.pdf"
