#!/usr/bin/env bash
# Extracts XVGs for all NVT/NPT EDRs into xvg/, naming as:
#   nvt_<TEMP>K_<FC>_T.xvg
#   npt_<TEMP>K_<FC>_P.xvg
#   npt_<TEMP>K_<FC>_rho.xvg
# Final unrestrained NPT → npt_303K_0_{P,rho}.xvg
set -euo pipefail

ROOT="$(cd "$(dirname "$0")/.."; pwd)"
EDIR="${1:-.}"
OUTDIR="${2:-xvg}"
mkdir -p "$OUTDIR"

det_idx_py="$ROOT/scripts/detect_energy_index.py"

# Map FC -> temperature K from config.yaml
temp_for_fc() {
  local fc="$1"
  python3 - <<'PY' "$fc" 2>/dev/null
import sys, yaml
fc = sys.argv[1]
cfg = yaml.safe_load(open("config.yaml"))
default = {"1000":100.0,"500":200.0,"200":303.15,"100":303.15,"50":303.15}
tmap = cfg.get("globals",{}).get("build",{}).get("temp_schedule_by_fc", default)
print(int(round(float(tmap.get(str(fc), default.get(str(fc), 303.15))))))
PY
}

# Helper to print status
ok(){ printf "✓ %s\n" "$*"; }

shopt -s nullglob
cd "$EDIR"

# NVT & NPT with FC in filename
for edr in nvt_fc*.edr npt_fc*.edr; do
  [[ -e "$edr" ]] || continue
  stem="${edr%.edr}"         # e.g., nvt_fc1000
  if [[ "$stem" =~ ^(nvt|npt)_fc([0-9]+)$ ]]; then
    typ="${BASH_REMATCH[1]}"
    fc="${BASH_REMATCH[2]}"
  else
    continue
  fi

  T="$(temp_for_fc "$fc")"   # integer K label
  if [[ "$typ" == "nvt" ]]; then
    Tidx="$("$det_idx_py" "$edr" "Temperature")"
    printf "%s\n0\n" "$Tidx" | gmx energy -f "$edr" -o "${OUTDIR}/nvt_${T}K_${fc}_T.xvg" >/dev/null
    ok "nvt_${T}K_${fc}_T.xvg"
  else
    Pidx="$("$det_idx_py" "$edr" "Pressure")"
    Ridx="$("$det_idx_py" "$edr" "Density")"
    printf "%s\n0\n" "$Pidx" | gmx energy -f "$edr" -o "${OUTDIR}/npt_${T}K_${fc}_P.xvg"   >/dev/null
    printf "%s\n0\n" "$Ridx" | gmx energy -f "$edr" -o "${OUTDIR}/npt_${T}K_${fc}_rho.xvg" >/dev/null
    ok "npt_${T}K_${fc}_{P,rho}.xvg"
  fi
done

# Final NPT (no FC in name)
if [[ -f npt_final.edr ]]; then
  Pidx="$("$det_idx_py" "npt_final.edr" "Pressure")"
  Ridx="$("$det_idx_py" "npt_final.edr" "Density")"
  printf "%s\n0\n" "$Pidx" | gmx energy -f npt_final.edr -o "${OUTDIR}/npt_303K_0_P.xvg"   >/dev/null || true
  printf "%s\n0\n" "$Ridx" | gmx energy -f npt_final.edr -o "${OUTDIR}/npt_303K_0_rho.xvg" >/dev/null || true
  ok "npt_303K_0_{P,rho}.xvg"
fi

# Optional: EM potential (not consumed by your grid, but handy)
if [[ -f em.edr ]]; then
  if TIDX="$("$det_idx_py" em.edr "Potential")"; then
    printf "%s\n0\n" "$TIDX" | gmx energy -f em.edr -o "${OUTDIR}/em_potential.xvg" >/dev/null || true
    ok "em_potential.xvg"
  fi
fi
