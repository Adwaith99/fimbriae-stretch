#!/usr/bin/env bash
# Extract XVGs for all NVT/NPT EDRs into xvg/, naming as:
#   nvt_<TEMP>K_<FC>_T.xvg
#   npt_<TEMP>K_<FC>_P.xvg
#   npt_<TEMP>K_<FC>_rho.xvg
# Also: em_potential.xvg if em.edr exists.
# Continues on errors; warns and summarizes at end.
set -euo pipefail
[[ "${DEBUG:-0}" == "1" ]] && set -x

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

ok(){ printf "✓ %s\n" "$*"; }
warn(){ printf "⚠ %s\n" "$*\n" >&2; }

shopt -s nullglob
cd "$EDIR"

n_ok=0
n_warn=0

extract_nvt() {
  local edr="$1" fc="$2" T="$3"
  local Tidx
  if ! Tidx="$("$det_idx_py" "$edr" "Temperature" 2>/dev/null)"; then
    warn "Failed to detect 'Temperature' in $edr"; ((n_warn++)); return
  fi
  if printf "%s\n0\n" "$Tidx" | gmx energy -f "$edr" -o "${OUTDIR}/nvt_${T}K_${fc}_T.xvg" >/dev/null 2>&1; then
    ok "nvt_${T}K_${fc}_T.xvg"; ((n_ok++))
  else
    warn "gmx energy failed for $edr (Temperature)"; ((n_warn++))
  fi
}

extract_npt() {
  local edr="$1" fc="$2" T="$3"
  local Pidx Ridx
  if ! Pidx="$("$det_idx_py" "$edr" "Pressure" 2>/dev/null)"; then
    warn "Failed to detect 'Pressure' in $edr"; ((n_warn++)); return
  fi
  if ! Ridx="$("$det_idx_py" "$edr" "Density" 2>/dev/null)"; then
    warn "Failed to detect 'Density' in $edr"; ((n_warn++)); return
  fi
  if printf "%s\n0\n" "$Pidx" | gmx energy -f "$edr" -o "${OUTDIR}/npt_${T}K_${fc}_P.xvg" >/dev/null 2>&1; then
    ok "npt_${T}K_${fc}_P.xvg"; ((n_ok++))
  else
    warn "gmx energy failed for $edr (Pressure)"; ((n_warn++))
  fi
  if printf "%s\n0\n" "$Ridx" | gmx energy -f "$edr" -o "${OUTDIR}/npt_${T}K_${fc}_rho.xvg" >/dev/null 2>&1; then
    ok "npt_${T}K_${fc}_rho.xvg"; ((n_ok++))
  else
    warn "gmx energy failed for $edr (Density)"; ((n_warn++))
  fi
}

# NVT & NPT with FC in filename
for edr in nvt_fc*.edr npt_fc*.edr; do
  [[ -e "$edr" ]] || continue
  stem="${edr%.edr}"  # e.g., nvt_fc1000
  if [[ "$stem" =~ ^(nvt|npt)_fc([0-9]+)$ ]]; then
    typ="${BASH_REMATCH[1]}"
    fc="${BASH_REMATCH[2]}"
  else
    warn "Skipping unexpected EDR filename '$edr'"; ((n_warn++)); continue
  fi
  T="$(temp_for_fc "$fc" || echo 303)"
  if [[ "$typ" == "nvt" ]]; then extract_nvt "$edr" "$fc" "$T"; else extract_npt "$edr" "$fc" "$T"; fi
done

# Final NPT (unrestrained)
if [[ -f npt_final.edr ]]; then
  if Pidx="$("$det_idx_py" "npt_final.edr" "Pressure" 2>/dev/null)"; then
    if printf "%s\n0\n" "$Pidx" | gmx energy -f npt_final.edr -o "${OUTDIR}/npt_303K_0_P.xvg" >/dev/null 2>&1; then
      ok "npt_303K_0_P.xvg"; ((n_ok++))
    else
      warn "gmx energy failed for npt_final.edr (Pressure)"; ((n_warn++))
    fi
  else
    warn "Failed to detect 'Pressure' in npt_final.edr"; ((n_warn++))
  fi
  if Ridx="$("$det_idx_py" "npt_final.edr" "Density" 2>/dev/null)"; then
    if printf "%s\n0\n" "$Ridx" | gmx energy -f npt_final.edr -o "${OUTDIR}/npt_303K_0_rho.xvg" >/dev/null 2>&1; then
      ok "npt_303K_0_rho.xvg"; ((n_ok++))
    else
      warn "gmx energy failed for npt_final.edr (Density)"; ((n_warn++))
    fi
  else
    warn "Failed to detect 'Density' in npt_final.edr"; ((n_warn++))
  fi
fi

# EM potential (optional)
if [[ -f em.edr ]]; then
  if TIDX="$("$det_idx_py" em.edr "Potential" 2>/dev/null)"; then
    if printf "%s\n0\n" "$TIDX" | gmx energy -f em.edr -o "${OUTDIR}/em_potential.xvg" >/dev/null 2>&1; then
      ok "em_potential.xvg"; ((n_ok++))
    else
      warn "gmx energy failed for em.edr (Potential)"; ((n_warn++))
    fi
  else
    warn "Failed to detect 'Potential' in em.edr"; ((n_warn++))
  fi
fi

echo "---"
echo "Extracted curves: ${n_ok}, warnings: ${n_warn}"
# If nothing extracted at all, return non-zero to hint something is wrong
if (( n_ok == 0 )); then
  echo "[ERR] No curves extracted. Check your *.edr files exist in $(pwd)" >&2
  exit 3
fi
