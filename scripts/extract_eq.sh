#!/usr/bin/env bash
# Extract XVGs for all NVT/NPT EDRs into xvg/, naming as your plot expects.
set -euo pipefail
[[ "${DEBUG:-0}" == "1" ]] && set -x

ROOT="$(cd "$(dirname "$0")/.."; pwd)"
EDIR="${1:-.}"
OUTDIR="${2:-xvg}"
mkdir -p "$OUTDIR"

det_idx_py="$ROOT/scripts/detect_energy_index.py"

temp_for_fc() {
  local fc="$1"
  python3 - <<'PY' "$fc" 2>/dev/null
import sys, yaml
fc=sys.argv[1]
cfg=yaml.safe_load(open("config.yaml"))
default={"1000":100.0,"500":200.0,"200":303.15,"100":303.15,"50":303.15}
tmap=cfg.get("globals",{}).get("build",{}).get("temp_schedule_by_fc",default)
print(int(round(float(tmap.get(str(fc), default.get(str(fc),303.15))))))
PY
}

ok(){ printf "✓ %s\n" "$*"; }
warn(){ printf "⚠ %s\n" "$*\n" >&2; }

shopt -s nullglob
cd "$EDIR"

n_ok=0
n_warn=0

is_nonempty() { [[ -s "$1" ]]; }

extract_T() {
  local edr="$1" out="$2"
  local idx
  if ! is_nonempty "$edr"; then warn "$edr is empty — skipping"; ((n_warn++)); return; fi
  if ! idx="$("$det_idx_py" "$edr" "temperature" 2>/dev/null)"; then warn "No Temperature in $edr"; ((n_warn++)); return; fi
  if printf "%s\n0\n" "$idx" | gmx energy -f "$edr" -o "$out" >/dev/null 2>&1; then ok "$(basename "$out")"; ((n_ok++)); else warn "gmx energy failed: $edr → $out"; ((n_warn++)); fi
}

extract_P() {
  local edr="$1" out="$2"
  local idx
  if ! is_nonempty "$edr"; then warn "$edr is empty — skipping"; ((n_warn++)); return; fi
  if ! idx="$("$det_idx_py" "$edr" "pressure" 2>/dev/null)"; then warn "No Pressure in $edr"; ((n_warn++)); return; fi
  if printf "%s\n0\n" "$idx" | gmx energy -f "$edr" -o "$out" >/dev/null 2>&1; then ok "$(basename "$out")"; ((n_ok++)); else warn "gmx energy failed: $edr → $out"; ((n_warn++)); fi
}

extract_rho() {
  local edr="$1" out="$2"
  local idx
  if ! is_nonempty "$edr"; then warn "$edr is empty — skipping"; ((n_warn++)); return; fi
  if ! idx="$("$det_idx_py" "$edr" "density" 2>/dev/null)"; then warn "No Density in $edr"; ((n_warn++)); return; fi
  if printf "%s\n0\n" "$idx" | gmx energy -f "$edr" -o "$out" >/dev/null 2>&1; then ok "$(basename "$out")"; ((n_ok++)); else warn "gmx energy failed: $edr → $out"; ((n_warn++)); fi
}

# NVT/NPT by FC
for edr in nvt_fc*.edr npt_fc*.edr; do
  [[ -e "$edr" ]] || continue
  stem="${edr%.edr}"
  if [[ "$stem" =~ ^(nvt|npt)_fc([0-9]+)$ ]]; then
    typ="${BASH_REMATCH[1]}"; fc="${BASH_REMATCH[2]}"
  else
    warn "Skipping unexpected EDR name '$edr'"; ((n_warn++)); continue
  fi
  T="$(temp_for_fc "$fc" || echo 303)"
  if [[ "$typ" == "nvt" ]]; then
    extract_T "$edr" "${OUTDIR}/nvt_${T}K_${fc}_T.xvg"
  else
    extract_P   "$edr" "${OUTDIR}/npt_${T}K_${fc}_P.xvg"
    extract_rho "$edr" "${OUTDIR}/npt_${T}K_${fc}_rho.xvg"
  fi
done

# Final NPT (unrestrained)
if [[ -f npt_final.edr ]]; then
  extract_P   "npt_final.edr" "${OUTDIR}/npt_303K_0_P.xvg"
  extract_rho "npt_final.edr" "${OUTDIR}/npt_303K_0_rho.xvg"
fi

# EM potential (optional)
if [[ -f em.edr ]]; then
  # Try both labels
  local idx
  if idx="$("$det_idx_py" "em.edr" "potential" 2>/dev/null)" || idx="$("$det_idx_py" "em.edr" "potential energy" 2>/dev/null)"; then
    if printf "%s\n0\n" "$idx" | gmx energy -f em.edr -o "${OUTDIR}/em_potential.xvg" >/dev/null 2>&1; then
      ok "em_potential.xvg"; ((n_ok++))
    else
      warn "gmx energy failed: em.edr → em_potential.xvg"; ((n_warn++))
    fi
  else
    warn "No Potential in em.edr"
  fi
fi

echo "---"
echo "Extracted curves: ${n_ok}, warnings: ${n_warn}"
# Only fail if nothing at all was extracted
(( n_ok == 0 )) && { echo "[ERR] No curves extracted — check your .edr files."; exit 3; }
