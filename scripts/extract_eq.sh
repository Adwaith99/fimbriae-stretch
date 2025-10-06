#!/usr/bin/env bash
# Extract XVGs for all NVT/NPT EDRs into xvg/, naming to match your plot script:
#   nvt_<TEMP>K_<FC>_T.xvg
#   npt_<TEMP>K_<FC>_P.xvg
#   npt_<TEMP>K_<FC>_rho.xvg
#   em_potential.xvg (optional)
set -euo pipefail
[[ "${DEBUG:-0}" == "1" ]] && set -x

ROOT="$(cd "$(dirname "$0")/.."; pwd)"
EDIR="${1:-.}"
OUTDIR="${2:-xvg}"
mkdir -p "$OUTDIR"

det_idx_py="$ROOT/scripts/detect_energy_index.py"  # kept as fallback

# FC -> integer temperature label (K) from config.yaml
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
is_nonempty(){ [[ -s "$1" ]]; }

# Try by NAME first (with aliases). Falls back to index if needed.
energy_by_name_or_index() {
  local edr="$1" out="$2" want="$3"  # want ∈ temperature|pressure|density|potential
  local -a names=()
  case "$want" in
    temperature) names=("Temperature" "Temp." "Temp");;
    pressure)    names=("Pressure" "Pres." "Press");;
    density)     names=("Density" "Mass density" "Rho" "ρ");;
    potential)   names=("Potential" "Potential Energy");;
    *)           names=("$want");;
  esac

  # 1) Try name-based selection
  for label in "${names[@]}"; do
    if printf "%s\n0\n" "$label" | gmx energy -f "$edr" -o "$out" >/dev/null 2>&1; then
      ok "$(basename "$out") (by name: $label)"
      return 0
    fi
  done

  # 2) Fallback to index detection (hybrid script), then use that 1-based index
  if idx="$("$det_idx_py" "$edr" "$want" 2>/dev/null)"; then
    if printf "%s\n0\n" "$idx" | gmx energy -f "$edr" -o "$out" >/dev/null 2>&1; then
      ok "$(basename "$out") (by index: $idx)"
      return 0
    fi
  fi

  warn "Failed to extract '$want' from $(basename "$edr")"
  return 1
}

shopt -s nullglob
cd "$EDIR"

n_ok=0
n_warn=0

# NVT/NPT with FC in filename
for edr in nvt_fc*.edr npt_fc*.edr; do
  [[ -e "$edr" ]] || continue
  stem="${edr%.edr}"                   # e.g., nvt_fc1000
  if [[ "$stem" =~ ^(nvt|npt)_fc([0-9]+)$ ]]; then
    typ="${BASH_REMATCH[1]}"; fc="${BASH_REMATCH[2]}"
  else
    warn "Skipping unexpected EDR name '$edr'"; ((n_warn++)); continue
  fi
  if ! is_nonempty "$edr"; then warn "$edr is empty — skipping"; ((n_warn++)); continue; fi
  T="$(temp_for_fc "$fc" || echo 303)"

  if [[ "$typ" == "nvt" ]]; then
    energy_by_name_or_index "$edr" "${OUTDIR}/nvt_${T}K_${fc}_T.xvg" temperature && ((n_ok++)) || ((n_warn++))
  else
    energy_by_name_or_index "$edr" "${OUTDIR}/npt_${T}K_${fc}_P.xvg"   pressure && ((n_ok++)) || ((n_warn++))
    energy_by_name_or_index "$edr" "${OUTDIR}/npt_${T}K_${fc}_rho.xvg" density  && ((n_ok++)) || ((n_warn++))
  fi
done

# Final NPT (unrestrained)
if [[ -f npt_final.edr ]] && is_nonempty npt_final.edr; then
  energy_by_name_or_index "npt_final.edr" "${OUTDIR}/npt_303K_0_P.xvg"   pressure && ((n_ok++)) || ((n_warn++))
  energy_by_name_or_index "npt_final.edr" "${OUTDIR}/npt_303K_0_rho.xvg" density  && ((n_ok++)) || ((n_warn++))
fi

# EM potential (optional)
if [[ -f em.edr ]] && is_nonempty em.edr; then
  energy_by_name_or_index "em.edr" "${OUTDIR}/em_potential.xvg" potential && ((n_ok++)) || ((n_warn++))
fi

echo "---"
echo "Extracted curves: ${n_ok}, warnings: ${n_warn}"
# Always succeed; plotting will use whatever exists
exit 0
