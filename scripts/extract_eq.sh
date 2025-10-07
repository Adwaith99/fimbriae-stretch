#!/usr/bin/env bash
set -euo pipefail
[[ "${DEBUG:-0}" == "1" ]] && set -x

ROOT="$(cd "$(dirname "$0")/.."; pwd)"
EDIR="${1:-.}"
OUTDIR="${2:-xvg}"
mkdir -p "$OUTDIR"

# Map FC -> temperature label (K) from config.yaml (fallbacks preserved)
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
is_nonempty(){ [[ -s "$1" ]]; }

# generic helper (by fixed index)
extract_idx() {
  # $1=edr $2=idx $3=out
  local edr="$1" idx="$2" out="$3"
  [[ -s "$edr" ]] || { warn "$(basename "$edr") is empty — skipping"; return 1; }
  rm -f -- "$out"   # avoid overwrite prompt
  if printf "%s\n0\n" "$idx" | gmx energy -f "$edr" -o "$out" >/dev/null 2>&1; then
    ok "$(basename "$out")"
    return 0
  else
    warn "gmx energy failed: $(basename "$edr") idx=$idx → $(basename "$out")"
    return 1
  fi
}

shopt -s nullglob
cd "$EDIR"

# -------- Staged runs (nvt_fc*.edr, npt_fc*.edr) --------
IDX_T=16
IDX_P=17
IDX_RHO=23

for edr in nvt_fc*.edr npt_fc*.edr; do
  [[ -e "$edr" ]] || continue
  stem="${edr%.edr}"
  if [[ "$stem" =~ ^(nvt|npt)_fc([0-9]+)$ ]]; then
    typ="${BASH_REMATCH[1]}"; fc="${BASH_REMATCH[2]}"
  else
    warn "Skipping unexpected EDR name '$edr'"; continue
  fi
  T="$(temp_for_fc "$fc" || echo 303)"
  if [[ "$typ" == "nvt" ]]; then
    extract_idx "$edr" "$IDX_T"   "${OUTDIR}/nvt_${T}K_${fc}_T.xvg"
  else
    extract_idx "$edr" "$IDX_P"   "${OUTDIR}/npt_${T}K_${fc}_P.xvg"
    extract_idx "$edr" "$IDX_RHO" "${OUTDIR}/npt_${T}K_${fc}_rho.xvg"
  fi
done

# -------- Final NPT (unrestrained): file may be npt_final.edr or npt.edr --------
FINAL_EDR=""
if   [[ -f npt_final.edr ]]; then FINAL_EDR="npt_final.edr"
elif [[ -f npt.edr       ]]; then FINAL_EDR="npt.edr"
fi

if [[ -n "$FINAL_EDR" ]]; then
  # Per your observation: Pressure=16, Density=22 in this file
  extract_idx "$FINAL_EDR" 16 "${OUTDIR}/npt_303K_0_P.xvg"
  extract_idx "$FINAL_EDR" 22 "${OUTDIR}/npt_303K_0_rho.xvg"
fi

# -------- EM potential (Potential = 11 here) --------
if [[ -f em.edr ]]; then
  extract_idx "em.edr" 11 "${OUTDIR}/em_potential.xvg"
fi

echo "Done."
