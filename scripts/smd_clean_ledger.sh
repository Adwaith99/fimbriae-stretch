#!/usr/bin/env bash
set -euo pipefail

# Clean smd_submitted.csv by removing rows that are truly completed on disk
LEDGER="manifests/smd_submitted.csv"

[[ -f "${LEDGER}" ]] || { echo "[clean-ledger] No ledger found: ${LEDGER}"; exit 0; }

TMP=$(mktemp)
trap 'rm -f "${TMP}"' EXIT

REMOVED=0
KEPT=0

while IFS=, read -r sys var rep spd sid || [[ -n "${sys:-}" ]]; do
  [[ -z "${sys:-}" ]] && continue
  
  RUN="smd/${sys}/${var}/v$(printf "%0.3f" "${spd}")/rep${rep}/${sid}"
  IS_DONE=0
  
  if [[ -f "${RUN}/pull.xtc" && -f "${RUN}/pull.log" ]]; then
    if grep -q "Finished mdrun" "${RUN}/pull.log" 2>/dev/null; then
      IS_DONE=1
    fi
  fi
  
  if [[ "${IS_DONE}" -eq 1 ]]; then
    echo "[clean-ledger] Removing completed: ${sys},${var},${rep},${spd},${sid}"
    REMOVED=$((REMOVED+1))
  else
    echo "${sys},${var},${rep},${spd},${sid}" >> "${TMP}"
    KEPT=$((KEPT+1))
  fi
done < "${LEDGER}"

mv "${TMP}" "${LEDGER}"
echo "[clean-ledger] Removed ${REMOVED} completed runs, kept ${KEPT} pending/incomplete."
