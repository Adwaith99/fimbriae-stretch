#!/usr/bin/env bash
set -euo pipefail

# Batch preprocess trajectories using GNU parallel
# Find all SMD run directories and process each in parallel (4 jobs)
# Usage:
#   ./scripts/batch_preproc_traj.sh                    # process all
#   ./scripts/batch_preproc_traj.sh fimA_WT            # process fimA_WT only
#   ./scripts/batch_preproc_traj.sh fimA_WT AtoD       # process fimA_WT/AtoD only

# Find repo root
find_root() {
  local p="$PWD"
  while [[ "$p" != "/" ]]; do
    [[ -f "$p/config.yaml" ]] && { echo "$p"; return 0; }
    p="$(dirname "$p")"
  done
  return 1
}
ROOT="$(find_root)"
if [[ -z "$ROOT" ]]; then
  echo "[batch-preproc] ERROR: cannot locate repo root (config.yaml)" >&2
  exit 2
fi

cd "${ROOT}"

# Check for parallel command
if ! command -v parallel >/dev/null 2>&1; then
  echo "[batch-preproc] ERROR: 'parallel' not found. Install GNU parallel or run preproc_traj.sh manually." >&2
  exit 2
fi

# Parse filters from arguments
FILTER_SYSTEM="${1:-}"
FILTER_VARIANT="${2:-}"

# Find all start_id directories (lowest level: smd/<system>/<variant>/v<speed>/rep<rep>/<start_id>)
echo "[batch-preproc] Scanning for SMD runs..." >&2

tmpfile=$(mktemp)
trap "rm -f ${tmpfile}" EXIT

find smd -mindepth 5 -maxdepth 5 -type d -name "s*" 2>/dev/null | sort > "${tmpfile}" || true

total=$(wc -l < "${tmpfile}")
if [[ "${total}" -eq 0 ]]; then
  echo "[batch-preproc] No SMD runs found" >&2
  exit 0
fi

echo "[batch-preproc] Found ${total} SMD run(s)" >&2

# Filter by system/variant if specified
if [[ -n "${FILTER_SYSTEM}" ]]; then
  grep "smd/${FILTER_SYSTEM}/" "${tmpfile}" > "${tmpfile}.filtered" || true
  mv "${tmpfile}.filtered" "${tmpfile}"
  filtered=$(wc -l < "${tmpfile}")
  echo "[batch-preproc] Filtered to ${filtered} run(s) for system=${FILTER_SYSTEM}" >&2
fi

if [[ -n "${FILTER_VARIANT}" ]]; then
  grep "/${FILTER_VARIANT}/" "${tmpfile}" > "${tmpfile}.filtered" || true
  mv "${tmpfile}.filtered" "${tmpfile}"
  filtered=$(wc -l < "${tmpfile}")
  echo "[batch-preproc] Filtered to ${filtered} run(s) for variant=${FILTER_VARIANT}" >&2
fi

# Check if anything left to process
if [[ ! -s "${tmpfile}" ]]; then
  echo "[batch-preproc] No runs match filters" >&2
  exit 0
fi

# Process in parallel (configurable jobs, default 4)
JOBS="${PREPROC_JOBS:-${3:-4}}"
if ! [[ "${JOBS}" =~ ^[0-9]+$ ]] || [[ "${JOBS}" -le 0 ]]; then
  JOBS=4
fi
echo "[batch-preproc] Starting parallel processing (${JOBS} jobs)..." >&2
cat "${tmpfile}" | parallel -j "${JOBS}" bash scripts/preproc_traj.sh {}

echo "[batch-preproc] Done" >&2

