#!/usr/bin/env bash
set -euo pipefail

ROOT="$(cd "$(dirname "$0")/.."; pwd)"
cd "$ROOT"

# Read GROMACS module from config
GMX_MOD="$(python3 - <<'PY'
import yaml
c=yaml.safe_load(open("config.yaml"))
print(c["globals"]["slurm"]["gromacs_module"])
PY
)"

declare -a SYS_LIST=()

if [[ -f manifests/systems.csv ]]; then
  # Skip header line robustly (use tail -n +2)
  while IFS=, read -r sys pdb bx by bz; do
    # defensive: skip empty or header-looking rows
    [[ -z "${sys:-}" ]] && continue
    [[ "${sys,,}" == "system" || "${sys,,}" == "name" ]] && continue
    SYS_LIST+=("$sys|$pdb|$bx|$by|$bz")
  done < <(tail -n +2 manifests/systems.csv)
else
  # Fallback: read from config.yaml
  while IFS= read -r line; do
    SYS_LIST+=("$line")
  done < <(python3 - <<'PY'
import yaml
c=yaml.safe_load(open("config.yaml"))
for s in c["systems"]:
    name=s["name"]; pdb=s["pdb"]; bx,by,bz=s["box"]
    print(f"{name}|{pdb}|{bx}|{by}|{bz}")
PY
)
fi

if (( ${#SYS_LIST[@]} == 0 )); then
  echo "No systems found to submit."
  exit 0
fi

for rec in "${SYS_LIST[@]}"; do
  IFS='|' read -r SYS PDB BX BY BZ <<< "$rec"
  echo "Submitting staged build for ${SYS}..."
  bash "$ROOT/pipelines/fimA_build_submit_staged.sh" "$SYS" "$PDB" "$BX" "$BY" "$BZ" "$GMX_MOD" ${BUILD_FORCE:+--force}
done
