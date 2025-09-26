#!/usr/bin/env bash
set -euo pipefail
ROOT="$(cd "$(dirname "$0")/.."; pwd)"
python3 "$ROOT/scripts/generate_manifest.py" >/dev/null

PARTITION=$(python3 - <<PY
import yaml,sys; c=yaml.safe_load(open("$ROOT/config.yaml"))
print(c["globals"]["slurm"]["partition"])
PY
)
GMXMOD=$(python3 - <<PY
import yaml; c=yaml.safe_load(open("$ROOT/config.yaml"))
print(c["globals"]["slurm"]["gromacs_module"])
PY
)

while IFS=, read -r system pdb box build_status; do
  [[ "$system" == "system" ]] && continue
  if [[ "$build_status" != "DONE" ]]; then
    echo "Submitting staged build for $system..."
    read -r BX BY BZ <<<"$box"
    bash "$ROOT/pipelines/fimA_build_submit_staged.sh" "${system}" "$ROOT/${pdb}" ${BX:-100} ${BY:-20} ${BZ:-20} "${GMXMOD}"
  else
    echo "Build already DONE for $system"
  fi
done < <(tail -n +2 "$ROOT/manifests/systems.csv")
