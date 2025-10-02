#!/usr/bin/env bash
set -euo pipefail

ROOT="$(cd "$(dirname "$0")/.."; pwd)"
cd "$ROOT"

# Get GROMACS module from config
GMX_MOD="$(python3 - <<'PY'
import yaml
c=yaml.safe_load(open("config.yaml"))
print(c["globals"]["slurm"]["gromacs_module"])
PY
)"

declare -a SYS_LIST=()

if [[ -f manifests/systems.csv ]]; then
  # Read WITHOUT a pipeline to avoid subshell array loss
  while IFS=, read -r name pdb bx by bz; do
    [[ "$name" == "name" ]] && continue         # skip header
    [[ -z "${name:-}" ]] && continue
    SYS_LIST+=("$name|$pdb|$bx|$by|$bz")
  done < manifests/systems.csv
else
  # Fallback: read straight from config.yaml
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
  # Call the staged submitter (it does all sbatch calls internally + resume/skip)
  bash "$ROOT/pipelines/fimA_build_submit_staged.sh" "$SYS" "$PDB" "$BX" "$BY" "$BZ" "$GMX_MOD" ${BUILD_FORCE:+--force}
done
