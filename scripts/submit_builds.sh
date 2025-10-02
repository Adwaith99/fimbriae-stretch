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

# Prefer manifests/systems.csv if present; otherwise read config.yaml
SYS_LIST=()
if [[ -f manifests/systems.csv ]]; then
  # Expected columns: name,pdb,box_x,box_y,box_z
  tail -n +2 manifests/systems.csv | while IFS=, read -r name pdb bx by bz; do
    [[ -z "$name" ]] && continue
    SYS_LIST+=("$name|$pdb|$bx|$by|$bz")
  done
else
  # Fallback: pull from config.yaml
  python3 - <<'PY'
import yaml
c=yaml.safe_load(open("config.yaml"))
for s in c["systems"]:
    name=s["name"]; pdb=s["pdb"]; bx,by,bz=s["box"]
    print(f"{name}|{pdb}|{bx}|{by}|{bz}")
PY
fi | while IFS= read -r line; do
  SYS_LIST+=("$line")
done

if [[ ${#SYS_LIST[@]} -eq 0 ]]; then
  echo "No systems found to submit."
  exit 0
fi

for rec in "${SYS_LIST[@]}"; do
  IFS='|' read -r SYS PDB BX BY BZ <<< "$rec"
  echo "Submitting staged build for ${SYS}..."
  # Call the *submitter* (it does all sbatch calls internally)
  bash "$ROOT/pipelines/fimA_build_submit_staged.sh" "$SYS" "$PDB" "$BX" "$BY" "$BZ" "$GMX_MOD" ${BUILD_FORCE:+--force}
done
