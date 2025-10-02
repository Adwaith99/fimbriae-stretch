#!/usr/bin/env bash
set -euo pipefail

ROOT="$(cd "$(dirname "$0")/.."; pwd)"
cd "$ROOT"

# GROMACS module from config
GMX_MOD="$(python3 - <<'PY'
import yaml
c=yaml.safe_load(open("config.yaml"))
print(c["globals"]["slurm"]["gromacs_module"])
PY
)"

# Emit "SYS|PDB|BX|BY|BZ" robustly from systems.csv (if present) else from config.yaml
python3 - <<'PY' | while IFS='|' read -r SYS PDB BX BY BZ; do
import csv, yaml, sys, os

def from_config():
    cfg=yaml.safe_load(open("config.yaml"))
    for s in cfg["systems"]:
        name=s.get("name") or s.get("system")
        pdb=s["pdb"]
        bx,by,bz = s["box"]
        print(f"{name}|{pdb}|{bx}|{by}|{bz}")

def from_csv():
    with open("manifests/systems.csv", newline='') as f:
        r=csv.DictReader(f)
        rows=list(r)
    if not rows:
        return False
    for row in rows:
        # Normalize keys / strip whitespace
        sysname=(row.get("system") or row.get("name") or "").strip()
        pdb=(row.get("pdb") or "").strip()
        bx=(row.get("box_x") or row.get("bx") or "").strip()
        by=(row.get("box_y") or row.get("by") or "").strip()
        bz=(row.get("box_z") or row.get("bz") or "").strip()
        # If anything essential missing, try to fill from config.yaml by name
        if not (sysname and pdb and bx and by and bz):
            cfg=yaml.safe_load(open("config.yaml"))
            match = next((s for s in cfg["systems"] if (s.get("name") or s.get("system"))==sysname), None)
            if match:
                if not pdb: pdb = match["pdb"]
                if not (bx and by and bz):
                    bx,by,bz = map(str, match["box"])
        # Only print valid rows
        if sysname and pdb and bx and by and bz:
            print(f"{sysname}|{pdb}|{bx}|{by}|{bz}")
    return True

if not os.path.exists("manifests/systems.csv") or not from_csv():
    from_config()
PY
do
  echo "Submitting staged build for ${SYS}..."
  bash "$ROOT/pipelines/fimA_build_submit_staged.sh" "$SYS" "$PDB" "$BX" "$BY" "$BZ" "$GMX_MOD" ${BUILD_FORCE:+--force}
done
