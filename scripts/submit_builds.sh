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

# Emit records as: SYS|PDB|BX|BY|BZ
emit_systems() {
python3 - <<'PY'
import os, csv, yaml

def from_config():
    cfg=yaml.safe_load(open("config.yaml"))
    for s in cfg["systems"]:
        name = s.get("name") or s.get("system")
        pdb  = s["pdb"]
        bx,by,bz = s["box"]
        print(f"{name}|{pdb}|{bx}|{by}|{bz}")

def from_csv():
    with open("manifests/systems.csv", newline='') as f:
        r = csv.DictReader(f)
        rows = list(r)
    if not rows:
        return False
    # map a row to normalized fields with fallback to config when needed
    cfg = yaml.safe_load(open("config.yaml"))
    cfg_by_name = { (s.get("name") or s.get("system")): s for s in cfg["systems"] }
    for row in rows:
        sysname = (row.get("system") or row.get("name") or "").strip()
        if not sysname:
            continue
        pdb = (row.get("pdb") or "").strip()
        bx  = (row.get("box_x") or row.get("bx") or "").strip()
        by  = (row.get("box_y") or row.get("by") or "").strip()
        bz  = (row.get("box_z") or row.get("bz") or "").strip()
        if not (pdb and bx and by and bz):
            # top up from config if missing
            s = cfg_by_name.get(sysname)
            if s:
                if not pdb: pdb = s["pdb"]
                if not (bx and by and bz):
                    bx,by,bz = map(str, s["box"])
        if sysname and pdb and bx and by and bz:
            print(f"{sysname}|{pdb}|{bx}|{by}|{bz}")
    return True

if not os.path.exists("manifests/systems.csv") or not from_csv():
    from_config()
PY
}

# Iterate over systems using process substitution (robust, no subshell pitfalls)
found_any=0
while IFS='|' read -r SYS PDB BX BY BZ; do
  [[ -z "${SYS:-}" ]] && continue
  found_any=1
  echo "Submitting staged build for ${SYS}..."
  bash "$ROOT/pipelines/fimA_build_submit_staged.sh" "$SYS" "$PDB" "$BX" "$BY" "$BZ" "$GMX_MOD" ${BUILD_FORCE:+--force}
done < <(emit_systems)

if [[ "$found_any" == "0" ]]; then
  echo "No systems found to submit."
fi
