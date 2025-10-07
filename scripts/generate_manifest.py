#!/usr/bin/env python3
import csv, yaml, os
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
CFG  = yaml.safe_load(open(ROOT/"config.yaml"))
MAN  = ROOT/"manifests/manifest.csv"
SYS  = ROOT/"manifests/systems.csv"
os.makedirs(ROOT/"manifests", exist_ok=True)

def read_csv_dict(p):
    if not p.exists(): return {}
    with open(p) as f:
        r = csv.DictReader(f)
        return {row["uid"]: row for row in r}

def uid(system, variant_id, speed, k, rep):
    return f"{system}|{variant_id}|v{speed:.2f}|k{k}|rep{rep:03d}"

globals = CFG["globals"]
speeds_global = globals["pull_speeds_nm_per_ns"]
nreps_default = globals["n_reps_default"]
kval = globals["k_kj_mol_nm2"]

existing = read_csv_dict(MAN)
rows = {}

# systems table (for build status)
srows = []
for sys in CFG["systems"]:
    sname = sys["name"]
    build_dir = ROOT/f"systems/{sname}/00_build"
    build_status = "DONE" if (build_dir/"BUILD_DONE").exists() else "PENDING"
    srows.append({"system": sname, "pdb": sys["pdb"], "box": " ".join(map(str, sys.get("box",[]))), "build_status": build_status})

with open(SYS, "w", newline="") as f:
    w = csv.DictWriter(f, fieldnames=["system","pdb","box","build_status"])
    w.writeheader(); w.writerows(srows)

for sys in CFG["systems"]:
    sname = sys["name"]
    sys_speeds = sys.get("speeds")  # may be None; fallback handled below
    sys_nreps  = sys.get("n_reps")
    for var in sys["variants"]:
        vid = var["id"]
        # Variant-level overrides take highest precedence
        v_speeds = var.get("speeds") or sys_speeds or speeds_global
        v_nreps  = var.get("n_reps") or sys_nreps or nreps_default
        for sp in v_speeds:
            for rep in range(1, v_nreps+1):
                U = uid(sname, vid, sp, kval, rep)
                row = {
                    "uid": U,
                    "system": sname,
                    "variant": vid,
                    "speed_nm_per_ns": f"{sp:.2f}",
                    "k_kj_mol_nm2": str(kval),
                    "rep": f"{rep:03d}",
                    "status": existing.get(U, {}).get("status", "PENDING"),
                }
                rows[U] = row

# mark retired
for U in existing:
    if U not in rows:
        e = existing[U]; e["status"] = "RETIRED"; rows[U] = e

with open(MAN, "w", newline="") as f:
    w = csv.DictWriter(f, fieldnames=["uid","system","variant","speed_nm_per_ns","k_kj_mol_nm2","rep","status"])
    w.writeheader(); w.writerows(rows.values())

print(f"Wrote {MAN} ({len(rows)} rows). Wrote {SYS}.")

