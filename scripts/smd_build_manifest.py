#!/usr/bin/env python3
import os, csv, yaml, sys

CFG = "config.yaml"
STARTS_DIR = "manifests/starts"
OUT = "manifests/smd_manifest.csv"

os.makedirs("manifests", exist_ok=True)

cfg = yaml.safe_load(open(CFG))
G   = cfg["globals"]
dt_ps       = float(G["dt_ps"])
k_kj        = float(G["k_kj_mol_nm2"])
target_ext  = float(G["target_extension_nm"])
axis        = cfg["globals"]["smd"]["axis"]
array_cap   = int(G["slurm"]["array_cap"])
default_speeds    = G.get("pull_speeds_nm_per_ns", [])
n_reps_default    = int(G.get("n_reps_default", 1))

FIELDNAMES = [
    "system","variant","replicate","speed_nm_per_ns",
    "k_kj_mol_nm2","dt_ps","target_extension_nm","axis",
    "perf_ns_per_day","start_time_ps","final_tpr","final_xtc",
    "start_id","anchor_chain","array_cap"
]

rows = []
for sysent in cfg["systems"]:
    system = sysent["name"]
    perf   = float(sysent.get("perf_ns_per_day", 50.0))
    for var in sysent.get("variants", []):
        variant = var["id"]
        anchor_chain = var["anchor"]["chain"]
        speeds = var.get("speeds", default_speeds)
        n_reps = int(var.get("n_reps", n_reps_default))

        starts_csv = os.path.join(STARTS_DIR, f"{system}__{variant}.csv")
        if not os.path.isfile(starts_csv):
            print(f"[smd-manifest] ERROR: missing starts CSV {starts_csv}. Run `make posteq-sample` first.", file=sys.stderr)
            sys.exit(2)

        starts = []
        with open(starts_csv) as f:
            r = csv.DictReader(f)
            starts = list(r)
        if not starts:
            print(f"[smd-manifest] ERROR: no starts in {starts_csv}", file=sys.stderr)
            sys.exit(2)

        # Option A: one start per replicate (cycle if fewer starts)
        assigned = [starts[(rep-1) % len(starts)] for rep in range(1, n_reps+1)]

        for speed in speeds:
            for rep in range(1, n_reps+1):
                st = assigned[rep-1]
                rows.append({
                    "system": system,
                    "variant": variant,
                    "replicate": rep,
                    "speed_nm_per_ns": speed,
                    "k_kj_mol_nm2": k_kj,
                    "dt_ps": dt_ps,
                    "target_extension_nm": target_ext,
                    "axis": axis,
                    "perf_ns_per_day": perf,
                    "start_time_ps": st["start_time_ps"],
                    "final_tpr": st["final_tpr"],
                    "final_xtc": st["final_xtc"],
                    "start_id": st["start_id"],
                    "anchor_chain": anchor_chain,
                    "array_cap": array_cap
                })

with open(OUT, "w", newline="") as f:
    w = csv.DictWriter(f, fieldnames=FIELDNAMES)
    w.writeheader()
    w.writerows(rows)

print(f"[smd-manifest] Wrote {OUT} ({len(rows)} rows)")
